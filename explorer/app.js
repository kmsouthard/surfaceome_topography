/* The explorer reads what code/explorer/build_explorer.py wrote under data/ and computes only what
   follows from a choice made on the page: an epitope's height from its residues, a bridge's gap
   from that height, and which proteins a gap excludes. The rules are the pipeline's
   (surfaceomeTopography.epitope and .interface); the pairs' shares of a contact are read, not
   recomputed. */
"use strict";

const S = { shards: {}, contacts: {}, cartoons: {}, sample: {}, draw: 0, sceneMode: "varied", extra: [] };
const $ = (id) => document.getElementById(id);
const esc = (s) => String(s ?? "").replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const fmt = (x, d = 1) => (x == null || Number.isNaN(x) ? "–" : Number(x).toFixed(d));
const pct = (x) => (x >= 99.95 ? "100" : x < 0.05 ? "0" : x < 10 ? x.toFixed(1) : x.toFixed(0)) + "%";
const KIND = { alphafold: ["AlphaFold model", "--model"], disorder: ["Disordered chain", "--disorder"],
  domain: ["Counted domain", "--domain"], sequence: ["No structure", "--sequence"] };

async function load(name) {
  const r = await fetch("data/" + name);
  if (!r.ok) throw new Error(name + ": " + r.status);
  return r.json();
}

/* ---------- the rules ---------- */

function parseResidues(text) {
  const out = new Set();
  for (const part of text.split(",").map((p) => p.trim()).filter(Boolean)) {
    const m = part.match(/^(\d+)(?:\s*[-–]\s*(\d+))?$/);
    if (!m) return null;
    const a = +m[1], b = +(m[2] ?? m[1]);
    if (b < a || b - a > 5000) return null;
    for (let r = a; r <= b; r++) out.add(r);
  }
  return out.size ? [...out].sort((x, y) => x - y) : null;
}

async function residueHeights(acc) {
  const i = [...acc].reduce((n, c) => n + c.charCodeAt(0), 0) % S.meta.shards;
  S.shards[i] ??= load(`residues/${String(i).padStart(2, "0")}.json`);
  return (await S.shards[i])[acc];
}

// An epitope's height: the mean of its residues' heights, with their range (epitope.epitope_height).
async function epitopeOn(acc, residues) {
  const rec = await residueHeights(acc);
  const hs = [];
  if (rec) for (const r of residues) { const h = rec[1][r - rec[0]]; if (h != null) hs.push(h / 100); }
  if (!hs.length) return null;
  return { height: hs.reduce((a, b) => a + b, 0) / hs.length, lowest: Math.min(...hs), highest: Math.max(...hs),
    outside: residues.length - hs.length };
}

// gap = max(epitope height + receptor + antibody, antigen height) (interface.antibody_bridge_heights)
const bridgeGap = (epitope, receptor, antigen) => Math.max(epitope + receptor + S.meta.antibody_length, antigen);

// fits beneath the gap, excluded when the margin taller, partly excluded between (interface.contact_exclusion)
const fit = (height, gap) => (height <= gap ? "fits" : height >= gap + S.meta.exclusion_margin ? "excluded" : "partly");
const FIT = { fits: "fits", partly: "partly excluded", excluded: "excluded" };
const status = (f) => `<span class="status ${f}">${FIT[f]}</span>`;

/* ---------- small pieces ---------- */

function tile(label, value, sub = "", hero = false) {
  return `<div class="tile"><div class="label">${esc(label)}</div><div class="value${hero ? " hero" : ""}">${value}</div>` +
    (sub ? `<div class="sub">${sub}</div>` : "") + "</div>";
}

function niceTicks(max, n = 5) {
  const raw = max / n, mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw);
  const out = [];
  for (let v = 0; v <= max * 1.0001; v += step) out.push(+v.toPrecision(12));
  return out;
}

const tip = $("tip");
document.addEventListener("mousemove", (e) => {
  const t = e.target.closest?.("[data-tip]");
  if (!t) { tip.hidden = true; return; }
  tip.textContent = t.dataset.tip;
  tip.hidden = false;
  const w = tip.offsetWidth, x = Math.min(e.clientX + 14, window.innerWidth - w - 8);
  tip.style.left = Math.max(8, x) + "px";
  tip.style.top = e.clientY + 16 + "px";
});
document.addEventListener("focusin", (e) => {
  const t = e.target.closest?.("[data-tip]");
  if (!t) return;
  const r = t.getBoundingClientRect();
  tip.textContent = t.dataset.tip; tip.hidden = false;
  tip.style.left = Math.max(8, r.left) + "px"; tip.style.top = r.bottom + 6 + "px";
});
document.addEventListener("focusout", () => { tip.hidden = true; });

/* ---------- charts ---------- */

// A histogram over the pipeline's log-spaced height bins. series: [{name, color, items: [[height, weight]]}]
function histogram(series, { ylabel, marker, unit = "%" }) {
  const edges = S.meta.height_bins, nb = edges.length - 1;
  const W = 640, H = 240, L = 44, R = 12, T = 30, B = 40, pw = W - L - R, ph = H - T - B;
  const lo = Math.log10(edges[0]), hi = Math.log10(edges[nb]);
  const x = (v) => L + ((Math.log10(Math.min(Math.max(v, edges[0]), edges[nb])) - lo) / (hi - lo)) * pw;
  const bins = series.map((s) => {
    const b = new Array(nb).fill(0);
    for (const [h, w] of s.items) {
      let i = edges.findIndex((e, k) => k < nb && h >= e && h < edges[k + 1]);
      if (i < 0) i = h < edges[0] ? 0 : nb - 1;
      b[i] += w;
    }
    return b;
  });
  const ticks = niceTicks(Math.max(...bins.flat(), 1e-9) * 1.05, 4), ymax = ticks[ticks.length - 1] || 1;
  const y = (v) => T + ph - (v / ymax) * ph;
  let g = "";
  for (const t of ticks) g += `<line x1="${L}" x2="${W - R}" y1="${y(t)}" y2="${y(t)}" stroke="var(--grid)"/>` +
    `<text class="axis-text" x="${L - 6}" y="${y(t) + 4}" text-anchor="end">${+t.toPrecision(3)}</text>`;
  for (const t of [1, 10, 100, 1000]) g += `<line x1="${x(t)}" x2="${x(t)}" y1="${T + ph}" y2="${T + ph + 4}" stroke="var(--axis)"/>` +
    `<text class="axis-text" x="${x(t)}" y="${T + ph + 16}" text-anchor="middle">${t}</text>`;
  g += `<line x1="${L}" x2="${W - R}" y1="${T + ph}" y2="${T + ph}" stroke="var(--axis)"/>` +
    `<text x="${L + pw / 2}" y="${H - 6}" text-anchor="middle">height, nm (log scale)</text>` +
    `<text x="${L - 34}" y="${T - 16}">${esc(ylabel)}</text>`;
  const n = series.length;
  for (let i = 0; i < nb; i++) {
    const x0 = x(edges[i]), x1 = x(edges[i + 1]), slot = x1 - x0, bw = Math.min(24, (slot - 2) / n - (n > 1 ? 1 : 0));
    const start = x0 + (slot - bw * n - (n - 1)) / 2;
    series.forEach((s, k) => {
      const v = bins[k][i];
      if (v <= 0) return;
      const top = y(v), h = T + ph - top, r = Math.min(3, h, bw / 2), bx = start + k * (bw + 1);
      g += `<path class="mark" tabindex="0" fill="var(${s.color})" d="M${bx},${T + ph}V${top + r}q0,${-r} ${r},${-r}H${bx + bw - r}q${r},0 ${r},${r}V${T + ph}Z" ` +
        `data-tip="${esc(s.name)}\n${fmt(edges[i], edges[i] < 10 ? 1 : 0)}–${fmt(edges[i + 1], edges[i + 1] < 10 ? 1 : 0)} nm: ${+v.toPrecision(3)}${unit}"/>`;
    });
  }
  if (marker) {
    const mx = x(marker.value), anchor = mx > W - 150 ? "end" : "start", dx = anchor === "end" ? -6 : 6;
    g += `<line x1="${mx}" x2="${mx}" y1="${T - 12}" y2="${T + ph}" stroke="var(--ink)" stroke-width="1.5"/>` +
      `<text class="strong" x="${mx + dx}" y="${T - 4}" text-anchor="${anchor}">${esc(marker.label)}</text>`;
  }
  const legend = n > 1 ? `<div class="legend">${series.map((s) => `<span><span class="swatch" style="background:var(${s.color})"></span>${esc(s.name)}</span>`).join("")}</div>` : "";
  return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(ylabel)} by height">${g}</svg>${legend}`;
}

// A protein's ectodomain drawn to scale on the membrane, its segments stacked as the pipeline stacks them.
function stackFigure(p, epitope) {
  const order = p.anchor === "first" ? [...p.stack].reverse() : p.stack;
  const total = Math.max(p.height, 0.01);
  const W = 320, H = 330, T = 16, base = 288, ph = base - T, cx = 62, cw = 52;
  const y = (h) => base - (h / total) * ph;
  const ticks = niceTicks(total, 5);
  let g = "";
  for (const t of ticks) g += `<line x1="${cx - 22}" x2="${cx - 16}" y1="${y(t)}" y2="${y(t)}" stroke="var(--axis)"/>` +
    `<text class="axis-text" x="${cx - 26}" y="${y(t) + 4}" text-anchor="end">${+t.toPrecision(3)}</text>`;
  g += `<line x1="${cx - 16}" x2="${cx - 16}" y1="${T}" y2="${base}" stroke="var(--axis)"/>` +
    `<text x="${cx - 46}" y="${T - 4}">nm</text>`;
  g += `<rect x="${cx - 16}" y="${base + 2}" width="${W - cx - 24}" height="14" rx="3" fill="var(--grid)"/>` +
    `<text class="axis-text" x="${W - 44}" y="${base + 30}" text-anchor="end">membrane</text>`;
  let below = 0, lastLabel = Infinity;
  for (const [kind, a, b, h] of order) {
    if (!h) continue;
    const y1 = y(below), y0 = y(below + h), hh = Math.max(y1 - y0 - 2, 1), [name, color] = KIND[kind];
    g += `<rect class="mark" tabindex="0" x="${cx}" y="${y0 + 1}" width="${cw}" height="${hh}" rx="${Math.min(4, hh / 2)}" fill="var(${color})" ` +
      `data-tip="${name}\nresidues ${a}–${b}\n${fmt(h, 2)} nm, from ${fmt(below, 1)} to ${fmt(below + h, 1)} nm"/>`;
    const mid = (y0 + y1) / 2;
    if (hh >= 12 && lastLabel - mid >= 13) {
      g += `<text x="${cx + cw + 10}" y="${mid + 4}">${name}, ${fmt(h, 1)} nm</text>`;
      lastLabel = mid;
    }
    below += h;
  }
  if (epitope) {
    const ex = cx - 8, yl = y(epitope.lowest), yh = y(epitope.highest), ym = y(epitope.height);
    g += `<line x1="${cx}" x2="${cx + cw}" y1="${ym}" y2="${ym}" stroke="var(--ink)" stroke-width="1.5"/>` +
      `<line x1="${ex}" x2="${ex}" y1="${yl}" y2="${yh}" stroke="var(--ink)" stroke-width="2" stroke-linecap="round"/>` +
      `<circle cx="${ex}" cy="${ym}" r="4.5" fill="var(--ink)" stroke="var(--surface)" stroke-width="2" ` +
      `data-tip="epitope at ${fmt(epitope.height, 1)} nm\nits residues span ${fmt(epitope.lowest, 1)} to ${fmt(epitope.highest, 1)} nm"/>`;
  }
  return `<svg class="stack-figure" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(p.gene)} ectodomain, ${fmt(p.height, 1)} nm, drawn to scale">${g}</svg>`;
}

// A contact as nested circles: each ring's area is the share of the contact held in one band of gaps.
function bullseye(rows) {
  const edges = [0, ...S.meta.bands, Infinity], total = rows.reduce((n, r) => n + r.share, 0) || 1;
  const bands = edges.slice(0, -1).map((lo, i) => {
    const hi = edges[i + 1], inside = rows.filter((r) => r.gap >= lo && r.gap < hi);
    const share = inside.reduce((n, r) => n + r.share, 0) / total;
    const top = inside.reduce((best, r) => (!best || r.share > best.share ? r : best), null);
    return { lo, hi, share, top, label: hi === Infinity ? `${lo}+ nm` : lo === 0 ? `under ${hi} nm` : `${lo}–${hi} nm`,
      topShare: top && share ? top.share / total / share : 0 };
  });
  const R = 120, C = 130;
  let g = "", cum = 1;
  for (let i = bands.length - 1; i >= 0; i--) {
    const b = bands[i];
    if (cum > 1e-9) g += `<circle class="mark" tabindex="0" cx="${C}" cy="${C}" r="${Math.sqrt(cum) * R}" fill="var(--band-${i})" stroke="var(${i ? "--surface" : "--muted"})" stroke-width="${i ? 2 : 1}" ` +
      `data-tip="${b.label}: ${pct(b.share * 100)} of the contact${b.top ? `\nmost of it: ${b.top.name} (${pct(b.topShare * 100)})` : ""}"/>`;
    cum -= b.share;
  }
  const legend = bands.map((b, i) => `<span><span class="swatch${i ? "" : " open"}" style="background:var(--band-${i})"></span><strong>${b.label}</strong> · ` +
    `${pct(b.share * 100)} of the contact${b.top ? ` · ${esc(b.top.name)} (${pct(b.topShare * 100)})` : ""}</span>`).reverse().join("");
  return `<div class="grid" style="margin-top:0;grid-template-columns:minmax(200px,260px) 1fr;align-items:center">` +
    `<svg viewBox="0 0 260 260" style="max-width:260px;margin:0 auto" role="img" aria-label="Share of the contact by gap height">${g}</svg>` +
    `<div class="legend column">${legend}</div></div>`;
}

// Where a protein of a given height fits across a contact: percent of it fitting, partly excluded, excluded.
function placed(height, rows) {
  const total = rows.reduce((n, r) => n + r.share, 0) || 1, out = { fits: 0, partly: 0, excluded: 0 };
  for (const r of rows) out[fit(height, r.gap)] += (r.share / total) * 100;
  return out;
}

function placedBar(name, height, rows) {
  const p = placed(height, rows), W = 300, colors = { fits: "--good", partly: "--warning", excluded: "--critical" };
  let x = 0, g = "";
  for (const k of ["fits", "partly", "excluded"]) {
    const w = (p[k] / 100) * W;
    if (w > 0.5) g += `<rect class="mark" tabindex="0" x="${x + 1}" y="2" width="${Math.max(w - 2, 1)}" height="12" rx="3" fill="var(${colors[k]})" ` +
      `data-tip="${esc(name)} ${FIT[k]} over ${pct(p[k])} of the contact"/>`;
    x += w;
  }
  return `<tr><td><strong>${esc(name)}</strong><br><span class="muted">${fmt(height, 1)} nm</span></td>` +
    `<td style="min-width:140px"><svg viewBox="0 0 ${W} 16" preserveAspectRatio="none" style="width:100%;height:16px" role="img" ` +
    `aria-label="${esc(name)}: fits ${pct(p.fits)}, partly excluded ${pct(p.partly)}, excluded ${pct(p.excluded)}">${g}</svg></td>` +
    `<td class="num">${pct(p.fits)}</td><td class="num">${pct(p.partly)}</td><td class="num">${pct(p.excluded)}</td></tr>`;
}

/* ---------- CellScape cartoons ---------- */

const RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281"];
const CAPSULE_R = 12;   // angstroms, as in the README's illustrations

async function cartoonOf(acc) {
  const i = [...acc].reduce((n, c) => (n * 31 + c.charCodeAt(0)) % 65521, 7) % 256;
  S.cartoons[i] ??= load(`cartoons/${String(i).padStart(3, "0")}.json`).catch(() => ({}));
  return (await S.cartoons[i])[acc] || null;
}

// A protein's colour by its height, on the README's ramp (2 to 60 nm, log scale).
function heightColor(nm) {
  const f = (Math.log(Math.max(nm, 0.1)) - Math.log(2)) / (Math.log(60) - Math.log(2));
  return RAMP[Math.round(Math.min(Math.max(f, 0), 1) * (RAMP.length - 1))];
}

// One protein as a drawable: its width and height in angstroms, and svg(x, base, k), which draws it with
// its left edge at x and its foot at base, k pixels to the angstrom. Without a cartoon it is a capsule.
function shape(p, c, color) {
  const cap = (cx, y0, len, x, base, k) => {
    const r = Math.min(CAPSULE_R, len / 2);
    return `<rect x="${x + (cx - r) * k}" y="${base - (y0 + len) * k}" width="${2 * r * k}" height="${len * k}" rx="${r * k}" ` +
      `fill="var(--sequence)" stroke="var(--ink)" stroke-opacity="0.45" stroke-width="0.6"/>`;
  };
  if (!c) {
    const h = Math.max(p.height * 10, 2);
    return { w: 2 * CAPSULE_R, h, svg: (x, base, k) => cap(CAPSULE_R, 0, h, x, base, k) };
  }
  return { w: c.w, h: c.h, svg: (x, base, k) => {
    let g = "";
    if (c.before > 0) g += cap(c.bx, 0, c.before + CAPSULE_R, x, base, k);
    if (c.after > 0) g += cap(c.tx, c.before + c.structure - CAPSULE_R, c.after + CAPSULE_R, x, base, k);
    for (const [shade, rings] of c.polys) {
      let d = "";
      for (const ring of rings) {
        for (let i = 0; i < ring.length; i += 2) d += (i ? "L" : "M") + (x + ring[i] * k).toFixed(1) + " " + (base - ring[i + 1] * k).toFixed(1);
        d += "Z";
      }
      g += `<path d="${d}" fill-rule="evenodd" fill="color-mix(in srgb, ${color} ${Math.round(58 + 42 * shade)}%, #000)" ` +
        `stroke="var(--ink)" stroke-opacity="0.5" stroke-width="0.5"/>`;
    }
    return g;
  } };
}

function cartoonFigure(p, c) {
  const sh = shape(p, c, "var(--model)"), W = 300, H = 330, base = 290, k = Math.min(262 / sh.h, 220 / sh.w, 6);
  const x = (W - sh.w * k) / 2;
  const bar = [50, 20, 10, 5, 2, 1, 0.5].find((nm) => nm * 10 * k <= 90) || 0.5;
  return `<svg class="stack-figure" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(p.gene)} drawn with CellScape">` +
    `<rect x="16" y="${base + 2}" width="${W - 32}" height="14" rx="3" fill="var(--grid)"/>` + sh.svg(x, base, k) +
    `<line x1="20" x2="${20 + bar * 10 * k}" y1="${H - 8}" y2="${H - 8}" stroke="var(--ink-2)" stroke-width="2" stroke-linecap="round"/>` +
    `<text x="${26 + bar * 10 * k}" y="${H - 4}">${bar} nm</text></svg>` +
    `<p class="muted small" style="margin:8px 0 0">${c ? "Its AlphaFold model outlined with CellScape and oriented by its topology" +
      (c.before + c.after > 0 ? "; grey capsules are height without a structure" : "") + ". The drawn extent is illustrative: the height is the number above."
      : "No AlphaFold model covers this ectodomain, so its whole height is drawn as a capsule."}</p>`;
}

// How a scene's proteins are coloured: each its own hue, by functional class, or by height on one ramp.
const HUES = ["--c1", "--c2", "--c3", "--c4", "--c5", "--c6", "--c7", "--c8"];
const CLASSES = [["Receptors", "--c1"], ["Transporters", "--c2"], ["Enzymes", "--c3"], ["Miscellaneous", "--c5"], ["Unclassified", "--sequence"]];
const classOf = (p) => (CLASSES.some(([name]) => name === p.family) ? p.family : "Unclassified");

// Proteins drawn from a surface by abundance, side by side on a membrane at one scale.
async function sceneFigure(items, seed, mode) {
  let state = [...String(seed)].reduce((n, ch) => (n * 31 + ch.charCodeAt(0)) >>> 0, 7) || 1;
  const random = () => { state ^= state << 13; state >>>= 0; state ^= state >>> 17; state ^= state << 5; state >>>= 0; return state / 4294967296; };
  const total = items.reduce((n, [, w]) => n + w, 0), picks = [];
  for (let i = 0; i < 24; i++) {
    let r = random() * total, j = 0;
    while (j < items.length - 1 && (r -= items[j][1]) > 0) j++;
    picks.push(items[j][0]);
  }
  const distinct = [...new Set(picks)];
  const cartoons = Object.fromEntries(await Promise.all(distinct.map(async (q) => [q.acc, await cartoonOf(q.acc)])));
  const color = (q) => mode === "height" ? heightColor(q.height)
    : mode === "class" ? `var(${CLASSES.find(([name]) => name === classOf(q))[1]})`
    : `var(${HUES[distinct.indexOf(q) % HUES.length]})`;
  const objs = picks.map((q) => ({ p: q, ...shape(q, cartoons[q.acc], color(q)) }));
  const tallest = [...new Set(objs.map((o) => o.h))].sort((a, b) => b - a);
  const limit = tallest.length > 1 && tallest[0] > 2 * tallest[1] ? tallest[1] + 140 : Infinity;   // break a protein twice the next tallest
  // proteins stand 3 angstroms apart; a narrow one keeps room for its label
  const pad = 3, VW = 900, slot = (o, k) => Math.max(o.w, 13 / k);
  let k = VW / (objs.reduce((n, o) => n + o.w, 0) + pad * (objs.length + 1));
  for (let i = 0; i < 4; i++) k = VW / (objs.reduce((n, o) => n + slot(o, k), 0) + pad * (objs.length + 1));
  const base = Math.min(tallest[0], limit) * k + 24, label = Math.max(...picks.map((q) => q.gene.length)) * 6.6 + 10;
  const H = base + 40 * k + label + 4;
  let g = `<rect x="0" y="${base}" width="${VW}" height="${40 * k}" rx="4" fill="var(--grid)"/>`, x = pad;
  objs.forEach((o, i) => {
    const w = slot(o, k), cut = o.h > limit, px = (x + (w - o.w) / 2) * k, cx = (x + w / 2) * k;
    g += `<g data-tip="${esc(o.p.gene)}\n${fmt(o.p.height, 1)} nm · ${esc(classOf(o.p).toLowerCase())}" ${cut ? `clip-path="url(#cut${i})"` : ""}>` +
      (cut ? `<clipPath id="cut${i}"><rect x="${px - 4}" y="${base - limit * k}" width="${o.w * k + 8}" height="${limit * k + 2}"/></clipPath>` : "") +
      `<rect class="hit" x="${px}" y="${base - Math.min(o.h, limit) * k}" width="${o.w * k}" height="${Math.min(o.h, limit) * k}"/>` + o.svg(px, base, k) + `</g>`;
    if (cut) g += `<text class="strong" x="${cx}" y="${base - limit * k - 6}" text-anchor="middle">${fmt(o.p.height, 0)} nm</text>`;
    g += `<text transform="translate(${cx + 4} ${base + 40 * k + 7}) rotate(90)">${esc(o.p.gene)}</text>`;
    x += w + pad;
  });
  const bar = [50, 20, 10, 5].find((nm) => nm * 10 * k <= 120) || 5;
  g += `<line x1="${VW - bar * 10 * k - 4}" x2="${VW - 4}" y1="8" y2="8" stroke="var(--ink-2)" stroke-width="2" stroke-linecap="round"/>` +
    `<text x="${VW - bar * 10 * k - 10}" y="12" text-anchor="end">${bar} nm</text>`;
  const legend = mode === "class" ? `<div class="legend">${CLASSES.filter(([name]) => picks.some((q) => classOf(q) === name))
      .map(([name, c]) => `<span><span class="swatch" style="background:var(${c})"></span>${name}</span>`).join("")}` +
      `<span class="muted">Functional class from Almén et al. 2009, as carried in the SURFY table</span></div>`
    : mode === "height" ? `<div class="legend"><span>shorter</span><span>${RAMP.map((c) => `<span class="swatch" style="background:${c};margin-right:1px"></span>`).join("")}</span><span>taller (2 to 60 nm)</span></div>` : "";
  return `<div class="scroll"><svg viewBox="0 0 ${VW} ${H}" style="min-width:560px" role="img" aria-label="Proteins sampled from the surface by abundance, drawn to scale">${g}</svg></div>${legend}`;
}

/* ---------- protein ---------- */

function findProtein(text) {
  const q = text.trim().split(/\s+/)[0].toUpperCase();
  return S.byGene[q] || S.byAcc[q] || S.byEntry[q] || null;
}

function bridgeTable(p, epitope) {
  const probes = S.meta.probes;
  const rows = S.meta.receptors.map((r) => {
    const gap = bridgeGap(epitope.height, r.height, p.height);
    const lo = bridgeGap(epitope.lowest, r.height, p.height), hi = bridgeGap(epitope.highest, r.height, p.height);
    return `<tr><td>${esc(r.gene)}</td><td class="num">${fmt(r.height, 1)}</td><td class="num"><strong>${fmt(gap, 1)}</strong></td>` +
      `<td class="num muted">${fmt(lo, 1)}–${fmt(hi, 1)}</td><td>${gap > p.height ? "epitope" : "antigen height"}</td>` +
      probes.map((q) => `<td>${status(fit(q.height, gap))}</td>`).join("") + "</tr>";
  }).join("");
  return `<div class="scroll"><table><thead><tr><th>Fc receptor</th><th class="num">height, nm</th><th class="num">gap, nm</th>` +
    `<th class="num">range</th><th>set by</th>${probes.map((q) => `<th>${esc(q.name)} <span class="muted">${fmt(q.height, 0)} nm</span></th>`).join("")}</tr></thead>` +
    `<tbody>${rows}</tbody></table></div>`;
}

async function renderProtein(p) {
  const body = $("protein-body");
  if (!p) { body.innerHTML = `<p class="muted">No protein of the height table matches.</p>`; return; }
  $("protein-search").value = p.gene;
  const cartoon = await cartoonOf(p.acc);
  const curated = S.meta.antibodies.find((a) => a.acc === p.acc);
  const taller = S.proteins.filter((q) => q.height < p.height).length / S.proteins.length * 100;
  const partners = S.pairs.filter((x) => x[0] === p.acc || x[1] === p.acc)
    .map((x) => ({ acc: x[0] === p.acc ? x[1] : x[0], gap: x[4], source: x[5] })).sort((a, b) => a.gap - b.gap);
  const methods = (p.methods || "").split(",").map((m) => ({ alphafold: "AlphaFold model", disorder: "disorder model",
    domain: "domain counting", sequence: "sequence length" }[m] || m)).join(", ");
  const held = { first: "anchored at its C-terminal end", last: "anchored at its N-terminal end", loop: "a loop anchored at both ends" }[p.anchor] || "";

  body.innerHTML =
    `<div class="title-row"><h2>${esc(p.gene)}</h2><span class="ink-2">${esc((p.name || "").split(" (")[0])}</span>` +
    `<a class="small" href="https://www.uniprot.org/uniprotkb/${esc(p.acc)}" target="_blank" rel="noopener">${esc(p.acc)}</a></div>` +
    `<div class="chips">${[p.location, p.tm != null ? `${p.tm} transmembrane segment${p.tm === 1 ? "" : "s"}` : null, p.cd,
      p.surfy ? `SURFY: ${p.surfy}` : null, p.family && p.family !== "Unclassified" ? p.family : null]
      .filter(Boolean).map((c) => `<span class="chip">${esc(c)}</span>`).join("")}</div>` +
    `<p class="muted small" style="margin:6px 0 0">Location and topology from UniProt; surface call from SURFY (Bausch-Fluck et al. 2018); functional class from Almén et al. 2009.</p>` +
    `<div class="tiles">${tile("Height", fmt(p.height, 1) + " nm", "an upper bound: fully extended", true)}` +
    tile("Percentile", `taller than ${pct(taller)}`, `of ${S.proteins.length.toLocaleString()} surfaceome proteins`) +
    tile("Ectodomain", `${p.ecd[0]}–${p.ecd[1]}`, `${p.ecd[1] - p.ecd[0]} residues${held ? ", " + held : ""}`) +
    tile("Method", `<span style="font-size:16px">${esc(methods || "–")}</span>`) + `</div>` +
    p.notes.map((n) => `<p class="note">${esc(n[0].toUpperCase() + n.slice(1))}.</p>`).join("") +
    `<div class="grid"><div class="card"><h3>${esc(p.gene)} on the membrane</h3>${cartoonFigure(p, cartoon)}</div>` +
    `<div class="card"><h3>Height by segment</h3>${stackFigure(p, null)}${stackLegend(p)}` +
    `<p class="small" style="margin:10px 0 0"><a href="#epitope/${esc(p.gene)}">Predict an antibody epitope on ${esc(p.gene)}</a>` +
    (curated ? ` <span class="muted">(${esc(curated.antibody.split(" (")[0])} targets it)</span>` : "") + `</p></div></div>` +
    `<div class="card"><h3>Height relative to the surfaceome</h3>` +
    histogram([{ name: "Surfaceome proteins", color: "--series-1", items: S.proteins.map((q) => [q.height, 1]) }],
      { ylabel: "proteins", unit: "", marker: { value: p.height, label: `${p.gene} ${fmt(p.height, 1)} nm` } }) + `</div>` +
    `<div class="card"><h3>Trans binding partners (${partners.length})</h3>` + (partners.length
      ? `<div class="scroll"><table><thead><tr><th>Partner</th><th class="num">partner height, nm</th><th class="num">gap, nm</th><th>source</th></tr></thead><tbody>` +
        partners.map((x) => { const q = S.byAcc[x.acc]; return `<tr><td>${q ? `<a href="#protein/${esc(q.gene)}">${esc(q.gene)}</a>` : esc(x.acc)}</td>` +
          `<td class="num">${q ? fmt(q.height, 1) : "–"}</td><td class="num">${fmt(x.gap, 1)}</td><td>${esc(x.source)}</td></tr>`; }).join("") + `</tbody></table></div>`
      : `<p class="muted small">No trans interaction is recorded for this protein.</p>`) + `</div>`;
}

function stackLegend(p) {
  return `<div class="legend">${Object.entries(KIND).filter(([k]) => p.stack.some((s) => s[0] === k && s[3]))
    .map(([, [n, c]]) => `<span><span class="swatch" style="background:var(${c})"></span>${n}</span>`).join("")}</div>`;
}

/* ---------- epitope ---------- */

const pubmed = (pmid) => `<a href="https://pubmed.ncbi.nlm.nih.gov/${esc(pmid)}/" target="_blank" rel="noopener">PMID ${esc(pmid)}</a>`;

// Where a curated epitope comes from, as a sentence: a structure, peptide mapping, or a domain.
function epitopeSource(ab) {
  const drug = esc(ab.antibody.split(" (")[0]), gene = esc(ab.gene);
  if (ab.basis === "structure") return `<span class="tag">curated · structure</span> Residues within 4 Å of ${drug} in its structure with ${gene} (PDB ${esc(ab.pdb)}, ${pubmed(ab.pmid)}).`;
  if (ab.basis === "peptide mapping") return `<span class="tag">curated · peptide mapping</span> The ${gene} peptides ${drug} binds (${pubmed(ab.pmid)}). No structure shows the contact residues, so the epitope is no finer than the peptides.`;
  if (ab.basis === "domain") return `<span class="tag">curated · domain</span> The ${gene} domain ${drug} is reported to bind (${pubmed(ab.pmid)}). The whole domain stands in for the epitope, so its height is the domain's middle.`;
  return `<span class="tag warn">no epitope known</span> Nothing shows where ${drug} binds ${gene}. The pipeline places its epitope at the membrane, the smallest gap its bridge can hold.`;
}

// The curated antibodies at a glance: where each binds, and the gap its bridge to CD16 would hold.
function antibodyTable(current) {
  const cd16 = S.meta.receptors.find((r) => r.gene === "FCGR3A") || S.meta.receptors[0];
  const ro = S.meta.probes.find((q) => q.name === "CD45RO");
  const rows = S.meta.antibodies.map((a) => {
    const antigen = S.byAcc[a.acc], known = a.residues != null, gap = bridgeGap(a.height, cd16.height, antigen.height);
    return `<tr${a.acc === current ? ' class="current"' : ""}><td><a href="#epitope/${esc(antigen.gene)}">${esc(a.antibody.split(" (")[0])}</a></td>` +
      `<td><a href="#protein/${esc(antigen.gene)}">${esc(antigen.gene)}</a></td><td class="num">${fmt(antigen.height, 1)}</td>` +
      `<td>${a.basis === "structure" ? `structure, PDB ${esc(a.pdb)}` : known ? esc(a.basis) : `<span class="muted">none</span>`}</td>` +
      `<td class="residues">${known ? esc(a.residues).replace(/,/g, ", ") : `<span class="muted">–</span>`}</td>` +
      `<td class="num">${known ? fmt(a.height, 1) : `<span class="muted">0, assumed</span>`}</td>` +
      `<td>${a.height <= S.meta.phagocytosis_range ? "within" : "beyond"}</td><td class="num">${fmt(gap, 1)}</td>` +
      `<td>${ro ? status(fit(ro.height, gap)) : ""}</td></tr>`;
  }).join("");
  return `<div class="scroll"><table><thead><tr><th>Antibody</th><th>Antigen</th><th class="num">antigen height, nm</th><th>Basis</th><th>Epitope residues</th>` +
    `<th class="num">epitope height, nm</th><th>10 nm range</th><th class="num">gap at ${esc(cd16.gene)}, nm</th><th>CD45RO</th></tr></thead><tbody>${rows}</tbody></table></div>`;
}

async function renderEpitope(p, residuesText) {
  const body = $("epitope-body");
  if (!p) { body.innerHTML = `<p class="muted">No protein of the height table matches.</p>`; return; }
  const ab = S.meta.antibodies.find((a) => a.acc === p.acc);
  const text = (residuesText ?? ab?.residues ?? "").replace(/\s/g, "");
  $("epitope-antigen").value = p.gene;
  $("epitope-antibody").value = ab ? p.acc : "";
  $("epitope-residues").value = text;
  $("epitope-residues").placeholder = `e.g. ${p.ecd[0] + 10}-${Math.min(p.ecd[0] + 40, p.ecd[1])}`;
  const residues = text ? parseResidues(text) : null;
  const epitope = residues ? await epitopeOn(p.acc, residues) : null;
  const isCurated = ab?.residues && text === ab.residues;

  let source = "", result = "";
  if (isCurated) source = epitopeSource(ab);
  else if (text) source = `<span class="tag warn">uncurated</span> Residues as entered. Nothing checks that an antibody binds them.`;
  else if (ab) source = epitopeSource(ab) + " Enter residues to test another position.";

  if (!text) result = `<p class="muted">Enter the residues an antibody binds on ${esc(p.gene)} (its ectodomain is ${p.ecd[0]}–${p.ecd[1]}) to predict the epitope's height, ` +
    `the membrane gap of the antibody bridge to each Fc receptor, and whether the phosphatases are excluded from that gap.</p>`;
  else if (!residues) result = `<p class="note">Write residues as ranges and single positions, like 579-625,630.</p>`;
  else if (!epitope) result = `<p class="note">None of those residues lies in the ectodomain of ${esc(p.gene)} (${p.ecd[0]}–${p.ecd[1]}).</p>`;
  else {
    const within = epitope.height <= S.meta.phagocytosis_range;
    result = `<div class="tiles">${tile("Epitope height", fmt(epitope.height, 1) + " nm", `its residues span ${fmt(epitope.lowest, 1)}–${fmt(epitope.highest, 1)} nm`, true)}` +
      tile("Antigen height", fmt(p.height, 1) + " nm", `<a href="#protein/${esc(p.gene)}">${esc(p.gene)}</a>, fully extended`) +
      tile("Distance from the target membrane", within ? "within 10 nm" : "beyond 10 nm", "phagocytosis falls beyond 10 nm (Bakalar et al. 2018)") + `</div>` +
      (epitope.outside ? `<p class="note">${epitope.outside} of the residues lie outside the ectodomain and were left out.</p>` : "") +
      p.notes.map((n) => `<p class="note">${esc(n[0].toUpperCase() + n.slice(1))}.</p>`).join("") +
      `<div class="grid" style="grid-template-columns:minmax(240px,320px) minmax(300px,1fr)"><div class="card"><h3>Position on ${esc(p.gene)}</h3>${stackFigure(p, epitope)}${stackLegend(p)}` +
      `<p class="muted small" style="margin:8px 0 0">The line marks the epitope's mean height; the bar beside it, the span of its residues.</p></div>` +
      `<div class="card"><h3>Predicted gap of the antibody bridge, by Fc receptor</h3>` + bridgeTable(p, epitope) +
      `<p class="muted small">Gap = the taller of (epitope height + Fc receptor + ${S.meta.antibody_length} nm of antibody) and the antigen itself, which must fit beneath. ` +
      `A phosphatase fits when no taller than the gap and is excluded when ${S.meta.exclusion_margin} nm or more taller. The range is the gap at the epitope's lowest and highest residue.</p></div></div>`;
  }
  body.innerHTML = (source ? `<p class="small ink-2" style="margin:10px 0 0">${source}</p>` : "") + result +
    `<div class="card"><h3>Curated antibody epitopes</h3>${antibodyTable(p.acc)}` +
    `<p class="muted small" style="margin:8px 0 0">Ten epitopes are read from a structure of the antibody on its antigen. Catumaxomab's is from peptide mapping and elotuzumab's is the domain it is reported to bind, ` +
    `so both are coarser. Tafasitamab and olaratumab have no epitope recorded; the pipeline places them at the membrane, the smallest gap their bridge can hold.</p></div>`;
}

/* ---------- cell surface ---------- */

function surfaceItems(list) {
  return list.map(([acc, w]) => [S.byAcc[acc], w]).filter(([p]) => p && p.height != null);
}

function surfaceStats(items) {
  const total = items.reduce((n, [, w]) => n + w, 0) || 1;
  return { n: items.length, mean: items.reduce((n, [p, w]) => n + p.height * w, 0) / total,
    tall: items.reduce((n, [p, w]) => n + (p.height >= 30 ? w : 0), 0) / total * 100, total };
}

function parseUpload(text) {
  const out = new Map();
  for (const line of text.split(/\r?\n/)) {
    const cells = line.split(/[,\t;]/).map((c) => c.trim().replace(/^"|"$/g, ""));
    const p = findProtein(cells[0] || ""), v = parseFloat(cells[1]);
    if (p && v > 0) out.set(p.acc, (out.get(p.acc) || 0) + v);
  }
  return [...out.entries()];
}

function renderSurface() {
  const label = (id) => (id === "upload" ? S.upload.name : S.meta.surfaces.find((s) => s.id === id)?.label);
  const get = (id) => (id === "upload" ? S.upload.items : S.surfaces[id]);
  const ids = [$("surface-a").value, $("surface-b").value].filter((id) => id && get(id));
  if (ids[1] === ids[0]) ids.pop();
  const sets = ids.map((id, k) => {
    const items = surfaceItems(get(id)), stats = surfaceStats(items);
    return { id, name: label(id), color: `--series-${k + 1}`, items, stats };
  });
  if (!sets.length) { $("surface-body").innerHTML = ""; return; }
  const first = sets[0];
  const tiles = sets.length === 1
    ? tile("Abundance-weighted mean height", fmt(first.stats.mean, 1) + " nm", esc(first.name), true) +
      tile("Surface proteins detected", first.stats.n.toLocaleString(), "with a height estimate") +
      tile("Abundance at 30 nm or taller", pct(first.stats.tall))
    : sets.map((s) => tile("Mean height, " + s.name, fmt(s.stats.mean, 1) + " nm",
        `${s.stats.n} proteins · ${pct(s.stats.tall)} of abundance at 30 nm or taller`)).join("");
  const top = [...first.items].sort((a, b) => b[1] - a[1]).slice(0, 25);
  let cum = 0;
  const draw = ++S.draw;
  $("surface-body").innerHTML = `<div class="tiles">${tiles}</div>` +
    `<div class="card"><div class="title-row" style="margin:0 0 8px"><h3 style="margin:0">The ${esc(first.name)} surface, to scale</h3>` +
    `<span class="segmented" role="group" aria-label="Colour">${[["varied", "Colourful"], ["class", "By class"], ["height", "By height"]]
      .map(([m, t]) => `<button data-mode="${m}" aria-pressed="${S.sceneMode === m}">${t}</button>`).join("")}</span>` +
    `<button id="resample" class="quiet">Resample</button></div><div id="scene"><p class="muted small">Drawing…</p></div>` +
    `<p class="muted small" style="margin:8px 0 0">24 proteins sampled in proportion to abundance, so an abundant protein appears more than once. ` +
    `Structures are AlphaFold models outlined with CellScape; grey capsules are height with no structure. A protein more than twice as tall as the next is cut off, with its height printed.</p></div>` +
    `<div class="card"><h3>Height distribution, weighted by abundance</h3>` +
    histogram(sets.map((s) => ({ name: s.name, color: s.color, items: s.items.map(([p, w]) => [p.height, w / s.stats.total * 100]) })),
      { ylabel: "% of abundance" }) + `</div>` +
    `<div class="card"><h3>Most abundant surface proteins, ${esc(first.name)}</h3><div class="scroll"><table><thead><tr><th>Protein</th>` +
    `<th class="num">height, nm</th><th class="num">share of abundance</th><th class="num">cumulative</th></tr></thead><tbody>` +
    top.map(([p, w]) => { const s = w / first.stats.total * 100; cum += s; return `<tr><td><a href="#protein/${esc(p.gene)}">${esc(p.gene)}</a></td>` +
      `<td class="num">${fmt(p.height, 1)}</td><td class="num">${pct(s)}</td><td class="num">${pct(cum)}</td></tr>`; }).join("") +
    `</tbody></table></div></div>`;
  const scene = async () => {
    const svg = await sceneFigure(first.items, `${first.id}:${S.sample[first.id] ?? 0}`, S.sceneMode);
    if (draw === S.draw) $("scene").innerHTML = svg;
  };
  $("resample").addEventListener("click", () => { S.sample[first.id] = (S.sample[first.id] ?? 0) + 1; scene(); });
  for (const b of document.querySelectorAll(".segmented button")) b.addEventListener("click", () => {
    S.sceneMode = b.dataset.mode;
    for (const o of document.querySelectorAll(".segmented button")) o.setAttribute("aria-pressed", o === b);
    scene();
  });
  scene();
}

/* ---------- contact ---------- */

async function renderContact() {
  const a = $("contact-a").value, b = $("contact-b").value, id = `${a}__${b}`;
  const entry = S.meta.contacts.find((c) => c.id === id);
  const abSel = $("contact-ab"), wanted = abSel.value;
  abSel.innerHTML = `<option value="">none</option>` + S.meta.antibodies.filter((x) => entry.antigens.includes(x.acc))
    .map((x) => `<option value="${x.acc}">${esc(x.antibody.split(" (")[0])} (${esc(x.gene)})</option>`).join("");
  abSel.value = entry.antigens.includes(wanted) ? wanted : "";
  const antigen = abSel.value, ab = S.meta.antibodies.find((x) => x.acc === antigen);
  const input = $("contact-epitope");
  $("contact-epitope-label").hidden = !antigen;
  if (antigen && input.dataset.for !== antigen) { input.value = ab.residues || ""; input.dataset.for = antigen; }

  S.contacts[id] ??= load(`contacts/${id}.json`);
  const c = await S.contacts[id];
  let epitope = 0, source = "";
  if (antigen) {
    const text = input.value.trim(), residues = text ? parseResidues(text) : null;
    if (ab.residues && text.replace(/\s/g, "") === ab.residues) {
      epitope = ab.height;
      source = epitopeSource(ab) + ` The epitope stands at ${fmt(epitope, 1)} nm.` + (ab.notes ? ` ${esc(ab.notes[0].toUpperCase() + ab.notes.slice(1))}.` : "");
    } else if (residues && (await epitopeOn(antigen, residues))) {
      epitope = (await epitopeOn(antigen, residues)).height;
      source = `<span class="tag warn">uncurated</span> Residues as typed, standing at ${fmt(epitope, 1)} nm on ${esc(ab.gene)}; nothing checks that an antibody binds them.`;
    } else {
      source = `<span class="tag warn">at the membrane</span> ` + (text ? "Those residues are not in the ectodomain; " : `Nothing shows where ${esc(ab.antibody)} binds; `) +
        `the epitope is taken at the membrane, the least gap the bridge can hold.`;
    }
  }
  const rows = c.variants[antigen].map(([ia, ib, g, share, bridge, antigenH]) => ({
    a: c.a[ia], b: c.b[ib], share, bridge, gap: bridge ? Math.max(g - antigenH + epitope, antigenH) : g,
    name: `${c.a[ia][0]}–${c.b[ib][0]}`.replace(/_/g, " ") })).sort((x, y) => y.share - x.share);
  const total = rows.reduce((n, r) => n + r.share, 0) || 1;
  const mean = rows.reduce((n, r) => n + r.gap * r.share, 0) / total;
  const close = rows.reduce((n, r) => n + (r.gap < 20 ? r.share : 0), 0) / total * 100;
  const labels = Object.fromEntries(S.meta.surfaces.map((s) => [s.id, s.label]));
  const ro = S.meta.probes.find((q) => q.name === "CD45RO"), roOut = ro ? placed(ro.height, rows).excluded : 0;

  const unitTable = (units, side) => {
    const seen = [...units].filter((u) => u[3] > 0).sort((x, y) => y[3] - x[3]).slice(0, 14);
    return `<div class="scroll"><table><thead><tr><th>${esc(side)}</th><th class="num">height, nm</th><th class="num">share of surface</th>` +
      `<th class="num">fits in</th><th class="num">excluded from</th></tr></thead><tbody>` +
      seen.map((u) => { const p = placed(u[2], rows), q = S.byAcc[u[1]]; return `<tr><td>${q && q.gene === u[0] ? `<a href="#protein/${esc(q.gene)}">${esc(u[0])}</a>` : esc(u[0].replace(/_/g, " "))}</td>` +
        `<td class="num">${fmt(u[2], 1)}</td><td class="num">${pct(u[3])}</td><td class="num">${pct(p.fits)}</td><td class="num">${pct(p.excluded)}</td></tr>`; }).join("") +
      `</tbody></table></div>`;
  };

  $("contact-body").innerHTML =
    (source ? `<p class="small ink-2" style="margin:10px 0 0">${source}</p>` : "") +
    `<div class="tiles">${tile("Mean membrane gap", fmt(mean, 1) + " nm", "weighted by each interaction's share of the contact", true)}` +
    tile("Trans interactions", rows.length, antigen ? `${rows.filter((r) => r.bridge).length} through the antibody` : "no antibody") +
    tile("Contact with a gap under 20 nm", pct(close)) + tile("CD45RO excluded from", pct(roOut), "of the contact, by size") + `</div>` +
    `<div class="card"><h3>Membrane gap across the contact</h3>${bullseye(rows)}` +
    `<p class="muted small">Ring area is the share of the contact at that gap height. Each band lists its dominant interaction and that interaction's share of the band.</p></div>` +
    `<div class="card"><h3>Size-based exclusion from the contact</h3>` +
    `<p class="small ink-2" style="margin:0 0 10px;max-width:78ch">A protein taller than the local gap is pushed out of it; this is how kinetic segregation clears the phosphatases CD45 and CD148 ` +
    `from close contacts. Each bar divides the contact area by whether a protein of that height fits beneath the gap there, is within ${S.meta.exclusion_margin} nm of fitting, or is excluded.</p>` +
    `<div class="scroll"><table><thead><tr><th>Protein</th><th>Share of the contact area</th>` +
    `<th class="num">fits</th><th class="num">partly excluded</th><th class="num">excluded</th></tr></thead><tbody>` +
    [...S.meta.probes, ...S.extra].map((q) => placedBar(q.name, q.height, rows)).join("") + `</tbody></table></div>` +
    `<div class="legend"><span><span class="swatch" style="background:var(--good)"></span>fits: no taller than the gap</span>` +
    `<span><span class="swatch" style="background:var(--warning)"></span>partly excluded: up to ${S.meta.exclusion_margin} nm taller</span>` +
    `<span><span class="swatch" style="background:var(--critical)"></span>excluded: ${S.meta.exclusion_margin} nm or more taller</span></div>` +
    `<div class="inline-form" style="margin-top:12px"><label>Test another protein, or a height in nm` +
    `<input id="probe-input" list="protein-list" placeholder="e.g. SPN or 25" autocomplete="off"></label>` +
    (S.extra.length ? `<button id="probe-clear" class="quiet">Clear added</button>` : "") + `</div>` +
    `<p class="muted small" style="margin:8px 0 0">Heights are the extended estimates, applied whether or not the cell expresses the protein. CD45 is its longest isoform (RABC) and CD45RO its shortest. ` +
    `The ${S.meta.exclusion_margin} nm margin is from Schmid et al. 2016.</p></div>` +
    `<div class="card"><h3>Trans interactions, ranked by share of the contact</h3><div class="scroll"><table><thead><tr><th>${esc(labels[a])}</th><th>${esc(labels[b])}</th>` +
    `<th class="num">gap, nm</th><th class="num">share of the contact</th><th></th></tr></thead><tbody>` +
    rows.slice(0, 20).map((r) => `<tr><td>${esc(r.a[0].replace(/_/g, " "))}</td><td>${esc(r.b[0].replace(/_/g, " "))}</td><td class="num">${fmt(r.gap, 1)}</td>` +
      `<td class="num">${pct(r.share / total * 100)}</td><td>${r.bridge ? `<span class="tag">antibody bridge</span>` : ""}</td></tr>`).join("") +
    `</tbody></table></div></div>` +
    `<div class="grid"><div class="card"><h3>Exclusion of ${esc(labels[a])} surface proteins</h3>${unitTable(c.a, "most abundant first")}</div>` +
    `<div class="card"><h3>Exclusion of ${esc(labels[b])} surface proteins</h3>${unitTable(c.b, "most abundant first")}</div></div>` +
    `<p class="muted small" style="margin:8px 0 0">"Fits in" and "excluded from" are shares of the contact area, by each protein's height.</p>`;

  $("probe-input").addEventListener("change", (e) => {
    const text = e.target.value.trim(), height = Number(text), found = findProtein(text);
    if (text && Number.isFinite(height) && height > 0) S.extra.push({ name: `${height} nm`, height });
    else if (found && ![...S.meta.probes, ...S.extra].some((q) => q.name === found.gene)) S.extra.push({ name: found.gene, height: found.height });
    else return;
    renderContact();
  });
  $("probe-clear")?.addEventListener("click", () => { S.extra = []; renderContact(); });
}

/* ---------- navigation ---------- */

function show(view) {
  for (const b of document.querySelectorAll(".tabs button")) b.setAttribute("aria-selected", b.dataset.view === view);
  for (const v of document.querySelectorAll(".view")) v.hidden = v.id !== "view-" + view;
}

async function route() {
  const [view, ...rest] = decodeURIComponent(location.hash.slice(1)).split("/");
  if ((view || "protein") !== S.view) window.scrollTo(0, 0);
  S.view = view || "protein";
  if (view === "surface") {
    show("surface");
    if (rest[0] && S.surfaces[rest[0]]) $("surface-a").value = rest[0];
    if (rest[1] != null) $("surface-b").value = S.surfaces[rest[1]] ? rest[1] : "";
    renderSurface();
  } else if (view === "contact") {
    show("contact");
    if (rest[0]) $("contact-a").value = rest[0];
    if (rest[1]) $("contact-b").value = rest[1];
    if (rest[2] != null) { $("contact-ab").innerHTML = `<option value="${esc(rest[2])}"></option>`; $("contact-ab").value = rest[2]; }
    await renderContact();
  } else if (view === "epitope") {
    show("epitope");
    await renderEpitope(findProtein(rest[0] || "ERBB2"), rest[1] === "-" ? "" : rest[1]);
  } else {
    show("protein");
    await renderProtein(findProtein(rest[0] || "PTPRC"));
  }
}

async function start() {
  [S.meta, S.proteins, S.pairs, S.surfaces] = await Promise.all(["meta.json", "proteins.json", "pairs.json", "surfaces.json"].map(load));
  S.byAcc = {}; S.byGene = {}; S.byEntry = {};
  for (const p of S.proteins) { S.byAcc[p.acc] = p; S.byEntry[p.entry] = p; S.byGene[p.gene.toUpperCase()] ??= p; }
  $("protein-list").innerHTML = S.proteins.map((p) => `<option value="${esc(p.gene)}">${esc(p.acc)}</option>`).join("");
  $("protein-count").textContent = `${S.proteins.length.toLocaleString()} proteins`;

  const groups = {};
  for (const s of S.meta.surfaces) (groups[s.dataset] ??= []).push(s);
  const options = Object.entries(groups).map(([d, list]) => `<optgroup label="${esc(d)}">` +
    list.map((s) => `<option value="${s.id}">${esc(s.label)}</option>`).join("") + "</optgroup>").join("");
  $("surface-a").innerHTML = options;
  $("surface-b").innerHTML = `<option value="">nothing</option>` + options;
  $("surface-a").value = "adult--monocyte";
  const label = Object.fromEntries(S.meta.surfaces.map((s) => [s.id, s.label]));
  const sideA = [...new Set(S.meta.contacts.map((c) => c.a))], sideB = [...new Set(S.meta.contacts.map((c) => c.b))];
  $("contact-a").innerHTML = sideA.map((id) => `<option value="${id}">${esc(label[id])}</option>`).join("");
  $("contact-b").innerHTML = sideB.map((id) => `<option value="${id}">${esc(label[id])}</option>`).join("");
  $("contact-a").value = "adult--natural-killer-cell";
  $("contact-b").value = "her2-positive-breast-carcinoma";

  for (const b of document.querySelectorAll(".tabs button")) b.addEventListener("click", () => {
    location.hash = b.dataset.view === "protein" ? `protein/${$("protein-search").value || "PTPRC"}`
      : b.dataset.view === "epitope" ? `epitope/${$("epitope-antigen").value || "ERBB2"}`
      : b.dataset.view === "surface" ? `surface/${$("surface-a").value}/${$("surface-b").value}`
      : `contact/${$("contact-a").value}/${$("contact-b").value}/${$("contact-ab").value}`;
  });
  $("protein-search").addEventListener("change", (e) => { const p = findProtein(e.target.value); if (p) location.hash = `protein/${p.gene}`; else renderProtein(null); });
  $("epitope-antibody").innerHTML = `<option value="">none chosen</option>` + S.meta.antibodies.map((a) =>
    `<option value="${a.acc}">${esc(a.antibody.split(" (")[0])} (${esc(a.gene)})${a.residues ? "" : ", no epitope known"}</option>`).join("");
  $("epitope-antigen").addEventListener("change", (e) => { const q = findProtein(e.target.value); if (q) location.hash = `epitope/${q.gene}`; else renderEpitope(null); });
  $("epitope-antibody").addEventListener("change", (e) => { if (e.target.value) location.hash = `epitope/${S.byAcc[e.target.value].gene}`; });
  $("epitope-residues").addEventListener("change", (e) => {
    const q = findProtein($("epitope-antigen").value), v = e.target.value.trim();
    if (q) location.hash = `epitope/${q.gene}/` + (v ? encodeURIComponent(v) : "-");
  });
  for (const id of ["surface-a", "surface-b"]) $(id).addEventListener("change", () => {
    if ($("surface-a").value === "upload" || $("surface-b").value === "upload") renderSurface();
    else location.hash = `surface/${$("surface-a").value}/${$("surface-b").value}`;
  });
  $("surface-file").addEventListener("change", async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const items = parseUpload(await file.text());
    if (!items.length) { $("surface-body").innerHTML = `<p class="note">No row of that file names a protein of the height table with an abundance.</p>`; return; }
    S.upload = { name: file.name, items };
    for (const id of ["surface-a", "surface-b"]) { $(id).querySelector('option[value="upload"]')?.remove(); $(id).insertAdjacentHTML("beforeend", `<option value="upload">${esc(file.name)}</option>`); }
    $("surface-a").value = "upload";
    renderSurface();
  });
  for (const id of ["contact-a", "contact-b", "contact-ab"]) $(id).addEventListener("change", () => {
    location.hash = `contact/${$("contact-a").value}/${$("contact-b").value}/${$("contact-ab").value}`;
  });
  $("contact-epitope").addEventListener("change", renderContact);

  const d = S.meta.databases;
  $("vintage").innerHTML = `UniProt ${esc(d.uniprot)}, AlphaFold ${esc(d.alphafold)}, Pfam ${esc(d.pfam)}, STRING ${esc((d.string || "").split(" ")[0])}, ` +
    `CellphoneDB ${esc((d.cellphonedb || "").split(" ")[0])}. Built from pipeline commit ${esc(S.meta.built_from)}. ` +
    `Code and methods: <a href="https://github.com/kmsouthard/surfaceome_topography">github.com/kmsouthard/surfaceome_topography</a>.`;
  $("loading").hidden = true;
  window.addEventListener("hashchange", route);
  await route();
}

start().catch((e) => { $("loading").textContent = "The tables could not be loaded: " + e.message; });
