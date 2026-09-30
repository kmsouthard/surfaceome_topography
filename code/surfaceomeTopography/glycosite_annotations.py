"""Locate glycosylation sites within extracellular domains.

Sites come from three sources, unioned: UniProt `CARBOHYD` features parsed from the proteome
GFF, GlycoMine predictions (`parse_glycoMine`), and -- since 2026 -- GlyGen experimental sites
(see `download_glygen`). `surface_glycosites` restricts a site list to the extracellular spans
of the surfaceome, which is what makes a per-segment glycan density meaningful.

That density sets the persistence length of every disordered segment, so this module is further
upstream of the height estimates than its size suggests. See `disorder_height`.
"""
def parse_glycosites(gff, surfaceome_accessions):
    glycosylation = gff[(gff['td'] == 'Glycosylation') & (gff['ID link'].isin(surfaceome_accessions))]

    glycosylation = glycosylation.assign(glycosite = glycosylation['start'])

    glycosylation.drop(columns = ['start', 'end'], inplace = True)

    return glycosylation

def surface_glycosites(glycosylation, surfaceome):

    surfaceome_gly = surfaceome[['ID link', 'start', 'end']].merge(glycosylation, on ='ID link', how ='left')

    ecd_glycosites = surfaceome_gly[
                    (surfaceome_gly['glycosite'] <= surfaceome_gly['end'])
                    &(surfaceome_gly['glycosite'] >= surfaceome_gly['start'])]

    return ecd_glycosites

def parse_glycoMine(glycoMine, surfaceome_accessions):

    glycoMine['ID link'] = glycoMine['UniProtID']

    glycosylation = glycoMine[glycoMine['ID link'].isin(surfaceome_accessions)]

    glycosylation = glycosylation.assign(glycosite = glycosylation['Site'])
    glycosylation.drop(columns = ['UniProtID', 'Site'], inplace = True)

    return glycosylation
