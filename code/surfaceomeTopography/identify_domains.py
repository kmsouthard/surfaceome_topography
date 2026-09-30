#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 14:25:46 2019

@author: southk

extracellular domain identification
    steps: load proteome gff and search for ECD annotation in 'notes' column
    how to deal with gpi anchored protiens: gpi/lypidation identified in gff
        some cases maybe not properly annotated eg CEACAM5 in mouse

for mouse proteome more annotations for gpi (150 vs 126) if tab file and
subcellular location [CC] column (GO term annotation) is used for ID
may need retrive annotations from both tab and gff files for full coverage
as proteomes may not all be annotated in the same way


unannotated proteins (uniprot) maybe need another identification method
cell membrane go terms and transmembrane domain annotation? Calculate ECD
from TMD annotation?


starting with just proteins that are annotated with an extracellular domain
    and gpi
    then compare to published surfaceome datasets, what is missing/how are
    they annotated
    Results: generally good human surfaceome coverage, most agreement with the
    machine learning surfaceome (2542 common/2886), both missed proteins and
    unique proteins their localization is uncertain

"""
import pandas as pd
from .database_parsing import fasta_subset
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC

def seq_len(df_gff):
    return df_gff.assign(seq_len = df_gff['end']-df_gff['start'])


def subset_ecds(df_gff, df_tab):

    #returns subset of sequences annotated as Extracellular in gff file
    df_ecds = df_gff[df_gff['note'].str.contains('Extracellular', na=False)].reset_index(drop=True)

    # df_loc = df_tab.loc[
    #             df_tab['Subcellular location [CC]'].str.contains('Cell membrane|Membrane|cell membrane|membrane',na=False)
    #               &
    #             df_tab['Subcellular location [CC]'].str.contains('Single-pass membrane protein|Multi-pass membrane protien',na=False)
    #                 ]
    # mem_prots = df_gff.merge(df_loc['Entry'], right_on = 'Entry', left_on = 'ID link', how = 'right')
    #
    # df_mem = mem_prots.loc[(mem_prots['td'] == 'Chain') | (mem_prots['td'] == 'Peptide')]
    #
    # df = pd.concat([df_mem, df_ecds])

    return seq_len(df_ecds)

def subset_icds(df_gff):

    #returns subset of sequences annotated as cytoplasmic in gff file
    df = df_gff[df_gff['note'].str.contains('Cytoplasmic', na=False)].reset_index(drop=True)
    return seq_len(df)


def subset_gpi(df_tab, df_gff):

    #returns subset of sequences with cell membrane and gpi annotations (tab)
    #retrieves sequences for modeling from gff file
    gpi = df_tab.loc[
                df_tab['Subcellular location [CC]'].str.contains('Cell membrane|Membrane|cell membrane|membrane',na=False)
              &
                df_tab['Subcellular location [CC]'].str.contains('GPI',na=False)
                ]
    gpi_prots = df_gff.merge(gpi['Entry'], right_on = 'Entry', left_on = 'ID link', how = 'right')

    df = gpi_prots.loc[(gpi_prots['td'] == 'Chain') | (gpi_prots['td'] == 'Peptide')]
    return seq_len(df)

def subset_tmd(df_gff):
    #return identified transmembrane domains
    df = df_gff[df_gff['td'].str.contains('Transmembrane', na =False)].reset_index(drop=True)

    return seq_len(df)

def order_domains(ecds, icds, tmds, gpi):
    #takes subseted_domain dataframes
    #concats, orders and lables
    ecd = ecds.assign(domain = 'ecd', source = 'ecd')
    icd = icds.assign(domain = 'icd', source = 'icd')
    tm = tmds.assign(domain = 'tmd', source = 'tmd')
    gpis = gpi.assign(domain = 'ecd', source = 'gpi').drop(columns = 'Entry')

    df = pd.concat([ecd,icd, tm, gpis]).reset_index(drop =True)

    #assign domain order
    df = df.assign(domain_order = df.sort_values('start').groupby('ID link')['start'].rank('dense'))

    return df

def split_seq(seq_len):
    #split sequence in the middle
    first_end = round(seq_len / 2)
    second_start = first_end + 1
    return first_end, second_start

def shorten_end(seq_len, start, end, c):
    first_end, second_start = split_seq(seq_len)
    if c == 0:
        return first_end + start
    if c == 1:
        return end

def truncate_start(seq_len, start, end, c):
    first_end, second_start = split_seq(seq_len)
    if c == 0:
        return start
    if c == 1:
        return (second_start+start)

def split_long(df):
    long = df[df['seq_len'] >= 3000]
    dup = pd.concat([long]*2, ignore_index=True)
    c = dup.groupby(['ID link']).cumcount()
    grouped = dup.assign(g = c)
    grouped['end'] = grouped.apply(lambda row: shorten_end(row['seq_len'], row['start'], row['end'], row['g']), axis=1)
    grouped['start'] = grouped.apply(lambda row: truncate_start(row['seq_len'], row['start'], row['end'], row['g']), axis=1)
    return grouped.dropna(subset=['start','end'])

def save_long_splits(ecds, fasta, path):
    if (ecds['seq_len']>= 3000).any() == True:
        long = split_long(ecds)
        long_seq =fasta_subset(long, fasta)
    #save long data
        SeqIO.write(long_seq, path+'_long_splits.fasta', "fasta")

    else:
        print('no long sequences')

def phyre_filter(df):
    phyre = df[(df['seq_len'] >= 30) & (df['seq_len'] < 3000)]
    return phyre.reset_index(drop=True)


def surfaceome_df(gpi, ecd, short = None, long = pd.DataFrame()):
    #combine dataframes for creating total ECDs dataframe
    #use
    df1 = gpi.assign(source = 'gpi')
    df2 = ecd.assign(source = 'ecd')

    if short is not None:
        df4 =short.assign(source = 'short')

    else:
        df2 = df2[~(df2['seq_len'] < 30)]
        short = ecd[(ecd['seq_len'] < 30)]
        df4 =short.assign(source = 'short')

    if long.empty == True:
        total = pd.concat([df1, df2, df4], sort = False)
        total.drop(['Entry'], axis =1, inplace = True)

    else:
        df3 =long.assign(source = 'long')
        total = pd.concat([df1, df2, df3, df4], sort = False)
        total.drop(['Entry', 'g'], axis =1, inplace = True)

    return total
