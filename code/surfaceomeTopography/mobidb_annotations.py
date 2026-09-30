#!/usr/bin/env python2
'''
Mobidb-lite disorder annotations
file could be processed into gff like file? and consensus disorder regions

takes proteome ID and finds disordered protein consensus annotations

Proteome ID
*UP000000803
*UP000005640
*UP000000589
*UP000002195
*UP000001940
*UP000000437

* = done


'''


import pandas as pd
import numpy as np

# Define request

#proteome = 'UP000005640'


#surfaceome = pd.read_csv(surfaceome_path+proteome+'.csv')

#accessions = surfaceome['ID link'].unique().tolist()


#expand list of lists into their own rows
def explode(df, lst_cols, fill_value='', preserve_index=False):
    # make sure `lst_cols` is list-alike
    if (lst_cols is not None
        and len(lst_cols) > 0
        and not isinstance(lst_cols, (list, tuple, np.ndarray, pd.Series))):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)
    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()
    # preserve original index values
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    res = (pd.DataFrame({
                col:np.repeat(df[col].values, lens)
                for col in idx_cols},
                index=idx)
             .assign(**{col:np.concatenate(df.loc[lens>0, col].values)
                            for col in lst_cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        res = (pd.concat([res, df.loc[lens==0, idx_cols]], sort=False)
                  .fillna(fill_value))
    # revert the original index order
    res = res.sort_index()
    # reset index if requested
    if not preserve_index:
        res = res.reset_index(drop=True)
    return res

def import_disorder(output_path, proteome, method = 'full'):
    disorder = pd.read_csv(output_path+proteome+'_disorder_'+method+'.txt', sep = '\t',
                      header = None, names = ['ID link', 'mobi', 'disordered_regions'])

    disorder['disordered_regions'] = disorder.disordered_regions.apply(lambda x: x[2:-2].split('],['))

    #Keep proteins that carry at least one region.  The test used to be `len(d) > 2`, and
    #`len(d)` counts *regions*, not characters -- a protein with no annotation is the string
    #"[]", which splits to one empty element, so the guard was plainly there to drop those.
    #But it also dropped every protein with one or two regions.  MobiDB-lite lists only the
    #disordered stretches and most proteins have one or two, so 2,396 of 2,645 were being
    #discarded -- 91% of that track.  On the surfaceome that cost 262 proteins their MobiDB
    #annotation and 17,265 disordered residues, 29% of the MobiDB-derived budget, and it is
    #a large part of why disorder looked so dependent on AlphaFold pLDDT.  The `full` track
    #annotates every residue with a structural state, so it nearly always has three or more
    #regions and is barely affected either way.
    #
    #A region needs two commas to be parsed (start, end, type), which is what distinguishes
    #a real single region from the empty-annotation placeholder.
    disorder = disorder[disorder['disordered_regions'].map(
        lambda d: any(str(r).count(',') >= 2 for r in d))]

    return explode(disorder, 'disordered_regions')

def parse_disorder(disorder):
    #assign start col
    disorder2 = disorder.assign(disorder_start = disorder['disordered_regions'].apply(lambda x: x.split(',')[0]))
    #assign end col
    disorder3 = disorder2.assign(disorder_end = disorder2['disordered_regions'].apply(lambda x: x.split(',')[1]))
    #assign type col
    disorder4 = disorder3.assign(disorder_type = disorder3['disordered_regions'].apply(lambda x: x.split(',')[2]))
    #clean type col
    disorder5 = disorder4.assign(disorder_type = disorder4.disorder_type.apply(lambda x: x.strip('"')))

    disorder6 = disorder5[disorder5.disorder_type.str.contains('D|d')]

    return disorder6

def surface_disorder(disorder, surfaceome):

    surfaceome_disorder = surfaceome.merge(disorder, on ='ID link', how ='left')
    surfaceome_disorder['disorder_start'] = pd.to_numeric(surfaceome_disorder['disorder_start'])
    surfaceome_disorder['disorder_end'] = pd.to_numeric(surfaceome_disorder['disorder_end'])

    ecd_disorder = surfaceome_disorder[
                        (surfaceome_disorder['end'] >= surfaceome_disorder['disorder_end'] - 2 )
                        &(surfaceome_disorder['start'] <= surfaceome_disorder['disorder_start'] + 2)] #allow some overlap with start/end

    ecd_disorder = ecd_disorder.assign(disorder_len = ecd_disorder['disorder_end'] - ecd_disorder['disorder_start'])

    return ecd_disorder[ecd_disorder.disorder_len > 0].drop_duplicates(keep = 'first')

def summarize_disorder(ecd_disorder, surfaceome):

    disorder_grouped = ecd_disorder[['ID link','disorder_start', 'disorder_end', 'start', 'end','disorder_len']].groupby(['ID link', 'start', 'end'], as_index =False)

    disorder_sum = ecd_disorder[['ID link','disorder_start', 'disorder_end', 'start', 'end','disorder_len']].groupby(['ID link', 'start', 'end'], as_index =False)['disorder_len'].sum()
    disorder_sum = disorder_sum.assign(disorder_first = disorder_grouped['disorder_start'].min()['disorder_start'],
                   disorder_last = disorder_grouped['disorder_end'].max()['disorder_end'])

    surfaceome_disorder = pd.merge(surfaceome, disorder_sum, on = ['ID link', 'start', 'end'], how = 'left')

    surfaceome_disorder = surfaceome_disorder.assign(percent_disorder = (surfaceome_disorder.disorder_len / surfaceome_disorder.seq_len) * 100)

    return surfaceome_disorder
