#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 17:39:21 2019

@author: southk
"""
import gzip
from pathlib import Path

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC


def _open_text(path):
    """Open a file for reading, transparently un-gzipping a ``.gz``.

    The proteome FASTA is 87 MB as text and is committed gzipped, so ``resolve()`` can hand
    back either form.  ``pd.read_csv`` decompresses by extension on its own; ``SeqIO.parse``
    needs the handle opened for it.
    """
    return gzip.open(path, "rt") if Path(path).suffix == ".gz" else open(path)


def import_gff(path):

    #read gff file and assign columns
    full = pd.read_csv(path, sep ='\t', comment ='#', names = ['ID link',
                       'up', 'td', 'start','end', '.', '-', '+', 'note', 'NaN'])

    #drop columns with no infromation
    full = full.drop(['up','.','-','+','NaN'], axis =1)

    #returns dataframe with uniprot IDs in first column
    return full


def import_tab(path):

    #read tab file
    full = pd.read_csv(path, sep = "\t")

    return full


def import_fasta(path):

    with _open_text(path) as handle:
        Full_seq = list(SeqIO.parse(handle, "fasta"))


    for seq_record in Full_seq:
        # Take the current sequence ID and split
        ID_list = str(seq_record.id).split('|')
        #specifiy ID as ID link
        ID = ID_list[1]
        #clean up name
        Name = ID_list[2].split('/')
        #add back to seq_record
        seq_record.id = ID
        seq_record.name = Name[0]

    record_dict = SeqIO.to_dict(Full_seq)
    return record_dict

def fasta_subset(df, record_dict):
    #create list of sequencences previously selected from gff file

    df.reset_index(inplace = True, drop = True)

    sub_seq = [] # Setup an empty list

    #loop over
    for idx in range(0, len(df['ID link'])):
        r_id = record_dict[df['ID link'][idx]].id+'[' +str(df['start'][idx].astype('int')) +':' + str(df['end'][idx].astype('int'))+']'
        r_seq = record_dict[df['ID link'][idx]].seq[int(df['start'][idx]): int(df['end'][idx])]
        # Add this record to our list
        r = SeqRecord(r_seq, id =r_id)
        sub_seq.append(r)

    print("Found %i sequences" % len(sub_seq))
    return sub_seq


#save ECD data
#SeqIO.write(ECD_seq, "mouse_proteome_minus_modeled_ECD_sequences.fasta", "fasta")
