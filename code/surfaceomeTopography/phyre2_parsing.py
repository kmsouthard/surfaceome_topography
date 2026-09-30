#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 14:45:17 2019

@author: southk

Processing Phyre2 output Data:

combining protein measurements and Phyre2 results annotation
"""

import numpy as np


         
            
    
    
    
    
def find_dims(df):
    #find max and min dimensions, need dataframe with measured dims
    df['max_dim'] = df[["x", "y", "z"]].max(axis=1)
    df['min_dim'] = df[["x", "y", "z"]].min(axis=1)
    
    return df


def largest_unique_ecd(df):
    #select for the largest ecd when there are multiple measured ecds
    
    idx = df.groupby(['ID link'])['max_dim'].transform(max) == df['max_dim']
    
    return df[idx]


def well_modeled(df, confidence = 90, coverage = 70):
    #identify the models above a confidence of sequence converage threshold
    df = df[(df['Confidence (%)'] >= confidence) & (df['Alignment coverage (%)'] >= coverage)]
    
    return df
    


#Testing code


#set path for protien measurements
            
#file pattern
#pattern = 'summaryinfo'
#pattern = 'summaryinfo_extended'
#pattern = 'some random string'

#test!
#info = read_summaryinfo(path, pattern)
