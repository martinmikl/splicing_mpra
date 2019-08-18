#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 17:14:57 2017

@author: martinm
"""

import pandas as pd
import numpy as np
import os
from Bio.Seq import Seq
#%%


library=pd.read_pickle('./design/LIBRARY/alllibraries210.pkl')

lib300=library[library.library=='cassette']
lib300 = lib300.dropna(axis=1, how='all')

rnareads=pd.Series('',index=lib300.index)

for filename in os.listdir('../rawdata/cas/'):
    if ('coveragePYTHON-' in filename):
        splitcov = pd.read_pickle('../rawdata/cas/' + filename)
        rnareads = rnareads.add(splitcov)

rnareads.to_pickle('../rawdata/cas/rnareads_cas_unbiased.pkl')


#%%


rnak562=pd.DataFrame()
exon2='AGAACTCCACAAACCCATC'
exon1='CTTGGAAGGCCGTCTCGTGG'
cryptic=pd.Series()
downstream1='CTCTCTAAAAAAAATCCTTC'

for var in rnareads.index:
    rr=pd.Series(rnareads.loc[var, 'K562'].split(' '))    
    rr.drop(0, inplace=True)
    for i in rr.index:
        if (exon2 not in rr.loc[i]):
            rr.drop(i, inplace=True)
    
    rrup=rr.apply(lambda x: x[114:134])
    casexon=str(Seq(lib300.varseq[var][142:162]).reverse_complement()).upper()
    if (casexon in rrup.values):
        rnak562.loc[var,'incl']=rrup.value_counts()[casexon]
    else:
        rnak562.loc[var,'incl']=0 
    if (exon1 in rrup.values):
        rnak562.loc[var,'excl']=rrup.value_counts()[exon1]
    else:
        rnak562.loc[var,'excl']=0
    if (downstream1 in rrup.values):
        rnak562.loc[var,'cryptic']=rrup.value_counts()[downstream1]
    else:
        rnak562.loc[var,'cryptic']=0
    rnak562.loc[var,'rnareads']=len(rr)
    rnak562.loc[var,'rawreads']=' '.join(rr)  

             
rnak562['fraction_incl']=rnak562.index.map(lambda x: rnak562.incl[x]/rnak562.rnareads[x])
rnak562['fraction_excl']=rnak562.index.map(lambda x: rnak562.excl[x]/rnak562.rnareads[x])
rnak562['fraction_canonical']=rnak562.index.map(lambda x: (rnak562.incl[x]+rnak562.excl[x])/rnak562.rnareads[x])
rnak562['fraction_cryptic']=rnak562.index.map(lambda x: rnak562.cryptic[x]/rnak562.rnareads[x])
rnak562['psi']=rnak562.index.map(lambda x: rnak562.incl[x]/(rnak562.incl[x]+rnak562.excl[x]))
rnak562['levelratiospl']=rnak562.index.map(lambda x: np.log2(rnak562.incl[x]/rnak562.excl[x]))
rnak562.to_pickle('../rawdata/cas/rnacas.pkl')


