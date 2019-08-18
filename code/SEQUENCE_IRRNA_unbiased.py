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

lib300=library[library.library=='ir']
lib300 = lib300.dropna(axis=1, how='all')

rnareads=pd.Series([''],index=lib300.index)

for filename in os.listdir('../rawdata/ir/'):
    if ('coveragePYTHON-' in filename):
        splitcov = pd.read_pickle('../rawdata/ir/' + filename)
        rnareads = rnareads.add(splitcov)

rnareads.to_pickle('../rawdata/ir/rnareads_ir_unbiased.pkl')


#%%

rnair=pd.DataFrame()
#exon2=
exon1='CTTGGAAGGCCGTCTCGTGG'
cryptic=pd.Series()
downstream1='CTCTCTAAAAAAAATCCTTC'

rna=pd.DataFrame()

for var in lib300[lib300.subset!='irmelanomavariants'].index:
    varsq=lib300.varseq162[var]
    ind=1
    intronstarts=[]
    intronends=[]
    reads=[]
    rr=pd.Series(rnareads.loc[var].split(' '))
    for test in rr:
    #    print(test)
        vspos=0
        test=str(Seq(test).reverse_complement())
        startpos=80
        testsub=test[startpos:90]
        while True:
    #        print(startpos)
            vspos=varsq.find(testsub)
            startpos-=1
            testsub=test[startpos:90]
            if (varsq.find(testsub)==-1)|(startpos<0):
                if (startpos<0):
                    vspos=0
                break
        intronends.append(vspos)
        if (vspos==0)|(varsq.find(test[startpos-9:startpos+1])==-1):
            intronstarts.append(0)
        elif (vspos==-1):
            intronstarts.append(-150)
        else:
            intronstarts.append(varsq.find(test[startpos-9:startpos+1])+10)
        reads.append(test)
        ind+=1    
    
    introns=pd.DataFrame([intronstarts,intronends, reads]).transpose()
    introns.columns=['intronstart','intronend', 'read']
    
    introns['intlen']=introns.index.map(lambda x: introns.intronend[x]-introns.intronstart[x])
    
    rna.loc[var,'rnareads']=len(introns)
    rna.loc[var,'unspliced']=len(introns[introns.intlen==0])
    rna.loc[var,'spliced']=len(introns[introns.intlen==float(lib300.intronlength[var])])
    rna.loc[var,'fraction_canonical']=(rna.loc[var,'spliced']+rna.loc[var,'unspliced'])/float(rna.loc[var,'rnareads'])
    rna.loc[var,'fraction_unknowndonor']=len(introns[introns.intlen==-1])/float(rna.loc[var,'rnareads'])
    introns.drop(introns[(introns.intlen<11)|(introns.intlen>140)].index, inplace=True)
    vals=introns.intlen.value_counts()
    if (float(lib300.intronlength[var]) in vals.index):
        vals.drop(float(lib300.intronlength[var]), inplace=True)
    rna.loc[var,'numberofcryptic']=len(vals[vals>rna.loc[var,'rnareads']*0.1])
    crypticreads=[]
    crypticintronlengths=[]
    for i in range(int(rna.loc[var,'numberofcryptic'])):
        crypticreads.append(introns[introns.intlen==vals[vals>rna.loc[var,'rnareads']*0.1].index[i]].read.value_counts().index[0])
        crypticintronlengths.append(str(vals.index[i]))
    rna.loc[var,'cryptic_read']=' '.join(crypticreads)
    rna.loc[var,'cryptic_intronlength']=' '.join(crypticintronlengths)
           


rna['fraction_spliced_total']=rna.index.map(lambda x: rna.spliced[x]/float(rna.rnareads[x]))
rna['fraction_unspliced_total']=rna.index.map(lambda x: rna.unspliced[x]/float(rna.rnareads[x]))
rna['fraction_spliced_canonical']=rna.index.map(lambda x: rna.spliced[x]/float(rna.unspliced[x] + rna.spliced[x]))
rna['levelratiospl_unbiased']=rna.index.map(lambda x: np.log2(rna.spliced[x]/float(rna.unspliced[x])))
rna['levelratiospl_againstallothers']=rna.index.map(lambda x: np.log2(rna.spliced[x]/float(rna.rnareads[x] - rna.spliced[x])))

rna.to_pickle('../rawdata/ir/rna_from_unbiased_mapping.pkl')





