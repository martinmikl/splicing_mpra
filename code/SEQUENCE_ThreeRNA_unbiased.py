#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 17:40:14 2017

@author: martinm
"""


import pandas as pd
import numpy as np
import os

#%%

library=pd.read_pickle('./design/LIBRARY/alllibraries210.pkl')

lib300=library[library.library=='three']
lib300 = lib300.dropna(axis=1, how='all')

rnareads=pd.Series([''],index=lib300.index)

for filename in os.listdir('../rawdata/three'):
    if ('coveragePYTHON-' in filename):
        splitcov = pd.read_pickle('../rawdata/three' + filename)
        rnareads = rnareads.add(splitcov)

rnareads.to_pickle('../rawdata/three/rnareads_three_unbiased.pkl')


#%%

upstream='GCCCCATACCTGAAGACCAAGTTTATCTGTGTGACACCGTAAGTGGCTTCCTTTCCCCGTTTTGCCTTCATTTCTAATATCCTCAGTTATCCCTGGGAATGGGACACTGGGTGAGAGTTAATCTGCCAAAGGTTGGAAGCCCCTGGGCTATGTTTAGTACTCAAAGTGACGACGTC'
rna=pd.DataFrame()

for var in lib300[lib300.subset!='threemelanomavariants'].index:
    varsq=(upstream + lib300.varseq162[var]).upper()
    ind=1
    intronstarts=[]
    intronends=[]
    reads=[]
    rr=pd.Series(rnareads.loc[var].split(' '))
    for test in rr:
    #    print(test)
        if (len(test)>100):
            vspos=0
            test=str(Seq(test).reverse_complement())
            startpos=90
            testsub=test[startpos:100]
            while True:
        #        print(startpos)
                vspos=varsq.find(testsub)
                startpos-=1
                testsub=test[startpos:100]
                if (varsq.find(testsub)==-1)|(startpos<0):
                    if (startpos<0):
                        vspos=0
                    break
            intronends.append(vspos)
            if (vspos==0)|(varsq.find(test[startpos-19:startpos+1])==-1):
                intronstarts.append(0)
            elif (vspos==-1):
                intronstarts.append(-150)
            else:
                intronstarts.append(varsq.find(test[startpos-19:startpos+1])+20)
            reads.append(test)
            ind+=1    
    
    introns=pd.DataFrame([intronstarts,intronends, reads]).transpose()
    introns.columns=['intronstart','intronend', 'read']
    
    introns['intlen']=introns.index.map(lambda x: introns.intronend[x]-introns.intronstart[x])

    if (len(introns)>0):    
        rna.loc[var,'rnareads']=len(introns)
        rna.loc[var,'unspliced']=len(introns[introns.intronend==0])
        rna.loc[var,'spliced1']=len(introns[introns.intlen==214])
        rna.loc[var,'spliced2']=len(introns[introns.intlen==214+int(lib300.diff_nt[var])])
        rna.loc[var,'fraction_canonical']=(rna.loc[var,'spliced1']+rna.loc[var,'spliced2'])/float(rna.loc[var,'rnareads'])
        rna.loc[var,'fraction_unknowndonor']=len(introns[introns.intlen==303])/float(rna.loc[var,'rnareads'])
        rna.loc[var,'deletion']=len(introns[introns.intlen==1])    
        introns.drop(introns[(introns.intlen<11)|(introns.intlen==214)].index, inplace=True)
        vals=introns.intlen.value_counts()
        if (214+int(lib300.diff_nt[var]) in vals.index):
            vals.drop(214+int(lib300.diff_nt[var]), inplace=True)
        rna.loc[var,'numberofcryptic']=len(vals[vals>rna.loc[var,'rnareads']*0.05])
        crypticreads=[]
        crypticintronlengths=[]
        for i in range(int(rna.loc[var,'numberofcryptic'])):
            crypticreads.append(introns[introns.intlen==vals[vals>rna.loc[var,'rnareads']*0.05].index[i]].read.value_counts().index[0])
            crypticintronlengths.append(str(vals.index[i]-214))
        rna.loc[var,'cryptic_read']=' '.join(crypticreads)
        rna.loc[var,'cryptic_diffnt']=' '.join(crypticintronlengths)
           

rna['fraction_spliced2_total']=rna.index.map(lambda x: rna.spliced2[x]/float(rna.rnareads[x]))
rna['fraction_spliced1_total']=rna.index.map(lambda x: rna.spliced1[x]/float(rna.rnareads[x]))
rna['fraction_spliced_canonical']=rna.index.map(lambda x: rna.spliced2[x]/float(rna.spliced1[x] + rna.spliced2[x]))
rna['levelratiospl_unbiased']=rna.index.map(lambda x: np.log2(rna.spliced2[x]/float(rna.spliced1[x])))
rna['levelratiospl_againstallothers']=rna.index.map(lambda x: np.log2(rna.spliced2[x]/float(rna.rnareads[x] - rna.spliced2[x])))
rna['fraction_deletion']=rna.index.map(lambda x: rna.deletion[x]/float(rna.rnareads[x]))
rna['fraction_unspliced']=rna.index.map(lambda x: rna.unspliced[x]/float(rna.rnareads[x]))
rna['fraction_canonicalunspliced']=rna.index.map(lambda x: (rna.unspliced[x] + rna.spliced1[x] + rna.spliced2[x])/float(rna.rnareads[x]))

rna.to_pickle('../rawdata/three/rna_from_unbiased_mapping.pkl')
