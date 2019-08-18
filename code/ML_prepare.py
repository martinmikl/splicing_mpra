#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 14:36:33 2017

@author: martinm
"""


import pandas as pd
from sklearn.cross_validation import train_test_split
import random
import os
#%%
#os.mkdir('./dataframes/ml')
#os.mkdir('./dataframes/ml/ir')
#os.mkdir('./dataframes/ml/cas')
#os.mkdir('./dataframes/ml/five')
#os.mkdir('./dataframes/ml/three')

#%% Load datasets from Supplementary Tables

irdf=pd.read_excel('./tables/TableS5.xlsx')

irdf.columns=['libindex','subset','name2','intronlength','intronstart_varseq','intronend_varseqnew','maxent5','maxent3',
      'exon1','donor','intron5','introncenter','intron3','acceptor','exon2',
'levelratiospl','wtratio','normratio','wav_stats','normwav','rnaperdna','normrnaperdna',
'changes','combination','first_ss','second_ss','firstss','secondss',
'rnareads_effective','number_reads','fraction_canonical','smoothednumberofpeaks','rawnumberofpeaks','noisestrengthlogwstd','noiseres_linear',
'noiseresgam','normnoise_linear','normnoise_gam','varseq']

irdf.set_index('libindex', inplace=True)
irdf['varseq162']=irdf.varseq.apply(lambda x: str(x[30:-18]).upper())

###

casdf=pd.read_excel('./tables/TableS6.xlsx')

casdf.columns=['libindex','subset','commonname','exonlength','exonstart_varseq','exonend_varseq','maxent5','maxent3',
       'intron1','acceptor','exon5','exoncenter','exon3','donor','intron2',
'levelratiospl','psi','wtratio','normratio','normpsi',
'changes','combination','first_ss','second_ss','firstss','secondss','first_SF','second_SF',
'rnareads_effective','fraction_canonical','varseq']

casdf.set_index('libindex', inplace=True)
casdf['varseq162']=casdf.varseq.apply(lambda x: str(x[30+15:-18]).upper())

###

fivedf=pd.read_excel('./tables/TableS7.xlsx')

fivedf.columns=['libindex','subset','commonname','donorpos1','donorpos2','diff_nt','maxent5first','maxent5second',
        'exon','donor1','alt5','altcenter','alt3','donor2','intron',
'levelratiospl','wtratio','normratio','wav123','normwav','rnaperdna','normrnaperdna',
'changes','first_ss','second_ss','firstss','secondss',
'rnareads_effective','number_reads','fraction_canonical','smoothednumberofpeaks','rawnumberofpeaks',
'noisestrengthlog','noiseresgam','normnoise_gam','varseq']

fivedf.set_index('libindex', inplace=True)
fivedf['varseq162']=fivedf.varseq.apply(lambda x: str(x[30+15:-18]).upper())

###

threedf=pd.read_excel('./tables/TableS8.xlsx')

threedf.columns=['libindex','subset','commonname','acceptorpos1','acceptorpos2','diff_nt','maxent3first','maxent3second',
                     'intron','acceptor1','alt5','altcenter','alt3','acceptor2','exon',
'levelratiospl','wtratio','normratio',
'changes','first_ss','second_ss','firstss','secondss',
'rnareads_effective','fraction_canonical','varseq']

threedf.set_index('libindex', inplace=True)
threedf['varseq162']=threedf.varseq.apply(lambda x: str(x[30+15:-18]).upper())

#%% IR [running this again will lead to different variants with the same sequence to be dropped and therefore models will need to be retrained with the new training sets]
'''
irl300=irdf[(irdf.subset!='irmelanomavariants')]

# drop indexes for 
indexes=[]
samevarseq=[]
for var, group in irl300.groupby(by='varseq162'):
	if (len(group)>1):
		samevarseq.append(var)
		indexes.append(list(irl300[irl300.varseq162==var].index))
		
indexestodrop=[]
for igroup in indexes:
	indexestodrop=indexestodrop + random.sample(igroup, len(igroup)-1)
		
irl300forml=irl300.drop(indexestodrop) 

#####

irlml,irleval, irlmly,irlevaly=train_test_split(irl300forml, irl300forml['levelratiospl'], \
												test_size=0.1, random_state=0)

irlml.to_pickle('./dataframes/ml/ir/irml.pkl')
irleval.to_pickle('./dataframes/ml/ir/ireval.pkl')
irlmly.to_pickle('./dataframes/ml/ir/irmly.pkl')
irlevaly.to_pickle('./dataframes/ml/ir/irevaly.pkl')
'''
#%% cassette [running this again will lead to different variants with the same sequence to be dropped and therefore models will need to be retrained with the new training sets]
'''
irl300=casdf[(casdf.subset!='cassettemelanomavariants')]

# drop indexes for 
indexes=[]
samevarseq=[]
for var, group in irl300.groupby(by='varseq162'):
	if (len(group)>1):
		samevarseq.append(var)
		indexes.append(list(irl300[irl300.varseq162==var].index))
		
indexestodrop=[]
for igroup in indexes:
	indexestodrop=indexestodrop + random.sample(igroup, len(igroup)-1)
		
irl300forml=irl300.drop(indexestodrop) 

#####

irlml,irleval, irlmly,irlevaly=train_test_split(irl300forml, irl300forml['levelratiospl'], \
												test_size=0.1, random_state=0)

irlml.to_pickle('./dataframes/ml/cas/casml.pkl')
irleval.to_pickle('./dataframes/ml/cas/caseval.pkl')
irlmly.to_pickle('./dataframes/ml/cas/casmly.pkl')
irlevaly.to_pickle('./dataframes/ml/cas/casevaly.pkl')
''' 
#%% five [running this again will lead to different variants with the same sequence to be dropped and therefore models will need to be retrained with the new training sets]
'''
irl300=fivedf[(fivedf.subset!='fivemelanomavariants')]

# drop indexes for 
indexes=[]
samevarseq=[]
for var, group in irl300.groupby(by='varseq162'):
    if (len(group)>1):
        samevarseq.append(var)
        indexes.append(list(irl300[irl300.varseq162==var].index))
        
indexestodrop=[]
for igroup in indexes:
    indexestodrop=indexestodrop + random.sample(igroup, len(igroup)-1)
        
irl300forml=irl300.drop(indexestodrop) 

#####

irlml,irleval, irlmly,irlevaly=train_test_split(irl300forml, irl300forml['levelratiospl'], \
                                                test_size=0.1, random_state=0)

irlml.to_pickle('./dataframes/ml/five/fiveml.pkl')
irleval.to_pickle('./dataframes/ml/five/fiveeval.pkl')
irlmly.to_pickle('./dataframes/ml/five/fivemly.pkl')
irlevaly.to_pickle('./dataframes/ml/five/fiveevaly.pkl')
'''
#%% three [running this again will lead to different variants with the same sequence to be dropped and therefore models will need to be retrained with the new training sets]
'''
irl300=threedf[(threedf.subset!='threemelanomavariants')]

# drop indexes for 
indexes=[]
samevarseq=[]
for var, group in irl300.groupby(by='varseq162'):
    if (len(group)>1):
        samevarseq.append(var)
        indexes.append(list(irl300[irl300.varseq162==var].index))
        
indexestodrop=[]
for igroup in indexes:
    indexestodrop=indexestodrop + random.sample(igroup, len(igroup)-1)
        
irl300forml=irl300.drop(indexestodrop) 

#####

irlml,irleval, irlmly,irlevaly=train_test_split(irl300forml, irl300forml['levelratiospl'], \
                                                test_size=0.1, random_state=0)

irlml.to_pickle('./dataframes/ml/three/threeml.pkl')
irleval.to_pickle('./dataframes/ml/three/threeeval.pkl')
irlmly.to_pickle('./dataframes/ml/three/threemly.pkl')
irlevaly.to_pickle('./dataframes/ml/three/threeevaly.pkl')
'''