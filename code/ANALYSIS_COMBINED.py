#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 12:28:21 2017

@author: martinm
"""

import pandas as pd
import numpy as np
import analysis_functions
from maxentpy import maxent
import RNA
import scipy.stats
import os
import pygam

#%%
library=pd.read_pickle('./design/LIBRARY/alllibraries210.pkl')


#%%
######## IR #############

lib300=pd.read_pickle('./dataframes/ir_lib300.pkl')

rnareads=pd.Series([''],index=lib300.index).astype(str)

for filename in os.listdir('../rawdata/ir/'):
    if ('coveragePYTHON-' in filename):
        splitcov = pd.read_pickle('../rawdata/ir/' + filename)
        rnareads = rnareads.add(splitcov)

rnareads.to_pickle('../rawdata/ir/rnareads.pkl')

rna_condition=analysis_functions.unbiased_mapping_ir(rnareads)
rna_condition_final=analysis_functions.prepare_rnadf_ir(rna_condition)
rna_condition_final.to_pickle('../rawdata/ir/rna_from_unbiased_mapping.pkl')
rna_condition_final.to_csv('../rawdata/ir/rna_from_unbiased_mapping.csv')

irdf=pd.read_pickle('../rawdata/ir/rna_from_unbiased_mapping.pkl')
irdf['maxent5']=irdf.index.map(lambda x: maxent.score5(irdf.varseq162[x][int(irdf.intronstart_varseq[x])-3:int(irdf.intronstart_varseq[x])+6]))
irdf['maxent3']=irdf.index.map(lambda x: maxent.score3(irdf.varseq162[x][int(irdf.intronend_varseqnew[x])-20:int(irdf.intronend_varseqnew[x])+3]))
irdf['maxentadd']=irdf.index.map(lambda x: irdf.maxent5[x]+irdf.maxent3[x])

irdf['exon1']=irdf.index.map(lambda x: RNA.fold(irdf[(irdf.intronstart_varseq>24)].varseq162[x][int(irdf[(irdf.intronstart_varseq>24)].intronstart_varseq[x])-24:int(irdf[(irdf.intronstart_varseq>24)].intronstart_varseq[x])])[1] if (irdf.intronstart_varseq[x]>24) else np.nan)
irdf['donor']=irdf.index.map(lambda x: RNA.fold(irdf[(irdf.intronstart_varseq>24)].varseq162[x][int(irdf[(irdf.intronstart_varseq>24)].intronstart_varseq[x])-12:int(irdf[(irdf.intronstart_varseq>24)].intronstart_varseq[x])+12])[1] if (irdf.intronstart_varseq[x]>24) else np.nan)
irdf['intron5']=irdf.index.map(lambda x: RNA.fold(irdf[(irdf.intronstart_varseq>24)].varseq162[x][int(irdf[(irdf.intronstart_varseq>24)].intronstart_varseq[x]):int(irdf[(irdf.intronstart_varseq>24)].intronstart_varseq[x])+24])[1] if (irdf.intronstart_varseq[x]>24) else np.nan)
irdf['introncenter']=irdf.index.map(lambda x: RNA.fold(irdf[(irdf.intronstart_varseq>24)].varseq162[x][(int(irdf[(irdf.intronstart_varseq>24)].intronend_varseqnew[x])-int(irdf[(irdf.intronstart_varseq>24)].intronstart_varseq[x]))/2-12:(int(irdf[(irdf.intronstart_varseq>24)].intronend_varseqnew[x])-int(irdf[(irdf.intronstart_varseq>24)].intronstart_varseq[x]))/2 + 12])[1] if (irdf.intronstart_varseq[x]>24) else np.nan)
irdf['intron3']=irdf.index.map(lambda x: RNA.fold(irdf[(irdf.intronstart_varseq>24)].varseq162[x][int(irdf[(irdf.intronstart_varseq>24)].intronend_varseqnew[x])-24:int(irdf[(irdf.intronstart_varseq>24)].intronend_varseqnew[x])])[1] if (irdf.intronstart_varseq[x]>24) else np.nan)
irdf['acceptor']=irdf.index.map(lambda x: RNA.fold(irdf[(irdf.intronstart_varseq>24)].varseq162[x][int(irdf[(irdf.intronstart_varseq>24)].intronend_varseqnew[x])-12:int(irdf[(irdf.intronstart_varseq>24)].intronend_varseqnew[x])+12])[1] if (irdf.intronstart_varseq[x]>24) else np.nan)
irdf['exon2']=irdf.index.map(lambda x: RNA.fold(irdf[(irdf.intronstart_varseq>24)].varseq[x][int(irdf[(irdf.intronstart_varseq>24)].intronend_varseqnew[x])+30:int(irdf[(irdf.intronstart_varseq>24)].intronend_varseqnew[x])+30+24])[1] if (irdf.intronstart_varseq[x]>24) else np.nan)

#### calculate noise residuals using a linear fit

[slope, const, a, b, c]=scipy.stats.linregress(irdf[(irdf.wav_stats>0)& \
    (irdf.wav_stats<8)&(irdf.number_reads>100)&(irdf.smoothednumberofpeaks==1)&(irdf.fraction_canonical>0.3)].wav_stats, irdf[(irdf.wav_stats>0)& \
    (irdf.wav_stats<8)&(irdf.number_reads>100)&(irdf.smoothednumberofpeaks==1)&(irdf.fraction_canonical>0.3)].noisestrengthlogwstd)

irdf['noiseres_linear']=irdf.index.map(lambda x: irdf.noisestrengthlogwstd[x]- \
       slope*irdf.wav_stats[x] + const if (irdf.wav_stats[x]>0)& \
    (irdf.wav_stats[x]<8)&(irdf.number_reads[x]>100)&(irdf.smoothednumberofpeaks[x]==1)&(irdf.fraction_canonical[x]>0.3) else np.nan)

#########
#irdfold=pd.read_pickle(martin + 'combined_analysis/dataframes/irdf_corrected_fromunbiasedmappingwithoutumis_July2018.pkl')

#### calculate noise residuals using a generalized additive model

meanplusnoise=irdf[(irdf.smoothednumberofpeaks==1)&(irdf.number_reads>100)&(irdf.fraction_canonical>0.3)][['wav_stats','rnaperdna','noisestrengthlogwstd']].dropna()

randomgaml=pygam.LinearGAM(pygam.s(0)+pygam.l(1)).gridsearch(meanplusnoise[['wav_stats','rnaperdna']].values, meanplusnoise.noisestrengthlogwstd.values, lam=[0.01,0.1,1,5,10])

gaml=pygam.LinearGAM(pygam.s(0,lam=1, n_splines=10)+pygam.l(1,lam=1)).fit(meanplusnoise[['wav_stats','rnaperdna']], meanplusnoise.noisestrengthlogwstd)

pred=gaml.predict(meanplusnoise[['wav_stats','rnaperdna']])

meanplusnoise['noisegampred']=pd.Series(pred, index=meanplusnoise.index)
meanplusnoise['noiseresgam']=meanplusnoise['noisestrengthlogwstd']-meanplusnoise['noisegampred']

irdf['noiseresgam']=meanplusnoise['noiseresgam']

#### Calculate expression levels based on RNA/DNA reads

irdf['rnaperdna'] = irdf.index.map(lambda x: np.log2(irdf.loc[x, 'rnareads']/irdf.loc[x, 'number_reads']) if (irdf.number_reads[x]>0) else np.nan)
irdf.loc[irdf.rnaperdna.dropna().index,'rnaperdna']=scipy.stats.zscore(irdf.rnaperdna.dropna())

##### calculate values for the corresponding wt sequence and normalized values for mean, noise and expression level

irdf['wtwav']=irdf.name2.apply(lambda x: irdf[(irdf.subset=='ir_filtered') & \
        (irdf.name2 == x)].wav_stats.dropna().mean())
irdf['normwav'] = irdf.index.map(lambda x: irdf.loc[x, 'wav_stats'] - irdf.loc[x, 'wtwav'])

irdf['wtnoise_linear']=irdf.name2.apply(lambda x: irdf[(irdf.subset=='ir_filtered') & \
        (irdf.name2 == x)].noiseres_linear.dropna().mean())
irdf['normnoise_linear'] = irdf.index.map(lambda x: irdf.loc[x, 'noiseres_linear'] - irdf.loc[x, 'wtnoise_linear'])

irdf['wtnoise_gam']=irdf.name2.apply(lambda x: irdf[(irdf.subset=='ir_filtered') & \
        (irdf.name2 == x)&(irdf.number_reads>200)&(irdf.fraction_canonical>0.5)].noiseresgam.dropna().mean())
irdf['normnoise_gam'] = irdf.index.map(lambda x: irdf.loc[x, 'noiseresgam'] - irdf.loc[x, 'wtnoise_gam'])

irdf['wtrnaperdna']=irdf.name2.apply(lambda x: irdf[(irdf.subset=='ir_filtered') & \
        (irdf.name2 == x)].rnaperdna.dropna().mean())
irdf['normrnaperdna'] = irdf.index.map(lambda x: irdf.loc[x, 'rnaperdna'] - irdf.loc[x, 'wtrnaperdna'])

irdf['intronlength']=irdf.intronlength.astype(int)


irdf.to_pickle('./dataframes/irdf.pkl')


irdf[['subset','name2','intronlength','intronstart_varseq','intronend_varseqnew','maxent5','maxent3',
      'exon1','donor','intron5','introncenter','intron3','acceptor','exon2',
'levelratiospl','wtratio','normratio','wav_stats','normwav','rnaperdna','normrnaperdna',
'changes','combination','first_ss','second_ss','firstss','secondss',
'rnareads_effective','number_reads','fraction_canonical','smoothednumberofpeaks','rawnumberofpeaks','noisestrengthlogwstd','noiseres_linear',
'noiseresgam','normnoise_linear','normnoise_gam','varseq']].to_csv('./tables/Table_irdf.csv')

irdf[['intronstart_varseq','intronend_varseqnew',
'levelratiospl','wav_stats','varseq']].to_csv('../rawdata_FORREPOSITORY/irdf.csv')


########## cassette ##############

caslib=library[library.library=='cassette']
caslib = caslib.dropna(axis=1, how='all')

rnareads=pd.Series([''],index=caslib.index).astype(str)

for filename in os.listdir('../rawdata/cas/'):
    if ('coveragePYTHON-' in filename):
        splitcov = pd.read_pickle('../rawdata/cas/' + filename)
        rnareads = rnareads.add(splitcov)

rnareads.to_pickle('../rawdata/cas/rnareads.pkl')

rna_condition=analysis_functions.unbiased_mapping_cas(rnareads)
rna_condition.to_pickle('../rawdata/cas/mappedreads.pkl')
rna_condition_final=analysis_functions.prepare_rnadf_cas(rna_condition)
rna_condition_final.to_pickle('../rawdata/cas/rna_from_unbiased_mapping.pkl')
rna_condition_final.to_csv('../rawdata/cas/rna_from_unbiased_mapping.csv')

casdf=pd.read_pickle('../rawdata/cas/rna_from_unbiased_mapping.pkl')
casdf['maxent5']=casdf.index.map(lambda x: maxent.score5(casdf.varseq162[x][-33:-24]))
casdf['maxent3']=casdf.index.map(lambda x: maxent.score3(casdf.varseq162[x][-50-int(casdf.exonlength[x]):-27-int(casdf.exonlength[x])]))
casdf['maxentadd']=casdf.index.map(lambda x: casdf.maxent5[x]+casdf.maxent3[x])

casdf['intron1']=casdf.index.map(lambda x: RNA.fold(casdf.varseq162[x][-27-int(casdf.exonlength[x])-24:-27-int(casdf.exonlength[x])])[1])
casdf['acceptor']=casdf.index.map(lambda x: RNA.fold(casdf.varseq162[x][-27-int(casdf.exonlength[x])-12:-27-int(casdf.exonlength[x])+12])[1])
casdf['exon5']=casdf.index.map(lambda x: RNA.fold(casdf.varseq162[x][-27-int(casdf.exonlength[x]):-27-int(casdf.exonlength[x])+24])[1])
casdf['exoncenter']=casdf.index.map(lambda x: RNA.fold(casdf.varseq162[x][-27-int(casdf.exonlength[x])/2-12:-27-int(casdf.exonlength[x])/2+12])[1])
casdf['exon3']=casdf.index.map(lambda x: RNA.fold(casdf.varseq162[x][-27-24:-27])[1])
casdf['donor']=casdf.index.map(lambda x: RNA.fold(casdf.varseq162[x][-27-12:-15])[1])
casdf['intron2']=casdf.index.map(lambda x: RNA.fold(casdf.varseq162[x][-27:-3])[1])

casdf['exonstart_varseq']=162-casdf.exonlength.astype(int)
casdf['exonend_varseq']=162

casdf.to_pickle('./dataframes/casdf.pkl')

casdf[['subset','commonname','exonlength','exonstart_varseq','exonend_varseq','maxent5','maxent3',
       'intron1','acceptor','exon5','exoncenter','exon3','donor','intron2',
'levelratiospl','psi','wtratio','normratio','normpsi',
'changes','combination','first_ss','second_ss','firstss','secondss','first_SF','second_SF',
'rnareads_effective','fraction_canonical','varseq']].to_csv('./tables/Table_casdf.csv')

casdf[['exonstart_varseq','exonend_varseq',
'levelratiospl','varseq']].to_csv('../rawdata_FORREPOSITORY/casdf.csv')


########## five #################


rna_condition=pd.read_pickle('../rawdata/five/rna_from_unbiased_mapping.pkl')
rna_condition=rna_condition.drop('levelratiospl', axis=1)
rna_condition_final=analysis_functions.prepare_rnadf_five(rna_condition)
rna_condition_final['maxent5first']=rna_condition_final.index.map(lambda x: maxent.score5(rna_condition_final.varseq162[x][48:57]))
rna_condition_final['maxent5second']=rna_condition_final.index.map(lambda x: maxent.score5(rna_condition_final.varseq162[x][48+int(rna_condition_final.diff_nt[x]):57+int(rna_condition_final.diff_nt[x])]))
rna_condition_final['maxentdiff']=rna_condition_final.index.map(lambda x: rna_condition_final.maxent5second[x]-rna_condition_final.maxent5first[x])

rna_condition_final['exon']=rna_condition_final.index.map(lambda x: RNA.fold(rna_condition_final.varseq162[x][51-24:51])[1])
rna_condition_final['donor1']=rna_condition_final.index.map(lambda x: RNA.fold(rna_condition_final.varseq162[x][51-12:51+12])[1])
rna_condition_final['alt5']=rna_condition_final.index.map(lambda x: RNA.fold(rna_condition_final.varseq162[x][51:51+24])[1])
rna_condition_final['altcenter']=rna_condition_final.index.map(lambda x: RNA.fold(rna_condition_final.varseq162[x][51+int(rna_condition_final.diff_nt[x])/2-12:51+int(rna_condition_final.diff_nt[x])/2+12])[1])
rna_condition_final['alt3']=rna_condition_final.index.map(lambda x: RNA.fold(rna_condition_final.varseq162[x][51+int(rna_condition_final.diff_nt[x])-24:51+int(rna_condition_final.diff_nt[x])])[1])
rna_condition_final['donor2']=rna_condition_final.index.map(lambda x: RNA.fold(rna_condition_final.varseq162[x][51+int(rna_condition_final.diff_nt[x])-12:51+int(rna_condition_final.diff_nt[x])+12])[1])
rna_condition_final['intron']=rna_condition_final.index.map(lambda x: RNA.fold(rna_condition_final.varseq162[x][51+int(rna_condition_final.diff_nt[x]):51+int(rna_condition_final.diff_nt[x])+24])[1])

rna_condition_final.to_pickle('../rawdata/five/five_from_unbiased_mapping.pkl')

fivedf=pd.read_pickle('../rawdata/five/five_from_unbiased_mapping.pkl')
fivedf['rnaperdna'] = fivedf.index.map(lambda x: np.log2(fivedf.loc[x, 'rnareads']/fivedf.loc[x, 'number_reads']) if (fivedf.number_reads[x]>0) else np.nan)
fivedf.loc[fivedf.rnaperdna.dropna().index,'rnaperdna']=scipy.stats.zscore(fivedf.rnaperdna.dropna())
fivedf['wtrnaperdna']=fivedf.commonname.apply(lambda x: fivedf[(fivedf.subset=='five_filtered') & \
        (fivedf.commonname == x)].rnaperdna.dropna().mean())
fivedf['normrnaperdna'] = fivedf.index.map(lambda x: fivedf.loc[x, 'rnaperdna'] - fivedf.loc[x, 'wtrnaperdna'])
fivedf['wtwav']=fivedf.commonname.apply(lambda x: fivedf[(fivedf.subset=='five_filtered') & \
        (fivedf.commonname == x)&(fivedf.smoothednumberofpeaks==1)].wav123.dropna().mean())
fivedf['normwav'] = fivedf.index.map(lambda x: fivedf.loc[x, 'wav123'] - fivedf.loc[x, 'wtwav'])


meanplusnoise=fivedf[(fivedf.smoothednumberofpeaks==1)&(fivedf.number_reads>100)&\
                 (fivedf.fraction_canonical>0.3)][['wav123','rnaperdna','noisestrengthlog']].dropna()

gam=pygam.LinearGAM(pygam.s(0,lam=1, n_splines=10)+pygam.l(1,lam=1)).fit(meanplusnoise[['wav123','rnaperdna']], meanplusnoise.noisestrengthlog)

pred=gam.predict(meanplusnoise[['wav123','rnaperdna']])

meanplusnoise['noisegampred']=pd.Series(pred, index=meanplusnoise.index)
meanplusnoise['noiseresgam']=meanplusnoise['noisestrengthlog']-meanplusnoise['noisegampred']

fivedf['noiseresgam']=meanplusnoise['noiseresgam']

fivedf['wtnoise_gam']=fivedf.commonname.apply(lambda x: fivedf[(fivedf.subset=='five_filtered') & \
        (fivedf.commonname == x)&(fivedf.number_reads>200)&(fivedf.fraction_canonical>0.5)].noiseresgam.dropna().mean())
fivedf['normnoise_gam'] = fivedf.index.map(lambda x: fivedf.loc[x, 'noiseresgam'] - fivedf.loc[x, 'wtnoise_gam'])

fivedf.to_pickle('./dataframes/fivedf.pkl')

fivedf[['subset','commonname','diff_nt','maxent5first','maxent5second',
        'exon','donor1','alt5','altcenter','alt3','donor2','intron',
'levelratiospl','wtratio','normratio','wav123','normwav','rnaperdna','normrnaperdna',
'changes','first_ss','second_ss','firstss','secondss',
'rnareads_effective','number_reads','fraction_canonical','smoothednumberofpeaks','rawnumberofpeaks',
'noisestrengthlog','noiseresgam','normnoise_gam','varseq']].to_csv('./tables/Table_fivedf.csv')

fivedf['donor1start']=51
fivedf['donor2start']=51+fivedf.diff_nt
fivedf[['donor1start','donor2start',
'levelratiospl','wav123','varseq']].to_csv('../rawdata_FORREPOSITORY/fivedf.csv')

########## three ################

rna_condition=pd.read_pickle('../rawdata/three/rna_from_unbiased_mapping.pkl')
rna_condition_final=analysis_functions.prepare_rnadf_three(rna_condition)
rna_condition_final['maxent3first']=rna_condition_final.index.map(lambda x: maxent.score3(rna_condition_final.varseq162[x][56:79]))
rna_condition_final['maxent3second']=rna_condition_final.index.map(lambda x: maxent.score3(rna_condition_final.varseq162[x][56+int(rna_condition_final.diff_nt[x]):79+int(rna_condition_final.diff_nt[x])]))
rna_condition_final['maxentdiff']=rna_condition_final.index.map(lambda x: rna_condition_final.maxent3second[x]-rna_condition_final.maxent3first[x])

rna_condition_final['intron']=rna_condition_final.index.map(lambda x: RNA.fold(rna_condition_final.varseq162[x][76-24:76])[1])
rna_condition_final['acceptor1']=rna_condition_final.index.map(lambda x: RNA.fold(rna_condition_final.varseq162[x][76-12:76+12])[1])
rna_condition_final['alt5']=rna_condition_final.index.map(lambda x: RNA.fold(rna_condition_final.varseq162[x][76:76+24])[1])
rna_condition_final['altcenter']=rna_condition_final.index.map(lambda x: RNA.fold(rna_condition_final.varseq162[x][76+int(rna_condition_final.diff_nt[x])/2-12:76+int(rna_condition_final.diff_nt[x])/2+12])[1])
rna_condition_final['alt3']=rna_condition_final.index.map(lambda x: RNA.fold(rna_condition_final.varseq162[x][76+int(rna_condition_final.diff_nt[x])-24:76+int(rna_condition_final.diff_nt[x])])[1])
rna_condition_final['acceptor2']=rna_condition_final.index.map(lambda x: RNA.fold(rna_condition_final.varseq162[x][76+int(rna_condition_final.diff_nt[x])-12:76+int(rna_condition_final.diff_nt[x])+12])[1])
rna_condition_final['exon']=rna_condition_final.index.map(lambda x: RNA.fold(rna_condition_final.varseq162[x][76+int(rna_condition_final.diff_nt[x]):76+int(rna_condition_final.diff_nt[x])+24])[1] if (int(rna_condition_final.diff_nt[x])<47) else np.nan)


rna_condition_final.to_pickle('./dataframes/threedf.pkl')

rna_condition_final[['subset','commonname','diff_nt','maxent3first','maxent3second',
                     'intron','acceptor1','alt5','altcenter','alt3','acceptor2','exon',
'levelratiospl','wtratio','normratio',
'changes','first_ss','second_ss','firstss','secondss',
'rnareads_effective','fraction_canonical','varseq']].to_csv('./tables/Table_threedf.csv')

threedf['acc1start']=76
threedf['acc2start']=76+threedf.diff_nt
threedf[['acc1start','acc2start',
'levelratiospl','varseq']].to_csv('../rawdata_FORREPOSITORY/threedf.csv')



#%%
'''
get all n numbers
'''

with open(martin + 'combined_analysis/nnumbers_IR.csv', 'a') as f:
    for subset in irdf.subset.unique():
        pd.DataFrame([' ',subset]).to_csv(f)
        irdf[(irdf.subset==subset)&(irdf.levelratiospl.isnull()==False)].changes.value_counts().to_csv(f)
        irdf[(irdf.subset==subset)&(irdf.levelratiospl.isnull()==False)].firstss.value_counts().to_csv(f)
        irdf[(irdf.subset==subset)&(irdf.levelratiospl.isnull()==False)].secondss.value_counts().to_csv(f)
        irdf[(irdf.subset==subset)&(irdf.levelratiospl.isnull()==False)].first_ss.value_counts().to_csv(f)
        irdf[(irdf.subset==subset)&(irdf.levelratiospl.isnull()==False)].second_ss.value_counts().to_csv(f)
        irdf[(irdf.subset==subset)&(irdf.levelratiospl.isnull()==False)].first_SF.value_counts().to_csv(f)
        irdf[(irdf.subset==subset)&(irdf.levelratiospl.isnull()==False)].second_SF.value_counts().to_csv(f)


with open(martin + 'combined_analysis/nnumbers_cas.csv', 'a') as f:
    for subset in casdf.subset.unique():
        pd.DataFrame([' ',subset]).to_csv(f)
        casdf[(casdf.subset==subset)&(casdf.levelratiospl_from_control.isnull()==False)].changes.value_counts().to_csv(f)
        casdf[(casdf.subset==subset)&(casdf.levelratiospl_from_control.isnull()==False)].firstss.value_counts().to_csv(f)
        casdf[(casdf.subset==subset)&(casdf.levelratiospl_from_control.isnull()==False)].secondss.value_counts().to_csv(f)
        casdf[(casdf.subset==subset)&(casdf.levelratiospl_from_control.isnull()==False)].first_ss.value_counts().to_csv(f)
        casdf[(casdf.subset==subset)&(casdf.levelratiospl_from_control.isnull()==False)].second_ss.value_counts().to_csv(f)
        casdf[(casdf.subset==subset)&(casdf.levelratiospl_from_control.isnull()==False)].first_SF.value_counts().to_csv(f)
        casdf[(casdf.subset==subset)&(casdf.levelratiospl_from_control.isnull()==False)].second_SF.value_counts().to_csv(f)


with open(martin + 'combined_analysis/nnumbers_five.csv', 'a') as f:
    for subset in fivedf.subset.unique():
        pd.DataFrame([' ',subset]).to_csv(f)
        fivedf[(fivedf.subset==subset)&(fivedf.levelratiospl.isnull()==False)].changes.value_counts().to_csv(f)
        fivedf[(fivedf.subset==subset)&(fivedf.levelratiospl.isnull()==False)].firstss.value_counts().to_csv(f)
        fivedf[(fivedf.subset==subset)&(fivedf.levelratiospl.isnull()==False)].secondss.value_counts().to_csv(f)
        fivedf[(fivedf.subset==subset)&(fivedf.levelratiospl.isnull()==False)].first_ss.value_counts().to_csv(f)
        fivedf[(fivedf.subset==subset)&(fivedf.levelratiospl.isnull()==False)].second_ss.value_counts().to_csv(f)
        fivedf[(fivedf.subset==subset)&(fivedf.levelratiospl.isnull()==False)].first_SF.value_counts().to_csv(f)
        fivedf[(fivedf.subset==subset)&(fivedf.levelratiospl.isnull()==False)].second_SF.value_counts().to_csv(f)


with open(martin + 'combined_analysis/nnumbers_three.csv', 'a') as f:
    for subset in threedf.subset.unique():
        pd.DataFrame([' ',subset]).to_csv(f)
        threedf[(threedf.subset==subset)&(threedf.levelratiospl.isnull()==False)].changes.value_counts().to_csv(f)
        threedf[(threedf.subset==subset)&(threedf.levelratiospl.isnull()==False)].firstss.value_counts().to_csv(f)
        threedf[(threedf.subset==subset)&(threedf.levelratiospl.isnull()==False)].secondss.value_counts().to_csv(f)
        threedf[(threedf.subset==subset)&(threedf.levelratiospl.isnull()==False)].first_ss.value_counts().to_csv(f)
        threedf[(threedf.subset==subset)&(threedf.levelratiospl.isnull()==False)].second_ss.value_counts().to_csv(f)
        threedf[(threedf.subset==subset)&(threedf.levelratiospl.isnull()==False)].first_SF.value_counts().to_csv(f)
        threedf[(threedf.subset==subset)&(threedf.levelratiospl.isnull()==False)].second_SF.value_counts().to_csv(f)





