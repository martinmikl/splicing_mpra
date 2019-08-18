#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 19:32:11 2019

@author: martinm
"""


import pandas as pd
import forprediction_irnoise
import random
import os
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_predict
from sklearn import feature_selection
from sklearn.ensemble import GradientBoostingRegressor
from sklearn import metrics
import scipy.stats

#%%
os.mkdir('./dataframes/ml/noise')
#%%

irdf=pd.read_excel('./tables/TableS5.xlsx')

irdf.columns=['libindex','subset','name2','intronlength','intronstart_varseq','intronend_varseqnew','maxent5','maxent3',
      'exon1','donor','intron5','introncenter','intron3','acceptor','exon2',
'levelratiospl','wtratio','normratio','wav_stats','normwav','rnaperdna','normrnaperdna',
'changes','combination','first_ss','second_ss','firstss','secondss',
'rnareads_effective','number_reads','fraction_canonical','smoothednumberofpeaks','rawnumberofpeaks','noisestrengthlogwstd','noiseres_linear',
'noiseresgam','normnoise_linear','normnoise_gam','varseq']

irdf.set_index('libindex', inplace=True)
irdf['varseq162']=irdf.varseq.apply(lambda x: str(x[30:-18]).upper())



# drop indexes for 
indexes=[]
samevarseq=[]
for var, group in irdf[(irdf.smoothednumberofpeaks==1)&(irdf.number_reads>100)&(irdf.fraction_canonical>0.3)].groupby(by='varseq162'):
	if (len(group)>1):
		samevarseq.append(var)
		indexes.append(list(irdf[(irdf.smoothednumberofpeaks==1)&(irdf.number_reads>100)&(irdf.fraction_canonical>0.3)&(irdf.varseq162==var)].index))
		
indexestodrop=[]
for igroup in indexes:
	indexestodrop=indexestodrop + random.sample(igroup, len(igroup)-1)
		
irdfforml=irdf[(irdf.smoothednumberofpeaks==1)&(irdf.number_reads>100)&(irdf.fraction_canonical>0.3)].drop(indexestodrop) 

#####

irnoiseml,irnoiseeval, irnoisemly,irnoiseevaly=train_test_split(irdfforml, irdfforml['noiseresgam'], \
												test_size=0.25, random_state=0)

irnoiseml.to_pickle('./dataframes/ml/noise/irnoiseml.pkl')
irnoiseeval.to_pickle('./dataframes/ml/noise/irnoiseeval.pkl')
irnoisemly.to_pickle('./dataframes/ml/noise/irnoisemly.pkl')
irnoiseevaly.to_pickle('./dataframes/ml/noise/irnoiseevaly.pkl')


#%%
irnoiseml=pd.read_pickle('./dataframes/ml/noise/irnoiseml.pkl')
irnoiseeval=pd.read_pickle('./dataframes/ml/noise/irnoiseeval.pkl')
irnoisemly=pd.read_pickle('./dataframes/ml/noise/irnoisemly.pkl')
irnoiseevaly=pd.read_pickle('./dataframes/ml/noise/irnoiseevaly.pkl')


#%% Intron retention

### make features with entire varseq
irdfforml['intronstart_wholevarseq']=irdfforml.intronstart_varseq.apply(lambda x: int(x)+30)
irdfforml['intronend_wholevarseq']=irdfforml.intronend_varseqnew.apply(lambda x: int(x)+30)

irdfforml[['varseq','intronstart_wholevarseq','intronend_wholevarseq']].to_csv('./data/ir_noise_allfeatures.csv', header=None)

# Make features

#df_features_ireval10=forprediction_ir.make_features_ir('./data/ireval_fortest_first10.csv')
df_features_irnoiseall=forprediction_irnoise.make_features_ir('./data/ir_noise_allfeatures.csv')

df_features_irnoiseall.to_pickle('./dataframes/ml/ir/Xs_allfeatures_ir_fornoise.pkl')
df_features_irnoiseall=pd.read_pickle('./dataframes/ml/ir/Xs_allfeatures_ir_fornoise.pkl')

####
# hyperparameter optimization for the different models

cols={'maxent':['maxent5','maxent3'],
      'secondary':['exon1','donor','intron5','introncenter','intron3','acceptor','exon2'],
      'hexamers':[x for x in df_features_irnoiseall.columns if (x.split('_')[-1] in ['exon1','intron','exon2'])&(len(x.split('_'))==2)],
      'motifs':[x for x in df_features_irnoiseall.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                'all features':':'}

def optimize_parameters_gbr_noise(Xs, Ys):
    cvpredscores=pd.DataFrame()
    for bestlr in [0.001,0.01, 0.05]:
        for bestnest in [100,300,400]:
            for bestdepth in [3,4,5,6,7]:
                clf=GradientBoostingRegressor(learning_rate=bestlr, \
                    n_estimators=bestnest, max_depth=bestdepth)
                cvpredgb = cross_val_predict(clf, Xs, \
                                             Ys, cv=5)
                cvpredscores.loc[str(str(bestlr) + "_" + str(bestnest) + '_' + str(bestdepth)),'r2score']=\
                                 metrics.r2_score(Ys,cvpredgb)
                cvpredscores.loc[str(str(bestlr) + "_" + str(bestnest) + '_' + str(bestdepth)),'pearsonr']=\
                                 scipy.stats.pearsonr(Ys,cvpredgb)[0]
                cvpredscores.loc[str(str(bestlr) + "_" + str(bestnest) + '_' + str(bestdepth)),'pearsonp']=\
                                 scipy.stats.pearsonr(Ys,cvpredgb)[1]
    return cvpredscores

Xs_train=df_features_irnoiseall.loc[irnoiseml.index].dropna()
Ys_train=irdf.loc[Xs_train.index, 'noiseresgam']
maxent_cvpredscores=optimize_parameters_gbr_noise(Xs_train[cols['maxent']], Ys_train)
#bestlr=0.01
#bestnest=100
#bestdepth=6
maxent_cvpredscores.to_pickle('./ml/noiseoptimization/maxent_cvpredscores.pkl')

secondary_cvpredscores=optimize_parameters_gbr_noise(Xs_train[cols['secondary']], Ys_train)
#bestlr=0.01
#bestnest=300
#bestdepth=4
secondary_cvpredscores.to_pickle('./ml/noiseoptimization/secondary_cvpredscores.pkl')

motifs_cvpredscores=optimize_parameters_gbr_noise(Xs_train[cols['motifs']], Ys_train)
#bestlr=0.01
#bestnest=300
#bestdepth=5
motifs_cvpredscores.to_pickle('./ml/noiseoptimization/motifs_cvpredscores.pkl')

clf=GradientBoostingRegressor(learning_rate=0.01, n_estimators=300,max_depth=5)
clf_sel=feature_selection.SelectFromModel(clf).fit(Xs_train[cols['motifs']], Ys_train)
Xs=Xs_train[cols['motifs']]
Xsmotifs_sel=Xs[Xs.columns[clf_sel.get_support(indices=True)]]
Xsmotifs_sel.to_pickle('./dataframes/ml/noise/Xs_motifs_IR_motifs_sel_noiseresgam.pkl')

Xsmotifs_sel=pd.read_pickle('./dataframes/ml/noise/Xs_motifs_IR_motifs_sel_noiseresgam.pkl')

hexamers_cvpredscores=optimize_parameters_gbr_noise(Xs_train[cols['hexamers']], Ys_train)
#bestlr=0.01
#bestnest=100
#bestdepth=6
hexamers_cvpredscores.to_pickle('./ml/noiseoptimization/hexamers_cvpredscores.pkl')

clf=GradientBoostingRegressor(learning_rate=0.01, n_estimators=100,max_depth=6)
clf_sel=feature_selection.SelectFromModel(clf).fit(Xs_train[cols['hexamers']], Ys_train)
Xs=Xs_train[cols['hexamers']]
Xshexamers_sel=Xs[Xs.columns[clf_sel.get_support(indices=True)]]
Xshexamers_sel.to_pickle('./dataframes/ml/noise/Xs_hexamers_IR_motifs_sel_noiseresgam.pkl')

Xshexamers_sel=pd.read_pickle('./dataframes/ml/noise/Xs_hexamers_IR_motifs_sel_noiseresgam.pkl')



clf=GradientBoostingRegressor(learning_rate=0.01, n_estimators=300,max_depth=5)
clf_sel=feature_selection.SelectFromModel(clf).fit(Xs_train, Ys_train)
Xsall_sel=Xs_train[Xs_train.columns[clf_sel.get_support(indices=True)]]
Xsall_sel.to_pickle('./dataframes/ml/noise/Xs_allfeatures_IR_sel_noiseresgam.pkl')

Xsall_sel=pd.read_pickle('./dataframes/ml/noise/Xs_allfeatures_IR_sel_noiseresgam.pkl')

allsel=Xsall_sel.columns


allfeatures_sel_cvpredscores=optimize_parameters_gbr_noise(Xs_train[allsel], Ys_train)
#bestlr=0.05
#bestnest=100
#bestdepth=6
allfeatures_sel_cvpredscores.to_pickle('./ml/noiseoptimization/allfeatures_sel_cvpredscores.pkl')


Zs=pd.DataFrame(irdf.loc[Xs_train.index, 'wav_stats'], index=Xs_train.index)
Zs['measured']=pd.Series(Ys_train)
Zs.dropna(inplace=True)

onlysplicingvalue_sel_cvpredscores=optimize_parameters_gbr_noise(Zs['wav_stats'].to_frame(), Zs['measured'])
#bestlr=0.05
#bestnest=100
#bestdepth=2
onlysplicingvalue_sel_cvpredscores.to_pickle('./ml/noiseoptimization/onlysplicingvalue_sel_cvpredscores.pkl')

#%%
### Train models

### Train models - with selected features
Xs_train=df_features_irnoiseall.loc[irnoiseml.index].dropna()
Ys_train=irdf.loc[Xs_train.index, 'noiseresgam']
forprediction_irnoise.train_all_models(Xs_train, irdf.loc[Xs_train.index, 'wav_stats'], Ys_train, 'irmodel_allsel_noiseresgam_trained')

### Test models - with selected features

Xs_test=df_features_irnoiseall.loc[irnoiseeval.index].dropna()
Ys_test=irdf.loc[Xs_test.index, 'noiseresgam']
forprediction_irnoise.make_prediction_all_models(Xs_test, irdf.loc[Xs_test.index, 'wav_stats'], Ys_test,\
                mode='irmodel_allsel_noiseresgam_trained', filename='irnoisegam_allsel_prediction_from_model')


