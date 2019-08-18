#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 14:36:15 2019

@author: martinm
"""

import pandas as pd
import forprediction_ir
import forprediction_cas
import forprediction_five
import forprediction_three
import seaborn as sns
from pylab import *
import scipy.io as sio
import scipy
import random
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.cross_validation import train_test_split

sns.set_context('poster')

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

#%% intronretention
### Run the predictor on the 10% test set put aside before building a model

irevaly=pd.read_pickle('./dataframes/ml/ir/irevaly.pkl')
ireval=irdf.loc[irevaly.index]
        
ireval['intronstart_wholevarseq']=ireval.intronstart_varseq.apply(lambda x: int(x)+30)
ireval['intronend_wholevarseq']=ireval.intronend_varseqnew.apply(lambda x: int(x)+30)

ireval[['varseq','intronstart_wholevarseq','intronend_wholevarseq']].to_csv('./data/ireval_fortest.csv', header=None)
# save separate file with only 10 sequences for testing
ireval[['varseq','intronstart_wholevarseq','intronend_wholevarseq']].head(n=10).to_csv('./data/ireval_fortest_first10.csv', header=None)
ireval[['varseq','intronstart_wholevarseq','intronend_wholevarseq']].head(n=3).to_csv('../../code_NCOMMSrev/code/data/ireval_first3.csv', header=None)
ireval[['varseq','intronstart_wholevarseq','intronend_wholevarseq']].head(n=1).to_csv('../../code_NCOMMSrev/code/data/ireval_first1.csv', header=None)

# Make features

#df_features_ireval10=forprediction_ir.make_features_ir('./data/ireval_fortest_first10.csv')
df_features_ireval=forprediction_ir.make_features_ir('./data/ireval_fortest.csv')

df_features_ireval.to_pickle('./dataframes/ml/ir/Xs_ireval.pkl')

# Run predictions

df_features_ireval=pd.read_pickle('./dataframes/ml/ir/Xs_ireval.pkl')

forprediction_ir.make_prediction(df_features_ireval, ireval.levelratiospl)
forprediction_ir.make_prediction_all_models(df_features_ireval, ireval.levelratiospl, \
                                              filename='ireval_prediction_from_model')


#%% cassette exons
### Run the predictor on the 10% test set put aside before building a model

casevaly=pd.read_pickle('./dataframes/ml/cas/casevaly.pkl')
caseval=casdf.loc[casevaly.index]

caseval['acceptor']=caseval.exonlength.apply(lambda x: 117-int(x))
caseval['donor']=117
        

caseval[['varseq162','acceptor','donor']].to_csv('./data/caseval_fortest.csv', header=None)
# save separate file with only 10 sequences for testing
caseval[['varseq162','acceptor','donor']].head(n=10).to_csv('./data/caseval_fortest_first10.csv', header=None)
caseval[['varseq162','acceptor','donor']].head(n=3).to_csv('./data/caseval_first3.csv', header=None)

# Make features

#df_features_caseval10=forprediction_cas.make_features_cas('./data/caseval_fortest_first10.csv')
df_features_caseval=forprediction_cas.make_features_cas('./data/caseval_fortest.csv')

df_features_caseval.to_pickle('./dataframes/ml/cas/Xs_caseval.pkl')

# Run predictions

df_features_caseval=pd.read_pickle('./dataframes/ml/cas/Xs_caseval.pkl')

forprediction_cas.make_prediction(df_features_caseval, caseval.levelratiospl)
forprediction_cas.make_prediction_all_models(df_features_caseval, caseval.levelratiospl, \
                                              filename='caseval_prediction_from_model')



############ Test on other data sets

### Test with vexseq data

# prepare file

# Data from Adamson et al., 2018, Genome Biol; formatted version used by Cheng et al., 2019, Genome Biol (MMSplice_test_pred.csv file obtained from github.com/gagneurlab/MMSplice_paper)

vexhep=pd.read_csv('../rawdata/vexseq/MMSplice_test_pred.csv')

vexhep['varpsi']=vexhep.index.map(lambda x: vexhep.HepG2_ref_psi[x] + vexhep.HepG2_delta_psi[x])
vexhep['wtlogratio']=vexhep.HepG2_ref_psi.apply(lambda x: np.log2(x/(100-x)) if x>0 else -10)
vexhep['varlogratio']=vexhep.varpsi.apply(lambda x: np.log2(x/(100-x)) if x>0 else -10)
vexhep['deltalogratio']=vexhep['varlogratio']-vexhep['wtlogratio']
vexhep['exonstart']=100
vexhep['exonend']=vexhep.width.apply(lambda x: int(x)+100)

vexhep[['ALT_SEQ','exonstart','exonend']].to_csv('./data/vexseq_test_forpredictor.csv', header=None)
vexhep[['ALT_SEQ','exonstart','exonend']].head(n=10).to_csv('./data/vexseq_test_forpredictor_first10.csv', header=None)
vexhep[['REF_SEQ','exonstart','exonend']].drop_duplicates().to_csv('./data/vexseq_test_forpredictor_wt.csv', header=None)

# Predict vexseq based on FACSseq model


df_features_cas_vexseq=forprediction_cas.make_features_cas('./data/vexseq_test_forpredictor.csv')
df_features_cas_vexseq.to_pickle('./dataframes/ml/cas/Xs_vexseq_test.pkl')

df_features_cas_vexseq_wt=forprediction_cas.make_features_cas('./data/vexseq_test_forpredictor_wt.csv')
df_features_cas_vexseq_wt.to_pickle('./dataframes/ml/cas/Xs_vexseq_test_wt.pkl')

#predict based on model trained on our data

df_features_cas_vexseq=pd.read_pickle('./dataframes/ml/cas/Xs_vexseq_test.pkl')

#make_prediction(df_features, vexhep.varpsi.head(n=10), mode='psi')
forprediction_cas.make_prediction(df_features_cas_vexseq, vexhep.varpsi/100, mode='psi')
forprediction_cas.make_prediction_all_models(df_features_cas_vexseq, vexhep.varpsi/100, 
                                             mode='casmodel_psi_trained',
                                             filename='vexseq_prediction_from_model_psi')

forprediction_cas.make_prediction(df_features_cas_vexseq, vexhep.varlogratio,
                                             filename='vexseq_prediction_from_model')

forprediction_cas.make_prediction_all_models(df_features_cas_vexseq, vexhep.varlogratio,
                                             filename='vexseq_prediction_from_model')

forprediction_cas.make_prediction_all_models(df_features_cas_vexseq_wt, vexhep.wtlogratio,
                                             filename='vexseq_prediction_wts_from_model')

forprediction_cas.make_gbr_cvpredict(df_features_cas_vexseq, vexhep.varlogratio,
                                             filename='vexseq_5foldCV')

### make figures

mdl=pickle.load(open('./ml/models/casmodel_logratio_trained.sav', 'rb'))

mutpred=pd.Series(mdl.predict(df_features_cas_vexseq), index=df_features_cas_vexseq.index, name='pred_mut').to_frame()
mutpred=mutpred.join(vexhep['varlogratio'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.varlogratio, mutpred.pred_mut, alpha=0.2)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['pred_mut','varlogratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['pred_mut','varlogratio']].corr(method='spearman').values[0][1]), fontsize=14)
plt.xlim(-11,5)
plt.ylim(-12,5)
f.savefig('./figures/FigS4/FigS4A_vexseq_predict_mutvalues_casmodel_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


mdl=pickle.load(open('./ml/models/casmodel_logratio_trained_maxentsecondary.sav', 'rb'))
cols=['maxent5','maxent3','intron1','acceptor','exon5','exoncenter','exon3','donor','intron2']

mutpred=pd.Series(mdl.predict(df_features_cas_vexseq[cols]), index=df_features_cas_vexseq.index, name='pred_mut').to_frame()
mutpred=mutpred.join(vexhep['varlogratio'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.varlogratio, mutpred.pred_mut, alpha=0.2)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['pred_mut','varlogratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['pred_mut','varlogratio']].corr(method='spearman').values[0][1]), fontsize=14)
plt.xlim(-11,5)
plt.ylim(-15,5)
f.savefig('./figures/FigS4/FigS4A_vexseq_predict_mutvalues_casmodel_logratio_trained_maxentsecondary.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()



### predict delta


mdl=pickle.load(open('./ml/models/casmodel_psi_trained.sav', 'rb'))

wtpred=pd.Series(mdl.predict(df_features_cas_vexseq_wt), index=df_features_cas_vexseq_wt.index, name='pred_wt').to_frame()
mutpred=pd.Series(mdl.predict(df_features_cas_vexseq), index=df_features_cas_vexseq.index, name='pred_mut').to_frame()

wtpred['REF_SEQ']=wtpred.index.map(lambda x: vexhep.REF_SEQ[x])
mutpred['REF_SEQ']=mutpred.index.map(lambda x: vexhep.REF_SEQ[x])
wtpred.set_index('REF_SEQ', inplace=True)

mutpred['wt_pred']=mutpred.REF_SEQ.apply(lambda x: wtpred.loc[x, 'pred_wt'])

mutpred['deltapred']=mutpred.pred_mut - mutpred.wt_pred
mutpred=mutpred.join(vexhep['HepG2_delta_psi'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.deltapred, mutpred.HepG2_delta_psi, alpha=0.2)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','HepG2_delta_psi']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','HepG2_delta_psi']].corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/FigS4/FigS4B_vexseq_predict_delta_casmodel_psi_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


mdl=pickle.load(open('./ml/models/casmodel_logratio_trained.sav', 'rb'))

wtpred=pd.Series(mdl.predict(df_features_cas_vexseq_wt), index=df_features_cas_vexseq_wt.index, name='pred_wt').to_frame()
mutpred=pd.Series(mdl.predict(df_features_cas_vexseq), index=df_features_cas_vexseq.index, name='pred_mut').to_frame()

wtpred['REF_SEQ']=wtpred.index.map(lambda x: vexhep.REF_SEQ[x])
mutpred['REF_SEQ']=mutpred.index.map(lambda x: vexhep.REF_SEQ[x])
wtpred.set_index('REF_SEQ', inplace=True)

mutpred['wt_pred']=mutpred.REF_SEQ.apply(lambda x: wtpred.loc[x, 'pred_wt'])

mutpred['deltapred']=mutpred.pred_mut - mutpred.wt_pred
mutpred=mutpred.join(vexhep['deltalogratio'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.deltapred, mutpred.deltalogratio, alpha=0.2)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','deltalogratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','deltalogratio']].corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/FigS4/FigS4B_vexseq_predict_delta_casmodel_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


mdl=pickle.load(open('./ml/models/casmodel_psi_trained_maxentsecondary.sav', 'rb'))
cols=['maxent5','maxent3','intron1','acceptor','exon5','exoncenter','exon3','donor','intron2']

wtpred=pd.Series(mdl.predict(df_features_cas_vexseq_wt[cols]), index=df_features_cas_vexseq_wt.index, name='pred_wt').to_frame()
mutpred=pd.Series(mdl.predict(df_features_cas_vexseq[cols]), index=df_features_cas_vexseq.index, name='pred_mut').to_frame()

wtpred['REF_SEQ']=wtpred.index.map(lambda x: vexhep.REF_SEQ[x])
mutpred['REF_SEQ']=mutpred.index.map(lambda x: vexhep.REF_SEQ[x])
wtpred.set_index('REF_SEQ', inplace=True)

mutpred['wt_pred']=mutpred.REF_SEQ.apply(lambda x: wtpred.loc[x, 'pred_wt'])

mutpred['deltapred']=mutpred.pred_mut - mutpred.wt_pred
mutpred=mutpred.join(vexhep['HepG2_delta_psi'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.deltapred, mutpred.HepG2_delta_psi, alpha=0.2)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','HepG2_delta_psi']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','HepG2_delta_psi']].corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/FigS4/FigS4B_vexseq_predict_delta_casmodel_psi_trained_maxentsecondary.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


mdl=pickle.load(open('./ml/models/casmodel_logratio_trained_maxentsecondary.sav', 'rb'))
cols=['maxent5','maxent3','intron1','acceptor','exon5','exoncenter','exon3','donor','intron2']

wtpred=pd.Series(mdl.predict(df_features_cas_vexseq_wt[cols]), index=df_features_cas_vexseq_wt.index, name='pred_wt').to_frame()
mutpred=pd.Series(mdl.predict(df_features_cas_vexseq[cols]), index=df_features_cas_vexseq.index, name='pred_mut').to_frame()

wtpred['REF_SEQ']=wtpred.index.map(lambda x: vexhep.REF_SEQ[x])
mutpred['REF_SEQ']=mutpred.index.map(lambda x: vexhep.REF_SEQ[x])
wtpred.set_index('REF_SEQ', inplace=True)

mutpred['wt_pred']=mutpred.REF_SEQ.apply(lambda x: wtpred.loc[x, 'pred_wt'])

mutpred['deltapred']=mutpred.pred_mut - mutpred.wt_pred
mutpred=mutpred.join(vexhep['deltalogratio'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.deltapred, mutpred.deltalogratio, alpha=0.2)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','deltalogratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','deltalogratio']].corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/FigS4/FigS4B_vexseq_predict_delta_casmodel_logratio_trained_maxentsecondary.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()




# build model based on VexSeq training data

vexhep_train=pd.read_csv(martin + 'combined_analysis/forrevision/MMSplice_paper/data/vexseq/HepG2_delta_PSI_CAGI_training_sided.csv')

vexhep_train['varpsi']=vexhep_train.index.map(lambda x: vexhep.HepG2_ref_psi[x] + vexhep.HepG2_delta_psi[x])
vexhep_train['varlogratio']=vexhep_train.varpsi.apply(lambda x: np.log2(x/(100-x)) if x>0 else -10)
vexhep_train['exonstart']=100
vexhep_train['exonend']=vexhep.width.apply(lambda x: int(x)+100)

vexhep_train[['ALT_SEQ','exonstart','exonend']].to_csv(martin + 'combined_analysis/forrevision/ml/vexseq_train_forpredictor.csv', header=None)
vexhep_train[['REF_SEQ','exonstart','exonend']].drop_duplicates().to_csv(martin + 'combined_analysis/forrevision/ml/vexseq_train_forpredictor_wt.csv', header=None)

df_train_features=make_features_cas(martin + 'combined_analysis/forrevision/ml/vexseq_train_forpredictor.csv')
df_train_features.to_pickle(martin + 'combined_analysis/forrevision/ml/Xs_vexseq_train.pkl')

df_train_features_wt=forprediction_cas.make_features_cas(martin + 'combined_analysis/forrevision/ml/vexseq_train_forpredictor_wt.csv')
df_train_features_wt.to_pickle(martin + 'combined_analysis/forrevision/ml/Xs_vexseq_train_wt.pkl')

df_train_features_vexseq=pd.read_pickle(martin + 'combined_analysis/forrevision/ml/Xs_vexseq_train.pkl')

forprediction_cas.train_all_models(df_train_features_vexseq, vexhep_train.varpsi, 'gbrmodel_trained_on_vexseq_chr1_8')
forprediction_cas.train_all_models(df_train_features_vexseq, vexhep_train.varlogratio, 'gbrmodel_logratio_trained_on_vexseq_chr1_8')

forprediction_cas.make_prediction(df_features_cas_vexseq, vexhep.varpsi, mode='gbrmodel_trained_on_vexseq_chr1_8')

forprediction_cas.make_prediction_all_models(df_features_cas_vexseq, vexhep.varpsi, mode='gbrmodel_trained_on_vexseq_chr1_8')
forprediction_cas.make_prediction_all_models(df_features_cas_vexseq, vexhep.varlogratio, 
                mode='gbrmodel_logratio_trained_on_vexseq_chr1_8', filename='vexseq_prediction_from_model_trained_on_vexseqdata')


# predict based on VexSeq-trained model


mdl=pickle.load(open('./ml/models/gbrmodel_logratio_trained_on_vexseq_chr1_8.sav', 'rb'))

wtpred=pd.Series(mdl.predict(df_features_cas_vexseq_wt), index=df_features_cas_vexseq_wt.index, name='pred_wt').to_frame()
mutpred=pd.Series(mdl.predict(df_features_cas_vexseq), index=df_features_cas_vexseq.index, name='pred_mut').to_frame()

wtpred['REF_SEQ']=wtpred.index.map(lambda x: vexhep.REF_SEQ[x])
mutpred['REF_SEQ']=mutpred.index.map(lambda x: vexhep.REF_SEQ[x])
wtpred.set_index('REF_SEQ', inplace=True)

mutpred['wt_pred']=mutpred.REF_SEQ.apply(lambda x: wtpred.loc[x, 'pred_wt'])

mutpred['deltapred']=mutpred.pred_mut - mutpred.wt_pred
mutpred=mutpred.join(vexhep['deltalogratio'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.deltapred, mutpred.deltalogratio, alpha=0.2)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','deltalogratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','deltalogratio']].corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/FigS4/FigS4X_vexseq_predict_delta_gbrmodel_logratio_trained_on_vexseq_chr1_8.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

#%%
### MAPSy
# Data from Soemedi et al., 2017, Nat Genetics; version used by Cheng et al., 2019, Genome Biol (mapsy.txt file obtained from github.com/gagneurlab/MMSplice_paper)

mapsy=pd.read_csv('../rawdata/mapsy/mapsy.txt')

mapsy['exonstart']=mapsy.index.map(lambda x: mapsy.WT_SEQ[x].find(mapsy.acceptor[x])+len(mapsy.acceptor[x]))
mapsy['exonend']=mapsy.index.map(lambda x: mapsy.WT_SEQ[x].find(mapsy.donor[x]))

mapsy['psi_wt']=mapsy.index.map(lambda x: mapsy.loc[x,'Vivo_WT_Spliced']/float(mapsy.loc[x,'Vivo_WT_Input'])*100)
mapsy['psi_mut']=mapsy.index.map(lambda x: mapsy.loc[x,'Vivo_MUT_Spliced']/float(mapsy.loc[x,'Vivo_MUT_Input'])*100)
mapsy['logpsiwt']=mapsy.psi_wt.apply(lambda x: np.log2(x) if x>0 else np.nan)
mapsy['logpsimut']=mapsy.psi_mut.apply(lambda x: np.log2(x) if x>0 else np.nan)

mapsy['psi_wt_vitro']=mapsy.index.map(lambda x: mapsy.loc[x,'Vitro_WT_Spliced']/float(mapsy.loc[x,'Vitro_WT_Input'])*100)
mapsy['psi_mut_vitro']=mapsy.index.map(lambda x: mapsy.loc[x,'Vitro_MUT_Spliced']/float(mapsy.loc[x,'Vitro_MUT_Input'])*100)
mapsy['logpsiwt_vitro']=mapsy.psi_wt_vitro.apply(lambda x: np.log2(x) if x>0 else np.nan)
mapsy['logpsimut_vitro']=mapsy.psi_mut_vitro.apply(lambda x: np.log2(x) if x>0 else np.nan)

mapsy['logpsimut_vitro']=mapsy.psi_mut_vitro.apply(lambda x: np.log2(x) if x>0 else np.nan)

mapsy[['MUT_SEQ','exonstart','exonend']].to_csv(martin + 'combined_analysis/forrevision/ml/mapsy_forpredictor.csv', header=None)

mapsy[['WT_SEQ','exonstart','exonend']].drop_duplicates().to_csv(martin + 'combined_analysis/forrevision/ml/mapsy_wts_forpredictor.csv', header=None)


df_features_mapsy=forprediction_cas.make_features_cas(martin + 'combined_analysis/forrevision/ml/mapsy_forpredictor.csv')
df_features_mapsy.to_pickle('./dataframes/ml/cas/Xs_mapsy.pkl')
df_features_mapsy=pd.read_pickle(martin + 'combined_analysis/forrevision/ml/Xs_mapsy.pkl')


df_features_mapsy_wt=forprediction_cas.make_features_cas(martin + 'combined_analysis/forrevision/ml/mapsy_wts_forpredictor.csv')
df_features_mapsy_wt.to_pickle('./dataframes/ml/cas/Xs_mapsy_wt.pkl')
df_features_mapsy_wt=pd.read_pickle(martin + 'combined_analysis/forrevision/ml/Xs_mapsy_wt.pkl')

#make_prediction(df_features, vexhep.varpsi.head(n=10), mode='psi')
forprediction_cas.make_prediction(df_features_mapsy, mapsy.logpsimut)
forprediction_cas.make_prediction_all_models(df_features_mapsy, mapsy.logpsimut, filename='mapsy_invivo_prediction_from_our_model')
#R=0.343
forprediction_cas.make_prediction(df_features_mapsy, mapsy.logpsimut_vitro)
forprediction_cas.make_prediction_all_models(df_features_mapsy, mapsy.logpsimut_vitro, filename='mapsy_invitro_prediction_from_our_model')
#R=0.334
forprediction_cas.make_prediction_all_models(df_features_mapsy_wt, mapsy.logpsiwt, filename='mapsy_wts_invivo_prediction_from_our_model')
#R=0.289
forprediction_cas.make_prediction_all_models(df_features_mapsy_wt, mapsy.logpsiwt_vitro, filename='mapsy_wts_invitro_prediction_from_our_model')
#R=0.296

### make figures

mdl=pickle.load(open('./ml/models/casmodel_logratio_trained.sav', 'rb'))


mutpred=pd.Series(mdl.predict(df_features_mapsy), index=df_features_mapsy.index, name='pred_mut').to_frame()
mutpred=mutpred.join(mapsy['logpsimut'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.logpsimut, mutpred.pred_mut, alpha=0.1)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['pred_mut','logpsimut']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['pred_mut','logpsimut']].corr(method='spearman').values[0][1]), fontsize=14)
plt.xlim(-4,10)
f.savefig('./figures/FigS4/FigS4A_mapsy_invivo_predict_mutvalues_casmodel_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


mutpred=pd.Series(mdl.predict(df_features_mapsy), index=df_features_mapsy.index, name='pred_mut').to_frame()
mutpred=mutpred.join(mapsy['logpsimut_vitro'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.logpsimut_vitro, mutpred.pred_mut, alpha=0.1)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['pred_mut','logpsimut_vitro']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['pred_mut','logpsimut_vitro']].corr(method='spearman').values[0][1]), fontsize=14)
plt.xlim(-4,12)
f.savefig('./figures/FigS4/FigS4A_mapsy_invitro_predict_mutvalues_casmodel_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


### predict delta

mdl=pickle.load(open('./ml/models/casmodel_logratio_trained.sav', 'rb'))

wtpred=pd.Series(mdl.predict(df_features_mapsy_wt), index=df_features_mapsy_wt.index, name='pred_wt').to_frame()
mutpred=pd.Series(mdl.predict(df_features_mapsy), index=df_features_mapsy.index, name='pred_mut').to_frame()

wtpred['REF_SEQ']=wtpred.index.map(lambda x: mapsy.WT_SEQ[x])
mutpred['REF_SEQ']=mutpred.index.map(lambda x: mapsy.WT_SEQ[x])
wtpred.set_index('REF_SEQ', inplace=True)

mutpred['wt_pred']=mutpred.REF_SEQ.apply(lambda x: wtpred.loc[x, 'pred_wt'])

mutpred['deltapred']=mutpred.pred_mut - mutpred.wt_pred
mutpred=mutpred.join(mapsy['log2_vivo_ratio'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.log2_vivo_ratio, mutpred.deltapred, alpha=0.1)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','log2_vivo_ratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','log2_vivo_ratio']].corr(method='spearman').values[0][1]), fontsize=14)
plt.xlim(-10,10)
plt.ylim(-10,10)
f.savefig('./figures/FigS4/FigSB_mapsy_predict_delta_casmodel_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


mutpred=mutpred.join(mapsy['log2_vitro_ratio'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.log2_vitro_ratio, mutpred.deltapred, alpha=0.1)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','log2_vitro_ratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','log2_vitro_ratio']].corr(method='spearman').values[0][1]), fontsize=14)
plt.xlim(-10,10)
plt.ylim(-10,10)
f.savefig('./figures/FigS4/FigS4B_mapsy_invitro_predict_delta_casmodel_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


mapsy_train, mapsy_test = train_test_split(wtpred, test_size=0.5, random_state=0)

gbr_model_mapsy_mut_vivo=GradientBoostingRegressor(learning_rate=0.2, n_estimators=300,\
                                    max_depth=5).fit(df_features_mapsy.loc[mapsy[[x in mapsy_train.index for x in mapsy.WT_SEQ]]['logpsimut'].dropna().index],
                                               mapsy[[x in mapsy_train.index for x in mapsy.WT_SEQ]].logpsimut.dropna()) 

gbr_model_mapsy_mut_vitro=GradientBoostingRegressor(learning_rate=0.2, n_estimators=300,\
                                    max_depth=5).fit(df_features_mapsy.loc[mapsy[[x in mapsy_train.index for x in mapsy.WT_SEQ]]['logpsimut_vitro'].dropna().index],
                                               mapsy[[x in mapsy_train.index for x in mapsy.WT_SEQ]].logpsimut_vitro.dropna()) 

gbr_model_mapsy_wt_vivo=GradientBoostingRegressor(learning_rate=0.2, n_estimators=300,\
                                    max_depth=5).fit(df_features_mapsy_wt.loc[mapsy[[x in mapsy_train.index for x in mapsy.WT_SEQ]]['logpsiwt'].dropna().index].dropna(),
                                               mapsy.loc[df_features_mapsy_wt.loc[mapsy[[x in mapsy_train.index for x in mapsy.WT_SEQ]]['logpsiwt'].dropna().index].dropna().index,'logpsiwt'].dropna()) 

gbr_model_mapsy_wt_vitro=GradientBoostingRegressor(learning_rate=0.2, n_estimators=300,\
                                    max_depth=5).fit(df_features_mapsy_wt.loc[mapsy[[x in mapsy_train.index for x in mapsy.WT_SEQ]]['logpsiwt_vitro'].dropna().index].dropna(),
                                               mapsy.loc[df_features_mapsy_wt.loc[mapsy[[x in mapsy_train.index for x in mapsy.WT_SEQ]]['logpsiwt_vitro'].dropna().index].dropna().index,'logpsiwt_vitro'].dropna()) 

df_features_mapsy_wt.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna()


### test on model trained on in vivo mapsy data - separate wt contexts for train and test

wtpred=pd.Series(gbr_model_mapsy_mut_vivo.predict(df_features_mapsy_wt.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna()), \
        index=df_features_mapsy_wt.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna().index, name='pred_wt').to_frame()
mutpred=pd.Series(gbr_model_mapsy_mut_vivo.predict(df_features_mapsy.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna()), \
        index=df_features_mapsy.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna().index, name='pred_mut').to_frame()

wtpred['REF_SEQ']=wtpred.index.map(lambda x: mapsy.WT_SEQ[x])
mutpred['REF_SEQ']=mutpred.index.map(lambda x: mapsy.WT_SEQ[x])
wtpred.set_index('REF_SEQ', inplace=True)

mutpred['wt_pred']=mutpred.REF_SEQ.apply(lambda x: wtpred.loc[x, 'pred_wt'])

mutpred['deltapred']=mutpred.pred_mut - mutpred.wt_pred
mutpred=mutpred.join(mapsy['log2_vivo_ratio'])
mutpred=mutpred.join(mapsy['log2_vitro_ratio'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.log2_vivo_ratio, mutpred.deltapred, alpha=0.1)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','log2_vivo_ratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','log2_vivo_ratio']].corr(method='spearman').values[0][1]), fontsize=14)
plt.xlim(-10,10)
plt.ylim(-10,10)
f.savefig('./figures/FigS4/FigS4B_mapsy_predict_delta_gbr_model_mapsy_mut_vivo_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.log2_vitro_ratio, mutpred.deltapred, alpha=0.1)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','log2_vitro_ratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','log2_vitro_ratio']].corr(method='spearman').values[0][1]), fontsize=14)
plt.xlim(-10,10)
plt.ylim(-10,10)
f.savefig('./figures/FigS4/FigS4B_mapsy_vitro_predict_delta_gbr_model_mapsy_mut_vivo_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()



### test on model trained on in vitro mapsy data - separate wt contexts for train and test

wtpred=pd.Series(gbr_model_mapsy_mut_vitro.predict(df_features_mapsy_wt.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna()), \
        index=df_features_mapsy_wt.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna().index, name='pred_wt').to_frame()
mutpred=pd.Series(gbr_model_mapsy_mut_vitro.predict(df_features_mapsy.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna()), \
        index=df_features_mapsy.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna().index, name='pred_mut').to_frame()

wtpred['REF_SEQ']=wtpred.index.map(lambda x: mapsy.WT_SEQ[x])
mutpred['REF_SEQ']=mutpred.index.map(lambda x: mapsy.WT_SEQ[x])
wtpred.set_index('REF_SEQ', inplace=True)

mutpred['wt_pred']=mutpred.REF_SEQ.apply(lambda x: wtpred.loc[x, 'pred_wt'])

mutpred['deltapred']=mutpred.pred_mut - mutpred.wt_pred
mutpred=mutpred.join(mapsy['log2_vivo_ratio'])
mutpred=mutpred.join(mapsy['log2_vitro_ratio'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.log2_vivo_ratio, mutpred.deltapred, alpha=0.1)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','log2_vivo_ratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','log2_vivo_ratio']].corr(method='spearman').values[0][1]), fontsize=14)
plt.xlim(-10,10)
plt.ylim(-10,10)
f.savefig('./figures/FigS4/FigS4B_mapsy_predict_delta_gbr_model_mapsy_mut_vitro_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.log2_vitro_ratio, mutpred.deltapred, alpha=0.1)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','log2_vitro_ratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','log2_vitro_ratio']].corr(method='spearman').values[0][1]), fontsize=14)
plt.xlim(-10,10)
plt.ylim(-10,10)
f.savefig('./figures/FigS4/FigS4B_mapsy_vitro_predict_delta_gbr_model_mapsy_mut_vitro_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()




### test on model trained on in vivo mapsy data, wt only - separate wt contexts for train and test

wtpred=pd.Series(gbr_model_mapsy_wt_vivo.predict(df_features_mapsy_wt.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna()), \
        index=df_features_mapsy_wt.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna().index, name='pred_wt').to_frame()
mutpred=pd.Series(gbr_model_mapsy_wt_vivo.predict(df_features_mapsy.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna()), \
        index=df_features_mapsy.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna().index, name='pred_mut').to_frame()

wtpred['REF_SEQ']=wtpred.index.map(lambda x: mapsy.WT_SEQ[x])
mutpred['REF_SEQ']=mutpred.index.map(lambda x: mapsy.WT_SEQ[x])
wtpred.set_index('REF_SEQ', inplace=True)

mutpred['wt_pred']=mutpred.REF_SEQ.apply(lambda x: wtpred.loc[x, 'pred_wt'])

mutpred['deltapred']=mutpred.pred_mut - mutpred.wt_pred
mutpred=mutpred.join(mapsy['log2_vivo_ratio'])
mutpred=mutpred.join(mapsy['log2_vitro_ratio'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.log2_vivo_ratio, mutpred.deltapred, alpha=0.1)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','log2_vivo_ratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','log2_vivo_ratio']].corr(method='spearman').values[0][1]), fontsize=14)
plt.xlim(-10,10)
plt.ylim(-10,10)
f.savefig('./figures/FigS4/FigS4B_mapsy_predict_delta_gbr_model_mapsy_wt_vivo_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.log2_vitro_ratio, mutpred.deltapred, alpha=0.1)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','log2_vitro_ratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','log2_vitro_ratio']].corr(method='spearman').values[0][1]), fontsize=14)
plt.xlim(-10,10)
plt.ylim(-10,10)
f.savefig('./figures/FigS4/FigS4B_mapsy_vitro_predict_delta_gbr_model_mapsy_wt_vivo_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()



### test on model trained on in vitro mapsy data, wt only - separate wt contexts for train and test

wtpred=pd.Series(gbr_model_mapsy_wt_vitro.predict(df_features_mapsy_wt.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna()), \
        index=df_features_mapsy_wt.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna().index, name='pred_wt').to_frame()
mutpred=pd.Series(gbr_model_mapsy_wt_vitro.predict(df_features_mapsy.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna()), \
        index=df_features_mapsy.loc[mapsy[[x in mapsy_test.index for x in mapsy.WT_SEQ]].index].dropna().index, name='pred_mut').to_frame()

wtpred['REF_SEQ']=wtpred.index.map(lambda x: mapsy.WT_SEQ[x])
mutpred['REF_SEQ']=mutpred.index.map(lambda x: mapsy.WT_SEQ[x])
wtpred.set_index('REF_SEQ', inplace=True)

mutpred['wt_pred']=mutpred.REF_SEQ.apply(lambda x: wtpred.loc[x, 'pred_wt'])

mutpred['deltapred']=mutpred.pred_mut - mutpred.wt_pred
mutpred=mutpred.join(mapsy['log2_vivo_ratio'])
mutpred=mutpred.join(mapsy['log2_vitro_ratio'])

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.log2_vivo_ratio, mutpred.deltapred, alpha=0.1)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','log2_vivo_ratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','log2_vivo_ratio']].corr(method='spearman').values[0][1]), fontsize=14)
plt.xlim(-10,10)
plt.ylim(-10,10)
f.savefig('./figures/FigS4/FigS4B_mapsy_predict_delta_gbr_model_mapsy_wt_vitro_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

f=plt.figure(figsize=(3,3))
plt.scatter(mutpred.log2_vitro_ratio, mutpred.deltapred, alpha=0.1)
plt.title('Pearson r = '+'{:.2g}'.format(mutpred[['deltapred','log2_vitro_ratio']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(mutpred[['deltapred','log2_vitro_ratio']].corr(method='spearman').values[0][1]), fontsize=14)
plt.xlim(-10,10)
plt.ylim(-10,10)
f.savefig('./figures/FigS4/FigS4B_mapsy_vitro_predict_delta_gbr_model_mapsy_wt_vitro_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()



### check performance by CV

from sklearn.cross_validation import cross_val_predict

gbr=GradientBoostingRegressor(learning_rate=0.2, n_estimators=300,\
                                        max_depth=5)

Zs=pd.DataFrame(df_features_mapsy.copy())
Zs['measured']=pd.Series(mapsy.logpsimut)
Zs.dropna(inplace=True)

ypred=cross_val_predict(gbr, Zs.iloc[:,:-1], Zs.iloc[:,-1])
pearsonr(ypred, Zs.iloc[:,-1])
(0.41775414377611253, 1.6131513696683902e-196)

f=plt.figure(figsize=(3,3))
plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
plt.title('$\mathregular{R^2}$ score=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
          '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=10)
f.savefig('./ml/plots/scatter_mapsy_invivo_prediction_from_5foldCV_ourmodelonMAPSyData.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

pd.Series(ypred).to_pickle('./ml/mapsy_prediction_from_5foldCV_ourmodelonMAPSyData.pkl')
ypred=pd.read_pickle('./ml/mapsy_prediction_from_5foldCV_ourmodelonMAPSyData.pkl')

Zs=pd.DataFrame(df_features_mapsy.copy())
Zs['measured']=pd.Series(mapsy.logpsimut_vitro)
Zs.dropna(inplace=True)

ypred=cross_val_predict(gbr, Zs.iloc[:,:-1], Zs.iloc[:,-1])
pearsonr(ypred, Zs.iloc[:,-1])
(0.54462830233975257, 0.0)

pd.Series(ypred).to_pickle('./ml/mapsy_invitro_prediction_from_5foldCV_ourmodelonMAPSyData.pkl')
ypred=pd.read_pickle('./ml/mapsy_invitro_prediction_from_5foldCV_ourmodelonMAPSyData.pkl')

f=plt.figure(figsize=(3,3))
plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
plt.title('$\mathregular{R^2}$ score=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
          '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=10)
f.savefig('./ml/plots/scatter_mapsy_invitro_prediction_from_5foldCV_ourmodelonMAPSyData.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)



gbr=GradientBoostingRegressor(learning_rate=0.2, n_estimators=300,\
                                        max_depth=5)

Zs=pd.DataFrame(df_features_mapsy.copy())
Zs['measured']=pd.Series(mapsy[(mapsy.log2_vivo_ratio>-20)&(mapsy.log2_vivo_ratio<20)].log2_vivo_ratio)
Zs.dropna(inplace=True)

ypred=cross_val_predict(gbr, Zs.iloc[:,:-1], Zs.iloc[:,-1])
pearsonr(ypred, Zs.iloc[:,-1])
(0.41775414377611253, 1.6131513696683902e-196)

f=plt.figure(figsize=(3,3))
plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
plt.title('$\mathregular{R^2}$ score=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
          '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=10)
f.savefig('./ml/plots/scatter_mapsy_invivo_prediction_from_5foldCV_ourmodelonMAPSyData.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

pd.Series(ypred).to_pickle('./ml/mapsy_prediction_from_5foldCV_ourmodelonMAPSyData.pkl')
ypred=pd.read_pickle('./ml/mapsy_prediction_from_5foldCV_ourmodelonMAPSyData.pkl')

Zs=pd.DataFrame(df_features_mapsy.copy())
Zs['measured']=pd.Series(mapsy.logpsimut_vitro)
Zs.dropna(inplace=True)

ypred=cross_val_predict(gbr, Zs.iloc[:,:-1], Zs.iloc[:,-1])
pearsonr(ypred, Zs.iloc[:,-1])
(0.54462830233975257, 0.0)

pd.Series(ypred).to_pickle('./ml/mapsy_invitro_prediction_from_5foldCV_ourmodelonMAPSyData.pkl')
ypred=pd.read_pickle('./ml/mapsy_invitro_prediction_from_5foldCV_ourmodelonMAPSyData.pkl')

f=plt.figure(figsize=(3,3))
plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
plt.title('$\mathregular{R^2}$ score=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
          '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=10)
f.savefig('./ml/plots/scatter_mapsy_invitro_prediction_from_5foldCV_ourmodelonMAPSyData.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


#%%
### Run the predictor on the 10% test set put aside before building a model

fiveevaly=pd.read_pickle('./dataframes/ml/five/fiveevaly.pkl')
fiveeval=fivedf.loc[fiveevaly.index]

fiveeval['donor1']=51
fiveeval['donor2']=fiveeval.diff_nt.apply(lambda x: 51+int(x))
        

fiveeval[['varseq162','donor1','donor2']].to_csv('./data/fiveeval_fortest.csv', header=None)
# save separate file with only 10 sequences for testing
fiveeval[['varseq162','donor1','donor2']].head(n=10).to_csv('./data/fiveeval_fortest_first10.csv', header=None)
fiveeval[['varseq162','donor1','donor2']].head(n=3).to_csv('./data/fiveeval_first3.csv', header=None)

# Make features

#df_features_fiveeval10=forprediction_five.make_features_five('./data/fiveeval_fortest_first10.csv')
df_features_fiveeval=forprediction_five.make_features_five('./data/fiveeval_fortest.csv')

df_features_fiveeval.to_pickle('./dataframes/ml/five/Xs_fiveeval.pkl')

# Run predictions

df_features_fiveeval=pd.read_pickle('./dataframes/ml/five/Xs_fiveeval.pkl')

forprediction_five.make_prediction(df_features_fiveeval, fiveeval.levelratiospl)
forprediction_five.make_prediction_all_models(df_features_fiveeval, fiveeval.levelratiospl, \
                                              filename='fiveeval_prediction_from_model')


### Run the predictor on data from Rosenberg et al., 2015, Cell

data=sio.loadmat('../additional/cell-2015/Reads.mat')
A5SS_data = data['A5SS']
A5SS_reads = np.array(A5SS_data.sum(1)).flatten()
A5SS_data = np.array(A5SS_data.todense())

A5SS_nn = find(A5SS_data.sum(axis=1))
A5SS_reads = A5SS_reads[A5SS_nn]
A5SS_data = A5SS_data[A5SS_nn]
A5SS_data_sum=A5SS_data.sum(axis=1)
A5SS_data = A5SS_data/A5SS_data.sum(axis=1)[:,newaxis]
A5SS_seqs = pd.read_csv('../additional/cell-2015/A5SS_Seqs.csv',index_col=0).Seq[A5SS_nn]

context='aagcagaagaacggcatcaaagtgaacttcaagatccgccacaacatcgaggtgcttggnnnnnnnnnnnnnnnnnnnnnnnnnggtcgacccaggttcgtgnnnnnnnnnnnnnnnnnnnnnnnnngaggtattcttatcaccttcgtggctacagagttt'.upper()

a5ss=pd.DataFrame(np.array([A5SS_data[testindexes,0],A5SS_data[testindexes,44],A5SS_data[testindexes,79]]).transpose(), index=testindexes, columns=['pos0','pos44','pos79'])

a5ss.replace(to_replace=0, value=0.00003, inplace=True)

a5ss['logratio']=np.log2(a5ss.pos44/a5ss.pos0)

make_prediction(df_features_A5SS1000random, a5ss['logratio'])
r2=0.24
pearsonr=0.49
make_prediction_all_models(df_features_A5SS1000random, a5ss['logratio'])

# take only variants with canonical site 

a5ss=pd.DataFrame(np.array([A5SS_data[:,0],A5SS_data[:,44],A5SS_data[:,79]]).transpose(), index=range(len(A5SS_data)), columns=['pos0','pos44','pos79'])
a5ss['rnareads']=pd.Series(A5SS_data_sum, index=a5ss.index)

testindexes0_44=random.sample(a5ss[(a5ss.pos0>0)&(a5ss.pos44>0)&(a5ss.rnareads>100)].index,1000)
testindexes0_79=random.sample(a5ss[(a5ss.pos0>0)&(a5ss.pos79>0)&(a5ss.rnareads>100)].index,1000)
testindexes44_79=random.sample(a5ss[(a5ss.pos79>0)&(a5ss.pos44>0)&(a5ss.rnareads>100)].index,1000)

#

with open('./data/Cell2015_A5SS_1000random_pos0pos44.csv', 'w') as f:
    for i in testindexes0_44:
        f.write(str(i)+','+ context[:52]+A5SS_seqs[i].upper()+\
                context[-9:]+',51,95\n')

df_features_A5SS1000random_pos0pos44=forprediction_five.make_features_five('./data/Cell2015_A5SS_1000random_pos0pos44.csv')
df_features_A5SS1000random_pos0pos44.to_pickle('./dataframes/ml/five/Xs_A5SS_1000random_pos0pos44.pkl')


df_features_A5SS1000random_pos0pos44=pd.read_pickle('./dataframes/ml/five/Xs_A5SS_1000random_pos0pos44.pkl')

a5ss['logratio']=np.log2(a5ss.pos44/a5ss.pos0)

forprediction_five.make_prediction(df_features_A5SS1000random_pos0pos44, a5ss.loc[df_features_A5SS1000random_pos0pos44.index,'logratio'])
forprediction_five.make_prediction_all_models(df_features_A5SS1000random_pos0pos44, a5ss.loc[df_features_A5SS1000random_pos0pos44.index,'logratio'], \
                           filename='Cell2015_prediction_from_model_1000random_pos0pos44')


#
### score cv

forprediction_five.make_gbr_cvpredict(df_features_A5SS1000random_pos0pos44, \
            a5ss.loc[df_features_A5SS1000random_pos0pos44.index,'logratio'],\
                    filename='A5SS1000random_pos0pos44')



#%%
### Run the predictor on the 10% test set put aside before building a model

threeevaly=pd.read_pickle('./dataframes/ml/three/threeevaly.pkl')
threeeval=threedf.loc[threeevaly.index]

threeeval['acceptor1']=76
threeeval['acceptor2']=threeeval.diff_nt.apply(lambda x: 76+int(x))
        

threeeval[['varseq162','acceptor1','acceptor2']].to_csv('./data/threeeval_fortest.csv', header=None)
# save separate file with only 10 sequences for testing
threeeval[['varseq162','acceptor1','acceptor2']].head(n=10).to_csv('./data/threeeval_fortest_first10.csv', header=None)
threeeval[['varseq162','acceptor1','acceptor2']].head(n=3).to_csv('./data/threeeval_first3.csv', header=None)

# Make features

df_features_threeeval10=forprediction_three.make_features_three('./data/threeeval_fortest_first10.csv')
df_features_threeeval=forprediction_three.make_features_three('./data/threeeval_fortest.csv')
 
df_features_threeeval.to_pickle('./dataframes/ml/three/Xs_threeeval.pkl')

# Run predictions

df_features_threeeval=pd.read_pickle('./dataframes/ml/three/Xs_threeeval.pkl')

forprediction_three.make_prediction(df_features_threeeval, threeeval.levelratiospl)
forprediction_three.make_prediction_all_models(df_features_threeeval, threeeval.levelratiospl, \
                                              filename='threeeval_prediction_from_model')


