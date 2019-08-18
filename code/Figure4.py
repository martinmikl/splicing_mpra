#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 10:54:53 2019

@author: martinm
"""


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats
import forprediction_ir
import forprediction_cas
import forprediction_five
import forprediction_three
import scipy.io as sio
import pickle

sns.set_context('poster')
import os

os.mkdir('./figures/Fig4')
os.mkdir('./figures/FigS4')
os.mkdir('./figures/FigS6')

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


#%%

ireval=pd.read_pickle('./dataframes/ml/ir/ireval.pkl')
        
df_features_ireval=pd.read_pickle('./dataframes/ml/ir/Xs_ireval.pkl')

forprediction_ir.make_prediction_all_models(df_features_ireval, ireval[ireval.fraction_canonical>0.3].levelratiospl, \
                                              filename='Fig4_ireval_prediction_from_model')


caseval=pd.read_pickle('./dataframes/ml/cas/caseval.pkl')

df_features_caseval=pd.read_pickle('./dataframes/ml/cas/Xs_caseval.pkl')

forprediction_cas.make_prediction_all_models(df_features_caseval, caseval[caseval.fraction_canonical>0.5].levelratiospl, \
                                              filename='Fig4_caseval_prediction_from_model')


fiveeval=pd.read_pickle('./dataframes/ml/five/fiveeval.pkl')

df_features_fiveeval=pd.read_pickle('./dataframes/ml/five/Xs_fiveeval.pkl')

forprediction_five.make_prediction_all_models(df_features_fiveeval, fiveeval[fiveeval.fraction_canonical>0.5].levelratiospl, \
                                              filename='Fig4_fiveeval_prediction_from_model')


threeeval=pd.read_pickle('./dataframes/ml/three/threeeval.pkl')

df_features_threeeval=pd.read_pickle('./dataframes/ml/three/Xs_threeeval.pkl')

forprediction_three.make_prediction_all_models(df_features_threeeval, threeeval[threeeval.fraction_canonical>0.3].levelratiospl, \
                                              filename='Fig4_threeeval_prediction_from_model')



#%% on data from Rosenberg et al., 2015

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


# take only variants with canonical site 

a5ss=pd.DataFrame(np.array([A5SS_data[:,0],A5SS_data[:,44],A5SS_data[:,79]]).transpose(), index=range(len(A5SS_data)), columns=['pos0','pos44','pos79'])
a5ss['rnareads']=pd.Series(A5SS_data_sum, index=a5ss.index)

#testindexes0_44=random.sample(a5ss[(a5ss.pos0>0)&(a5ss.pos44>0)&(a5ss.rnareads>100)].index,1000)
#testindexes0_79=random.sample(a5ss[(a5ss.pos0>0)&(a5ss.pos79>0)&(a5ss.rnareads>100)].index,1000)

#

with open('./data/Cell2015_A5SS_1000random_pos0pos44.csv', 'w') as f:
    for i in testindexes0_44:
        f.write(str(i)+','+ context[:52]+A5SS_seqs[i].upper()+\
                context[-9:]+',51,95\n')

df_features_A5SS1000random_pos0pos44=forprediction_five.make_features_five('./data/Cell2015_A5SS_1000random_pos0pos44.csv')
df_features_A5SS1000random_pos0pos44.to_pickle('./dataframes/ml/five/Xs_A5SS_1000random_pos0pos44.pkl')

with open('./data/Cell2015_A5SS_1000random_pos0pos79.csv', 'w') as f:
    for i in testindexes0_79:
        f.write(str(i)+','+ context[:52]+A5SS_seqs[i].upper()+\
                context[-9:]+',51,130\n')

df_features_A5SS1000random_pos0pos79=forprediction_five.make_features_five('./data/Cell2015_A5SS_1000random_pos0pos79.csv')
df_features_A5SS1000random_pos0pos79.to_pickle('./dataframes/ml/five/Xs_A5SS_1000random_pos0pos79.pkl')




df_features_A5SS1000random_pos0pos44=pd.read_pickle('./dataframes/ml/five/Xs_A5SS_1000random_pos0pos44.pkl')

a5ss['logratio']=np.log2(a5ss.pos44/a5ss.pos0)

forprediction_five.make_prediction_all_models(df_features_A5SS1000random_pos0pos44, a5ss.loc[df_features_A5SS1000random_pos0pos44.index,'logratio'], \
                           filename='Cell2015_prediction_from_model_1000random_pos0pos44')

mdl=pickle.load(open('./ml/models/fivemodel_logratio_trained.sav', 'rb'))

pred=mdl.predict(df_features_A5SS1000random_pos0pos44)

f=plt.figure(figsize=(3,3))
plt.scatter(a5ss.loc[df_features_A5SS1000random_pos0pos44.index,'logratio'], pred, alpha=0.2)
plt.title('Pearson r = '+'{:.2g}'.format(scipy.stats.pearsonr(a5ss.loc[df_features_A5SS1000random_pos0pos44.index,'logratio'], pred)[0])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(scipy.stats.spearmanr(a5ss.loc[df_features_A5SS1000random_pos0pos44.index,'logratio'], pred)[0]), fontsize=14)
plt.xlim(-10,10)
plt.ylim(-10,10)
f.savefig('./figures/FigS6/FigS6A_cell2015_predict_logratio0_44_fivemodel_logratio_trained.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


#
#
#df_features_A5SS1000random_pos0pos79=pd.read_pickle('./dataframes/ml/five/Xs_A5SS_1000random_pos0pos79.pkl')
#
#a5ss['logratio']=np.log2(a5ss.pos79/a5ss.pos0)
#
#forprediction_five.make_prediction(df_features_A5SS1000random_pos0pos79, a5ss.loc[df_features_A5SS1000random_pos0pos79.index,'logratio'])
##r2=
##pearsonr=
#forprediction_five.make_prediction_all_models(df_features_A5SS1000random_pos0pos79, a5ss.loc[df_features_A5SS1000random_pos0pos79.index,'logratio'], \
#                           filename='Cell2015_prediction_from_model_1000random_pos0pos79')
#
#
#df_features_A5SS1000random_pos44pos79=pd.read_pickle('./dataframes/ml/five/Xs_A5SS_1000random_pos44pos79.pkl')
#
#a5ss['logratio']=np.log2(a5ss.pos79/a5ss.pos44)
#
#forprediction_five.make_prediction(df_features_A5SS1000random_pos44pos79, a5ss.loc[df_features_A5SS1000random_pos44pos79.index,'logratio'])
#r2=
#pearsonr=
#forprediction_five.make_prediction_all_models(df_features_A5SS1000random_pos44pos79, a5ss.loc[df_features_A5SS1000random_pos44pos79.index,'logratio'], \
#                           filename='Cell2015_prediction_from_model_1000random_pos44pos79')
#
#### score cv
#
#forprediction_five.make_gbr_cvpredict(df_features_A5SS1000random_pos0pos44, \
#            a5ss.loc[df_features_A5SS1000random_pos0pos44.index,'logratio'],\
#                    filename='A5SS1000random_pos0pos44')
#forprediction_five.make_gbr_cvpredict(df_features_A5SS1000random_pos0pos79, \
#            a5ss.loc[df_features_A5SS1000random_pos0pos79.index,'logratio'],\
#                    filename='A5SS1000random_pos0pos79')
#forprediction_five.make_gbr_cvpredict(df_features_A5SS1000random_pos44pos79, \
#            a5ss.loc[df_features_A5SS1000random_pos44pos79.index,'logratio'],\
#                    filename='A5SS1000random_pos44pos79')
#

#%% Secondary structure

### IR

second = irdf[(irdf.subset == 'irsecvar')]

def splitsf(ss):
    if (ss=='endogenous'):
        return 'no change'
    else:
        return ss.split(' ')[0]


def splitpos(ss):
    if (ss=='endogenous'):
        return 'no change'
    else:
        return ss.split(' ')[-1]


second['firstsschange']=second.firstss.apply(lambda x: splitsf(x))
second['firstsspos']=second.firstss.apply(lambda x: splitpos(x))
second['secondsschange']=second.secondss.apply(lambda x: splitsf(x))
second['secondsspos']=second.secondss.apply(lambda x: splitpos(x))
second['location']=second.index.map(lambda x: second.secondsspos[x] if (second.firstsspos[x]=='no change') else second.firstsspos[x])
second['change']=second.index.map(lambda x: second.secondsschange[x] if (second.firstsspos[x]=='no change') else second.firstsschange[x])

f=plt.figure(figsize=(6.5,3))
ax=sns.pointplot(data=second, x='location', y='normratio', \
               hue='change', hue_order=['comp','rev','comp2','rev2'], \
            palette='Paired', join=False, dodge=0.4, errwidth=2, scale=0.5)
ax.set_ylim(-4.3,2)
ax.set_xticklabels(['donor-12', 'donor+15', 'acc-37', 'acc+3'])
ax.set_xlabel('location of sequence change')
ax.set_ylabel('normalized ratio [log2]')
plt.axhline(y=0, color='gray',linewidth=2)
ax.legend(ax.get_legend_handles_labels()[0],('comp 1','rev-comp 1','comp 2','rev-comp 2'), \
          bbox_to_anchor=(0.8,1.07),ncol=2, fontsize=12, columnspacing=0.5)
f.savefig('./figures/FigS4/FigS4B_secondary_ir_overview_pointplot.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

scipy.stats.mannwhitneyu(second[(second.location=='[start-12:start-3]')&(second.change=='comp')].normratio.dropna(),\
                         second[(second.location=='[start-12:start-3]')&(second.change=='rev')].normratio.dropna())  
# p=0.000466

scipy.stats.mannwhitneyu(second[(second.location=='[start+15:start+24]')&(second.change=='comp')].normratio.dropna(),\
                         second[(second.location=='[start+15:start+24]')&(second.change=='rev')].normratio.dropna())  
# p=0.0066

scipy.stats.mannwhitneyu(second[(second.location=='[end-37:end-28]')&(second.change=='comp')].normratio.dropna(),\
                         second[(second.location=='[end-37:end-28]')&(second.change=='rev')].normratio.dropna())  
# p=0.068

scipy.stats.mannwhitneyu(second[(second.location=='[end+3:end+12]')&(second.change=='comp')].normratio.dropna(),\
                         second[(second.location=='[end+3:end+12]')&(second.change=='rev')].normratio.dropna())  
# p=0.0164

scipy.stats.mannwhitneyu(second[(second.location=='[start-12:start-3]')&(second.change=='comp2')].normratio.dropna(),\
                         second[(second.location=='[start-12:start-3]')&(second.change=='rev2')].normratio.dropna())  
# p=0.00425

scipy.stats.mannwhitneyu(second[(second.location=='[start+15:start+24]')&(second.change=='comp2')].normratio.dropna(),\
                         second[(second.location=='[start+15:start+24]')&(second.change=='rev2')].normratio.dropna())  
# p=0.033

scipy.stats.mannwhitneyu(second[(second.location=='[end-37:end-28]')&(second.change=='comp2')].normratio.dropna(),\
                         second[(second.location=='[end-37:end-28]')&(second.change=='rev2')].normratio.dropna())  
# p=0.24

scipy.stats.mannwhitneyu(second[(second.location=='[end+3:end+12]')&(second.change=='comp2')].normratio.dropna(),\
                         second[(second.location=='[end+3:end+12]')&(second.change=='rev2')].normratio.dropna())  
# p=0.023


### FIVE


second = fivedf[(fivedf.subset == 'fivesecvar')]
'''
# fix problem in design file
for i in second[(second.firstss=='exonrev at [51-12:51-3]')].index:
    print(second.commonname[i])
    exonrev=str(Seq(fivedf.varseq162[i][51-24:51-15]).reverse_complement())
    exoncomp=str(Seq(fivedf.varseq162[i][51-24:51-15]).complement())
    if (second.varseq162[i][51-12:51-3]==exonrev):
        second.firstss[i]='exonrev at [51-12:51-3]'
    elif (second.varseq162[i][51-12:51-3]==exoncomp):
        second.firstss[i]='exoncomp at [51-12:51-3]'
    else:
        print('unknown')
        print(second.varseq162[i][25:51])

for i in second[(second.secondss=='altexoncomp at [51+diff-14:51+diff-5]')].index:
    print(second.commonname[i])
    altexonrev=str(Seq(fivedf.varseq162[i][51+int(fivedf.diff_nt[i])-26:51+int(fivedf.diff_nt[i])-17]).reverse_complement()) 
    altexoncomp=str(Seq(fivedf.varseq162[i][51+int(fivedf.diff_nt[i])-26:51+int(fivedf.diff_nt[i])-17]).complement())
    if (second.varseq162[i][51+int(fivedf.diff_nt[i])-14:51+int(fivedf.diff_nt[i])-5]==altexonrev):
        print('rev')
        print(second.varseq162[i][51:51+int(fivedf.diff_nt[i])])
        second.secondss[i]='altexonrev at [51+diff-14:51+diff-5]'
    elif (second.varseq162[i][51+int(fivedf.diff_nt[i])-14:51+int(fivedf.diff_nt[i])-5]==altexoncomp):
        print('comp')
        print(second.varseq162[i][51:51+int(fivedf.diff_nt[i])])
        second.secondss[i]='altexoncomp at [51+diff-14:51+diff-5]'
    else:
        print('unknown')
        print(second.varseq162[i][51:51+int(fivedf.diff_nt[i])])
'''
#analysis
def splitsf(ss):
    if (ss=='endogenous'):
        return 'no change'
    else:
        return ss.split(' ')[0]
    
def splitpos(ss):
    if (ss=='endogenous'):
        return 'no change'
    else:
        return ss.split(' ')[-1]

second['firstsschange']=second.firstss.apply(lambda x: splitsf(x))
second['firstsspos']=second.firstss.apply(lambda x: splitpos(x))
second['secondsschange']=second.secondss.apply(lambda x: splitsf(x))
second['secondsspos']=second.secondss.apply(lambda x: splitpos(x))

second['location']=second.index.map(lambda x: second.secondsspos[x] if (second.firstsspos[x]=='no change') else second.firstsspos[x])
second['change']=second.index.map(lambda x: second.secondsschange[x] if (second.firstsspos[x]=='no change') else second.firstsschange[x])

f=plt.figure(figsize=(6.5,3))
ax1=sns.pointplot(data=second[['location','normratio','change']].dropna(), x='location', y='normratio', order=['[51-12:51-3]','[51+15:51:51+24]','[51+diff-14:51+diff-5]','[51+diff+15:51+diff+24]'],\
               hue='change', hue_order=['comp','rev','comp2','rev2'], \
            palette='Paired', join=False, dodge=0.4, errwidth=2, scale=0.5)
ax1.set_ylim(-4,8)
ax1.set_xticklabels(['first-12', 'first+15','second-14', 'second+15'])
ax1.set_xlabel('location of sequence change')
ax1.set_ylabel('normalized ratio [log2]')
plt.axhline(y=0, color='gray',linewidth=2)
ax1.legend(ax1.get_legend_handles_labels()[0],('comp 1', 'rev-comp 1', 'comp 2', 'rev-comp 2'), \
          bbox_to_anchor=(1,1),ncol=2, fontsize=12, columnspacing=0.2,handletextpad=0.2)
f.savefig('./figures/FigS4/FigS4B_secondary_five_overview_pointplot.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


scipy.stats.mannwhitneyu(second[(second.location=='[51-12:51-3]')&(second.change=='comp')].normratio.dropna(),\
                         second[(second.location=='[51-12:51-3]')&(second.change=='rev')].normratio.dropna())  
# p=0.00032

scipy.stats.mannwhitneyu(second[(second.location=='[51+15:51:51+24]')&(second.change=='comp')].normratio.dropna(),\
                         second[(second.location=='[51+15:51:51+24]')&(second.change=='rev')].normratio.dropna())  
# p=0.00394

scipy.stats.mannwhitneyu(second[(second.location=='[51+diff-14:51+diff-5]')&(second.change=='comp')].normratio.dropna(),\
                         second[(second.location=='[51+diff-14:51+diff-5]')&(second.change=='rev')].normratio.dropna())  
# p=0.167

scipy.stats.mannwhitneyu(second[(second.location=='[51+diff+15:51+diff+24]')&(second.change=='comp')].normratio.dropna(),\
                         second[(second.location=='[51+diff+15:51+diff+24]')&(second.change=='rev')].normratio.dropna())  
# p=0.00176

scipy.stats.mannwhitneyu(second[(second.location=='[51-12:51-3]')&(second.change=='comp2')].normratio.dropna(),\
                         second[(second.location=='[51-12:51-3]')&(second.change=='rev2')].normratio.dropna())  
# p=0.0012

scipy.stats.mannwhitneyu(second[(second.location=='[51+15:51:51+24]')&(second.change=='comp2')].normratio.dropna(),\
                         second[(second.location=='[51+15:51:51+24]')&(second.change=='rev2')].normratio.dropna())  
# p=0.0087

scipy.stats.mannwhitneyu(second[(second.location=='[51+diff-14:51+diff-5]')&(second.change=='comp2')].normratio.dropna(),\
                         second[(second.location=='[51+diff-14:51+diff-5]')&(second.change=='rev2')].normratio.dropna())  
# p=0.287

scipy.stats.mannwhitneyu(second[(second.location=='[51+diff+15:51+diff+24]')&(second.change=='comp2')].normratio.dropna(),\
                         second[(second.location=='[51+diff+15:51+diff+24]')&(second.change=='rev2')].normratio.dropna())  
# p=0.0177


