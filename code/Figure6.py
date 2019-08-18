#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 14:25:31 2019

@author: martinm
"""


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats
import os
import random
import pygam

sns.set_context('poster')

os.mkdir('./figures/Fig6')
os.mkdir('./figures/FigS8')

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
data300=pd.read_pickle('../rawdata/ir_protein/data300.pkl')

ratiocalcbins=pd.read_excel('../rawdata/ir_protein/calculatedGFPmCherryratiosforbins.xlsx')
xval=np.log2(ratiocalcbins.loc[:,'median']*100)
xvalwidth=[]

for i in range(len(xval)):
    if (i==15):
        xvalwidth.append(0.8)
    else:
        xvalwidth.append(xval.loc[i+2] - xval.loc[i+1])


#%%
f=plt.figure(figsize=(4,3))
ax=f.add_subplot(111)
with sns.axes_style('dark'):
    plt.bar(xval, data300.loc[30692,'summedreads'], width=xvalwidth)
ax.xaxis.grid() 
plt.title('mean splicing value: ' + '{:.2f}'.format(irdf.wav_stats[30692]) + \
          '\nsplicing noise strength [log2]: ' + '{:.2f}'.format(irdf.noisestrengthlogwstd[30692]) + \
          '\nnoise residual [log2]: ' + '{:.2f}'.format(irdf.noiseresgam[30692]), fontsize=14)
plt.xlim(0,8)
f.savefig('./figures/FigS8/FigS8A_example_lownoise.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


f=plt.figure(figsize=(4,3))
ax=f.add_subplot(111)
with sns.axes_style('dark'):
    plt.bar(xval, data300.loc[29976,'summedreads'], width=xvalwidth)
ax.xaxis.grid() 
plt.title('mean splicing value: ' + '{:.2f}'.format(irdf.wav_stats[29976]) + \
          '\nsplicing noise strength [log2]: ' + '{:.2f}'.format(irdf.noisestrengthlogwstd[29976]) + \
          '\nnoise residual [log2]: ' + '{:.2f}'.format(irdf.noiseresgam[29976]), fontsize=14)
plt.xlim(0,8)
f.savefig('./figures/FigS8/FigS8A_example_highnoise.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()



f=plt.figure(figsize=(4,3))
ax=f.add_subplot(111)
with sns.axes_style('dark'):
    plt.bar(xval, data300.loc[34440,'summedreads'], width=xvalwidth)
ax.xaxis.grid() 
plt.title('mean splicing value: ' + '{:.2f}'.format(irdf.wav_stats[34440]) + \
          '\nsplicing noise strength [log2]: ' + '{:.2f}'.format(irdf.noisestrengthlogwstd[34440]) + \
          '\nnoise residual [log2]: ' + '{:.2f}'.format(irdf.noiseresgam[34440]), fontsize=14)
plt.xlim(0,8)
f.savefig('./figures/FigS8/FigS8A_example_mediumnoise.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()




#%% Check if noise readout for barcode control sets show less variability than expected by chance

df=irdf[(irdf.smoothednumberofpeaks==1)&(irdf.number_reads>100)&(irdf.fraction_canonical>0.3)]

samevarseqmin5=[]
for var, group in df.groupby(by='varseq162'):
    if (len(group)>4):
        samevarseqmin5.append(var)


dffive=fivedf[(fivedf.smoothednumberofpeaks==1)&(fivedf.number_reads>100)&\
                 (fivedf.fraction_canonical>0.3)]

samevarseqmin5five=[]
for var, group in dffive.groupby(by='varseq162'):
    if (len(group)>4):
        samevarseqmin5five.append(var)

os.mkdir('./figures/FigS6/noiseresofMultipleBarcodes')   

def compare_barcodecontrolnoise_to_random(df, varseqs, noise_measure, size_of_random_set, filename=None):
    '''
    varseqs: list of sequences of library variable regions (excluding the barcode and the primer, length 162)
    random set size: 10000
    '''
    pvals = pd.DataFrame()
    for i in varseqs:
        indexes=list(df[df.varseq162==i].index)
        print(indexes[0])
        print(np.log2(np.var(df.loc[indexes, noise_measure].dropna())))
            
    # test significance by comparing to a random distribution
    
        randvars = []
        for c in range(size_of_random_set):
            randvars.append(np.var(random.sample(df[(df.wav_stats>=df.loc[indexes,'wav_stats'].min()) \
                                         & (df.wav_stats<=df.loc[indexes,'wav_stats'].max())].drop(indexes)[noise_measure].dropna(), len(indexes))))
        
        
        pvals.loc[indexes[0],'genename']=df.loc[indexes[0], 'name2']
        pvals.loc[indexes[0],'pval']=float(sum(ivar <= np.var(df.loc[indexes,noise_measure].dropna()) for ivar in randvars))/float(len(randvars))  ### for NOP16: 0.0283
        
        
        f = plt.figure(figsize=(4,3))
        plt.hist(np.log2(randvars), bins=93, color=sns.xkcd_rgb['medium blue'], linewidth=0)
        plt.axvline(x = np.log2(np.var(df.loc[indexes, noise_measure].dropna())), color=sns.xkcd_rgb['dark blue'])
        plt.xlim(-6,3)
        plt.xticks([-6,-3,0,3])
        if filename!=None:
            f.savefig('./figures/FigS6/noiseresofMultipleBarcodes/' + filename + '_' + str(indexes[0]) + '_vs_randomdistribution.png', \
                  dpi = 300, format='png', bbox_inches='tight', frameon=True)
    return pvals


#

pvals_gam_min5vars30perc=compare_barcodecontrolnoise_to_random(df, samevarseqmin5, 'noiseresgam',10000, 'barcodecontrols_noiseresgam_min5vars30perc_ir')
pvals_gam_min5vars30percfive=compare_barcodecontrolnoise_to_random(dffive.rename({'wav123':'wav_stats','commonname':'name2'}, axis=1), samevarseqmin5five, 'noiseresgam',10000, 'barcodecontrols_noiseresgam_min5vars30perc_five')

pvals_gam_min5vars30percfive.to_pickle('./data/noiseresgamofMultipleBarcodes_min5vars30perc_five.pkl')
pvals_gam_min5vars30perc.to_pickle('./data/noiseresgamofMultipleBarcodes_min5vars30perc_ir.pkl')


###
pvals_gam_min5vars30percfive=pd.read_pickle('./data/noiseresgamofMultipleBarcodes_min5vars30perc_five.pkl')
pvals_gam_min5vars30perc=pd.read_pickle('./data/noiseresgamofMultipleBarcodes_min5vars30perc_ir.pkl')

f = plt.figure(figsize=(4,1))
sns.swarmplot(pvals_gam_min5vars30perc.pval)
plt.xlim(-0.02,1.02)
plt.xticks([0,0.5,1])
plt.xlabel('fraction of random sets with\nvariance<variance(barcode control set)')
f.savefig('./figures/Fig6/Fig6D_ir_bccontrols_vs_randomdistribution_overview.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


f = plt.figure(figsize=(4,1))
sns.swarmplot(pvals_gam_min5vars30percfive.pval)
plt.xlim(-0.02,1.02)
plt.xticks([0,0.5,1])
plt.xlabel('fraction of random sets with\nvariance<variance(barcode control set)')
f.savefig('./figures/Fig6/Fig6D_five_bccontrols_vs_randomdistribution_overview.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)




#%%

################
# Overview plots, relationship mean splicing values - noise
################
meanplusnoise=irdf[(irdf.smoothednumberofpeaks==1)&(irdf.number_reads>100)&(irdf.fraction_canonical>0.3)][['wav_stats','rnaperdna','noisestrengthlogwstd']].dropna()

gaml=pygam.LinearGAM(pygam.s(0,lam=1, n_splines=10)+pygam.l(1,lam=1)).fit(meanplusnoise[['wav_stats','rnaperdna']], meanplusnoise.noisestrengthlogwstd)

pred=gaml.predict(meanplusnoise[['wav_stats','rnaperdna']])

meanplusnoise['noisegampred']=pd.Series(pred, index=meanplusnoise.index)
meanplusnoise['noiseresgam']=meanplusnoise['noisestrengthlogwstd']-meanplusnoise['noisegampred']


f=plt.figure(figsize=(4,4))
plt.scatter(meanplusnoise.wav_stats,\
    meanplusnoise.noisestrengthlogwstd, s=10, alpha=0.2, color=sns.xkcd_rgb['medium blue'])
plt.plot(meanplusnoise.wav_stats, meanplusnoise.noisegampred, '.',  color=sns.xkcd_rgb['light green'], alpha=0.2, markersize=5)
plt.xlabel('splicing value')
plt.ylabel('splicing noise strength [log2]')
plt.ylim(-9,3)
plt.xlim(0,7.2)
plt.title('Pearson r = '+'{:.2g}'.format(irdf[['wav_stats','noisestrengthlogwstd']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(irdf[['wav_stats','noisestrengthlogwstd']].corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/Fig6/Fig6A_IR_wavstats_vs_noise_plusnoiseresgamprediction.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

f=plt.figure(figsize=(4,4))
plt.scatter(meanplusnoise.rnaperdna,\
    meanplusnoise.noisestrengthlogwstd, s=10, alpha=0.2, color=sns.xkcd_rgb['medium blue'])
plt.xlabel('RNA/DNA reads [log2]')
plt.ylabel('splicing noise strength [log2]')
plt.ylim(-10,3)
plt.xlim(-3,2)
plt.title('Pearson r = '+'{:.2g}'.format(irdf[['rnaperdna','noisestrengthlogwstd']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(irdf[['rnaperdna','noisestrengthlogwstd']].corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/FigS8/FigS8B_IR_rnaperdna_vs_noisestrength.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

f=plt.figure(figsize=(4,4))
plt.scatter(irdf.wav_stats,\
    irdf.noiseresgam, s=10, alpha=0.2, color=sns.xkcd_rgb['medium blue'])
plt.xlabel('splicing value')
plt.ylabel('noise residual')
plt.ylim(-6,6)
plt.xlim(0,7.2)
plt.title('Pearson r = '+'{:.2g}'.format(irdf[['wav_stats','noiseresgam']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(irdf[['wav_stats','noiseresgam']].corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/FigS8/FigS8B_ir_protein_vs_noiseres.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


f=plt.figure(figsize=(4,4))
plt.scatter(irdf.rnaperdna,\
    irdf.noiseresgam, s=10, alpha=0.2, color=sns.xkcd_rgb['medium blue'])
plt.xlabel('RNA/DNA reads [log2]')
plt.ylabel('noise residual')
plt.ylim(-6,6)
plt.xlim(-3,2.2)
plt.title('Pearson r = '+'{:.2g}'.format(irdf[['noiseresgam','rnaperdna']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2e}'.format(irdf[['noiseresgam','rnaperdna']].corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/FigS8/FigS8B_ir_noiseres_vs_rnaperdna.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()




meanplusnoise=fivedf[(fivedf.smoothednumberofpeaks==1)&(fivedf.number_reads>100)&\
                 (fivedf.fraction_canonical>0.3)][['wav123','rnaperdna','noisestrengthlog']].dropna()

gam=pygam.LinearGAM(pygam.s(0,lam=1, n_splines=10)+pygam.l(1,lam=1)).fit(meanplusnoise[['wav123','rnaperdna']], meanplusnoise.noisestrengthlog)

pred=gam.predict(meanplusnoise[['wav123','rnaperdna']])

meanplusnoise['noisegampred']=pd.Series(pred, index=meanplusnoise.index)
meanplusnoise['noiseresgam']=meanplusnoise['noisestrengthlog']-meanplusnoise['noisegampred']


f=plt.figure(figsize=(4,4))
plt.scatter(meanplusnoise.wav123,\
    meanplusnoise.noisestrengthlog, s=10, alpha=0.2, color=sns.xkcd_rgb['medium blue'])
plt.plot(meanplusnoise.wav123, meanplusnoise.noisegampred, '.',  color=sns.xkcd_rgb['light green'], alpha=0.3, markersize=5)
plt.xlabel('splicing value')
plt.ylabel('splicing noise strength [log2]')
plt.ylim(-10.5,6)
plt.xlim(0.5,7.2)
plt.title('Pearson r = '+'{:.2g}'.format(fivedf[['wav123','noisestrengthlog']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(fivedf[['wav123','noisestrengthlog']].corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/Fig6/Fig6B_five_wavstats_vs_noise_plusnoiseresgamprediction.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()



f=plt.figure(figsize=(4,4))
plt.scatter(fivedf.rnaperdna,\
    fivedf.noisestrengthlog, s=10, alpha=0.2, color=sns.xkcd_rgb['medium blue'])
plt.xlabel('RNA/DNA reads [log2]')
plt.ylabel('splicing noise strength [log2]')
plt.ylim(-17,7)
plt.title('Pearson r = '+'{:.2g}'.format(fivedf[['rnaperdna','noisestrengthlog']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(fivedf[['rnaperdna','noisestrengthlog']].corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/FigS8/FigS8C_five_expression_vs_noisestrength.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


f=plt.figure(figsize=(4,4))
plt.scatter(fivedf.wav123,\
    fivedf.noiseresgam, s=10, alpha=0.2, color=sns.xkcd_rgb['medium blue'])
plt.xlabel('splicing value')
plt.ylabel('noise residual')
plt.ylim(-10.10)
plt.xlim(0.5,7.2)
plt.title('Pearson r = '+'{:.2g}'.format(fivedf[['wav123','noiseresgam']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(fivedf[['wav123','noiseresgam']].corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/FigS8/FigS8C_five_protein_vs_noiseres.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()



f=plt.figure(figsize=(4,4))
plt.scatter(fivedf.rnaperdna,\
    fivedf.noiseresgam, s=10, alpha=0.2, color=sns.xkcd_rgb['medium blue'])
plt.ylabel('noise residual')
plt.xlabel('RNA/DNA reads [log2]')
plt.ylim(-10.10)
plt.xlim(-3,2.2)
plt.title('Pearson r = '+'{:.2g}'.format(fivedf[['noiseresgam','rnaperdna']].corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2e}'.format(fivedf[['noiseresgam','rnaperdna']].corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/FigS8/FigS8C_five_rnaperdna_vs_noiseres.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


################
# Noise prediction
################

'''
see file noiseprediction.py
'''


################
# Effect of splice site on noise residual
################


df = irdf[(irdf.subset=='irconstvar')&(irdf.number_reads>100)&(irdf.smoothednumberofpeaks==1)&(irdf.fraction_canonical>0.3)]
df['changes']=df.index.map(lambda x: df.first_ss[x][-2:] + '_' + df.second_ss[x][:2])

f=plt.figure(figsize=(3,3))
ax=sns.boxplot(data=df, x='changes', y='noiseresgam', \
                 order=['us_en','ve_co','GC_co','GT_U1','AT_U1'], \
                       palette=[sns.xkcd_rgb['light blue'],sns.xkcd_rgb['light green'],sns.xkcd_rgb['light green'],sns.xkcd_rgb['light red'],sns.xkcd_rgb['light red']])
ax=sns.swarmplot(data=df, x='changes', y='noiseresgam', \
                 order=['us_en','ve_co','GC_co','GT_U1','AT_U1'], \
                       palette=[sns.xkcd_rgb['medium blue'],'g','g',sns.xkcd_rgb['red'],sns.xkcd_rgb['red']])
plt.legend('')
plt.xlabel('')
plt.ylabel('noise residual')
ax.set_xticklabels(['endogenous','consensus GT', 'consensus GC', 'minor GT', 'minor AT'], rotation=45, horizontalalignment='right')
f.savefig('./figures/Fig6/Fig6G_ir_noise.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
plt.show()


'''
for i in df.changes.unique():
    print 'donor: ' + df[df.changes==i].first_ss.values[0] + ', acceptor: ' + df[df.changes==i].second_ss.values[0]
    print scipy.stats.wilcoxon(df[df.changes==i].normnoise_gam.dropna())

donor: U12_GT, acceptor: U12
WilcoxonResult(statistic=22.0, pvalue=0.32806478170083675)
donor: U12_AT, acceptor: U12
WilcoxonResult(statistic=44.0, pvalue=0.91651190786389403)
donor: endogenous, acceptor: constBP
WilcoxonResult(statistic=37.0, pvalue=0.19144640357050802)
donor: constitutive, acceptor: constitutive
WilcoxonResult(statistic=50.0, pvalue=0.87529118521097971)
donor: GC, acceptor: constitutive
WilcoxonResult(statistic=35.0, pvalue=0.027857098015306787)
donor: nonsplicable, acceptor: nonsplicable
WilcoxonResult(statistic=13.0, pvalue=0.86577237499262139)
donor: endogenous, acceptor: endogenous
WilcoxonResult(statistic=25.0, pvalue=0.47690655496758394)


for i in df.changes.unique():
    print 'donor: ' + df[df.changes==i].first_ss.values[0] + ', acceptor: ' + df[df.changes==i].second_ss.values[0]
    print scipy.stats.mannwhitneyu(df[df.changes=='us_en'].noiseresgam.dropna(), df[df.changes==i].noiseresgam.dropna())

donor: U12_GT, acceptor: U12
MannwhitneyuResult(statistic=127.0, pvalue=0.14075731625501137)
donor: U12_AT, acceptor: U12
MannwhitneyuResult(statistic=171.0, pvalue=0.45402197101306446)
donor: endogenous, acceptor: constBP
MannwhitneyuResult(statistic=122.0, pvalue=0.10920459330962085)
donor: constitutive, acceptor: constitutive
MannwhitneyuResult(statistic=129.0, pvalue=0.25504592171861662)
donor: GC, acceptor: constitutive
MannwhitneyuResult(statistic=90.0, pvalue=0.013366838199360343)
donor: nonsplicable, acceptor: nonsplicable
MannwhitneyuResult(statistic=57.0, pvalue=0.20865078000641107)
donor: endogenous, acceptor: endogenous
MannwhitneyuResult(statistic=84.5, pvalue=0.48974903124857661)


'''

####


df = fivedf[(fivedf.subset=='fiveconstvar')&(fivedf.number_reads>100)&(fivedf.smoothednumberofpeaks==1)&(fivedf.fraction_canonical>0.3)]
df['changes']=df.index.map(lambda x: df.first_ss[x][-2:] + '_' + df.second_ss[x][:2])

f=plt.figure(figsize=(4,3))
ax=sns.boxplot(data=df[['changes','noiseresgam']].dropna(), x='changes', y='noiseresgam', \
                 order=['us_en','ve_en','us_co','ve_co','GC_en','us_GC','le_en','us_un'], \
                 palette=[sns.xkcd_rgb['light blue'],sns.xkcd_rgb['light green'],sns.xkcd_rgb['light green'],sns.xkcd_rgb['light green'],sns.xkcd_rgb['light green'],sns.xkcd_rgb['light green'],sns.xkcd_rgb['light purple'],sns.xkcd_rgb['light purple']])
ax=sns.swarmplot(data=df[['changes','noiseresgam']].dropna(), x='changes', y='noiseresgam', \
                 order=['us_en','ve_en','us_co','ve_co','GC_en','us_GC','le_en','us_un'], \
                 palette=[sns.xkcd_rgb['medium blue'],'g','g','g','g','g','purple','purple'])
plt.legend('')
plt.xlabel('')
plt.ylabel('noise residual')
plt.ylim(-5,7)
ax.set_xticklabels(['both endogenous','first consensus GT', ' second consensus GT', 'both consensus GT','first consensus GC','second consensus GC', 'first nonspliceable', 'second nonspliceable'], rotation=45, horizontalalignment='right')
f.savefig('./figures/FigS8/FigS8D_five_noise.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
plt.show()

'''
for i in df.changes.unique():
    print 'first splice site: ' + df[df.changes==i].first_ss.values[0] + ', second splice site: ' + df[df.changes==i].second_ss.values[0]
    print scipy.stats.mannwhitneyu(df[df.changes=='us_en'].noiseresgam.dropna(), df[df.changes==i].noiseresgam.dropna())


first splice site: endogenous, second splice site: endogenous
MannwhitneyuResult(statistic=200.0, pvalue=0.49459938408241588)
first splice site: constitutive, second splice site: endogenous
MannwhitneyuResult(statistic=172.0, pvalue=0.1157985852265308)
first splice site: unsplicable, second splice site: endogenous
MannwhitneyuResult(statistic=113.0, pvalue=0.27161933712453745)
first splice site: endogenous, second splice site: GC
MannwhitneyuResult(statistic=151.0, pvalue=0.28644246490111286)
first splice site: constitutive, second splice site: constitutive
MannwhitneyuResult(statistic=218.0, pvalue=0.38973252526399682)
first splice site: endogenous, second splice site: constitutive
MannwhitneyuResult(statistic=151.0, pvalue=0.13968034516453481)
first splice site: endogenous, second splice site: unsplicable
MannwhitneyuResult(statistic=116.0, pvalue=0.083048538529121918)
first splice site: GC, second splice site: endogenous
MannwhitneyuResult(statistic=153.0, pvalue=0.21924825546065002)
'''

################
# Effect of splicing factor binding sites on noise (and mean splicing value)
################

df=irdf[(irdf.subset=='irSFvar')&(irdf.first_ss.isnull()==False)&(irdf.second_ss.isnull()==False)&(irdf.number_reads>100)&(irdf.smoothednumberofpeaks==1)&(irdf.fraction_canonical>0.3)]

def splitsf(ss):
    if (ss=='endogenous'):
        return 'no SF'
    else:
        return ss.split(' ')[0]
    
def splitpos(ss):
    if (ss=='endogenous'):
        return 'no SF'
    else:
        return ss.split(' ')[1]

df['firstssSF']=df.first_ss.dropna().apply(lambda x: splitsf(x))
df['firstsspos']=df.first_ss.dropna().apply(lambda x: splitpos(x))
df['secondssSF']=df.second_ss.dropna().apply(lambda x: splitsf(x))
df['secondsspos']=df.second_ss.dropna().apply(lambda x: splitpos(x))
       
def determine_location(var):
    if ('-' in df.firstsspos[var]):
        return 'exon up'
    elif ('+' in df.secondsspos[var]):
        return 'exon down'
    else:
        return 'intron'
    
df['location']=df.index.map(lambda x: determine_location(x))        
df['factor']=df.index.map(lambda x: df.firstssSF[x] if (df.secondssSF[x]=='no SF') else df.secondssSF[x]) 
df['absnormratio']=df.normratio.apply(lambda x: np.abs(x))

df['factor_loc']=df.index.map(lambda x: df.factor[x]+' ' +df.location[x])


factornames=['SRSF1','SRSF2','SRSF5','SRSF6','hnRNPA1','hnRNPG','hnRNPU']
f=plt.figure(figsize=(6,3))
ax=sns.pointplot(data=df, x='factor', y='normwav', order=factornames, \
                 hue='location', hue_order=['exon up','intron','exon down'], join=False, dodge=0.3, errwidth=3, scale=0.6)
ax.set_xticklabels(factornames, rotation=45, horizontalalignment='right')
ax.legend(title=None,ncol=3, bbox_to_anchor=[1,1.15], fontsize=14, handletextpad=0.1)
plt.ylabel('effect on splicing value')
plt.xlabel('motif inserted')
plt.axhline(y=0, color='gray')
f.savefig('./figures/FigS8/FigS8E_IRprotein_SFvar_bylocation.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


factornames=['SRSF1','SRSF2','SRSF5','SRSF6','hnRNPA1','hnRNPG','hnRNPU']
f=plt.figure(figsize=(6,3))
ax=sns.pointplot(data=df, x='factor', y='normnoise_gam', order=factornames, \
                 hue='location', hue_order=['exon up','intron','exon down'], join=False, dodge=0.3, errwidth=3, scale=0.6)
ax.set_xticklabels(factornames, rotation=45, horizontalalignment='right')
ax.legend(title=None,ncol=3, bbox_to_anchor=[1,1.15], fontsize=14, handletextpad=0.1)
plt.ylabel('effect on noise residual')
plt.xlabel('motif inserted')
plt.axhline(y=0, color='gray')
f.savefig('./figures/FigS8/FigS8E_IRnoise_SFvar_bylocation.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


factornames=['SRSF1','SRSF2','SRSF5','SRSF6','hnRNPA1','hnRNPG','hnRNPU']
f=plt.figure(figsize=(6,3))
ax=sns.pointplot(data=df, x='factor', y='normrnaperdna', order=factornames, \
                 hue='location', hue_order=['exon up','intron','exon down'], join=False, dodge=0.3, errwidth=3, scale=0.6)
ax.set_xticklabels(factornames, rotation=45, horizontalalignment='right')
ax.legend(title=None,ncol=3, bbox_to_anchor=[1,1.15], fontsize=14, handletextpad=0.1)
plt.ylabel('effect on RNA/DNA reads')
plt.xlabel('motif inserted')
plt.axhline(y=0, color='gray')
f.savefig('./figures/FigS8/FigS8E_IRexpression_SFvar_bylocation.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



f=plt.figure(figsize=(3,3))
plt.errorbar(x=df[['factor','normnoise_gam','normwav']].dropna().pivot_table(index='factor', values='normwav').values,
             y=df[['factor','normnoise_gam','normwav']].dropna().pivot_table(index='factor', values='normnoise_gam').values,
            xerr=df[['factor','normnoise_gam','normwav']].dropna().pivot_table(index='factor', values='normwav', aggfunc=scipy.stats.sem).values,
            yerr=df[['factor','normnoise_gam','normwav']].dropna().pivot_table(index='factor', values='normnoise_gam', aggfunc=scipy.stats.sem).values,
        fmt='bo ', markersize=8, color=sns.xkcd_rgb['medium blue'], alpha=0.5)
plt.axhline(y=0, color='gray', alpha=0.5)
plt.axvline(x=0, color='gray', alpha=0.5)
plt.xlim(-0.85,0.25)
f.savefig('./figures/Fig6/Fig6H_ir_sf_proteinvsnoise.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


'''
for i in df.factor.unique():
    print i
    print wilcoxon(df[['factor','normnoise_gam','normwav']].dropna()[(df.factor==i)].normnoise_gam.dropna())
    print np.mean(df[['factor','normnoise_gam','normwav']].dropna()[(df.factor==i)].normnoise_gam.dropna())
    print wilcoxon(df[['factor','normnoise_gam','normwav']].dropna()[(df.factor==i)].normwav.dropna())
    print np.mean(df[['factor','normnoise_gam','normwav']].dropna()[(df.factor==i)].normwav.dropna())
    

SRSF1
WilcoxonResult(statistic=142.0, pvalue=0.2588240374104106)
0.48214070682
WilcoxonResult(statistic=84.0, pvalue=0.01164825778943783)
-0.478875518067
SRSF2
WilcoxonResult(statistic=261.0, pvalue=0.72752278093102407)
0.184964199124
WilcoxonResult(statistic=260.0, pvalue=0.714148319549191)
0.0295809979694
SRSF5
WilcoxonResult(statistic=129.0, pvalue=0.23760036260183626)
0.574143643926
WilcoxonResult(statistic=91.0, pvalue=0.031862681910020826)
-0.382425310093
SRSF6
WilcoxonResult(statistic=159.0, pvalue=0.13059163494446679)
0.404586952314
WilcoxonResult(statistic=190.0, pvalue=0.38203416302245696)
-0.361272076011
KeEnh
WilcoxonResult(statistic=184.0, pvalue=0.13467575894495332)
0.472538851881
WilcoxonResult(statistic=191.0, pvalue=0.17224594522309822)
-0.191562826149
KeSil
WilcoxonResult(statistic=139.0, pvalue=0.22965425863748201)
0.0207133054634
WilcoxonResult(statistic=57.0, pvalue=0.0015175839700344678)
-0.648944594444
KeNeu
WilcoxonResult(statistic=174.0, pvalue=0.71856673070224608)
-0.0554844632414
WilcoxonResult(statistic=118.0, pvalue=0.088049992627778476)
-0.387840400694
hnRNPA1
WilcoxonResult(statistic=212.0, pvalue=0.90533271020781647)
0.0643462936042
WilcoxonResult(statistic=86.0, pvalue=0.0044627877475162236)
-0.58506142216
hnRNPG
WilcoxonResult(statistic=151.0, pvalue=0.034602462838034358)
0.47555047264
WilcoxonResult(statistic=233.0, pvalue=0.56213911107258718)
-0.218530578004
hnRNPU
WilcoxonResult(statistic=59.0, pvalue=0.0053549952564592319)
0.68054655084
WilcoxonResult(statistic=145.0, pvalue=0.63773289005018841)
-0.0738114243324
'''

