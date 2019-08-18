#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 11:34:46 2019

@author: martinm
"""


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats
import numpy as np
from statsmodels.stats.multitest import fdrcorrection
import os

sns.set_context('poster')

os.mkdir('./figures/Fig2')
os.mkdir('./figures/FigS2')

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

#%% Consensus/nonspliceable donor and acceptor sites

df = casdf[(casdf.subset=='cassetteconstvar')].drop('changes', axis=1)
df['changes']=df.index.map(lambda x: df.first_ss[x][-2:] + '_' + df.second_ss[x][:2])

f=plt.figure(figsize=(3,4))
ax=sns.swarmplot(data=df, x='changes', y='levelratiospl', \
                 order=['us_en','ve_co','ve_GC','le_un'], \
                       palette=[sns.xkcd_rgb['medium blue'],'g','g','purple'])
plt.legend('')
plt.xlabel('')
plt.ylabel('splicing ratio\nexon included/skipped [log2]')
plt.ylim(-15,10)
ax.set_xticklabels(['endogenous','consensus GT', 'consensus GC', 'nonspliceable'], rotation=45, horizontalalignment='right')
f.savefig('./figures/Fig2/Figure2B_cas.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
plt.show()


df = fivedf[(fivedf.subset=='fiveconstvar')].drop('changes', axis=1)
df['changes']=df.index.map(lambda x: df.first_ss[x][-2:] + '_' + df.second_ss[x][:2])

f=plt.figure(figsize=(4,4))
ax=sns.swarmplot(data=df, x='changes', y='levelratiospl', order=['us_en','ve_en','us_co','ve_co','le_en','us_un'], \
                 palette=[sns.xkcd_rgb['medium blue'],sns.xkcd_rgb['apple green'],sns.xkcd_rgb['apple green'],'g','purple','purple'])
plt.legend('')
plt.xlabel('')
plt.ylabel('splicing ratio 2$^{nd}$/1$^{st}$ site [log2]')
plt.ylim(-15,15)
ax.set_xticklabels(['endogenous','first consensus', ' second consensus', 'both consensus', 'first nonspliceable', 'second nonspliceable'], rotation=45, horizontalalignment='right')
f.savefig('./figures/Fig2/Figure2C_five.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
plt.show()



df = threedf[(threedf.subset=='threeconstvar')].drop('changes', axis=1)
df['changes']=df.index.map(lambda x: df.first_ss[x][-2:] + '_' + df.second_ss[x][:2])

f=plt.figure(figsize=(4,4))
ax=sns.swarmplot(data=df[df.second_ss!='constBP'], x='changes', y='levelratiospl', order=['us_en','ve_en','us_co','ve_co','le_en','us_un'],\
                 palette=[sns.xkcd_rgb['medium blue'],sns.xkcd_rgb['apple green'],sns.xkcd_rgb['apple green'],'g','purple','purple'])
plt.legend('')
plt.xlabel('')
plt.ylabel('splicing ratio 2$^{nd}$/1$^{st}$ site [log2]')
plt.ylim(-15,15)
ax.set_xticklabels(['endogenous','first consensus', ' second consensus', 'both consensus','first nonspliceable', 'second nonspliceable'], rotation=45, horizontalalignment='right')
f.savefig('./figures/Fig2/Figure2D_three.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
plt.show()




### Branchpoint 
### intronretention

df = irdf[(irdf.subset=='irconstvar')]
df['changes']=df.index.map(lambda x: df.first_ss[x][-2:] + '_' + df.second_ss[x][:2])
dfresh=df.pivot(index='name2',columns='changes',values='levelratiospl')

f=plt.figure(figsize=(2,3))
ax=sns.boxplot(data=df, x='changes', y='normratio', order=['us_en','us_co'], color=sns.xkcd_rgb['medium blue'])
plt.legend('')
plt.xlabel('')
plt.ylabel('normalized ratio [log2]')
plt.ylim(-10,15)
plt.title('p= {:.2e}'.format(scipy.stats.wilcoxon(dfresh[['us_en','us_co']].dropna().us_en, \
          dfresh[['us_en','us_co']].dropna().us_co)[1]), fontsize=12)
ax.set_xticklabels(['endogenous','canonical\nbranchpoint'], rotation=45, horizontalalignment='right', \
                   multialignment='center')
f.savefig('./figures/FigS2/FigureS2A_ir.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
plt.show()

### cassette
df = casdf[(casdf.subset=='cassetteconstvar')]
df['changes']=df.index.map(lambda x: df.first_ss[x][-2:] + '_' + df.second_ss[x][:2])

dfresh=df.pivot(index='commonname',columns='changes',values='levelratiospl')

# test for effect of introduction of a consensus branchpoint sequence

scipy.stats.wilcoxon(dfresh[['us_en','BP_en']].dropna().us_en, dfresh[['us_en','BP_en']].dropna().BP_en)


f=plt.figure(figsize=(2,3))
ax=sns.boxplot(data=df, x='changes', y='normratio', \
                 order=['us_en','BP_en'], color=sns.xkcd_rgb['medium blue'])
plt.legend('')
plt.xlabel('')
plt.ylabel('normalized ratio [log2]')
plt.ylim(-10,15)
plt.title('p= {:.2e}'.format(scipy.stats.wilcoxon(dfresh[['us_en','BP_en']].dropna().us_en, \
          dfresh[['us_en','BP_en']].dropna().BP_en)[1]), fontsize=12)
ax.set_xticklabels(['endogenous','canonical\nbranchpoint'], rotation=45, horizontalalignment='right', \
                   multialignment='center')
f.savefig('./figures/FigS2/FigureS2B_cas.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
plt.show()

### Three

df = threedf[(threedf.subset=='threeconstvar')]
df['changes']=df.index.map(lambda x: df.first_ss[x][-2:] + '_' + df.second_ss[x][:2])

dfresh=df[df.second_ss!='constitutive'].pivot(index='commonname',columns='changes',values='levelratiospl')

# test for effect of introduction of a consensus branchpoint sequence

scipy.stats.wilcoxon(dfresh[['us_en','BP_en']].dropna().us_en, dfresh[['us_en','BP_en']].dropna().BP_en)
# pvalue=0.75831237426103459
scipy.stats.wilcoxon(dfresh[['us_en','us_co']].dropna().us_en, dfresh[['us_en','us_co']].dropna().us_co)
# pvalue=0.039474350914502791
scipy.stats.wilcoxon(dfresh[['us_en','BP_co']].dropna().us_en, dfresh[['us_en','BP_co']].dropna().BP_co)
# pvalue=0.19812291175396723

f=plt.figure(figsize=(4,3))
ax=sns.boxplot(data=df[(df.second_ss!='constitutive')&(df.normratio>-20)][['changes','normratio']].dropna(), x='changes', y='normratio', order=['us_en','BP_en','us_co','BP_co'], color=sns.xkcd_rgb['medium blue'])
plt.legend('')
plt.xlabel('')
plt.ylabel('normalized ratio [log2]')
plt.ylim(-10,15)
ax.set_xticklabels(['endogenous','first branchpoint\ncanonical',  'second branchpoint\ncanonical','both branchpoints\ncanonical'], rotation=45, horizontalalignment='right')
f.savefig('./figures/FigS2/FigureS2C_three.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
plt.show()

#%%

df = fivedf[(fivedf.subset == 'fiveswitchvar')]

df['firstlen']=df.firstss.apply(lambda x: int(x.split(' ')[-1][:-2]))
df['secondlen']=df.secondss.apply(lambda x: int(x.split(' ')[-1][:-2]))
df['firstss']=df.firstss.apply(lambda x: x.split(' ')[0])
df['secondss']=df.secondss.apply(lambda x: x.split(' ')[0])

df['firststart']=df.index.map(lambda x: np.min([51, 51+int(df.firstlen[x])]))
df['firstend']=df.index.map(lambda x: np.max([51, 51+int(df.firstlen[x])]))

df['secondstart']=df.index.map(lambda x: np.min([51+int(df.diff_nt[x]), 51+int(df.diff_nt[x])+int(df.firstlen[x])]))
df['secondend']=df.index.map(lambda x: np.max([51+int(df.diff_nt[x]), 51+int(df.diff_nt[x])+int(df.firstlen[x])]))

df['firstseq']=df.index.map(lambda x: str(df.varseq162[x][df.firststart[x]:df.firstend[x]]).upper())
df['secondseq']=df.index.map(lambda x: str(df.varseq162[x][df.secondstart[x]:df.secondend[x]]).upper())

df['maxentfirstwt']=df.index.map(lambda x: fivedf[(fivedf.subset=='five_filtered')&(fivedf.commonname==df.commonname[x])].maxent5first.dropna().mean())
df['maxentsecondwt']=df.index.map(lambda x: fivedf[(fivedf.subset=='five_filtered')&(fivedf.commonname==df.commonname[x])].maxent5second.dropna().mean())

df['gene_len']=df.index.map(lambda x: '_'.join([df.commonname[x],str(df.firstlen[x])]))

dfcomp=df[(df.firstlen!=-15)&(df.firstss==df.secondss)].pivot(index='gene_len',columns='firstss',values='levelratiospl')

f=plt.figure(figsize=(3,3))
sns.regplot(data=dfcomp, x='firstss',y='secondss', fit_reg=False)
plt.xlim(-14,15)
plt.ylim(-14,15)
plt.xlabel('first ss sequence')
plt.ylabel('second ss sequence')
plt.annotate("r= {:.4f}".format(scipy.stats.pearsonr(dfcomp.dropna().firstss, dfcomp.dropna().secondss)[0]) +\
", p={:.2e}".format(scipy.stats.pearsonr(dfcomp.dropna().firstss, dfcomp.dropna().secondss)[1]), xy=(-13.5,15), fontsize=14)
f.savefig('./figures/Fig2/Figure2E_five.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


fivedf['maxentdiff']=fivedf.index.map(lambda x: fivedf.maxent5second[x]-fivedf.maxent5first[x])

fivedfs=fivedf.sort_values(by='maxentdiff')


def bootstrap(data, n=1000, func=np.mean, p=0.95):
    """
    Generate `n` bootstrap samples, evaluating `func`
    at each resampling. `bootstrap` returns a function,
    which can be called to obtain confidence intervals
    of interest.
    """
    simulations = list()
    sample_size = len(data)
    for c in range(n):
        itersample = np.random.choice(data, size=sample_size, replace=True)
        simulations.append(func(itersample))
    simulations.sort()
    u_pval = (1+p)/2.
    l_pval = (1-u_pval)
    l_indx = int(np.floor(n*l_pval))
    u_indx = int(np.floor(n*u_pval))
    return func(data), simulations[l_indx],simulations[u_indx]

five_runningav=pd.DataFrame(columns=['meanratio','lowerbound','upperbound'])
for i in range(150,len(fivedfs.levelratiospl.dropna())-150):
    five_runningav.loc[i,['meanratio','lowerbound','upperbound']]=pd.Series(bootstrap(fivedfs.levelratiospl.dropna().values[i-150:i+150]),index=['meanratio', 'lowerbound', 'upperbound'])

five_runningav.astype(float)

f, ax=plt.subplots(1,1,figsize=(6,3))
ax.scatter(fivedfs.maxentdiff, fivedfs.levelratiospl, color=sns.xkcd_rgb['light blue'], alpha=0.5)
ax.scatter(fivedfs[['maxentdiff','levelratiospl']].dropna().maxentdiff.values[150:-150], \
               five_runningav.meanratio,
                      color=sns.xkcd_rgb['medium blue'], alpha=0.4)
ax.fill_between(fivedfs[['maxentdiff','levelratiospl']].dropna().maxentdiff.values[150:-150],five_runningav.lowerbound.astype(float),five_runningav.upperbound.astype(float))
plt.xlim(-11,11)
plt.ylim(-15,15)
plt.axvline(x=0, linewidth=2, color='grey')
plt.axhline(y=0, linewidth=2, color='grey')
plt.xlabel('difference in splice site strength')
plt.ylabel('splicing ratio [log2]')
f.savefig('./figures/FigS2/FigureS2D_five_ci.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


### Three

threedf['maxentdiff']=threedf.index.map(lambda x: threedf.maxent3second[x]-threedf.maxent3first[x])

threedfs=threedf.sort_values(by='maxentdiff')

three_runningav=pd.DataFrame(columns=['meanratio','lowerbound','upperbound'])
for i in range(150,len(threedfs.levelratiospl.dropna())-150):
    three_runningav.loc[i,['meanratio','lowerbound','upperbound']]=pd.Series(bootstrap(threedfs.levelratiospl.dropna().values[i-150:i+150]),index=['meanratio', 'lowerbound', 'upperbound'])

three_runningav.astype(float)

f, ax=plt.subplots(1,1,figsize=(6,3))
ax.scatter(threedfs.maxentdiff, threedfs.levelratiospl, color=sns.xkcd_rgb['light blue'], alpha=0.5)
ax.scatter(threedfs[['maxentdiff','levelratiospl']].dropna().maxentdiff.values[150:-150], \
               three_runningav.meanratio,
                      color=sns.xkcd_rgb['medium blue'], alpha=0.4)
ax.fill_between(threedfs[['maxentdiff','levelratiospl']].dropna().maxentdiff.values[150:-150],three_runningav.lowerbound.astype(float),three_runningav.upperbound.astype(float))
plt.xlim(-11,11)
plt.ylim(-15,15)
plt.axvline(x=0, linewidth=2, color='grey')
plt.axhline(y=0, linewidth=2, color='grey')
plt.xlabel('difference in splice site strength')
plt.ylabel('splicing ratio [log2]')
f.savefig('./figures/FigS2/FigureS2D_three_ci.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


#%% Splicing factor binding sites

df=casdf[(casdf.subset=='cassetteSFvar')&(casdf.first_ss.isnull()==False)&(casdf.second_ss.isnull()==False)]

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
        return 'intron up'
    elif ('+' in df.secondsspos[var]):
        return 'intron down'
    else:
        return 'exon'
    
df['location']=df.index.map(lambda x: determine_location(x))        
df['factor']=df.index.map(lambda x: df.firstssSF[x] if (df.secondssSF[x]=='no SF') else df.secondssSF[x])        


factornames=['SRSF1','SRSF6','hnRNPA1','hnRNPU']
f=plt.figure(figsize=(4,3))
ax=sns.boxplot(data=df[(df.wtratio<-1)], x='factor', y='normratio', \
                 hue='location', order=factornames)
ax.set_xticklabels(factornames, rotation=45, horizontalalignment='right')
ax.legend_.remove()
plt.ylabel('normalized ratio [log2]')
plt.xlabel('motif inserted')
plt.axhline(y=0, color='gray')
plt.ylim(-10,12)
f.savefig('./figures/Fig2/Figure2E_cas_low_efficiency_contexts.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

factornames=['SRSF1','SRSF6','hnRNPA1','hnRNPU']
f=plt.figure(figsize=(4,3))
ax=sns.boxplot(data=df[(df.wtratio>1)], x='factor', y='normratio', \
                 hue='location', order=factornames)
ax.set_xticklabels(factornames, rotation=45, horizontalalignment='right')
ax.legend_.remove()
plt.ylabel('normalized ratio [log2]')
plt.xlabel('motif inserted')
plt.axhline(y=0, color='gray')
plt.ylim(-15,5)
f.savefig('./figures/Fig2/Figure2E_cas_high_efficiency_contexts.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


###

#%% Overview of all motifs inserted


allfactornames=['SRSF1','SRSF2','SRSF5','SRSF6','hnRNPA1','hnRNPG','hnRNPU','KeEnh','KeEnh1','KeEnh2','KeEnh3','KeSil1','KeSil2','KeSil3','KeNeu1','KeNeu2','KeNeu3','RosenbergGENsil',
'RosenbergE5enh','RosenbergE5sil','RosenbergI5enh','RosenbergI5sil','RosenbergE3enh','RosenbergE3sil','RosenbergI3enh','RosenbergI3sil']

litmotifs=pd.DataFrame(index=allfactornames)
litmotifsp=pd.DataFrame(index=allfactornames)

###############
df=irdf[((irdf.subset=='irSFRosenberg')|(irdf.subset=='irSFvar'))&(irdf.first_ss.isnull()==False)&(irdf.second_ss.isnull()==False)]

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
    elif ('+' in df.firstsspos[var]):
        return 'intron start'
    elif ('-' in df.secondsspos[var]):
        return 'intron end'
    
df['location']=df.index.map(lambda x: determine_location(x))        
df['factor']=df.index.map(lambda x: df.firstssSF[x] if (df.secondssSF[x]=='no SF') else df.secondssSF[x])        

dfke=df.pivot_table(index='factor', columns='location', values='normratio')
dfke=dfke[['exon up','intron start','intron end','exon down']]

dfkep=pd.DataFrame(index=dfke.index,columns=['exon up','intron start','intron end','exon down'])
for f in df.factor.unique():
    for l in df.location.unique():
        dfkep.loc[f,l]=scipy.stats.wilcoxon(df[(df.factor==f)&(df.location==l)].normratio.dropna())[1]

pcorrsdonor=fdrcorrection(dfkep[['exon up','intron start']].stack(), alpha=0.05)[1]
pcorrsacceptor=fdrcorrection(dfkep[['intron end','exon down']].stack(), alpha=0.05)[1]

corrpvaldonor=dict(zip(dfkep[['exon up','intron start']].stack().values,pcorrsdonor))
corrpvalacceptor=dict(zip(dfkep[['intron end','exon down']].stack().values,pcorrsacceptor))

dfkepcorrdonor=dfkep[['exon up','intron start']].applymap(lambda x: corrpvaldonor[x] if np.isnan(x)==False else np.nan)
dfkepcorracceptor=dfkep[['intron end','exon down']].applymap(lambda x: corrpvalacceptor[x] if np.isnan(x)==False else np.nan)

litmotifs=litmotifs.join(dfke)
litmotifsp=litmotifsp.join(dfkepcorrdonor.join(dfkepcorracceptor))




df=casdf[((casdf.subset=='cassetteSFRosenberg')|(casdf.subset=='cassetteSFvar'))&(casdf.first_ss.isnull()==False)&(casdf.second_ss.isnull()==False)]

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
        return 'intron up'
    elif ('+' in df.secondsspos[var]):
        return 'intron down'
    elif ('-' in df.secondsspos[var]):
        return 'exon end'
    else:
        return 'exon start'
    
df['location']=df.index.map(lambda x: determine_location(x))        
df['factor']=df.index.map(lambda x: df.firstssSF[x] if (df.secondssSF[x]=='no SF') else df.secondssSF[x])        


dfke=df.pivot_table(index='factor', columns='location', values='normratio')
dfke=dfke[['intron up','exon start','exon end','intron down']]

dfkep=pd.DataFrame(index=dfke.index,columns=['intron up','exon start','exon end','intron down'])
for f in df.factor.unique():
    for l in df.location.unique():
        dfkep.loc[f,l]=scipy.stats.wilcoxon(df[(df.factor==f)&(df.location==l)].normratio.dropna())[1]

pcorrsdonor=fdrcorrection(dfkep[['exon end','intron down']].stack(), alpha=0.05)[1]
pcorrsacceptor=fdrcorrection(dfkep[['intron up','exon start']].stack(), alpha=0.05)[1]

corrpvaldonor=dict(zip(dfkep[['exon end','intron down']].stack().values,pcorrsdonor))
corrpvalacceptor=dict(zip(dfkep[['intron up','exon start']].stack().values,pcorrsacceptor))

dfkepcorrdonor=dfkep[['exon end','intron down']].applymap(lambda x: corrpvaldonor[x] if np.isnan(x)==False else np.nan)
dfkepcorracceptor=dfkep[['intron up','exon start']].applymap(lambda x: corrpvalacceptor[x] if np.isnan(x)==False else np.nan)

litmotifs=litmotifs.join(dfke)
litmotifsp=litmotifsp.join(dfkepcorrdonor.join(dfkepcorracceptor))



df=fivedf[((fivedf.subset=='fiveSFRosenberg')|(fivedf.subset=='fiveSFvar'))&(fivedf.first_ss.isnull()==False)&(fivedf.second_ss.isnull()==False)&(fivedf.second_ss != 'YuE1')&(fivedf.second_ss != 'YuI1')&(fivedf.second_ss != 'YuI2')&(fivedf.first_ss != 'YuE1')&(fivedf.first_ss != 'YuI1')&(fivedf.first_ss != 'YuI2')]

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
        return 'exon'
    elif ('+' in df.secondsspos[var]):
        return 'intron'
    elif ('-' in df.secondsspos[var]):
        return 'alt end'
    else:
        return 'alt start'
    
df['location']=df.index.map(lambda x: determine_location(x))        
df['factor']=df.index.map(lambda x: df.firstssSF[x] if (df.secondssSF[x]=='no SF') else df.secondssSF[x])        


dfke=df[df.firstsspos!=df.secondsspos].pivot_table(index='factor', columns='location', values='normratio')
dfke=dfke[['exon','alt start','alt end','intron']]

dfkep=pd.DataFrame(index=dfke.index,columns=['exon','alt start','alt end','intron'])
for f in df.factor.unique():
    for l in df.location.unique():
        dfkep.loc[f,l]=scipy.stats.wilcoxon(df[(df.factor==f)&(df.location==l)].normratio.dropna())[1]

pcorrsfirst=fdrcorrection(dfkep[['exon','alt start']].stack(), alpha=0.05)[1]
pcorrssecond=fdrcorrection(dfkep[['alt end','intron']].stack(), alpha=0.05)[1]

corrpvalfirst=dict(zip(dfkep[['exon','alt start']].stack().values,pcorrsfirst))
corrpvalsecond=dict(zip(dfkep[['alt end','intron']].stack().values,pcorrssecond))

dfkepcorrfirst=dfkep[['exon','alt start']].applymap(lambda x: corrpvalfirst[x] if np.isnan(x)==False else np.nan)
dfkepcorrsecond=dfkep[['alt end','intron']].applymap(lambda x: corrpvalsecond[x] if np.isnan(x)==False else np.nan)

litmotifs=litmotifs.join(dfke)
litmotifsp=litmotifsp.join(dfkepcorrfirst.join(dfkepcorrsecond))



df=threedf[((threedf.subset=='threeSFRosenberg')|(threedf.subset=='threeSFvar'))&(threedf.first_ss.isnull()==False)&(threedf.second_ss.isnull()==False)]

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
        return 'intron'
    elif ('+' in df.secondsspos[var]):
        return 'exon'
    elif ('-' in df.secondsspos[var]):
        return 'alt end'
    else:
        return 'alt start'
    
df['location']=df.index.map(lambda x: determine_location(x))        
df['factor']=df.index.map(lambda x: df.firstssSF[x] if (df.secondssSF[x]=='no SF') else df.secondssSF[x])        


dfke=df[df.firstsspos!=df.secondsspos].pivot_table(index='factor', columns='location', values='normratio')
dfke=dfke[['intron','alt start','exon']]

dfkep=pd.DataFrame(index=dfke.index,columns=['intron','alt start','exon'])
for f in df.factor.unique():
    for l in df.location.unique():
        dfkep.loc[f,l]=scipy.stats.wilcoxon(df[(df.factor==f)&(df.location==l)].normratio.dropna())[1]

pcorrsfirst=fdrcorrection(dfkep[['intron','alt start']].stack(), alpha=0.05)[1]
pcorrssecond=fdrcorrection(dfkep[['alt start','exon']].stack(), alpha=0.05)[1]

corrpvalfirst=dict(zip(dfkep[['intron','alt start']].stack().values,pcorrsfirst))
corrpvalsecond=dict(zip(dfkep[['alt start','exon']].stack().values,pcorrssecond))

dfkepcorrfirst=dfkep[['intron','alt start']].applymap(lambda x: corrpvalfirst[x] if np.isnan(x)==False else np.nan)
dfkepcorrsecond=dfkep[['alt start','exon']].applymap(lambda x: corrpvalsecond[x] if np.isnan(x)==False else np.nan)

litmotifs=litmotifs.join(dfke, rsuffix='_three')
litmotifsp=litmotifsp.join(dfkepcorrfirst.join(dfkepcorrsecond,rsuffix='_threesecond'), rsuffix='_three')

allmotifs=['SRSF1','SRSF2','SRSF5','SRSF6','hnRNPA1','hnRNPG','hnRNPU','GACGTC','AGAAGA','GAAGAT','GCAAGA','CCAGCA','CTTTTA','CTAGTA','AAAGAG','AAACTT','AACCTT','GTGGGG','CACCGC','GGTGGG','TTGTTC','CGAACC','CGAAGA','GGGGGG','TCTAAC','CCAAGC']
litmotifs.index=allmotifs
litmotifs.to_pickle('./data/sfkerosforheatmap.pkl')
litmotifsp.index=allmotifs
litmotifsp.to_pickle('./data/sfkerosforheatmap_wilcoxon.pkl')




litmotifsreordered=pd.DataFrame([litmotifs['exon up'], litmotifs['intron start'], \
                                litmotifs['exon end'], litmotifs['intron down'], \
                                litmotifs['exon']*-1, litmotifs['alt start']*-1, litmotifs['alt end'], litmotifs['intron']])

litmotifspreordered=pd.DataFrame([litmotifsp['exon up'], litmotifsp['intron start'], \
                                litmotifsp['exon end'], litmotifsp['intron down'], \
                                litmotifsp['exon'], litmotifsp['alt start'], litmotifsp['alt start'], litmotifsp['intron']])

f=plt.figure(figsize=(12,9))
ax=sns.heatmap(litmotifsreordered.T, annot=litmotifspreordered.T, annot_kws={'fontsize':12})
ax.set_xticklabels([])
plt.axvline(x=2, color='black')
plt.axvline(x=4, color='black')
plt.axvline(x=6, color='grey')
f.savefig('./figures/FigS2/FigureS2E_OverviewHeatmap_donors.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

litmotifsreordered=pd.DataFrame([litmotifs['intron end'], litmotifs['exon down'], \
                                litmotifs['intron up'], litmotifs['exon start'], \
                                litmotifs['intron_three']*-1, litmotifs['alt start_three']*-1, litmotifs['alt start_three'], litmotifs['exon_three']])

litmotifspreordered=pd.DataFrame([litmotifsp['intron end'], litmotifsp['exon down'], \
                            litmotifsp['intron up'], litmotifsp['exon start'], \
                            litmotifsp['intron_three'], litmotifsp['alt start_three'], litmotifsp['alt start_three'], litmotifsp['exon_three']])


f=plt.figure(figsize=(12,9))
ax=sns.heatmap(litmotifsreordered.T, annot=litmotifspreordered.T, annot_kws={'fontsize':12})
ax.set_xticklabels([])
plt.axvline(x=2, color='black')
plt.axvline(x=4, color='black')
plt.axvline(x=6, color='grey')
f.savefig('./figures/FigS2/FigureS2E_OverviewHeatmap_acceptors.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)







f=plt.figure(figsize=(3,3))
plt.scatter(litmotifs['intron start'], litmotifs['intron end'])
plt.xlabel("motif effect in the intron\nclose to the donor")
plt.ylabel("motif effect in the intron\nclose to the acceptor")
plt.xticks([-2,-1,0,1,2,3,4])
plt.yticks([-2,-1,0,1,2,3,4])
plt.xlim(-2,4)
plt.ylim(-2,4)
plt.plot([-2,4],[-2,4], color='gray', linewidth=2, alpha=0.5)
plt.annotate("r = {:.2f}".format(scipy.stats.pearsonr(litmotifs['intron start'].dropna(), litmotifs['intron end'].dropna())[0]) +\
"\np = {:.1e}".format(scipy.stats.pearsonr(litmotifs['intron start'].dropna(), litmotifs['intron end'].dropna())[1]), xy=(-1.5,3), fontsize=14)
f.savefig('./figures/FigS2/FigureS2F_motifeffects_ir.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



f=plt.figure(figsize=(3,3))
plt.scatter(litmotifs['exon start'], litmotifs['exon end'])
plt.xlabel("motif effect in the exon\nclose to the acceptor")
plt.ylabel("motif effect in the exon\nclose to the donor")
plt.xticks([-6,-4,-2,0,2,4])
plt.yticks([-6,-4,-2,0,2,4])
plt.xlim(-6,4)
plt.ylim(-6,4)
plt.plot([-6,4],[-6,4], color='gray', linewidth=2, alpha=0.5)
plt.annotate("r = {:.2f}".format(scipy.stats.pearsonr(litmotifs['exon start'].dropna(), litmotifs['exon end'].dropna())[0]) +\
"\np = {:.1e}".format(scipy.stats.pearsonr(litmotifs['exon start'].dropna(), litmotifs['exon end'].dropna())[1]), xy=(-5,2), fontsize=14)
f.savefig('./figures/FigS2/FigureS2F_motifeffects_cas.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(3,3))
plt.scatter(litmotifs['alt start'], litmotifs['alt end'])
plt.xlabel("motif effect in the alt. exon\nclose to the first donor")
plt.ylabel("motif effect in the alt. exon\nclose to the second donor")
plt.xticks([-4,-2,0,2,4,6])
plt.yticks([-4,-2,0,2,4,6])
plt.xlim(-4,6)
plt.ylim(-4,6)
plt.plot([-4,6],[-4,6], color='gray', linewidth=2, alpha=0.5)
plt.annotate("r = {:.2f}".format(scipy.stats.pearsonr(litmotifs['alt start'].dropna(), litmotifs['alt end'].dropna())[0]) +\
"\np = {:.1e}".format(scipy.stats.pearsonr(litmotifs['alt start'].dropna(), litmotifs['alt end'].dropna())[1]), xy=(1,-3), fontsize=14)
f.savefig('./figures/FigS2/FigureS2F_motifeffects_five.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


#%% Combinatorial effects between splicing factor binding sites




df = casdf[(casdf.subset == 'cassetteSFcombvar')]

df['firstssSF']=df.first_SF.apply(lambda x: x.split(' at position ')[0])
df['firstsspos']=df.first_SF.apply(lambda x: x.split(' at position ')[1])
df['secondssSF']=df.second_SF.apply(lambda x: x.split(' at position ')[0])
df['secondsspos']=df.second_SF.apply(lambda x: x.split(' at position ')[1])

def sfloc1(ss):
    if (int(df.firstsspos[ss])<117-df.exonlength[ss]):
        return 'intron up'
    elif (int(df.firstsspos[ss])>117):
        return 'intron down'
    else:
        return 'exon'

def sfloc2(ss):
    if (int(df.secondsspos[ss])<117-df.exonlength[ss]):
        return 'intron up'
    elif (int(df.secondsspos[ss])>117):
        return 'intron down'
    else:
        return 'exon'

     
df['firstsfloc']=df.index.map(lambda x: sfloc1(x))
df['secondsfloc']=df.index.map(lambda x: sfloc2(x))


# heatmap including location - all

sfs= df.firstssSF.unique()

sflocmat = pd.DataFrame()

for xv in sfs:
    for xloc in ['intron up','exon','intron down']:
        for yv in sfs:
            for yloc in ['intron up','exon','intron down']:
                locav = np.average(df[(df.firstssSF==xv) & \
                (df.firstsfloc==xloc) & (df.secondssSF==yv) & \
                (df.secondsfloc==yloc)].normratio.dropna())
                if (np.isnan(locav)==False):
                    sflocmat.loc[str(xv + ' ' + xloc),str(yv + ' ' + yloc)] =  locav

sflocmat = sflocmat[sflocmat.columns].astype(float)
                                        
sflocmatp = pd.DataFrame("",index=sflocmat.index, columns=sflocmat.columns)

for xv in sfs:
    for xloc in ['intron up','exon','intron down']:
        for yv in sfs:
            for yloc in ['intron up','exon','intron down']:
                if len(df[(df.firstssSF==xv) & \
                (df.firstsfloc==xloc) & (df.secondssSF==yv) & \
                (df.secondsfloc==yloc)].normratio.dropna())>0:
                    locav = scipy.stats.wilcoxon(df[(df.firstssSF==xv) & \
                    (df.firstsfloc==xloc) & (df.secondssSF==yv) & \
                    (df.secondsfloc==yloc)].normratio.dropna())[1]
                    if (np.isnan(locav)==False)&(len(df[(df.firstssSF==xv) & \
                            (df.firstsfloc==xloc) & (df.secondssSF==yv) & \
                            (df.secondsfloc==yloc)].normratio.dropna())>2)&\
                            (str(xv + ' ' + xloc) in sflocmat.index)&(str(yv + ' ' + yloc) in sflocmat.columns):
                        if (locav<0.05):
                            sflocmatp.loc[str(xv + ' ' + xloc),str(yv + ' ' + yloc)] ="*"
                        else:
                            sflocmatp.loc[str(xv + ' ' + xloc),str(yv + ' ' + yloc)] =""

f=plt.figure(figsize=(6,3))
ax = sns.heatmap(data=sflocmat, center=0, linewidths=2, annot=sflocmatp, fmt='s', annot_kws={'fontsize':14, 'color':'white'})
plt.xlabel('second splicing factor')
plt.ylabel('first splicing factor')
ax.set_xticklabels(sflocmat.columns, rotation=45, horizontalalignment='right')
f.savefig('./figures/FigS2/FigureS2G_combinatorial_cas.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

#%% Recoding and methylation

### cassette

nuc = casdf[(casdf.subset == 'cassettenuc')]


# effect of recoding and GC content of the cassette exon on splicing ratios

nucrec = nuc[[('perf' not in x) and 'periodicity' not in x and 'control' not in x and 'GC' not in x and 'CG' not in x for x in nuc.changes]]

nucrec['changes']=nucrec.changes.apply(lambda x: 'randomly' if 'rand' in x else x)

f=plt.figure(figsize=(3,3))
ax=sns.boxplot(data=nucrec,\
                 x='changes', y='normratio', fliersize=3, \
                 color=sns.xkcd_rgb['medium blue'])
plt.axhline(y=0, color='gray',linewidth=2)
plt.xlabel('exon recoding')
plt.ylabel('normalized ratio [log2]')
ax.set_xticklabels(['not recoded','max GC','min GC','randomly'], rotation=45, horizontalalignment='right')
f.savefig('./figures/Fig2/Figure2F_recoding_cas.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

# effect of insertion of CG (and GC as a control) on splicing ratios of cassette exons

nucmeth = nuc[['GC' in x or 'CG' in x for x in nuc.changes]]

nucmeth['recod']=nucmeth.changes.apply(lambda x: x.split(', ')[0].lower() if (',' in x) else 'not recoded')
nucmeth['cpg']=nucmeth.changes.apply(lambda x: x.split(', ')[1].lower() if (',' in x) else x.lower())
    
nucmeth['cg'] = nucmeth.cpg.apply(lambda x: x.split(' ')[0])
nucmeth['modification'] = nucmeth.cpg.apply(lambda x: ' '.join(np.array(x.split(' '))[2:]))
nucmeth['normratio_cg'] = nucmeth.index.map(lambda x: nucmeth.levelratiospl[x] - \
       nucmeth[(nucmeth.modification==nucmeth.modification[x])&(nucmeth.cg=='gc')&\
    (nucmeth.commonname==nucmeth.commonname[x]) &(nucmeth.recod==nucmeth.recod[x])].levelratiospl.values[0] \
    if (nucmeth.cg[x]=='cg') else np.nan)


f=sns.factorplot(data=nucmeth[(nucmeth.cg=='cg')],\
                 x='modification', y='normratio_cg', order=['upstream intron every 6nt','last 20 nt of the exon every 6nt','exon every 9nt','exon every 6nt','upstream intron and in exon every 6nt'],\
                 kind='point', size=3, aspect=2,
                 join=False, dodge=0.3, errwidth=3, scale=0.6)
f.map(plt.axhline, y=0, color='gray',linewidth=2)
plt.xlabel('CG/GC (periodicity)')
plt.ylabel('effect of CG [log2]\nnormalized to GC')
f.set_xticklabels(['intron up (6 nt)','exon (6 nt, only 3 times)','exon (9 nt)','exon (6 nt)','intron+exon (6 nt)'], rotation=45, horizontalalignment='right')
f.savefig('./figures/Fig2/Figure2G_methylation_cas.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


### intronretention

nuc = irdf[(irdf.subset == 'irnuc')]

nucmeth = nuc[['GC' in x or 'CG' in x for x in nuc.changes]]


# effect of insertion of CG (and GC as a control) on splicing of introns

nucmeth['recod']=nucmeth.changes.apply(lambda x: x.split(', ')[0].lower() if (',' in x) else 'not recoded')
nucmeth['cpg']=nucmeth.changes.apply(lambda x: x.split(', ')[1].lower() if (',' in x) else x.lower())
    
nucmeth['cg'] = nucmeth.cpg.apply(lambda x: x.split(' ')[0])
nucmeth['modification'] = nucmeth.cpg.apply(lambda x: ' '.join(np.array(x.split(' '))[2:]))
nucmeth['normratio_cg'] = nucmeth.index.map(lambda x: nucmeth.levelratiospl[x] - \
       nucmeth[(nucmeth.modification==nucmeth.modification[x])&(nucmeth.cg=='gc')&\
    (nucmeth.name2==nucmeth.name2[x]) &(nucmeth.recod==nucmeth.recod[x])].levelratiospl.values[0] \
    if ((nucmeth.cg[x]=='cg')&(len(nucmeth[(nucmeth.modification==nucmeth.modification[x])&(nucmeth.cg=='gc')&\
    (nucmeth.name2==nucmeth.name2[x]) &(nucmeth.recod==nucmeth.recod[x])].levelratiospl.values)>0)) else np.nan)

f=sns.factorplot(data=nucmeth[(nucmeth.cg=='cg')],\
                 x='modification', y='normratio_cg', order=['intron every 9nt','intron every 6nt','exon every 9nt','exon every 6nt','exon and intron every 6nt'],\
                 kind='point', size=3, aspect=1.6,
                 join=False, dodge=0.3, errwidth=3, scale=0.6)
f.map(plt.axhline, y=0, color='gray',linewidth=2)
plt.xlabel('CG/GC (periodicity)')
plt.ylabel('effect of CG [log2]\nnormalized to GC')
f.set_xticklabels(['intron (9 nt)','intron (6 nt)','exon (9 nt)','exon (6 nt)','intron+exon (6 nt)'], rotation=45, horizontalalignment='right')
f.savefig('./figures/Fig2/Figure2H_methylation_ir.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

