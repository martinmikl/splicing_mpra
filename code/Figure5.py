#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 10:20:54 2019

@author: martinm
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats
import os

sns.set_context('poster')

os.mkdir('./figures/Fig5')
os.mkdir('./figures/FigS7')

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

# plot barcode control groups

def plot_barcode_controls(df, samplename, figname=None):
    samevarseq=[]
    for var, group in df.groupby(by='varseq162'):
        if (len(group[samplename].dropna())>5):
            samevarseq.append(var)
    barcodectrlratio= []
    if (len(samevarseq)>0):
        for i in samevarseq:
            indexes=list(df[df.varseq162==i].index)
            barcodectrlratio.append(df.loc[indexes,samplename].values)
    
        bcratio=pd.DataFrame(barcodectrlratio)
        bcratio['meanval']=bcratio.median(axis=1)
        bcratio.sort_values(by='meanval', inplace=True)
        bcratio=bcratio.transpose()
        bcratio.drop('meanval', inplace=True)
        bcratio.dropna(how='all',axis=1)
        
        f=plt.figure(figsize=(len(samevarseq)/2,3))
        ax=sns.boxplot(data=bcratio, palette='Blues',fliersize=3)
        ax.set_xticklabels('')
        plt.xlabel('groups of multiple barcodes for same variant')
        plt.ylabel('splicing value')
        if figname!=None:
            f.savefig('./figures/Fig5/'+figname+'.png', \
                      dpi = 300, format='png', bbox_inches='tight', frameon=True)
        f.show()

plot_barcode_controls(irdf[(irdf.number_reads>200)&(irdf.smoothednumberofpeaks==1)],'wav_stats', figname='FigS7A_ir')
plot_barcode_controls(fivedf[(fivedf.number_reads>200)&(fivedf.smoothednumberofpeaks==1)],'wav123', figname='FigS7A_five')

plot_barcode_controls(irdf[(irdf.number_reads>200)&(irdf.smoothednumberofpeaks==1)],'noiseresgam', figname='FigS7A_irnoise')
plot_barcode_controls(fivedf[(fivedf.number_reads>200)&(fivedf.smoothednumberofpeaks==1)],'noiseresgam', figname='FigS7A_fivenoise')


# distribution of splicing values

f=plt.figure(figsize=(4,2))
plt.hist(irdf[(irdf.smoothednumberofpeaks==1)&(irdf.number_reads>100)&(irdf.fraction_canonical>0.3)].wav_stats.dropna(), bins=100, linewidth=0)
plt.ylabel('# variants')
plt.xlim(0.5,7)
plt.xlabel('splicing value [a.u.]')
f.savefig('./figures/Fig5/Fig5A_ir.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,2))
plt.hist(fivedf[(fivedf.smoothednumberofpeaks==1)&(fivedf.number_reads>100)&(fivedf.fraction_canonical>0.3)].wav123.dropna(), bins=100, linewidth=0)
plt.ylabel('# variants')
plt.xlim(0.5,7)
plt.xlabel('splicing value [a.u.]')
f.savefig('./figures/Fig5/Fig5A_five.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


# splicing ratio vs protein based splicing value

f=plt.figure(figsize=(3.5,3))
plt.scatter(irdf[(irdf.fraction_canonical>0.3)&(irdf.number_reads>100)&(irdf.smoothednumberofpeaks==1)&(irdf.rnareads_effective>100)][['wav_stats','levelratiospl','rnaperdna']].dropna().levelratiospl, \
            irdf[(irdf.fraction_canonical>0.3)&(irdf.number_reads>100)&(irdf.smoothednumberofpeaks==1)&(irdf.rnareads_effective>100)][['wav_stats','levelratiospl','rnaperdna']].dropna().wav_stats, \
                c=irdf[(irdf.fraction_canonical>0.3)&(irdf.number_reads>100)&(irdf.smoothednumberofpeaks==1)&(irdf.rnareads_effective>100)][['wav_stats','levelratiospl','rnaperdna']].dropna().rnaperdna, \
                      s=10, cmap='Blues', linewidth=0, alpha=0.3)
plt.ylabel('splicing value - protein')
plt.xlabel('splicing ratio - RNA')
plt.xlim(-12,12)
plt.colorbar(fraction=1.0/6)
plt.title('Pearson r = '+'{:.2g}'.format(irdf[(irdf.fraction_canonical>0.3)&(irdf.number_reads>100)&(irdf.smoothednumberofpeaks==1)&(irdf.rnareads_effective>100)][['wav_stats','levelratiospl']].dropna().corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(irdf[(irdf.fraction_canonical>0.3)&(irdf.number_reads>100)&(irdf.smoothednumberofpeaks==1)&(irdf.rnareads_effective>100)][['wav_stats','levelratiospl']].dropna().corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/Fig5/Fig5B_ir.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(3.5,3))
plt.scatter(fivedf[(fivedf.fraction_canonical>0.5)&(fivedf.number_reads>100)&(fivedf.smoothednumberofpeaks==1)&(fivedf.rnareads_effective>100)][['wav123','levelratiospl','rnaperdna']].dropna().levelratiospl, \
            fivedf[(fivedf.fraction_canonical>0.5)&(fivedf.number_reads>100)&(fivedf.smoothednumberofpeaks==1)&(fivedf.rnareads_effective>100)][['wav123','levelratiospl','rnaperdna']].dropna().wav123, \
                c=fivedf[(fivedf.fraction_canonical>0.5)&(fivedf.number_reads>100)&(fivedf.smoothednumberofpeaks==1)&(fivedf.rnareads_effective>100)][['wav123','levelratiospl','rnaperdna']].dropna().rnaperdna, \
                      s=10, cmap='Blues', linewidth=0, alpha=0.3)
plt.ylabel('splicing value - protein')
plt.xlabel('splicing ratio - RNA')
plt.xlim(-15,15)
plt.colorbar(fraction=1.0/6, ticks=[-1.5,-1,-0.5,0,0.5,1,1.5])
plt.title('Pearson r = '+'{:.2g}'.format(fivedf[(fivedf.fraction_canonical>0.5)&(fivedf.number_reads>100)&(fivedf.smoothednumberofpeaks==1)&(fivedf.rnareads_effective>100)][['wav123','levelratiospl']].dropna().corr().values[0][1])
          + '\nSpearman $\\rho$ = '+'{:.2g}'.format(fivedf[(fivedf.fraction_canonical>0.5)&(fivedf.number_reads>100)&(fivedf.smoothednumberofpeaks==1)&(fivedf.rnareads_effective>100)][['wav123','levelratiospl']].dropna().corr(method='spearman').values[0][1]), fontsize=14)
f.savefig('./figures/Fig5/Fig5B_five.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

# splicing ratio vs expression level - RNA

irdfs=irdf.loc[irdf.levelratiospl.dropna().index].sort_values(by='levelratiospl')

means=irdfs[(irdfs.smoothednumberofpeaks==1)&(irdfs.number_reads>100)&(irdfs.fraction_canonical>0.3)&
      (irdfs.rnareads_effective>100)].rnaperdna.dropna().rolling(window=100, center=True).mean()
stds=irdfs[(irdfs.smoothednumberofpeaks==1)&(irdfs.number_reads>100)&(irdfs.fraction_canonical>0.3)&
      (irdfs.rnareads_effective>100)].rnaperdna.dropna().rolling(window=100, center=True).std()

z=scipy.stats.norm.ppf(0.975)

meanstd=means.to_frame().join(stds, rsuffix='_std')
meanstd['meanminusci']=meanstd.index.map(lambda x: meanstd.rnaperdna[x]-z*(meanstd.rnaperdna_std[x]/100.0**(1.0/2)))
meanstd['meanplusci']=meanstd.index.map(lambda x: meanstd.rnaperdna[x]+z*(meanstd.rnaperdna_std[x]/100.0**(1.0/2)))

f=plt.figure(figsize=(3,3))
plt.scatter(irdf[(irdf.smoothednumberofpeaks==1)&(irdf.number_reads>100)&(irdf.fraction_canonical>0.3)&(irdf.rnareads_effective>100)].levelratiospl,\
    irdf[(irdf.smoothednumberofpeaks==1)&(irdf.number_reads>100)&(irdf.fraction_canonical>0.3)&(irdf.rnareads_effective>100)].rnaperdna, s=10, alpha=0.3, color=sns.xkcd_rgb['light blue'])
plt.plot(irdfs[(irdfs.smoothednumberofpeaks==1)&(irdfs.number_reads>100)&(irdfs.fraction_canonical>0.3)&(irdfs.rnareads_effective>100)].levelratiospl, \
              meanstd.rnaperdna, '.', markersize=3)
plt.fill_between(irdfs[(irdfs.smoothednumberofpeaks==1)&(irdfs.number_reads>100)&(irdfs.fraction_canonical>0.3)&(irdfs.rnareads_effective>100)].levelratiospl, \
               meanstd.meanminusci, meanstd.meanplusci,
                      color=sns.xkcd_rgb['medium blue'], alpha=0.4)
plt.xlabel('splicing ratio [log2]')
plt.ylabel('RNA/DNA reads [log2]')
plt.xlim(-12,12)
plt.ylim(-2,2)
f.savefig('./figures/Fig5/Fig5D_ir.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()



fivedfs=fivedf.loc[fivedf.levelratiospl.dropna().index].sort_values(by='levelratiospl')

means=fivedfs[(fivedfs.smoothednumberofpeaks==1)&(fivedfs.number_reads>100)&(fivedfs.fraction_canonical>0.3)&
      (fivedfs.rnareads_effective>100)].rnaperdna.dropna().rolling(window=100, center=True).mean()
stds=fivedfs[(fivedfs.smoothednumberofpeaks==1)&(fivedfs.number_reads>100)&(fivedfs.fraction_canonical>0.3)&
      (fivedfs.rnareads_effective>100)].rnaperdna.dropna().rolling(window=100, center=True).std()

z=scipy.stats.norm.ppf(0.975)

meanstd=means.to_frame().join(stds, rsuffix='_std')
meanstd['meanminusci']=meanstd.index.map(lambda x: meanstd.rnaperdna[x]-z*(meanstd.rnaperdna_std[x]/100.0**(1.0/2)))
meanstd['meanplusci']=meanstd.index.map(lambda x: meanstd.rnaperdna[x]+z*(meanstd.rnaperdna_std[x]/100.0**(1.0/2)))

f=plt.figure(figsize=(3,3))
plt.scatter(fivedf[(fivedf.smoothednumberofpeaks==1)&(fivedf.number_reads>100)&(fivedf.fraction_canonical>0.3)&(fivedf.rnareads_effective>100)].levelratiospl,\
    fivedf[(fivedf.smoothednumberofpeaks==1)&(fivedf.number_reads>100)&(fivedf.fraction_canonical>0.3)&(fivedf.rnareads_effective>100)].rnaperdna, s=10, alpha=0.3, color=sns.xkcd_rgb['light blue'])
plt.plot(fivedfs[(fivedfs.smoothednumberofpeaks==1)&(fivedfs.number_reads>100)&(fivedfs.fraction_canonical>0.3)&(fivedfs.rnareads_effective>100)].levelratiospl, \
              meanstd.rnaperdna, '.', markersize=3)
plt.fill_between(fivedfs[(fivedfs.smoothednumberofpeaks==1)&(fivedfs.number_reads>100)&(fivedfs.fraction_canonical>0.3)&(fivedfs.rnareads_effective>100)].levelratiospl, \
               meanstd.meanminusci, meanstd.meanplusci,
                      color=sns.xkcd_rgb['medium blue'], alpha=0.4)
plt.xlabel('splicing ratio [log2]')
plt.ylabel('RNA/DNA reads [log2]')
plt.xlim(-14,14)
plt.ylim(-2,2)
f.savefig('./figures/Fig5/Fig5E_five.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


### mutation of donor and acceptor sites

df = irdf[(irdf.subset=='irconstvar')&(irdf.rnareads_effective>100)&(irdf.number_reads>100)&(irdf.smoothednumberofpeaks==1)&(irdf.fraction_canonical>0.5)]
df['changes']=df.index.map(lambda x: df.first_ss[x][-2:] + '_' + df.second_ss[x][:2])

f=plt.figure(figsize=(3,3))
ax=sns.swarmplot(data=df, x='changes', y='levelratiospl', \
                 order=['us_en','ve_co','GC_co','GT_U1','AT_U1'], \
                       palette=[sns.xkcd_rgb['medium blue'],'g','g','r','r'])
plt.legend('')
plt.ylim(-10,10)
plt.xlabel('')
plt.ylabel('splicing ratio [log2]')
ax.set_xticklabels(['endogenous','consensus GT', 'consensus GC', 'minor GT', 'minor AT'], rotation=45, horizontalalignment='right')
f.savefig('./figures/Fig5/Fig5H_ir_rna.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
plt.show()

f=plt.figure(figsize=(3,3))
ax=sns.swarmplot(data=df, x='changes', y='wav_stats', \
                 order=['us_en','ve_co','GC_co','GT_U1','AT_U1'], \
                       palette=[sns.xkcd_rgb['medium blue'],'g','g','r','r'])
plt.legend('')
plt.ylim(0,8)
plt.xlabel('')
plt.ylabel('splicing value [a.u.]')
ax.set_xticklabels(['endogenous','consensus GT', 'consensus GC', 'minor GT', 'minor AT'], rotation=45, horizontalalignment='right')
f.savefig('./figures/Fig5/Fig5H_ir_protein.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
plt.show()

f=plt.figure(figsize=(3,3))
ax=sns.swarmplot(data=df, x='changes', y='rnaperdna', \
                 order=['us_en','ve_co','GC_co','GT_U1','AT_U1'], \
                       palette=[sns.xkcd_rgb['medium blue'],'g','g','r','r'])
plt.legend('')
plt.xlabel('')
plt.ylabel('RNA/DNA reads [log2]')
ax.set_xticklabels(['endogenous','consensus GT', 'consensus GC', 'minor GT', 'minor AT'], rotation=45, horizontalalignment='right')
f.savefig('./figures/Fig5/Fig5H_ir_expr.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
plt.show()



#%% intron gc content

def calculateGCintron(idx):
    if (irdf.intronstart_varseq[idx]==np.nan)|(irdf.intronend_varseqnew[idx]==np.nan):
        return np.nan
    else:
        start=int(irdf.intronstart_varseq[idx])
        end=int(irdf.intronend_varseqnew[idx])
        lngth=end-start
        return (irdf.varseq162[idx][start:end].upper().count('G')+irdf.varseq162[idx][start:end].upper().count('C'))/float(lngth)
    
irdf['intronGC']=irdf.index.map(lambda x: calculateGCintron(x))


def calculateGCexonup(idx):
    if (irdf.intronstart_varseq[idx]==np.nan)|(irdf.intronend_varseqnew[idx]==np.nan):
        return np.nan
    else:
        start=int(irdf.intronstart_varseq[idx])
        return (irdf.varseq162[idx][:start].upper().count('G')+irdf.varseq162[idx][0:start].upper().count('C'))/float(start)
    
irdf['exonupGC']=irdf.index.map(lambda x: calculateGCexonup(x))

def calculateGCexondown(idx):
    if (irdf.intronstart_varseq[idx]==np.nan)|(irdf.intronend_varseqnew[idx]==np.nan):
        return np.nan
    else:
        end=int(irdf.intronend_varseqnew[idx])
        return (irdf.varseq162[idx][end:].upper().count('G')+irdf.varseq162[idx][end:].upper().count('C'))/float(162-end)
    
irdf['exondownGC']=irdf.index.map(lambda x: calculateGCexondown(x))

def calculateGCvarseq(idx):
    return (irdf.varseq162[idx].upper().count('G')+irdf.varseq162[idx].upper().count('C'))/float(162)
    
irdf['varseqGC']=irdf.index.map(lambda x: calculateGCvarseq(x))


irdf['intronGCrel']=irdf.index.map(lambda x: irdf.intronGC[x]-(irdf.exonupGC[x]+irdf.exondownGC[x])/2)


f=plt.figure(figsize=(6,2.5))
sns.heatmap(irdf[(irdf.wav_stats>0)& \
    (irdf.wav_stats<8)&(irdf.rnareads_effective>100)&(irdf.number_reads>100)&(irdf.smoothednumberofpeaks==1)&(irdf.fraction_canonical>0.3)]\
[['levelratiospl','wav_stats','rnaperdna','intronGC','intronGCrel']].dropna().corr(),\
annot=irdf[(irdf.wav_stats>0)& \
    (irdf.wav_stats<8)&(irdf.rnareads_effective>100)&(irdf.number_reads>100)&(irdf.smoothednumberofpeaks==1)&(irdf.fraction_canonical>0.3)]\
[['levelratiospl','wav_stats','rnaperdna','intronGC','intronGCrel']].dropna().corr(), \
annot_kws={'fontsize':12})
f.savefig('./figures/Fig5/Fig5C_ir_heatmap_corr_readouts_inclGC.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

#%%

### recoding


nuc = irdf[(irdf.subset == 'irnuc')]

nucmeth = nuc[['GC' in x or 'CG' in x for x in nuc.changes]]

nucrec = nuc[[('perf' not in x) and 'control' not in x and 'GC' not in x and 'CG' not in x for x in nuc.changes]]

# effect of recoding and GC content on splicing ratios

nucrec['exonchanges']=nucrec.changes.apply(lambda x: x.split('exon')[1] if 'exon' in x else ' not recoded')
nucrec['intronchanges']=nucrec.changes.apply(lambda x: x.split('intron')[1] if 'intron' in x else ' not recoded')
nucrec['exonchanges']=nucrec.exonchanges.apply(lambda x: x.split(' and')[0] if ' and' in x else x)
nucrec['intronchanges']=nucrec.intronchanges.apply(lambda x: x.split(' and')[0] if ' and' in x else x)
nucrec['exonchanges']=nucrec.index.map(lambda x: nucrec.intronchanges[x] if ('recoded' not in nucrec.exonchanges[x]) else nucrec.exonchanges[x])

nucrec['intronchanges']=nucrec.intronchanges.apply(lambda x: 'randomly' if 'rand' in x else x)
nucrec['exonchanges']=nucrec.exonchanges.apply(lambda x: 'randomly' if 'rand' in x else x)

f=plt.figure(figsize=(6,3))
plt.axhline(y=0, color='gray',linewidth=2)
ax = sns.boxplot(data=nucrec, x = 'exonchanges', y = 'normwav', hue='intronchanges',\
                 palette=[sns.xkcd_rgb['medium blue'],'r','g',sns.xkcd_rgb['dark yellow']], fliersize=3)
plt.axhline(y=0, color='black',linewidth=2,alpha=0.7)
plt.xlabel('exon recoding')
plt.ylabel('effect on\nsplicing value [a.u.]')
plt.ylim(-6,6)
ax.set_xticklabels(['not\nrecoded','max\nGC','min\nGC','randomly'])
ax.legend().set_visible(False)
f.savefig('./figures/FigS7/FigureS5D_recoding_ir_normwav.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(6,3))
plt.axhline(y=0, color='gray',linewidth=2)
ax = sns.boxplot(data=nucrec, x = 'exonchanges', y = 'normratio', hue='intronchanges',\
                 palette=[sns.xkcd_rgb['medium blue'],'r','g',sns.xkcd_rgb['dark yellow']], fliersize=3)
plt.axhline(y=0, color='black',linewidth=2,alpha=0.7)
plt.xlabel('exon recoding')
plt.ylabel('effect on\nsplicing ratio [log2]')
ax.set_xticklabels(['not\nrecoded','max\nGC','min\nGC','randomly'])
ax.legend().set_visible(False)
f.savefig('./figures/FigS7/FigureS5D_recoding_ir_normratio.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(6,3))
plt.axhline(y=0, color='gray',linewidth=2)
ax = sns.boxplot(data=nucrec, x = 'exonchanges', y = 'normrnaperdna', hue='intronchanges',\
                 palette=[sns.xkcd_rgb['medium blue'],'r','g',sns.xkcd_rgb['dark yellow']], fliersize=3)
plt.axhline(y=0, color='black',linewidth=2,alpha=0.7)
plt.xlabel('exon recoding')
plt.ylabel('effect on\nRNA/DNA reads [log2]')
ax.set_xticklabels(['not\nrecoded','max\nGC','min\nGC','randomly'])
ax.legend().set_visible(False)
f.savefig('./figures/FigS7/FigureS5D_recoding_ir_normrnaperdna.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)



#%% Hexamers correlations

sixmerfeatures=pd.read_pickle('../rawdata/ir/sixmerfeatures.pkl')

corrs=pd.DataFrame(index=sixmerfeatures.columns)

irhexamers=sixmerfeatures.join(irdf.levelratiospl)
irhexamers=irhexamers[(irhexamers.levelratiospl>-15)&(irhexamers.levelratiospl<15)]

for i in corrs.index:
    corrs.loc[i,'r_withrna']=scipy.stats.pearsonr(irhexamers[i], irhexamers['levelratiospl'])[0]
    corrs.loc[i,'p_withrna']=scipy.stats.pearsonr(irhexamers[i], irhexamers['levelratiospl'])[1]

irhexamers=sixmerfeatures.join(irdf.wav_stats)
irhexamers=irhexamers[(irhexamers.wav_stats>0)&(irhexamers.wav_stats<8)]

for i in corrs.index:
    corrs.loc[i,'r_withprotein']=scipy.stats.pearsonr(irhexamers[i], irhexamers['wav_stats'])[0]
    corrs.loc[i,'p_withprotein']=scipy.stats.pearsonr(irhexamers[i], irhexamers['wav_stats'])[1]

irhexamers=sixmerfeatures.join(irdf.rnaperdna)
irhexamers=irhexamers[(irhexamers.rnaperdna>-15)&(irhexamers.rnaperdna<15)]

for i in corrs.index:
    corrs.loc[i,'r_withexpression']=scipy.stats.pearsonr(irhexamers[i], irhexamers['rnaperdna'])[0]
    corrs.loc[i,'p_withexpression']=scipy.stats.pearsonr(irhexamers[i], irhexamers['rnaperdna'])[1]

corrs['seqmotif']=corrs.index.map(lambda x: x.split('_')[0])
corrs['location']=corrs.index.map(lambda x: x.split('_')[1])
corrs['deltar']=corrs.index.map(lambda x: corrs.r_withprotein[x] - corrs.r_withrna[x])

corrs['factorloc']=corrs.index.map(lambda x: corrs.seqmotif[x] + '_' + corrs.location[x])

corrs['GCcont']=corrs.seqmotif.apply(lambda x: (x.count('C')+x.count('G'))/6.0)

corrs.to_pickle('./data/corrs_IR_sixmerfeatures_rnaproteinexpression.pkl')
corrs=pd.read_pickle('./data/corrs_IR_sixmerfeatures_rnaproteinexpression.pkl')

corrs.sort_values(by='deltar', inplace=True)
corrs.drop(corrs[corrs.deltar.isnull()].index, inplace=True)

f=plt.figure(figsize=(12,3))
ax=sns.barplot(data=corrs[corrs.location=='intron'],x='seqmotif',y='deltar')
ax.set_xticklabels([])
plt.xlabel('6mers')
f.savefig('./figures/FigS7/FigS7E_intronmotifs_IR_deltar_barplot.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


f=plt.figure(figsize=(6,3))
ax=sns.barplot(data=corrs[corrs.location=='intron'].iloc[:15], x='seqmotif', y='deltar', color=sns.xkcd_rgb['medium blue'])
ax.set_xticklabels(corrs[corrs.location=='intron'].iloc[:15].seqmotif.values, rotation=45, horizontalalignment='right')
f.savefig('./figures/FigS7/FigS7E_intronmotifs_IR_deltar_barplot_lowest15.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(6,3))
ax=sns.barplot(data=corrs[corrs.location=='intron'].iloc[-15:], x='seqmotif', y='deltar', color=sns.xkcd_rgb['medium blue'])
ax.set_xticklabels(corrs[corrs.location=='intron'].iloc[-15:].seqmotif.values, rotation=45, horizontalalignment='right')
f.savefig('./figures/FigS7/FigS7E_intronmotifs_IR_deltar_barplot_highest15.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(5,4))
plt.axhline(y=0, linewidth=2, color='grey', alpha=0.5)
plt.axvline(x=0, linewidth=2, color='grey', alpha=0.5)
plt.scatter(corrs[(corrs.location=='intron')&((corrs.p_withrna<0.05/4096)|(corrs.p_withprotein<0.05/4096))].r_withrna, 
                  corrs[(corrs.location=='intron')&((corrs.p_withrna<0.05/4096)|(corrs.p_withprotein<0.05/4096))].r_withprotein, 
        c=corrs[(corrs.location=='intron')&((corrs.p_withrna<0.05/4096)|(corrs.p_withprotein<0.05/4096))].GCcont, 
                s=20, cmap='Blues', linewidth=0)
#plt.scatter(corrs.r_withrna[['GTAAGT_intron','TATTTA_intron','TTAACT_intron','AATATT_intron']], 
#                  corrs.r_withprotein[['GTAAGT_intron','TATTTA_intron','TTAACT_intron','AATATT_intron']],  
#                s=20, c='Red', linewidth=0)
#plt.scatter(corrs.r_withrna[['CCCTGC_intron','CCCCGC_intron','CCTCCA_intron','CTATCC_intron']], 
#                  corrs.r_withprotein[['CCCTGC_intron','CCCCGC_intron','CCTCCA_intron','CTATCC_intron']],  
#                s=20, c='magenta', linewidth=0)
plt.xlim(-0.32, 0.32)
plt.ylim(-0.27,0.22)
plt.colorbar(ticks=[0,0.5,1])
plt.xlabel('Pearson r - RNA\n(hexamer count to splicing ratio)')
plt.ylabel('Pearson r - protein\n(hexamer count to splicing value)')
f.savefig('./figures/FigS7/FigS7G_intronmotifs_rwithrna_vs_rwithprotein.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(5,3))
ax=sns.boxplot(data=corrs, x='GCcont', y='deltar', hue='location', hue_order=['exon1','intron','exon2'], fliersize=2)
ax.set_xticklabels(['0/6','1/6','2/6','3/6','4/6','5/6','6/6'])
plt.ylim(-0.25,0.2)
plt.legend(bbox_to_anchor=(1.5,1))
f.savefig('./figures/FigS7/FigS7F_motifgccontent_vs_deltar_boxplot_forexonsandintron.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)
