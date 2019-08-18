#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 17:31:33 2019

@author: martinm
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats
import os

sns.set_context('poster')

os.mkdir('./figures/Fig3')
os.mkdir('./figures/FigS3')

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

comb=fivedf[(fivedf.subset == 'five_comb_with_constitutive')]

comb['exonup']=comb.changes.apply(lambda x: x.split(' ')[1] if (x.split(' ')[0]=='exon') else 'alt')
comb['introndown']=comb.changes.apply(lambda x: str(x.split(' ')[1]) if (x.split(' ')[0]=='intron') else 'alt')
comb['changed']=comb.changes.apply(lambda x: x.split(' ')[0])
comb['constname']=comb.changes.apply(lambda x: str(x.split(' ')[1]))

comb['ratioadjacent']=comb.index.map(lambda x: comb.levelratiospl[x] if comb.changed[x]=='intron' else comb.levelratiospl[x]*-1)
comb['normratioadjacent']=comb.index.map(lambda x: comb.normratio[x] if comb.changed[x]=='intron' else comb.normratio[x]*-1)

comb['nnumber']=comb.index.map(lambda x: len(comb[(comb.introndown==comb.introndown[x])]))

endog=pd.DataFrame(comb[comb.wtratio>-25].wtratio.unique().T,index=comb[comb.wtratio>-25].commonname.unique(), columns=['levelratiospl'])
endog['constname']='none'
######################




genesorder=list(comb[(comb.nnumber>5)].introndown.unique())
labelsorder=np.copy(genesorder)
labelsorder[0]='endogenous'


f=plt.figure(figsize=(4.5,3))
ax=sns.boxplot(data=comb[(comb.nnumber>5)&(comb.changed=='exon')][['constname','levelratiospl']].append(endog).dropna(),x='constname', y='levelratiospl',\
            order=['none']+genesorder[1:], palette=[sns.xkcd_rgb['light grey']]+[sns.xkcd_rgb['light blue']]*8)
ax=sns.swarmplot(data=comb[(comb.nnumber>5)&(comb.changed=='exon')][['constname','levelratiospl']].append(endog).dropna(),x='constname', y='levelratiospl',\
            split=True,order=['none']+genesorder[1:], palette=[sns.xkcd_rgb['dark grey']]+[sns.xkcd_rgb['medium blue']]*8,s=4)
ax.set_xticklabels(['none']+genesorder[1:],rotation=45, horizontalalignment='right')    
plt.ylabel('splicing ratio 2$^{nd}$/1$^{st}$ site [log2]')
plt.xlabel('constitutive splice site')
plt.ylim(-15,15)
handles, labels=ax.get_legend_handles_labels()
f.savefig('./figures/Fig3/Figure3A_five_exon.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4.5,3))
ax=sns.boxplot(data=comb[(comb.nnumber>5)&(comb.changed=='intron')][['constname','levelratiospl']].append(endog).dropna(),x='constname', y='levelratiospl',\
            order=['none']+genesorder[1:], palette=[sns.xkcd_rgb['light grey']]+[sns.xkcd_rgb['light green']]*8)
ax=sns.swarmplot(data=comb[(comb.nnumber>5)&(comb.changed=='intron')][['constname','levelratiospl']].append(endog).dropna(),x='constname', y='levelratiospl',\
            split=True,order=['none']+genesorder[1:], palette=[sns.xkcd_rgb['dark grey']]+[sns.xkcd_rgb['dark green']]*8,s=4)
ax.set_xticklabels(['none']+genesorder[1:],rotation=45, horizontalalignment='right')    
plt.ylabel('splicing ratio 2$^{nd}$/1$^{st}$ site [log2]')
plt.xlabel('constitutive splice site')
plt.ylim(-15,15)
f.savefig('./figures/Fig3/Figure3B_five_intron.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


comb['name_comb']=comb.index.map(lambda x: str(comb.commonname[x]) + '_' + str(comb.constname[x]))
'''
comb['levelratiospl_adjacent']=comb.index.map(lambda x: comb.levelratiospl[x]*-1 if comb.changed[x]=='exon' else comb.levelratiospl[x])


scipy.stats.wilcoxon(comb.pivot_table(index='name_comb',columns='changed',values='levelratiospl_adjacent').dropna().exon,
                     comb.pivot_table(index='name_comb',columns='changed',values='levelratiospl_adjacent').dropna().intron)


WilcoxonResult(statistic=7051.0, pvalue=0.0041008158157462031)
mean
exon      3.079961
intron    0.500206
std
exon      7.282378
intron    6.979200
'''


rr,pp=scipy.stats.pearsonr(comb.pivot(index='name_comb',columns='changed',values='levelratiospl').dropna()['exon'],comb.pivot(index='name_comb',columns='changed',values='levelratiospl').dropna()['intron'])
#(0.42649200101783341, 6.9067908306314377e-10)

f=plt.figure(figsize=(3,3))
sns.regplot(data=comb.pivot(index='name_comb',columns='changed',values='levelratiospl'), x='exon',y='intron')
plt.xlabel('splicing ratio 2$^{nd}$/1$^{st}$ site [log2]\nupstream exon replaced')
plt.ylabel('splicing ratio 2$^{nd}$/1$^{st}$ site [log2]\ndownstream intron replaced')
plt.ylim(-15,15)
plt.xlim(-15,15)
plt.annotate("r={:.2f}".format(rr) + \
             "\np={:.0e}".format(pp),xy=(0,-14.5), fontsize=14)
f.savefig('./figures/Fig3/Figure3C_five_intronvsexoneffect.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


### three


comb=threedf[(threedf.subset == 'three_comb_with_constitutive')]

comb['exondown']=comb.changes.apply(lambda x: x.split(' ')[1] if (x.split(' ')[0]=='exon') else 'alt')
comb['intronup']=comb.changes.apply(lambda x: x.split(' ')[1] if (x.split(' ')[0]=='intron') else 'alt')
comb['changed']=comb.changes.apply(lambda x: x.split(' ')[0])
comb['constname']=comb.changes.apply(lambda x: x.split(' ')[1])

comb['ratioadjacent']=comb.index.map(lambda x: comb.levelratiospl[x] if comb.changed[x]=='exon' else comb.levelratiospl[x]*-1)
comb['normratioadjacent']=comb.index.map(lambda x: comb.normratio[x] if comb.changed[x]=='exon' else comb.normratio[x]*-1)


comb['nnumber']=comb.index.map(lambda x: len(comb[(comb.intronup==comb.intronup[x])]))

endog=pd.DataFrame(comb[comb.wtratio>-25].wtratio.unique().T,index=comb[comb.wtratio>-25].commonname.unique(), columns=['levelratiospl'])
endog['constname']='none'

######################

genesorder=list(comb[(comb.nnumber>5)].intronup.unique())
genesorder.remove('alt')


f=plt.figure(figsize=(4.5,3))
ax=sns.boxplot(data=comb[(comb.nnumber>5)&(comb.changed=='exon')][['constname','levelratiospl']].append(endog).dropna(),x='constname', y='levelratiospl',\
            order=['none']+genesorder, palette=[sns.xkcd_rgb['light grey']]+[sns.xkcd_rgb['light blue']]*8)
ax=sns.swarmplot(data=comb[(comb.nnumber>5)&(comb.changed=='exon')][['constname','levelratiospl']].append(endog).dropna(),x='constname', y='levelratiospl',\
            split=True,order=['none']+genesorder, palette=[sns.xkcd_rgb['dark grey']]+[sns.xkcd_rgb['medium blue']]*8,s=4)
ax.set_xticklabels(['none']+genesorder,rotation=45, horizontalalignment='right')    
plt.ylabel('splicing ratio 2$^{nd}$/1$^{st}$ site [log2]')
plt.xlabel('constitutive splice site')
plt.ylim(-17,15)
f.savefig('./figures/FigS3/FigureS3A_three_exon.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4.5,3))
ax=sns.boxplot(data=comb[(comb.nnumber>5)&(comb.changed=='intron')][['constname','levelratiospl']].append(endog).dropna(),x='constname', y='levelratiospl',\
            order=['none']+genesorder, palette=[sns.xkcd_rgb['light grey']]+[sns.xkcd_rgb['light green']]*8)
ax=sns.swarmplot(data=comb[(comb.nnumber>5)&(comb.changed=='intron')][['constname','levelratiospl']].append(endog).dropna(),x='constname', y='levelratiospl',\
            split=True,order=['none']+genesorder, palette=[sns.xkcd_rgb['dark grey']]+[sns.xkcd_rgb['dark green']]*8,s=4)
ax.set_xticklabels(['none']+genesorder,rotation=45, horizontalalignment='right')    
plt.ylabel('splicing ratio 2$^{nd}$/1$^{st}$ site [log2]')
plt.xlabel('constitutive splice site')
plt.ylim(-17,15)
f.savefig('./figures/FigS3/FigureS3B_three_intron.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



comb['name_comb']=comb.index.map(lambda x: str(comb.commonname[x]) + '_' + str(comb.constname[x]))

'''
comb['levelratiospl_adjacent']=comb.index.map(lambda x: comb.levelratiospl[x]*-1 if comb.changed[x]=='intron' else comb.levelratiospl[x])


scipy.stats.wilcoxon(comb.pivot_table(index='name_comb',columns='changed',values='levelratiospl_adjacent').dropna().exon,
                     comb.pivot_table(index='name_comb',columns='changed',values='levelratiospl_adjacent').dropna().intron)


WilcoxonResult(statistic=2167.0, pvalue=2.2077336201368034e-29)
mean
exon      0.148442
intron    7.135059
std
exon      4.312391
intron    5.160345
'''


rr,pp=scipy.stats.pearsonr(comb.pivot(index='name_comb',columns='changed',values='levelratiospl').dropna()['exon'],comb.pivot(index='name_comb',columns='changed',values='levelratiospl').dropna()['intron'])
#(0.0024551848229035279, 0.97007273864052046)

f=plt.figure(figsize=(3,3))
sns.regplot(data=comb.pivot(index='name_comb',columns='changed',values='levelratiospl'), x='exon',y='intron')
plt.xlabel('splicing ratio 2$^{nd}$/1$^{st}$ site [log2]\ndownstream exon replaced')
plt.ylabel('splicing ratio 2$^{nd}$/1$^{st}$ site [log2]\nupstream intron replaced')
plt.ylim(-17,15)
plt.xlim(-12,12)
plt.annotate("r={:.2f}".format(rr) + \
             "\np={:.2f}".format(pp),xy=(4,9), fontsize=14)
f.savefig('./figures/FigS3/FigureS3C_three_intronvsexoneffect.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


###


comb=irdf[(irdf.subset == 'ir_comb_with_constitutive')]

comb['exon1']=comb.combination.apply(lambda x: x.split('\t')[0])
comb['intron']=comb.combination.apply(lambda x: x.split('\t')[1])
comb['exon2']=comb.combination.apply(lambda x: x.split('\t')[2])
comb['constname']=comb.index.map(lambda x: 'endogenous' if (comb.exon1[x]=='alt')&(comb.intron[x]=='alt')&(comb.exon2[x]=='alt') else [int(s) for s in comb.combination[x].split('\t') if s.isdigit()][0])

comb['config']=comb.index.map(lambda x: 'endogenous' if (comb.exon1[x]=='alt')&(comb.intron[x]=='alt')&(comb.exon2[x]=='alt') else '')
comb['config']=comb.index.map(lambda x: 'exon1 ' if (comb.exon1[x]!='alt') else comb.config[x])
comb['config']=comb.index.map(lambda x: comb.config[x] + 'intron ' if (comb.intron[x]!='alt') else comb.config[x])
comb['config']=comb.index.map(lambda x: comb.config[x] + 'exon2' if (comb.exon2[x]!='alt') else comb.config[x])

comb['config']=comb.config.apply(lambda x: x[:-1] if (x[-1]==' ') else x)

comb['nnumber']=comb.index.map(lambda x: len(comb[(comb.intron==comb.intron[x])]))


comb['name_comb']=comb.index.map(lambda x: comb.name2[x] + '_' + str(comb.constname[x]))

combpivoted=comb[(comb.config!='endogenous')].pivot_table(index='name_comb',\
        columns='config',values='levelratiospl')
combpivoted.columns=['exon up','both exons','exon up+intron','full construct','exon down','intron','intron+exon down']

colorder=['exon up','intron','exon down','both exons','exon up+intron','intron+exon down','full construct']

f=plt.figure(figsize=(4,4))
ax=sns.boxplot(data=combpivoted, color=sns.xkcd_rgb['light blue'], order=colorder)
sns.swarmplot(data=combpivoted, color=sns.xkcd_rgb['medium blue'], order=colorder, size=3)
ax.set_xticklabels(colorder, rotation=45, horizontalalignment='right')
f.savefig('./figures/FigS3/FigureS3F_ir_distributionsofgroups.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

tests=pd.DataFrame('',index=colorder,columns=colorder)
for i in colorder:
    for j in colorder:
        tests.loc[i,j]=scipy.stats.mannwhitneyu(combpivoted[i],combpivoted[j])[1]

f=plt.figure(figsize=(6,4))
ax=sns.heatmap(tests.astype(float), annot=tests, fmt='.0e', annot_kws={'fontsize':12})
ax.set_xticklabels(colorder, rotation=45, horizontalalignment='right')
f.savefig('./figures/FigS3/FigureS3G_ir_heatmap_pvals_diffbetweengroups.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=sns.clustermap(combpivoted.corr(),center=0, figsize=(4,4), cbar_kws={'ticks':[-1,0,1]})
plt.setp(f.ax_heatmap.get_xticklabels(), rotation=45, horizontalalignment='right')
plt.setp(f.ax_heatmap.get_yticklabels(), rotation=0, horizontalalignment='left')
plt.setp(f.ax_heatmap.set_xlabel(''))
plt.setp(f.ax_heatmap.set_ylabel(''))
#        f.cax.set_visible(False)
f.savefig('./figures/Fig3/Figure3E_ir_heatmap.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)



f=plt.figure(figsize=(3,3))
ax=sns.regplot(data=combpivoted, x='exon up',y='exon down', fit_reg=False)
plt.xlim(-15,15)
plt.ylim(-15,15)
plt.xlabel('splicing ratio\nby upstream exon')
plt.ylabel('splicing ratio\nby downstream exon')
plt.annotate("r={:.2f}".format(scipy.stats.pearsonr(combpivoted[['exon up','exon down']].dropna()['exon up'], combpivoted[['exon up','exon down']].dropna()['exon down'])[0]) + \
             ", p={:.0e}".format(scipy.stats.pearsonr(combpivoted[['exon up','exon down']].dropna()['exon up'], combpivoted[['exon up','exon down']].dropna()['exon down'])[1]),xy=(-10,12), fontsize=14)
f.savefig('./figures/FigS3/FigureS3D_ir_exon1VSexon2.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(3,3))
ax=sns.regplot(data=combpivoted, x='exon up',y='intron', fit_reg=False)
plt.xlim(-15,15)
plt.ylim(-15,15)
plt.xlabel('splicing ratio\nby upstream exon')
plt.ylabel('splicing ratio\nby intron')
plt.annotate("r={:.2f}".format(scipy.stats.pearsonr(combpivoted[['exon up','intron']].dropna().intron, combpivoted[['exon up','intron']].dropna()['exon up'])[0]) + \
             ", p={:.0e}".format(scipy.stats.pearsonr(combpivoted[['exon up','intron']].dropna().intron, combpivoted[['exon up','intron']].dropna()['exon up'])[1]),xy=(-10,12), fontsize=14)
f.savefig('./figures/FigS3/FigureS3D_ir_exon1VSintron.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)
  
f=plt.figure(figsize=(3,3))
ax=sns.regplot(data=combpivoted, x='exon down',y='intron', fit_reg=False)
plt.xlim(-15,15)
plt.ylim(-15,15)
plt.xlabel('splicing ratio\nby downstream exon')
plt.ylabel('splicing ratio\nby intron')
plt.annotate("r={:.2f}".format(scipy.stats.pearsonr(combpivoted[['exon down','intron']].dropna().intron, combpivoted[['exon down','intron']].dropna()['exon down'])[0]) + \
             ", p={:.0e}".format(scipy.stats.pearsonr(combpivoted[['exon down','intron']].dropna().intron, combpivoted[['exon down','intron']].dropna()['exon down'])[1]),xy=(-10,12), fontsize=14)
f.savefig('./figures/FigS3/FigureS3D_ir_exon2VSintron.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(3,3))
ax=sns.regplot(data=combpivoted, x='both exons',y='intron', fit_reg=False)
plt.xlim(-15,15)
plt.ylim(-15,15)
plt.xlabel('splicing ratio\nby both exons')
plt.ylabel('splicing ratio\nby intron')
plt.annotate("r={:.2f}".format(scipy.stats.pearsonr(combpivoted[['both exons','intron']].dropna().intron, combpivoted[['both exons','intron']].dropna()['both exons'])[0]) + \
             ", p={:.0e}".format(scipy.stats.pearsonr(combpivoted[['both exons','intron']].dropna().intron, combpivoted[['both exons','intron']].dropna()['both exons'])[1]),xy=(-10,12), fontsize=14)
f.savefig('./figures/FigS3/FigureS3D_ir_bothexonsVSintron.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(3,3))
ax=sns.regplot(data=combpivoted, x='full construct',y='intron', fit_reg=False)
plt.xlim(-15,15)
plt.ylim(-15,15)
plt.xlabel('splicing ratio\nby full construct')
plt.ylabel('splicing ratio\nby intron')
plt.annotate("r={:.2f}".format(scipy.stats.pearsonr(combpivoted[['full construct','intron']].dropna().intron, combpivoted[['full construct','intron']].dropna()['full construct'])[0]) + \
             ", p={:.0e}".format(scipy.stats.pearsonr(combpivoted[['full construct','intron']].dropna().intron, combpivoted[['full construct','intron']].dropna()['full construct'])[1]),xy=(-10,12), fontsize=14)
f.savefig('./figures/FigS3/FigureS3D_ir_fullconstructVSintron.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


### cassette exons


comb=casdf[(casdf.subset == 'cassette_comb_with_constitutive')]

comb['exon1']=comb.combination.apply(lambda x: x.split('\t')[0])
comb['intron']=comb.combination.apply(lambda x: x.split('\t')[1])
comb['exon2']=comb.combination.apply(lambda x: x.split('\t')[2])
comb['constname']=comb.index.map(lambda x: 'endogenous' if (comb.exon1[x]=='alt')&(comb.intron[x]=='alt')&(comb.exon2[x]=='alt') else [int(s) for s in comb.combination[x].split('\t') if s.isdigit()][0])

comb['config']=comb.index.map(lambda x: 'endogenous' if (comb.exon1[x]=='alt')&(comb.intron[x]=='alt')&(comb.exon2[x]=='alt') else '')
comb['config']=comb.index.map(lambda x: 'intron1 ' if (comb.exon1[x]!='alt') else comb.config[x])
comb['config']=comb.index.map(lambda x: comb.config[x] + 'exon ' if (comb.intron[x]!='alt') else comb.config[x])
comb['config']=comb.index.map(lambda x: comb.config[x] + 'intron2' if (comb.exon2[x]!='alt') else comb.config[x])

comb['config']=comb.config.apply(lambda x: x[:-1] if (x[-1]==' ') else x)

comb['nnumber']=comb.index.map(lambda x: len(comb[(comb.intron==comb.intron[x])]))

comb['name_comb']=comb.index.map(lambda x: comb.commonname[x] + '_' + str(comb.constname[x]))

combpivoted=comb[(comb.config!='endogenous')].pivot_table(index='name_comb',\
        columns='config',values='levelratiospl')
combpivoted.columns=['exon','exon+intron down','intron up','intron up+exon','full construct','both introns','intron down']

colorder=['intron up','exon','intron down','both introns','intron up+exon','exon+intron down','full construct']

f=plt.figure(figsize=(4,4))
ax=sns.boxplot(data=combpivoted, color=sns.xkcd_rgb['light blue'], order=colorder)
sns.swarmplot(data=combpivoted, color=sns.xkcd_rgb['medium blue'], order=colorder, size=3)
ax.set_xticklabels(colorder, rotation=45, horizontalalignment='right')
f.savefig('./figures/FigS3/FigureS3H_cas_distributionsofgroups.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

tests=pd.DataFrame('',index=colorder,columns=colorder)
for i in colorder:
    for j in colorder:
        tests.loc[i,j]=scipy.stats.mannwhitneyu(combpivoted[i],combpivoted[j])[1]

f=plt.figure(figsize=(6,4))
ax=sns.heatmap(tests.astype(float), annot=tests, annot_kws={'fontsize':12})
ax.set_xticklabels(colorder, rotation=45, horizontalalignment='right')
f.savefig('./figures/FigS3/FigureS3I_cas_heatmap_pvals_diffbetweengroups.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=sns.clustermap(combpivoted.corr(),center=0, figsize=(4,4), cbar_kws={'ticks':[-1,0,1]})
plt.setp(f.ax_heatmap.get_xticklabels(), rotation=45, horizontalalignment='right')
plt.setp(f.ax_heatmap.get_yticklabels(), rotation=0, horizontalalignment='left')
plt.setp(f.ax_heatmap.set_xlabel(''))
plt.setp(f.ax_heatmap.set_ylabel(''))
f.savefig('./figures/Fig3/Figure3F_cas_heatmap.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


combpivoted=comb[([comb.config[x] in ['exon','intron1','intron2','intron1 intron2','intron1 exon intron2'] for x in comb.index])].pivot(index='name_comb',\
            columns='config',values='levelratiospl').dropna()

f=plt.figure(figsize=(3,3))
ax=sns.regplot(data=combpivoted, x='intron1',y='exon', fit_reg=False)
plt.xlim(-12,8)
plt.ylim(-12,8)
plt.xlabel('splicing ratio\nby upstream intron')
plt.ylabel('splicing ratio\nby exon')
plt.annotate("r={:.2f}".format(scipy.stats.pearsonr(combpivoted.exon, combpivoted.intron1)[0]) + \
             ", p={:.0e}".format(scipy.stats.pearsonr(combpivoted.exon, combpivoted.intron1)[1]),xy=(-11,6.5), fontsize=14)
f.savefig('./figures/FigS3/FigureS3E_cas_intron1VSexon.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)
  

f=plt.figure(figsize=(3,3))
ax=sns.regplot(data=combpivoted, x='intron2',y='exon', fit_reg=False)
plt.xlim(-12,8)
plt.ylim(-12,8)
plt.xlabel('splicing ratio\nby downstream intron')
plt.ylabel('splicing ratio\nby exon')
plt.annotate("r={:.2f}".format(scipy.stats.pearsonr(combpivoted.exon, combpivoted.intron2)[0]) + \
             ", p={:.0e}".format(scipy.stats.pearsonr(combpivoted.exon, combpivoted.intron2)[1]),xy=(-11,6.5), fontsize=14)
f.savefig('./figures/FigS3/FigureS3E_cas_intron2VSexon.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(3,3))
ax=sns.regplot(data=combpivoted, x='intron1 intron2',y='exon', fit_reg=False)
plt.xlim(-12,8)
plt.ylim(-12,8)
plt.xlabel('splicing ratio\nby both introns')
plt.ylabel('splicing ratio\nby exon')
plt.annotate("r={:.2f}".format(scipy.stats.pearsonr(combpivoted.exon, combpivoted['intron1 intron2'])[0]) + \
             ", p={:.0e}".format(scipy.stats.pearsonr(combpivoted.exon, combpivoted['intron1 intron2'])[1]),xy=(-11,6.5), fontsize=14)
f.savefig('./figures/FigS3/FigureS3E_cas_intron1intron2VSexon.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(3,3))
ax=sns.regplot(data=combpivoted, x='intron1 exon intron2',y='exon', fit_reg=False)
plt.xlim(-12,8)
plt.ylim(-12,8)
plt.xlabel('splicing ratio\nfull construct')
plt.ylabel('splicing ratio\nby exon')
plt.annotate("r={:.2f}".format(scipy.stats.pearsonr(combpivoted.exon, combpivoted['intron1 exon intron2'])[0]) + \
             ", p={:.0e}".format(scipy.stats.pearsonr(combpivoted.exon, combpivoted['intron1 exon intron2'])[1]),xy=(-11,6.5), fontsize=14)
f.savefig('./figures/FigS3/FigureS3E_cas_allVSexon.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)
 

