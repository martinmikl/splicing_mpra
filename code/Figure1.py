#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 17:31:33 2019

@author: martinm
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

sns.set_context('poster')

os.mkdir('./figures/Fig1')
os.mkdir('./figures/FigS1')

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

f=plt.figure(figsize=(3,3))
plt.hist(irdf[(irdf.rnareads_effective>100)].levelratiospl, bins=100, orientation=u'horizontal', linewidth=0)
plt.ylabel('splicing ratio\nspliced/unspliced [log2]')
plt.ylim(-15,15)
plt.xlabel('# variants')
f.savefig('./figures/Fig1/Figure1B_ir.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(3,3))
plt.hist(casdf[(casdf.rnareads_effective>100)].levelratiospl, bins=100, orientation=u'horizontal', linewidth=0)
plt.ylabel('splicing ratio\nexon included/skipped [log2]')
plt.ylim(-15,15)
plt.xlabel('# variants')
plt.xticks([0,50,100,150])
plt.xlim(0,150)
f.savefig('./figures/Fig1/Figure1B_cas.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(3,3))
plt.hist(fivedf[(fivedf.rnareads_effective>100)].levelratiospl, bins=100, orientation=u'horizontal', linewidth=0)
plt.ylabel('splicing ratio\n2$^{nd}$/1$^{st}$ site [log2]')
plt.ylim(-15,15)
plt.xlabel('# variants')
plt.xticks([0,50,100,150])
plt.xlim(0,150)
f.savefig('./figures/Fig1/Figure1B_five.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(3,3))
plt.hist(threedf[(threedf.rnareads_effective>100)].levelratiospl, bins=100, orientation=u'horizontal', linewidth=0)
plt.ylabel('splicing ratio\n2$^{nd}$/1$^{st}$ site [log2]')
plt.ylim(-15,15)
plt.xlabel('# variants')
plt.xticks([0,50,100,150])
plt.xlim(0,150)
f.savefig('./figures/Fig1/Figure1B_three.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)



f=plt.figure(figsize=(3,2))
irdf.wtratio.drop_duplicates().hist()
plt.xlim(-15,15)
plt.title('mean +/- s.d. : '+'{:.2f}'.format(irdf.wtratio.drop_duplicates().mean())\
          + '+/-'+'{:.2f}'.format(irdf.wtratio.drop_duplicates().std()),fontsize=14)
f.savefig('./figures/Fig1/FigureS1B_ir_wtratios_histogram.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(3,2))
casdf.wtratio.drop_duplicates().hist()
plt.xlim(-15,15)
plt.title('mean +/- s.d. : '+'{:.2f}'.format(casdf.wtratio.drop_duplicates().mean())\
          + '+/-'+'{:.2f}'.format(casdf.wtratio.drop_duplicates().std()),fontsize=14)
f.savefig('./figures/Fig1/FigureS1B_cas_wtratios_histogram.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(3,2))
fivedf.wtratio.drop_duplicates().hist()
plt.xlim(-15,15)
plt.title('mean +/- s.d. : '+'{:.2f}'.format(fivedf.wtratio.drop_duplicates().mean())\
          + '+/-'+'{:.2f}'.format(fivedf.wtratio.drop_duplicates().std()),fontsize=14)
f.savefig('./figures/Fig1/FigureS1B_five_wtratios_histogram.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(3,2))
threedf.wtratio.drop_duplicates().hist()
plt.xlim(-15,15)
plt.title('mean +/- s.d. : '+'{:.2f}'.format(threedf.wtratio.drop_duplicates().mean())\
          + '+/-'+'{:.2f}'.format(threedf.wtratio.drop_duplicates().std()),fontsize=14)
f.savefig('./figures/Fig1/FigureS1B_three_wtratios_histogram.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

####


def plot_barcode_controls(df, samplename):
    samevarseq=[]
    for var, group in df.groupby(by='varseq162'):
        if (len(group)>5):
            samevarseq.append(var)
    barcodectrlratio= []
    if (len(samevarseq)>0):
        for i in samevarseq:
            indexes=list(df[df.varseq162==i].index)
            barcodectrlratio.append(df.loc[indexes,'levelratiospl'].values)
    
        bcratio=pd.DataFrame(barcodectrlratio)
        bcratio['meanval']=bcratio.mean(axis=1)
        bcratio.sort_values(by='meanval', inplace=True)
        bcratio=bcratio.transpose()
        bcratio.drop('meanval', inplace=True)
        bcratio.dropna(how='all',axis=1)
        
        f=plt.figure(figsize=(12,3))
        ax=sns.boxplot(data=bcratio, palette='Blues',fliersize=3)
        plt.ylim([-15,15])
        ax.set_xticklabels('')
        plt.xlabel('groups of multiple barcodes for same variant')
        plt.ylabel('splicing ratio [log2]')
        f.savefig('./figures/'+ samplename + '_multiple_barcodes_RNA_logratio.png', \
                  dpi = 300, format='png', bbox_inches='tight', frameon=True)
        f.show()

plot_barcode_controls(irdf[(irdf.rnareads_effective>100)],'Fig1/Fig1C_ir')
plot_barcode_controls(casdf[(casdf.rnareads_effective>100)],'FigS1/FigS1C_cas')
plot_barcode_controls(fivedf[(fivedf.rnareads_effective>100)],'FigS1/FigS1C_five')
plot_barcode_controls(threedf[threedf.rnareads_effective>100],'FigS1/FigS1C_three')

#%% technical replicates

irdf_reps=pd.read_pickle('./dataframes/irdf_reps.pkl')

f=plt.figure(figsize=(3,3))
plt.scatter(irdf_reps.levelratiospl_from_HS_0hold_0h, irdf_reps.levelratiospl_from_HS_0hnew_0h, alpha=0.1)
plt.xlabel('replicate 1')
plt.ylabel('replicate 2')
plt.annotate("r={:.3f}".format(irdf_reps[['levelratiospl_from_HS_0hold_0h','levelratiospl_from_HS_0hnew_0h']].corr().values[0][1]), xy=(-13,10), fontsize=14)
f.savefig('./figures/FigS1/FigureS1A_ir.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


casdf_reps=pd.read_pickle('./dataframes/casdf_reps.pkl')

f=plt.figure(figsize=(3,3))
plt.scatter(casdf_reps[casdf_reps.rnareads_effective>100].levelratiospl4, casdf_reps[casdf_reps.rnareads_effective>100].levelratiospl5, alpha=0.1)
plt.xlabel('replicate 1')
plt.ylabel('replicate 2')
plt.annotate("r={:.3f}".format(casdf_reps[(casdf_reps.levelratiospl4>-20)&(casdf_reps.levelratiospl4<20)&(casdf_reps.levelratiospl5>-20)&(casdf_reps.levelratiospl5<20)][['levelratiospl4','levelratiospl5']].corr().iloc[0,1]), xy=(-13,5), fontsize=14)
f.savefig('./figures/FigS1/FigureS1A_cas.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

