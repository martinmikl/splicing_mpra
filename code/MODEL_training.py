#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 15:22:47 2019

@author: martinm
"""

import pandas as pd
import forprediction_five
import forprediction_three
import forprediction_cas
import forprediction_ir

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

#%% Intron retention
### Train model for predicting logratio

irmly=pd.read_pickle('./dataframes/ml/ir/irmly.pkl')
irml=irdf.loc[irmly.index]

### make features with entire varseq
irml['intronstart_wholevarseq']=irml.intronstart_varseq.apply(lambda x: int(x)+30)
irml['intronend_wholevarseq']=irml.intronend_varseqnew.apply(lambda x: int(x)+30)

irml[['varseq','intronstart_wholevarseq','intronend_wholevarseq']].to_csv('./data/irtrain.csv', header=None)

# Make features

#df_features_ireval10=forprediction_ir.make_features_ir('./data/ireval_fortest_first10.csv')
df_features_irtrain=forprediction_ir.make_features_ir('./data/irtrain.csv')

df_features_irtrain.to_pickle('./dataframes/ml/ir/Xs_trainnew_ir.pkl')

Xs=pd.read_pickle('./dataframes/ml/ir/Xs_trainnew_ir.pkl')
Ys=irml[irml.fraction_canonical>0.3].loc[Xs.index,'levelratiospl'].dropna()
Zs=pd.DataFrame(Xs.copy())
Zs['measured']=pd.Series(Ys)
Zs.dropna(inplace=True)

forprediction_ir.train_all_models(Zs.iloc[:,:-1], Zs.iloc[:,-1], 'irmodel_logratio_trained')
'' 


#%% Cassette exons
### Train model for predicting logratio

casmly=pd.read_pickle('./dataframes/ml/cas/casmly.pkl')
casml=casdf.loc[casmly.index]

motifs_train=pd.read_pickle('./dataframes/ml/cas/Xs_cas_motifs_rna_sel_.pkl')
sixmers_train=pd.read_pickle('./dataframes/ml/cas/Xs_cas_sixmer_rna_sel.pkl')

Xs=casml[['maxent5','maxent3','intron1','acceptor','exon5','exoncenter','exon3','donor','intron2']].dropna()
Ys=casml[casml.fraction_canonical>0.5].loc[Xs.index,'levelratiospl'].dropna()
Xs=Xs.loc[Ys.index]

Xs=Xs.join(motifs_train.loc[Xs.index])
Xs=Xs.join(sixmers_train.loc[Xs.index])
Xs.dropna(inplace=True)
Ys=Ys[Xs.index]

Xs.to_pickle('./dataframes/ml/cas/Xs_train_cas.pkl')
Ys.to_pickle('./dataframes/ml/cas/Ys_train_cas.pkl')


Xs=casml[['maxent5','maxent3','intron1','acceptor','exon5','exoncenter','exon3','donor','intron2']].dropna()
Ys=casml[casml.fraction_canonical>0.5].loc[Xs.index,'psi'].dropna()
Xs=Xs.loc[Ys.index]

Xs=Xs.join(motifs_train.loc[Xs.index])
Xs=Xs.join(sixmers_train.loc[Xs.index])
Xs.dropna(inplace=True)
Ys=Ys[Xs.index]

Xs.to_pickle('./dataframes/ml/cas/Xs_train_cas_psi.pkl')
Ys.to_pickle('./dataframes/ml/cas/Ys_train_cas_psi.pkl')



#
# Train model on full feature set and on all subsets separately - logratio
Xs=pd.read_pickle('./dataframes/ml/cas/Xs_train_cas.pkl')
Ys=pd.read_pickle('./dataframes/ml/cas/Ys_train_cas.pkl')

forprediction_cas.train_gbr_model(Xs, Ys, 'casmodel_logratio_trained')
forprediction_cas.train_all_models(Xs, Ys, 'casmodel_logratio_trained')

# Train model on full feature set and on all subsets separately - psi
Xs=pd.read_pickle('./dataframes/ml/cas/Xs_train_cas_psi.pkl')
Ys=pd.read_pickle('./dataframes/ml/cas/Ys_train_cas_psi.pkl')

forprediction_cas.train_gbr_model(Xs, Ys, 'casmodel_psi_trained')
forprediction_cas.train_all_models(Xs, Ys, 'casmodel_psi_trained')



#%% FIVE

# assemble features for training set
fivemly=pd.read_pickle('./dataframes/ml/five/fivemly.pkl')
fiveml=fivedf.loc[fivemly.index]

motifs_train=pd.read_pickle('./dataframes/ml/five/Xs_five_motifs_rna_sel_.pkl')
sixmers_train=pd.read_pickle('./dataframes/ml/five/Xs_five_sixmer_rna_sel.pkl')

Xs=fiveml[['maxent5first','maxent5second','exon','donor1','alt5','altcenter','alt3','donor2','intron']].dropna()
Ys=fiveml[fiveml.fraction_canonical>0.5].loc[Xs.index,'levelratiospl'].dropna()
Xs=Xs.loc[Ys.index]

Xs=Xs.join(motifs_train.loc[Xs.index])
Xs=Xs.join(sixmers_train.loc[Xs.index])
Xs.dropna(inplace=True)
Ys=Ys[Xs.index]

Xs.to_pickle('./dataframes/ml/five/Xs_train_five.pkl')
Ys.to_pickle('./dataframes/ml/five/Ys_train_five.pkl')
#
Xs=pd.read_pickle('./dataframes/ml/five/Xs_train_five.pkl')
Ys=pd.read_pickle('./dataframes/ml/five/Ys_train_five.pkl')

# Train model on full feature set and on all subsets separately
forprediction_five.train_all_models(Xs, Ys, 'fivemodel_logratio_trained')

#%% THREE

# assemble features for training set
threemly=pd.read_pickle('./dataframes/ml/three/threemly.pkl')
threeml=threedf.loc[threemly.index]

motifs_train=pd.read_pickle('./dataframes/ml/three/Xs_three_motifs_rna_sel_.pkl')
sixmers_train=pd.read_pickle('./dataframes/ml/three/Xs_three_sixmer_rna_sel.pkl')

Xs=threeml[['maxent3first','maxent3second','intron','acceptor1','alt5','altcenter','alt3','acceptor2','exon']].dropna()
Ys=threeml[threeml.fraction_canonical>0.3].loc[Xs.index,'levelratiospl'].dropna()
Xs=Xs.loc[Ys.index]

Xs=Xs.join(motifs_train.loc[Xs.index])
Xs=Xs.join(sixmers_train.loc[Xs.index])
Xs.dropna(inplace=True)
Ys=Ys[Xs.index]

Xs.to_pickle('./dataframes/ml/three/Xs_train_three.pkl')
Ys.to_pickle('./dataframes/ml/three/Ys_train_three.pkl')
#
Xs=pd.read_pickle('./dataframes/ml/three/Xs_train_three.pkl')
Ys=pd.read_pickle('./dataframes/ml/three/Ys_train_three.pkl')

# Train model on full feature set and on all subsets separately
forprediction_three.train_all_models(Xs, Ys, 'threemodel_logratio_trained')
