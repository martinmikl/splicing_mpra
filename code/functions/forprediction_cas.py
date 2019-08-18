#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import RNA
import pickle
from maxentpy import maxent
from Bio.motifs import matrix
from Bio.motifs import parse
from Bio.Seq import Seq
from Bio.Seq import IUPAC
from sklearn.ensemble import GradientBoostingRegressor
from scipy.stats import pearsonr
from sklearn import metrics
from sklearn.cross_validation import cross_val_predict



#%% for revision

def get_maxent(df):
    mxnt=pd.DataFrame(index=df.index)
    mxnt['maxent5']=df.index.map(lambda x: maxent.score5(str(df.sequence[x][int(df.exonend[x])-3:int(df.exonend[x])+6])))
    mxnt['maxent3']=df.index.map(lambda x: maxent.score3(str(df.sequence[x][int(df.exonstart[x])-20:int(df.exonstart[x])+3])))
    return mxnt

def get_secstruct(df):
    sec=pd.DataFrame(index=df.index)
    sec['intron1']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.exonstart[x])-24:int(df.exonstart[x])]))[1])
    sec['acceptor']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.exonstart[x])-12:int(df.exonstart[x])+12]))[1])
    sec['exon5']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.exonstart[x]):int(df.exonstart[x])+24]))[1])
    sec['exoncenter']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.exonstart[x])+(int(df.exonend[x])-int(df.exonstart[x]))/2-12:int(df.exonstart[x])+(int(df.exonend[x])-int(df.exonstart[x]))/2+12]))[1])
    sec['exon3']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.exonend[x])-24:int(df.exonend[x])]))[1])
    sec['donor']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.exonend[x])-12:int(df.exonend[x])+12]))[1])
    sec['intron2']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.exonend[x]):int(df.exonend[x])+24]))[1])
    return sec

def get_sixmers(df):
    Xs_sel=pd.read_pickle('./dataframes/ml/cas/Xs_cas_sixmer_rna_sel.pkl')
    mtfs=[]
    for i in Xs_sel.columns:
        if i.split('_')[0] not in mtfs:
            mtfs.append(i.split('_')[0])
    
    sixmerfeatures=pd.DataFrame(index=df.index)
    for i in mtfs:
        sixmerfeatures[i + '_intron1'] = df.index.map(lambda x: df.sequence[x][:df.exonstart[x]].count(i))
        sixmerfeatures[i + '_exon'] = df.index.map(lambda x: df.sequence[x][df.exonstart[x]:df.exonend[x]].count(i))
        sixmerfeatures[i + '_intron2'] = df.index.map(lambda x: df.sequence[x][int(df.exonend[x]):].count(i))
    return sixmerfeatures[Xs_sel.columns]
    
def get_motifscores(df):      
    with open('../additional/ATtRACT/pwm_transposed.txt', 'r') as f:
        records=parse(f, 'jaspar')
    Xs_sel=pd.read_pickle('./dataframes/ml/cas/Xs_cas_motifs_rna_sel_.pkl')
    mtfs=[]
    for i in Xs_sel.columns:
        if i.split('__')[0] not in mtfs:
            mtfs.append(i.split('__')[0])
    
    def find_motifs_cas(varid):
        motifs=pd.Series()
        for pos, seq in mm.counts.log_odds().search(Seq(df.sequence[varid], \
                        alphabet=IUPAC.IUPACUnambiguousDNA()), threshold=0, both=False):
            motifs.loc[pos]=seq
        motifs_up=motifs[motifs.index<df.exonstart[varid]]
        motifs_alt=motifs[(motifs.index>df.exonstart[varid])&(motifs.index<df.exonend[varid])]
        motifs_down=motifs[motifs.index>df.exonend[varid]]
        
        return list([motifs_up.sum(), motifs_alt.sum(), motifs_down.sum()])
    
    juncasmotifs=pd.DataFrame(index=df.index)
    
    database=pd.read_table('../additional/ATtRACT/ATtRACT_db.txt')
    database.drop_duplicates('Matrix_id', inplace=True)
    database.set_index('Matrix_id', inplace=True)
    
    for mm in records:
#        if (mm.name in database[database.Organism=='Homo_sapiens'].index):
        if (mm.name in mtfs):
            mm.counts.__class__=matrix.PositionWeightMatrix
            juncasmotifs['motifs']=df.index.map(lambda x: find_motifs_cas(x))
            juncasmotifs[[str(mm.name) + '__score_motifs_up',str(mm.name) + '__score_motifs_alt',\
                          str(mm.name) +  '__score_motifs_down']]=juncasmotifs.motifs.apply(lambda x: pd.Series(x))
            juncasmotifs.drop('motifs',axis=1, inplace=True)
    
    return juncasmotifs[Xs_sel.columns]


def make_features_cas(filename):
    '''
    input is a csv file with four columns: identifier, sequence, exonstart, exonend; 0-based, No header
    '''
    df=pd.read_csv(filename, header=None)
    df.columns=['identifier','sequence','exonstart','exonend']
    df['exonstart']=df['exonstart'].astype(int)
    df['exonend']=df['exonend'].astype(int)
    df.set_index('identifier', inplace=True)
    mxnt=get_maxent(df)
    sec=get_secstruct(df)
    motifs=get_motifscores(df)
    sixmers=get_sixmers(df)
    features=pd.concat([mxnt,sec,motifs,sixmers], axis=1)
    return features

def train_gbr_model(Xs, Ys, modelname):
    gbr_model=GradientBoostingRegressor(learning_rate=0.2, n_estimators=300,\
                                        max_depth=5).fit(Xs,Ys)    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'.sav','wb'))

def train_gbr_model_maxent(Xs, Ys, modelname):
    cols=['maxent5','maxent3']
    gbr_model=GradientBoostingRegressor(learning_rate=0.2, n_estimators=300,\
                                        max_depth=4).fit(Xs[cols],Ys)
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_maxent.sav','wb'))

def train_gbr_model_secondary(Xs, Ys, modelname):
    cols=['intron1','acceptor','exon5','exoncenter','exon3','donor','intron2']
    gbr_model=GradientBoostingRegressor(learning_rate=0.1, n_estimators=300,\
                                        max_depth=6).fit(Xs[cols],Ys)
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_secondary.sav','wb'))

def train_gbr_model_maxentsecondary(Xs, Ys, modelname):
    cols=['maxent5','maxent3','intron1','acceptor','exon5','exoncenter','exon3','donor','intron2']
    gbr_model=GradientBoostingRegressor(learning_rate=0.3, n_estimators=300,\
                                        max_depth=5).fit(Xs[cols],Ys)
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_maxentsecondary.sav','wb'))

def train_gbr_model_hexamers(Xs, Ys, modelname):
    cols=[x for x in Xs.columns if (x.split('_')[-1] in ['intron1','exon','intron2'])&(len(x.split('_'))==2)]
    gbr_model=GradientBoostingRegressor(learning_rate=0.2, n_estimators=200,\
                                        max_depth=5).fit(Xs[cols],Ys)
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_hexamers.sav','wb'))

def train_gbr_model_motifs(Xs, Ys, modelname):
    cols=[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)]
    gbr_model=GradientBoostingRegressor(learning_rate=0.1, n_estimators=200,\
                                        max_depth=6).fit(Xs[cols],Ys)
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_motifs.sav','wb'))

def train_gbr_model_hexamersmotifs(Xs, Ys, modelname):
    cols=[x for x in Xs.columns if (x.split('_')[-1] in ['intron1','exon','intron2'])&(len(x.split('_'))==2)]+[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)]
    gbr_model=GradientBoostingRegressor(learning_rate=0.1, n_estimators=200,\
                                        max_depth=6).fit(Xs[cols],Ys)
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_hexamersmotifs.sav','wb'))
    
def make_prediction(Xs, Ys=pd.Series([0]), wtratios=None, mode='logratio', filename=None):
    '''
    possible modes:
        logratio (default)
        psi
        deltapsi - wtratios must be provided
        [str] this model is loaded and used for prediction
    saves image to file if filename is provided
    '''
    if Ys.sum()!=0:
        if mode=='logratio':
            mdl=pickle.load(open('./ml/models/casmodel_logratio_trained.sav', 'rb'))
            Zs=pd.DataFrame(Xs.copy())
            Zs['measured']=pd.Series(Ys)
            Zs.dropna(inplace=True)
        elif mode=='psi':
            mdl=pickle.load(open('./ml/models/casmodel_psi_trained.sav', 'rb'))
            Zs=pd.DataFrame(Xs.copy())
            Zs['measured']=pd.Series(Ys)*100
            Zs.dropna(inplace=True)
        elif mode=='deltapsi':
            mdl=pickle.load(open('./ml/models/casmodel_normpsi_trained.sav', 'rb'))
            Zs=pd.DataFrame(data=wtratios, index=Xs.index, columns=['wtpsi']).join(pd.DataFrame(Xs.copy()))
            Zs['measured']=pd.Series(Ys)*100
            Zs.dropna(inplace=True)
        else:
            mdl=pickle.load(open('./ml/models/'+mode+'.sav', 'rb'))
            Zs=pd.DataFrame(Xs.copy())
            Zs['measured']=pd.Series(Ys)
            Zs.dropna(inplace=True)
        ypred=mdl.predict(Zs.iloc[:,:-1])
        if mode=='psi':
            ypred_corr=[0 if x<0 else x for x in ypred]
            ypred=[1 if x>1 else x for x in ypred_corr]
        f=plt.figure(figsize=(3,3))
        plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
        plt.title('$\mathregular{R^2}$ score=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
                  '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=14)
        if filename!=None:
            f.savefig('./ml/plots/'+filename +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    else:
        if mode=='logratio':
            mdl=pickle.load(open('./ml/models/casmodel_logratio_trained.sav', 'rb'))
            Zs=pd.DataFrame(Xs.copy())
            Zs.dropna(inplace=True)
        elif mode=='psi':
            mdl=pickle.load(open('./ml/models/casmodel_psi_trained.sav', 'rb'))
            Zs=pd.DataFrame(Xs.copy())
            Zs.dropna(inplace=True)
        elif mode=='deltapsi':
            mdl=pickle.load(open('./ml/models/casmodel_normpsi_trained.sav', 'rb'))
            Zs=pd.DataFrame(data=wtratios, index=Xs.index, columns=['wtpsi']).join(pd.DataFrame(Xs.copy()))
            Zs.dropna(inplace=True)
        else:
            mdl=pickle.load(open('./ml/models/'+mode+'.sav', 'rb'))
            Zs=pd.DataFrame(Xs.copy())
            Zs.dropna(inplace=True)
        ypred=mdl.predict(Zs.iloc[:,:-1])
        if mode=='psi':
            ypred_corr=[0 if x<0 else x for x in ypred]
            ypred=[1 if x>1 else x for x in ypred_corr]
        
        return pd.Series(ypred, index=Zs.index)

def train_all_models(Xs, Ys, modelname):
    '''
    trains and saves models separately based on maxentscan, 
    secondary structure, hexamer counts and motif scores, 
    as well as a combined model
    '''
    train_gbr_model_maxent(Xs, Ys, modelname)
    train_gbr_model_secondary(Xs, Ys, modelname)
    train_gbr_model_maxentsecondary(Xs, Ys, modelname)
    train_gbr_model_hexamers(Xs, Ys, modelname)
    train_gbr_model_motifs(Xs, Ys, modelname)
    train_gbr_model_hexamersmotifs(Xs, Ys, modelname)
    train_gbr_model(Xs, Ys, modelname)

def make_prediction_all_models(Xs, Ys, mode='logratio', filename=None):
    '''
    possible modes:
        logratio (default)
        [str] this model is loaded and used for prediction
    saves image to file if filename is provided
    '''
    predictionscores=pd.DataFrame(columns=['r2score','pearsonr'])
    cols={'maxent':['maxent5','maxent3'],
          'secondary':['intron1','acceptor','exon5','exoncenter','exon3','donor','intron2'],
           'maxentsecondary':['maxent5','maxent3','intron1','acceptor','exon5','exoncenter','exon3','donor','intron2'],
         'hexamers':[x for x in Xs.columns if (x.split('_')[-1] in ['intron1','exon','intron2'])&(len(x.split('_'))==2)],
          'motifs':[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
          'hexamersmotifs':[x for x in Xs.columns if (x.split('_')[-1] in ['intron1','exon','intron2'])&(len(x.split('_'))==2)]+[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                    'all features':':'}
    for feat in ['maxent','secondary','maxentsecondary','hexamers','motifs','hexamersmotifs']:
        if mode=='logratio':
            mdl=pickle.load(open('./ml/models/casmodel_logratio_trained_'+feat+'.sav', 'rb'))
            Zs=pd.DataFrame(Xs[cols[feat]])
            Zs['measured']=pd.Series(Ys)
            Zs.dropna(inplace=True)
        else:
            mdl=pickle.load(open('./ml/models/'+mode+'_'+feat+'.sav', 'rb'))
            Zs=pd.DataFrame(Xs[cols[feat]])
            Zs['measured']=pd.Series(Ys)
            Zs.dropna(inplace=True)
        ypred=mdl.predict(Zs.iloc[:,:-1])
        if mode=='casmodel_psi_trained':
            ypred_corr=[0 if x<0 else x for x in ypred]
            ypred=[1 if x>1 else x for x in ypred_corr]
        predictionscores.loc[feat,'r2score']= metrics.r2_score(Zs.iloc[:,-1], ypred)
        predictionscores.loc[feat,'pearsonr']= pearsonr(Zs.iloc[:,-1], ypred)[0]
        f=plt.figure(figsize=(3,3))
        plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
        plt.title('$\mathregular{R^2}$ score=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
                  '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=14)
        plt.xlim(-15,10)
        plt.ylim(-15,10)
        plt.plot([-15,10],[-15,10], '--',color='gray',linewidth=2,alpha=0.5)
        if filename!=None:
            f.savefig('./ml/plots/scatter_'+filename +'_'+feat+'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    if mode=='logratio':
        mdl=pickle.load(open('./ml/models/casmodel_logratio_trained.sav', 'rb'))
        Zs=pd.DataFrame(Xs[:])
        Zs['measured']=pd.Series(Ys)
        Zs.dropna(inplace=True)
    else:
        mdl=pickle.load(open('./ml/models/'+mode+'.sav', 'rb'))
        Zs=pd.DataFrame(Xs[:])
        Zs['measured']=pd.Series(Ys)
        Zs.dropna(inplace=True)
    ypred=mdl.predict(Zs.iloc[:,:-1])
    if mode=='casmodel_psi_trained':
        ypred_corr=[0 if x<0 else x for x in ypred]
        ypred=[1 if x>1 else x for x in ypred_corr]
    predictionscores.loc['all features','r2score']= metrics.r2_score(Zs.iloc[:,-1], ypred)
    predictionscores.loc['all features','pearsonr']= pearsonr(Zs.iloc[:,-1], ypred)[0]
    f=plt.figure(figsize=(3,3))
    plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
    plt.title('$\mathregular{R^2}$ score=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
              '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=14)
    plt.xlim(-15,10)
    plt.ylim(-15,10)
    plt.plot([-15,10],[-15,10], '--',color='gray',linewidth=2,alpha=0.5)
    if filename!=None:
        f.savefig('./ml/plots/scatter_'+filename +'_all_features.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
        
    f=plt.figure(figsize=(4,3))
    with sns.axes_style('dark'):
        ax=predictionscores.pearsonr.plot(kind='bar', color=sns.xkcd_rgb['medium blue'])
    ax.yaxis.grid()    
    ax.set_xticklabels('')
    plt.ylim(0,np.ceil(predictionscores.pearsonr.max()*10)/10.0)
    plt.ylabel('Pearson r')
    if filename!=None:
        f.savefig('./ml/plots/'+filename+'_modelbuilding_pearsonr_summary.png', \
                  dpi = 300, format='png', bbox_inches='tight', frameon=True)
        f.show()

    f=plt.figure(figsize=(4,3))
    with sns.axes_style('dark'):
        ax=predictionscores.r2score.plot(kind='bar', color=sns.xkcd_rgb['medium blue'])
    ax.yaxis.grid()    
    ax.set_xticklabels('')
    plt.ylim(0,1)
    plt.ylabel('$\mathregular{R^2}$ score')
    if filename!=None:
        f.savefig('./ml/plots/'+filename+'_modelbuilding_R2scores_summary.png', \
                  dpi = 300, format='png', bbox_inches='tight', frameon=True)
        f.show()
    

def predict_splicingratio_all_models(Xs, mode='logratio', filename=None):
    '''
    possible modes:
        logratio (default)
        [str] this model is loaded and used for prediction
    saves image to file if filename is provided
    '''
    predictionscores=pd.DataFrame(columns=Xs.index)
    cols={'maxent':['maxent5','maxent3'],
          'secondary':['intron1','acceptor','exon5','exoncenter','exon3','donor','intron2'],
           'maxentsecondary':['maxent5','maxent3','intron1','acceptor','exon5','exoncenter','exon3','donor','intron2'],
         'hexamers':[x for x in Xs.columns if (x.split('_')[-1] in ['intron1','exon','intron2'])&(len(x.split('_'))==2)],
          'motifs':[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
          'hexamersmotifs':[x for x in Xs.columns if (x.split('_')[-1] in ['intron1','exon','intron2'])&(len(x.split('_'))==2)]+[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                    'all features':':'}
    for feat in ['maxent','secondary','maxentsecondary','hexamers','motifs','hexamersmotifs']:
        if mode=='logratio':
            mdl=pickle.load(open('./ml/models/casmodel_logratio_trained_'+feat+'.sav', 'rb'))
            Zs=pd.DataFrame(Xs[cols[feat]])
            Zs.dropna(inplace=True)
        else:
            mdl=pickle.load(open('./ml/models/'+mode+'_'+feat+'.sav', 'rb'))
            Zs=pd.DataFrame(Xs[cols[feat]])
            Zs.dropna(inplace=True)
        ypred=mdl.predict(Zs)
        if mode=='casmodel_psi_trained':
            ypred_corr=[0 if x<0 else x for x in ypred]
            ypred=[1 if x>1 else x for x in ypred_corr]
        predictionscores.loc[feat]=pd.Series(ypred, index=Xs.index)
    if mode=='logratio':
        mdl=pickle.load(open('./ml/models/casmodel_logratio_trained.sav', 'rb'))
        Zs=pd.DataFrame(Xs[:])
        Zs.dropna(inplace=True)
    else:
        mdl=pickle.load(open('./ml/models/'+mode+'.sav', 'rb'))
        Zs=pd.DataFrame(Xs[:])
        Zs.dropna(inplace=True)
    ypred=mdl.predict(Zs)
    if mode=='casmodel_psi_trained':
        ypred_corr=[0 if x<0 else x for x in ypred]
        ypred=[1 if x>1 else x for x in ypred_corr]
    predictionscores.loc['all features']=pd.Series(ypred, index=Xs.index)
    predictionscores.T.to_csv(filename + '.csv')

            
def make_gbr_cvpredict(Xs, Ys, filename=None):

    Zs=pd.DataFrame(Xs)
    Zs['measured']=pd.Series(Ys)
    Zs.dropna(inplace=True)
    
    gbr=GradientBoostingRegressor(learning_rate=0.2, n_estimators=300,\
                                        max_depth=5)
    ypred=cross_val_predict(gbr, Zs.iloc[:,:-1], Zs.iloc[:,-1], cv=10)
    f=plt.figure(figsize=(3,3))
    plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
    plt.title('$\mathregular{R^2}$ score=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
              '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=14)
    if filename!=None:
        f.savefig('./ml/plots/'+filename +'_predict_by_10foldcv.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    