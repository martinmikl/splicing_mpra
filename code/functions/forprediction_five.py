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
    mxnt['maxent5first']=df.index.map(lambda x: maxent.score5(str(df.sequence[x][int(df.donor1[x])-3:int(df.donor1[x])+6])))
    mxnt['maxent5second']=df.index.map(lambda x: maxent.score5(str(df.sequence[x][int(df.donor2[x])-3:int(df.donor2[x])+6])))
    return mxnt

def get_secstruct(df):
    sec=pd.DataFrame(index=df.index)
    sec['exon']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.donor1[x])-24:int(df.donor1[x])]))[1])
    sec['donor1']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.donor1[x])-12:int(df.donor1[x])+12]))[1])
    sec['alt5']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.donor1[x]):int(df.donor1[x])+24]))[1])
    sec['altcenter']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.donor1[x])+(int(df.donor2[x])-int(df.donor1[x]))/2-12:int(df.donor1[x])+(int(df.donor2[x])-int(df.donor1[x]))/2+12]))[1])
    sec['alt3']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.donor2[x])-24:int(df.donor2[x])]))[1])
    sec['donor2']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.donor2[x])-12:int(df.donor2[x])+12]))[1])
    sec['intron']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.donor2[x]):int(df.donor2[x])+24]))[1])
    return sec

def get_sixmers(df):
    Xs_sel=pd.read_pickle('./dataframes/ml/five/Xs_five_sixmer_rna_sel.pkl')
    mtfs=[]
    for i in Xs_sel.columns:
        if i.split('_')[0] not in mtfs:
            mtfs.append(i.split('_')[0])
    
    sixmerfeatures=pd.DataFrame(index=df.index)
    for i in mtfs:
        sixmerfeatures[i + '_exon'] = df.index.map(lambda x: df.sequence[x][:df.donor1[x]].count(i))
        sixmerfeatures[i + '_alt'] = df.index.map(lambda x: df.sequence[x][df.donor1[x]:df.donor2[x]].count(i))
        sixmerfeatures[i + '_intron'] = df.index.map(lambda x: df.sequence[x][int(df.donor2[x]):].count(i))
    return sixmerfeatures[Xs_sel.columns]
    
def get_motifscores(df):      
    with open('../additional/ATtRACT/pwm_transposed.txt', 'r') as f:
        records=parse(f, 'jaspar')
    Xs_sel=pd.read_pickle('./dataframes/ml/five/Xs_five_motifs_rna_sel_.pkl')
    mtfs=[]
    for i in Xs_sel.columns:
        if i.split('__')[0] not in mtfs:
            mtfs.append(i.split('__')[0])
    
    def find_motifs(varid):
        motifs=pd.Series()
        for pos, seq in mm.counts.log_odds().search(Seq(df.sequence[varid], \
                        alphabet=IUPAC.IUPACUnambiguousDNA()), threshold=0, both=False):
            motifs.loc[pos]=seq
        motifs_up=motifs[motifs.index<df.donor1[varid]]
        motifs_alt=motifs[(motifs.index>df.donor1[varid])&(motifs.index<df.donor2[varid])]
        motifs_down=motifs[motifs.index>df.donor2[varid]]
        
        return list([motifs_up.sum(), motifs_alt.sum(), motifs_down.sum()])
    
    motifscores=pd.DataFrame(index=df.index)
    
    database=pd.read_table('../additional/ATtRACT/ATtRACT_db.txt')
    database.drop_duplicates('Matrix_id', inplace=True)
    database.set_index('Matrix_id', inplace=True)
    
    for mm in records:
        if (mm.name in mtfs):
            mm.counts.__class__=matrix.PositionWeightMatrix
            motifscores['motifs']=df.index.map(lambda x: find_motifs(x))
            motifscores[[str(mm.name) + '__score_motifs_up',str(mm.name) + '__score_motifs_alt',\
                          str(mm.name) +  '__score_motifs_down']]=motifscores.motifs.apply(lambda x: pd.Series(x))
            motifscores.drop('motifs',axis=1, inplace=True)
    
    return motifscores[Xs_sel.columns]

def make_features_five(filename):
    '''
    input is a csv file with four columns: identifier, sequence, donor1, donor2; No header; 0-based
    '''
    df=pd.read_csv(filename, header=None)
    df.columns=['identifier','sequence','donor1','donor2']
    df['donor1']=df['donor1'].astype(int)
    df['donor2']=df['donor2'].astype(int)
    df.set_index('identifier', inplace=True)
    mxnt=get_maxent(df)
    sec=get_secstruct(df)
#    paired=get_paired(df)
    motifs=get_motifscores(df)
    sixmers=get_sixmers(df)
    features=pd.concat([mxnt,sec,motifs,sixmers], axis=1)
    return features

def train_gbr_model(Xs, Ys, modelname):
    gbr_model=GradientBoostingRegressor(learning_rate=0.2, n_estimators=300,\
                                        max_depth=5).fit(Xs,Ys)
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'.sav','wb'))

def train_gbr_model_maxent(Xs, Ys, modelname):
    cols=['maxent5first','maxent5second']
    gbr_model=GradientBoostingRegressor(learning_rate=0.2, n_estimators=300,\
                                        max_depth=4).fit(Xs[cols],Ys)
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_maxent.sav','wb'))

def train_gbr_model_secondary(Xs, Ys, modelname):
    cols=['exon','donor1','alt5','altcenter','alt3','donor2','intron']
    gbr_model=GradientBoostingRegressor(learning_rate=0.2, n_estimators=300,\
                                        max_depth=5).fit(Xs[cols],Ys)
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_secondary.sav','wb'))

def train_gbr_model_maxentsecondary(Xs, Ys, modelname):
    cols=['maxent5first','maxent5second','exon','donor1','alt5','altcenter','alt3','donor2','intron']
    gbr_model=GradientBoostingRegressor(learning_rate=0.3, n_estimators=300,\
                                        max_depth=5).fit(Xs[cols],Ys)
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_maxentsecondary.sav','wb'))

def train_gbr_model_hexamers(Xs, Ys, modelname):
    cols=[x for x in Xs.columns if (x.split('_')[-1] in ['intron','alt','exon'])&(len(x.split('_'))==2)]
    gbr_model=GradientBoostingRegressor(learning_rate=0.2, n_estimators=300,\
                                        max_depth=6).fit(Xs[cols],Ys)
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_hexamers.sav','wb'))

def train_gbr_model_motifs(Xs, Ys, modelname):
    cols=[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)]
    gbr_model=GradientBoostingRegressor(learning_rate=0.1, n_estimators=200,\
                                        max_depth=6).fit(Xs[cols],Ys)
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_motifs.sav','wb'))

def train_gbr_model_hexamersmotifs(Xs, Ys, modelname):
    cols=[x for x in Xs.columns if (x.split('_')[-1] in ['intron','alt','exon'])&(len(x.split('_'))==2)]+[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)]
    gbr_model=GradientBoostingRegressor(learning_rate=0.1, n_estimators=200,\
                                        max_depth=6).fit(Xs[cols],Ys)
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_hexamersmotifs.sav','wb'))

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
    cols={'maxent':['maxent5first','maxent5second'],
          'secondary':['exon','donor1','alt5','altcenter','alt3','donor2','intron'],
            'maxentsecondary':['maxent5first','maxent5second','exon','donor1','alt5','altcenter','alt3','donor2','intron'],
          'hexamers':[x for x in Xs.columns if (x.split('_')[-1] in ['intron','alt','exon'])&(len(x.split('_'))==2)],
          'motifs':[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
            'hexamersmotifs':[x for x in Xs.columns if (x.split('_')[-1] in ['intron','alt','exon'])&(len(x.split('_'))==2)]+[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                    'all features':':'}
    for feat in ['maxent','secondary','maxentsecondary','hexamers','motifs','hexamersmotifs']:
        if mode=='logratio':
            mdl=pickle.load(open('./ml/models/fivemodel_logratio_trained_'+feat+'.sav', 'rb'))
            Zs=pd.DataFrame(Xs[cols[feat]])
            Zs['measured']=pd.Series(Ys)
            Zs.dropna(inplace=True)
        else:
            mdl=pickle.load(open('./ml/models/'+mode+'_'+feat+'.sav', 'rb'))
            Zs=pd.DataFrame(Xs[cols[feat]])
            Zs['measured']=pd.Series(Ys)
            Zs.dropna(inplace=True)
        ypred=mdl.predict(Zs.iloc[:,:-1])
        predictionscores.loc[feat,'r2score']= metrics.r2_score(Zs.iloc[:,-1], ypred)
        predictionscores.loc[feat,'pearsonr']= pearsonr(Zs.iloc[:,-1], ypred)[0]
        f=plt.figure(figsize=(3,3))
        plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
        plt.title('$\mathregular{R^2}$ score=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
                  '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=14)
        plt.xlim(-15,15)
        plt.ylim(-15,15)
        plt.plot([-15,15],[-15,15], '--',color='gray',linewidth=2,alpha=0.5)
        if filename!=None:
            f.savefig('./ml/plots/scatter_'+filename +'_'+feat+'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    if mode=='logratio':
        mdl=pickle.load(open('./ml/models/fivemodel_logratio_trained.sav', 'rb'))
        Zs=pd.DataFrame(Xs.values, index=Xs.index, columns=Xs.columns)
        Zs['measured']=pd.Series(Ys)
        Zs.dropna(inplace=True)
    else:
        mdl=pickle.load(open('./ml/models/'+mode+'.sav', 'rb'))
        Zs=pd.DataFrame(Xs.values, index=Xs.index, columns=Xs.columns)
        Zs['measured']=pd.Series(Ys)
        Zs.dropna(inplace=True)
    ypred=mdl.predict(Zs.iloc[:,:-1])
    predictionscores.loc['all features','r2score']= metrics.r2_score(Zs.iloc[:,-1], ypred)
    predictionscores.loc['all features','pearsonr']= pearsonr(Zs.iloc[:,-1], ypred)[0]
    f=plt.figure(figsize=(3,3))
    plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
    plt.title('$\mathregular{R^2}$ score=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
              '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=14)
    plt.xlim(-15,15)
    plt.ylim(-15,15)
    plt.plot([-15,15],[-15,15], '--',color='gray',linewidth=2,alpha=0.5)
    if filename!=None:
        f.savefig('./ml/plots/scatter_'+filename +'_all_features.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
        
    f=plt.figure(figsize=(4,3))
    with sns.axes_style('dark'):
        ax=predictionscores.pearsonr.plot(kind='bar', color=sns.xkcd_rgb['medium blue'])
    ax.yaxis.grid()    
    ax.set_xticklabels('')
    plt.ylim(0,1)
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
    cols={'maxent':['maxent5first','maxent5second'],
          'secondary':['exon','donor1','alt5','altcenter','alt3','donor2','intron'],
            'maxentsecondary':['maxent5first','maxent5second','exon','donor1','alt5','altcenter','alt3','donor2','intron'],
          'hexamers':[x for x in Xs.columns if (x.split('_')[-1] in ['intron','alt','exon'])&(len(x.split('_'))==2)],
          'motifs':[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
            'hexamersmotifs':[x for x in Xs.columns if (x.split('_')[-1] in ['intron','alt','exon'])&(len(x.split('_'))==2)]+[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                    'all features':':'}
    for feat in ['maxent','secondary','maxentsecondary','hexamers','motifs','hexamersmotifs']:
        if mode=='logratio':
            mdl=pickle.load(open('./ml/models/fivemodel_logratio_trained_'+feat+'.sav', 'rb'))
            Zs=pd.DataFrame(Xs[cols[feat]])
            Zs.dropna(inplace=True)
        else:
            mdl=pickle.load(open('./ml/models/'+mode+'_'+feat+'.sav', 'rb'))
            Zs=pd.DataFrame(Xs[cols[feat]])
            Zs.dropna(inplace=True)
        predictionscores.loc[feat]= pd.Series(mdl.predict(Zs), index=Xs.index)
    if mode=='logratio':
        mdl=pickle.load(open('./ml/models/fivemodel_logratio_trained.sav', 'rb'))
        Zs=pd.DataFrame(Xs.values, index=Xs.index, columns=Xs.columns)
        Zs.dropna(inplace=True)
    else:
        mdl=pickle.load(open('./ml/models/'+mode+'.sav', 'rb'))
        Zs=pd.DataFrame(Xs.values, index=Xs.index, columns=Xs.columns)
        Zs.dropna(inplace=True)
    predictionscores.loc['all features']= pd.Series(mdl.predict(Zs), index=Xs.index)
    predictionscores.T.to_csv(filename + '.csv')
        

def make_prediction(Xs, Ys, mode='logratio', filename=None):
    '''
    possible modes:
        logratio (default)
        [str] this model is loaded and used for prediction
    saves image to file if filename is provided
    '''
    if mode=='logratio':
        mdl=pickle.load(open('./ml/models/fivemodel_logratio_trained.sav', 'rb'))
        Zs=pd.DataFrame(Xs)
        Zs['measured']=pd.Series(Ys)
        Zs.dropna(inplace=True)
    else:
        mdl=pickle.load(open('./ml/models/'+mode+'.sav', 'rb'))
        Zs=pd.DataFrame(Xs)
        Zs['measured']=pd.Series(Ys)
        Zs.dropna(inplace=True)
    ypred=mdl.predict(Zs.iloc[:,:-1])
    f=plt.figure(figsize=(3,3))
    plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
    plt.title('r2=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
              '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=14)
    if filename!=None:
        f.savefig('./ml/plots/'+filename +'.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
        
def cvpredict_by_gene(Xs, Ys, df, genelist=None):
    Zs=pd.DataFrame(Xs.copy())
    Zs['measured']=pd.Series(Ys)
    Zs.dropna(inplace=True)
    dfml=df.loc[Zs.index].dropna(how='all')
    predictionscores=pd.DataFrame()
    if genelist==None:
        genelist=dfml.commonname.dropna().unique()
    for i in genelist:
        if len(dfml[dfml.commonname==i])>50:
            gbr=GradientBoostingRegressor(learning_rate=0.2, n_estimators=300,\
                    max_depth=5).fit(Zs.loc[dfml[dfml.commonname!=i].index].iloc[:,:-1],Zs.loc[dfml[dfml.commonname!=i].index].iloc[:,-1])
            ypred=gbr.predict(Zs.loc[dfml[dfml.commonname==i].index].iloc[:,:-1])
            predictionscores.loc[i,'r2score']= metrics.r2_score(Zs.loc[dfml[dfml.commonname==i].index].iloc[:,-1], ypred)
            predictionscores.loc[i,'pearsonr']= pearsonr(Zs.loc[dfml[dfml.commonname==i].index].iloc[:,-1], ypred)[0]
    return predictionscores
        

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
    
    