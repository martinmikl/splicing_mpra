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
import itertools



#%% for revision

def get_maxent(df):
    mxnt=pd.DataFrame(index=df.index)
    mxnt['maxent5']=df.index.map(lambda x: maxent.score5(str(df.sequence[x][int(df.intronstart[x])-3:int(df.intronstart[x])+6])))
    mxnt['maxent3']=df.index.map(lambda x: maxent.score3(str(df.sequence[x][int(df.intronend[x])-20:int(df.intronend[x])+3])))
    return mxnt

def get_secstruct(df):
    sec=pd.DataFrame(index=df.index)
    sec['exon1']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.intronstart[x])-24:int(df.intronstart[x])]))[1])
    sec['donor']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.intronstart[x])-12:int(df.intronstart[x])+12]))[1])
    sec['intron5']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.intronstart[x]):int(df.intronstart[x])+24]))[1])
    sec['introncenter']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.intronstart[x])+(int(df.intronend[x])-int(df.intronstart[x]))/2-12:int(df.intronstart[x])+(int(df.intronend[x])-int(df.intronstart[x]))/2+12]))[1])
    sec['intron3']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.intronend[x])-24:int(df.intronend[x])]))[1])
    sec['acceptor']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.intronend[x])-12:int(df.intronend[x])+12]))[1])
    sec['exon2']=df.index.map(lambda x: RNA.fold(str(df.sequence[x][int(df.intronend[x]):int(df.intronend[x])+24]))[1])
    return sec

def get_sixmers(df):
    sixmers=[]
    for i in itertools.product(['G','C','A','T'], repeat=6):
        sixmers.append(''.join(list(i)))
    
    sixmerfeatures=pd.DataFrame(index=df.index)
    for i in sixmers:
        sixmerfeatures[i + '_exon1'] = df.index.map(lambda x: df.sequence[x][:df.intronstart[x]].count(i))
        sixmerfeatures[i + '_intron'] = df.index.map(lambda x: df.sequence[x][df.intronstart[x]:df.intronend[x]].count(i))
        sixmerfeatures[i + '_exon2'] = df.index.map(lambda x: df.sequence[x][int(df.intronend[x]):].count(i))
    return sixmerfeatures
    
def get_motifscores(df):      
    with open('../additional/ATtRACT/pwm_transposed.txt', 'r') as f:
        records=parse(f, 'jaspar')
    Xs_sel=pd.read_pickle('/net/mraid08/export/genie/Runs/Martin/ATtRACT/irmotifs_scores.pkl')
    mtfs=[]
    for i in Xs_sel.columns:
        if i.split('__')[0] not in mtfs:
            mtfs.append(i.split('__')[0])
    
    def find_motifs_ir(varid):
        motifs=pd.Series()
        for pos, seq in mm.counts.log_odds().search(Seq(df.sequence[varid], \
                        alphabet=IUPAC.IUPACUnambiguousDNA()), threshold=0, both=False):
            motifs.loc[pos]=seq
        motifs_up=motifs[motifs.index<df.intronstart[varid]]
        motifs_alt=motifs[(motifs.index>df.intronstart[varid])&(motifs.index<df.intronend[varid])]
        motifs_down=motifs[motifs.index>df.intronend[varid]]
        
        return list([motifs_up.sum(), motifs_alt.sum(), motifs_down.sum()])
    
    junirmotifs=pd.DataFrame(index=df.index)
    
    database=pd.read_table('../additional/ATtRACT/ATtRACT_db.txt')
    database.drop_duplicates('Matrix_id', inplace=True)
    database.set_index('Matrix_id', inplace=True)
    
    for mm in records:
#        if (mm.name in database[database.Organism=='Homo_sapiens'].index):
        if (mm.name in mtfs):
            mm.counts.__class__=matrix.PositionWeightMatrix
            junirmotifs['motifs']=df.index.map(lambda x: find_motifs_ir(x))
            junirmotifs[[str(mm.name) + '__score_motifs_up',str(mm.name) + '__score_motifs_alt',\
                          str(mm.name) +  '__score_motifs_down']]=junirmotifs.motifs.apply(lambda x: pd.Series(x))
            junirmotifs.drop('motifs',axis=1, inplace=True)
    
    return junirmotifs[Xs_sel.columns]


def make_features_ir(filename):
    '''
    input is a csv file with four columns: identifier, sequence, intronstart, intronend; 0-based, No header
    '''
    df=pd.read_csv(filename, header=None)
    df.columns=['identifier','sequence','intronstart','intronend']
    df['intronstart']=df['intronstart'].astype(int)
    df['intronend']=df['intronend'].astype(int)
    df.set_index('identifier', inplace=True)

    mxnt=get_maxent(df)
    sec=get_secstruct(df)
    motifs=get_motifscores(df)
    sixmers=get_sixmers(df)
    features=pd.concat([mxnt,sec,motifs,sixmers], axis=1)
    return features

def train_gbr_model(Xs, Ys, modelname):
    Xscols=pd.read_pickle('./dataframes/ml/noise/Xs_allfeatures_IR_sel_noiseresgam.pkl')
    Zs=pd.DataFrame(Xs[Xscols.columns].values, index=Xs.index, columns=Xscols.columns)
    Zs['measured']=pd.Series(Ys)
    Zs.dropna(inplace=True)
    gbr_model=GradientBoostingRegressor(learning_rate=0.05, n_estimators=100,\
                                        max_depth=6).fit(Zs.iloc[:,:-1],Zs.iloc[:,-1])    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'.sav','wb'))

def train_gbr_model_inclsplicingvalue(Xs, meansplicingvalues, Ys, modelname):
    Xscols=pd.read_pickle('./dataframes/ml/noise/Xs_allfeatures_IR_sel_noiseresgam.pkl')
    Zs=pd.DataFrame(Xs[Xscols.columns].values, index=Xs.index, columns=Xscols.columns)
    Zs=Zs.join(meansplicingvalues.loc[Zs.index])
    Zs['measured']=pd.Series(Ys)
    Zs.dropna(inplace=True)
    gbr_model=GradientBoostingRegressor(learning_rate=0.05, n_estimators=100,\
                                        max_depth=6).fit(Zs.iloc[:,:-1],Zs.iloc[:,-1])    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_inclsplicingvalue.sav','wb'))

def train_gbr_model_onlysplicingvalue(meansplicingvalues, Ys, modelname):
    Zs=meansplicingvalues.to_frame()
    Zs['measured']=pd.Series(Ys)
    Zs.dropna(inplace=True)
    gbr_model=GradientBoostingRegressor(learning_rate=0.05, n_estimators=100,\
                                        max_depth=2).fit(Zs.iloc[:,:-1],Zs.iloc[:,-1])    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_onlysplicingvalue.sav','wb'))

def train_gbr_model_maxent(Xs, Ys, modelname):
    cols=['maxent5','maxent3']
    Zs=pd.DataFrame(Xs[cols].values, index=Xs.index, columns=cols)
    Zs['measured']=pd.Series(Ys)
    Zs.dropna(inplace=True)
    gbr_model=GradientBoostingRegressor(learning_rate=0.01, n_estimators=100,\
                                        max_depth=5).fit(Zs.iloc[:,:-1],Zs.iloc[:,-1])    
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_maxent.sav','wb'))

def train_gbr_model_secondary(Xs, Ys, modelname):
    cols=['exon1','donor','intron5','introncenter','intron3','acceptor','exon2']
    Zs=pd.DataFrame(Xs[cols].values, index=Xs.index, columns=cols)
    Zs['measured']=pd.Series(Ys)
    Zs.dropna(inplace=True)
    gbr_model=GradientBoostingRegressor(learning_rate=0.01, n_estimators=300,\
                                        max_depth=4).fit(Zs.iloc[:,:-1],Zs.iloc[:,-1])    
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_secondary.sav','wb'))

def train_gbr_model_hexamers(Xs, Ys, modelname):
    Xscols=pd.read_pickle('./dataframes/ml/noise/Xs_hexamers_IR_motifs_sel_noiseresgam.pkl')
    Zs=pd.DataFrame(Xs[Xscols.columns].values, index=Xs.index, columns=Xscols.columns)
    Zs['measured']=pd.Series(Ys)
    Zs.dropna(inplace=True)
    gbr_model=GradientBoostingRegressor(learning_rate=0.01, n_estimators=100,\
                                        max_depth=6).fit(Zs.iloc[:,:-1],Zs.iloc[:,-1])    
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_hexamers.sav','wb'))

def train_gbr_model_motifs(Xs, Ys, modelname):
    Xscols=pd.read_pickle('./dataframes/ml/noise/Xs_motifs_IR_motifs_sel_noiseresgam.pkl')
    Zs=pd.DataFrame(Xs[Xscols.columns].values, index=Xs.index, columns=Xscols.columns)
    Zs['measured']=pd.Series(Ys)
    Zs.dropna(inplace=True)
    gbr_model=GradientBoostingRegressor(learning_rate=0.01, n_estimators=300,\
                                        max_depth=5).fit(Zs.iloc[:,:-1],Zs.iloc[:,-1])    
    
    pickle.dump(gbr_model, open('./ml/models/'+modelname+'_motifs.sav','wb'))


def train_all_models(Xs, meansplicingvalues, Ys, modelname):
    '''
    trains and saves models separately based on maxentscan, 
    secondary structure, hexamer counts and motif scores, 
    as well as a combined model
    '''
    train_gbr_model_maxent(Xs, Ys, modelname)
    train_gbr_model_secondary(Xs, Ys, modelname)
    train_gbr_model_hexamers(Xs, Ys, modelname)
    train_gbr_model_motifs(Xs, Ys, modelname)
    train_gbr_model(Xs, Ys, modelname)
    train_gbr_model_onlysplicingvalue(meansplicingvalues, Ys, modelname)
    train_gbr_model_inclsplicingvalue(Xs, meansplicingvalues, Ys, modelname)


def make_prediction_all_models(Xs, meansplicingvalues, Ys, mode='noiseres_gam', filename=None):
    '''
    possible modes:
        noiseres_gam (default)
        [str] this model is loaded and used for prediction
    saves image to file if filename is provided
    '''
    predictionscores=pd.DataFrame(columns=['r2score','pearsonr'])
    Xshex=pd.read_pickle('./dataframes/ml/noise/Xs_hexamers_IR_motifs_sel_noiseresgam.pkl')
    Xsmot=pd.read_pickle('./dataframes/ml/noise/Xs_motifs_IR_motifs_sel_noiseresgam.pkl')
    Xsall_sel=pd.read_pickle('./dataframes/ml/noise/Xs_allfeatures_IR_sel_noiseresgam.pkl')

    cols={'maxent':['maxent5','maxent3'],
          'secondary':['exon1','donor','intron5','introncenter','intron3','acceptor','exon2'],
          'hexamers':Xshex.columns,
          'motifs':Xsmot.columns,
                    'all features':':'}
    for feat in ['maxent','secondary','hexamers','motifs']:
        if mode=='noiseres_gam':
            mdl=pickle.load(open('./ml/models/irmodel_noiseres_gam_trained_'+feat+'.sav', 'rb'))
            Zs=pd.DataFrame(Xs[cols[feat]].values, index=Xs.index, columns=cols[feat])
            Zs['measured']=pd.Series(Ys)
            Zs.dropna(inplace=True)
        else:
            mdl=pickle.load(open('./ml/models/'+mode+'_'+feat+'.sav', 'rb'))
            Zs=pd.DataFrame(Xs[cols[feat]].values, index=Xs.index, columns=cols[feat])
            Zs['measured']=pd.Series(Ys)
            Zs.dropna(inplace=True)
        ypred=mdl.predict(Zs.iloc[:,:-1])
        predictionscores.loc[feat,'r2score']= metrics.r2_score(Zs.iloc[:,-1], ypred)
        predictionscores.loc[feat,'pearsonr']= pearsonr(Zs.iloc[:,-1], ypred)[0]
        f=plt.figure(figsize=(3,3))
        plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
        plt.title('$\mathregular{R^2}$ score=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
                  '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=14)
        if filename!=None:
            f.savefig('./ml/plots/scatter_'+filename +'_'+feat+'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    if mode=='noiseres_gam':
        mdl=pickle.load(open('./ml/models/irmodel_noiseres_gam_trained.sav', 'rb'))
        Zs=pd.DataFrame(Xs[Xsall_sel.columns].values, index=Xs.index, columns=Xsall_sel.columns)
        Zs['measured']=pd.Series(Ys)
        Zs.dropna(inplace=True)
    else:
        mdl=pickle.load(open('./ml/models/'+mode+'.sav', 'rb'))
        Zs=pd.DataFrame(Xs[Xsall_sel.columns].values, index=Xs.index, columns=Xsall_sel.columns)
        Zs['measured']=pd.Series(Ys)
        Zs.dropna(inplace=True)
    ypred=mdl.predict(Zs.iloc[:,:-1])
    predictionscores.loc['all features','r2score']= metrics.r2_score(Zs.iloc[:,-1], ypred)
    predictionscores.loc['all features','pearsonr']= pearsonr(Zs.iloc[:,-1], ypred)[0]
    f=plt.figure(figsize=(3,3))
    plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
    plt.title('$\mathregular{R^2}$ score=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
              '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=14)
    if filename!=None:
        f.savefig('./ml/plots/scatter_'+filename +'_all_features.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

    if mode=='noiseres_gam':
        mdl=pickle.load(open('./ml/models/irmodel_noiseres_gam_trained_onlysplicingvalue.sav', 'rb'))
        Zs=meansplicingvalues.to_frame()
        Zs['measured']=pd.Series(Ys)
        Zs.dropna(inplace=True)
    else:
        mdl=pickle.load(open('./ml/models/'+mode+'_onlysplicingvalue.sav', 'rb'))
        Zs=meansplicingvalues.to_frame()
        Zs['measured']=pd.Series(Ys)
        Zs.dropna(inplace=True)
    ypred=mdl.predict(Zs.iloc[:,:-1])
    predictionscores.loc['only splicing value','r2score']= metrics.r2_score(Zs.iloc[:,-1], ypred)
    predictionscores.loc['only splicing value','pearsonr']= pearsonr(Zs.iloc[:,-1], ypred)[0]
    f=plt.figure(figsize=(3,3))
    plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
    plt.title('$\mathregular{R^2}$ score=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
              '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=14)
    if filename!=None:
        f.savefig('./ml/plots/scatter_'+filename +'_only_splicingvalues.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    if mode=='noiseres_gam':
        mdl=pickle.load(open('./ml/models/irmodel_noiseres_gam_trained_inclsplicingvalue.sav', 'rb'))
        Zs=pd.DataFrame(Xs[Xsall_sel.columns].values, index=Xs.index, columns=Xsall_sel.columns)
        Zs=Zs.join(meansplicingvalues.loc[Zs.index])
        Zs['measured']=pd.Series(Ys)
        Zs.dropna(inplace=True)
    else:
        mdl=pickle.load(open('./ml/models/'+mode+'_inclsplicingvalue.sav', 'rb'))
        Zs=pd.DataFrame(Xs[Xsall_sel.columns].values, index=Xs.index, columns=Xsall_sel.columns)
        Zs=Zs.join(meansplicingvalues.loc[Zs.index])
        Zs['measured']=pd.Series(Ys)
        Zs.dropna(inplace=True)
    ypred=mdl.predict(Zs.iloc[:,:-1])
    predictionscores.loc['incl. splicing value','r2score']= metrics.r2_score(Zs.iloc[:,-1], ypred)
    predictionscores.loc['incl. splicing value','pearsonr']= pearsonr(Zs.iloc[:,-1], ypred)[0]
    f=plt.figure(figsize=(3,3))
    plt.scatter(Zs.iloc[:,-1], ypred, alpha=0.1)
    plt.title('$\mathregular{R^2}$ score=' + '{:.3f}'.format(metrics.r2_score(Zs.iloc[:,-1], ypred)) + '\nPearson r='+ \
              '{:.3f}'.format(pearsonr(Zs.iloc[:,-1], ypred)[0]), fontsize=14)
    if filename!=None:
        f.savefig('./ml/plots/scatter_'+filename +'_all_features_incl_splicingvalues.png', \
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
    plt.ylim(0,np.ceil(predictionscores.r2score.max()*10)/10.0)
    plt.ylabel('$\mathregular{R^2}$ score')
    if filename!=None:
        f.savefig('./ml/plots/'+filename+'_modelbuilding_R2scores_summary.png', \
                  dpi = 300, format='png', bbox_inches='tight', frameon=True)
        f.show()
        