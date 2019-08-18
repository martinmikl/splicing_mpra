#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 16:49:20 2017

@author: martinm
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# folder
# 
def send_cross_val_score_toq(folder, Xs, Ys, bestlr, bestnest, bestdepth):
    Xsfile=folder + 'Xs.pkl'
    Ysfile=folder + 'Ys.pkl'
    Xs.to_pickle(Xsfile)
    Ys.to_pickle(Ysfile)

    with open('/net/mraid08/export/genie/Runs/Martin/forml/runsth','w') as f:
        f.write('python2.7 /net/mraid08/export/genie/Runs/Martin/forml/testml.py ' \
        + folder + ' ' + Xsfile + ' ' + Ysfile + ' ' + str(bestlr) + ' ' +\
        str(bestnest) + ' ' + str(bestdepth))

    os.system('qp.pl -f /net/mraid08/export/genie/Runs/Martin/forml/runsth -max_u 100 -q himem7.q ')

def send_cvpredict_featimp_toq(folder, Xs, Ys, prefix, bestlr, bestnest, bestdepth):
    Xsfile=folder + 'Xs.pkl'
    Ysfile=folder + 'Ys.pkl'
    Xs.to_pickle(Xsfile)
    Ys.to_pickle(Ysfile)

    with open('/net/mraid08/export/genie/Runs/Martin/forml/runsth','w') as f:
        f.write('~/anaconda2/bin/python2.7 /net/mraid08/export/genie/Runs/Martin/forml/cvpredict_featimp.py ' \
        + folder + ' ' + Xsfile + ' ' + Ysfile + ' ' + prefix + ' ' + str(bestlr) + ' ' +\
        str(bestnest) + ' ' + str(bestdepth))

    os.system('qp.pl -f /net/mraid08/export/genie/Runs/Martin/forml/runsth -max_u 100 -q himem7.q ')

def send_cvpredict_featimp_sixmers_toq(folder, Xs, Ys, prefix, bestlr, bestnest, bestdepth):
    Xsfile=folder + 'Xs_'+ prefix + '.pkl'
    Ysfile=folder + 'Ys_'+ prefix + '.pkl'
    Xs.to_pickle(Xsfile)
    Ys.to_pickle(Ysfile)

    with open('/net/mraid08/export/genie/Runs/Martin/forml/runsth','w') as f:
        f.write('~/anaconda2/bin/python2.7 /net/mraid08/export/genie/Runs/Martin/forml/cvpredict_featimp_sixmers.py ' \
        + folder + ' ' + Xsfile + ' ' + Ysfile + ' ' + prefix + ' ' + str(bestlr) + ' ' +\
        str(bestnest) + ' ' + str(bestdepth))

    os.system('qp.pl -f /net/mraid08/export/genie/Runs/Martin/forml/runsth -max_u 100 -q himem7.q ')

def send_cvpredict_featimp_motifs_toq(folder, Xs, Ys, prefix, bestlr, bestnest, bestdepth):
    Xsfile=folder + 'Xs_'+ prefix + '.pkl'
    Ysfile=folder + 'Ys_'+ prefix + '.pkl'
    Xs.to_pickle(Xsfile)
    Ys.to_pickle(Ysfile)

    with open('/net/mraid08/export/genie/Runs/Martin/forml/runsth','w') as f:
        f.write('~/anaconda2/bin/python2.7 /net/mraid08/export/genie/Runs/Martin/forml/cvpredict_featimp_motifs.py ' \
        + folder + ' ' + Xsfile + ' ' + Ysfile + ' ' + prefix + ' ' + str(bestlr) + ' ' +\
        str(bestnest) + ' ' + str(bestdepth))

    os.system('qp.pl -f /net/mraid08/export/genie/Runs/Martin/forml/runsth -max_u 100 -q himem7.q ')

def send_permutation_test_gbr_toq(Xs, Ys, cv, n_permutations, outfile):
    Xsfile=folder + 'Xs.pkl'
    Ysfile=folder + 'Ys.pkl'
    Xs.to_pickle(Xsfile)
    Ys.to_pickle(Ysfile)

    with open('/net/mraid08/export/genie/Runs/Martin/forml/runsth','w') as f:
        f.write('python2.7 /net/mraid08/export/genie/Runs/Martin/forml/permutation_test_gbr.py ' \
        + ' ' + Xsfile + ' ' + Ysfile + ' ' + str(cv) + ' ' +\
        str(n_permutations) + ' ' + str(outfile))

    os.system('qp.pl -f /net/mraid08/export/genie/Runs/Martin/forml/runsth -max_u 100 -q himem7.q ')

def send_permutation_test_randomforest_toq(Xs, Ys, cv, n_permutations, outfile):
    Xsfile=folder + 'Xs.pkl'
    Ysfile=folder + 'Ys.pkl'
    Xs.to_pickle(Xsfile)
    Ys.to_pickle(Ysfile)

    with open('/net/mraid08/export/genie/Runs/Martin/forml/runsth','w') as f:
        f.write('python2.7 /net/mraid08/export/genie/Runs/Martin/forml/permutation_test_randomforest.py ' \
        + ' ' + Xsfile + ' ' + Ysfile + ' ' + str(cv) + ' ' +\
        str(n_permutations) + ' ' + str(outfile))

    os.system('qp.pl -f /net/mraid08/export/genie/Runs/Martin/forml/runsth -max_u 100 -q himem7.q ')


def send_optimize_parameters_gbr_toq(folder, Xs, Ys, prefix):
    Xsfile=folder + 'Xs.pkl'
    Ysfile=folder + 'Ys.pkl'
    Xs.to_pickle(Xsfile)
    Ys.to_pickle(Ysfile)

    with open('/net/mraid08/export/genie/Runs/Martin/forml/runsth','w') as f:
        f.write('~/anaconda2/bin/python2.7 /net/mraid08/export/genie/Runs/Martin/forml/optimize_parameters_gbr.py ' \
        + folder + ' ' + Xsfile + ' ' + Ysfile + ' ' + prefix)

    os.system('qp.pl -f /net/mraid08/export/genie/Runs/Martin/forml/runsth -max_u 100 -q himem7.q ')

def send_optimize_parameters_xgb_toq(folder, Xs, Ys, prefix):
    Xsfile=folder + 'Xs.pkl'
    Ysfile=folder + 'Ys.pkl'
    Xs.to_pickle(Xsfile)
    Ys.to_pickle(Ysfile)

    with open('/net/mraid08/export/genie/Runs/Martin/forml/runsth','w') as f:
        f.write('/usr/wisdom/python/bin/python /net/mraid08/export/genie/Runs/Martin/forml/optimize_parameters_xgb.py ' \
        + folder + ' ' + Xsfile + ' ' + Ysfile + ' ' + prefix)

    os.system('qp.pl -f /net/mraid08/export/genie/Runs/Martin/forml/runsth -max_u 100 -q himem7.q ')

def send_optimize_parameters_xgb_toqpy(folder, Xs, Ys, prefix):
    Xsfile=folder + 'Xs.pkl'
    Ysfile=folder + 'Ys.pkl'
    Xs.to_pickle(Xsfile)
    Ys.to_pickle(Ysfile)

    with open('/net/mraid08/export/genie/Runs/Martin/forml/runsth','w') as f:
        f.write('/usr/wisdom/python/bin/python /net/mraid08/export/genie/Runs/Martin/forml/optimize_parameters_xgb.py ' \
        + folder + ' ' + Xsfile + ' ' + Ysfile + ' ' + prefix)

    os.system('qp.py -q himem7.q -trds_def=2 /net/mraid08/export/genie/Runs/Martin/forml/runsth')


def send_optimize_parameters_gbr_noise_toq(folder, Xs, Ys, prefix):
    Xsfile=folder + 'Xs.pkl'
    Ysfile=folder + 'Ys.pkl'
    Xs.to_pickle(Xsfile)
    Ys.to_pickle(Ysfile)

    with open('/net/mraid08/export/genie/Runs/Martin/forml/runsth','w') as f:
        f.write('~/anaconda2/bin/python2.7 /net/mraid08/export/genie/Runs/Martin/forml/optimize_parameters_gbr_noise.py ' \
        + folder + ' ' + Xsfile + ' ' + Ysfile + ' ' + prefix)

    os.system('qp.pl -f /net/mraid08/export/genie/Runs/Martin/forml/runsth -max_u 100 -q himem7.q ')

def send_optimize_parameters_randomforest_toq(folder, Xs, Ys, prefix):
    Xsfile=folder + 'Xs.pkl'
    Ysfile=folder + 'Ys.pkl'
    Xs.to_pickle(Xsfile)
    Ys.to_pickle(Ysfile)

    with open('/net/mraid08/export/genie/Runs/Martin/forml/runsth','w') as f:
        f.write('python2.7 /net/mraid08/export/genie/Runs/Martin/forml/optimize_parameters_randomforest.py ' \
        + folder + ' ' + Xsfile + ' ' + Ysfile + ' ' + prefix)

    os.system('qp.pl -f /net/mraid08/export/genie/Runs/Martin/forml/runsth -max_u 100 -q himem7.q ')


def process_optimization_file_gbr(filepath):
    import math
    cvscs=pd.read_csv(filepath, header=None)
    cvscs.columns=['settings','r2']
    sett=cvscs.settings.apply(lambda x: pd.Series(x.split('_'), index=['lr','nest','depth']))
    cvscs=cvscs.join(sett)
    cvscs['lr']= cvscs['lr'].astype(float)
    cvscs['r2']=cvscs['r2'].astype(float)
    cvscs['nest']=cvscs['nest'].astype(int)
    cvscs['depth']=cvscs['depth'].astype(int)
    f=sns.lmplot(data=cvscs, x='lr', y='r2', hue='nest', col='depth',fit_reg=False, size=4)
    f.savefig(filepath[:-4] + '.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

def process_optimization_file_randomforest(filepath):
    cvscs=pd.read_csv(filepath, header=None)
    cvscs.columns=['settings','r2']
    sett=cvscs.settings.apply(lambda x: pd.Series(x.split('_'), index=['nest','depth']))
    cvscs=cvscs.join(sett)
    cvscs['r2']=cvscs['r2'].astype(float)
    cvscs['nest']=cvscs['nest'].astype(int)
    cvscs['depth']=cvscs['depth'].astype(int)
    f=sns.lmplot(data=cvscs, x='depth', y='r2', hue='nest',fit_reg=False, size=4)
    plt.ylim(0,cvscs.r2.max()*1.1)
    f.savefig(filepath[:-4]  + '.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

def make_binary_features(df, col):
    dfout=pd.DataFrame(False,index=df.index, columns=df[col].unique())
    for i in df[df[col].isnull()==False].index:
        dfout.loc[i,df.loc[i, col]]=True
    return dfout
    



'''
send_cross_val_score_toq(folderml, Xs, Ys, bestlr, bestnest, bestdepth)
send_optimize_parameters_gbr_toq(folderml, Xs, Ys, 'fourmer_')
send_optimize_parameters_gbr_toq(folderml, Xs, Ys, 'paired_')
send_optimize_parameters_randomforest_toq(folderml, Xs, Ys, 'fourmer_')
send_optimize_parameters_randomforest_toq(folderml, Xs, Ys, 'paired_')


# test different learning rates - with cv

cvpredscores=pd.Series()
for a in [0.01,0.1,0.2,0.3, 0.4]:
    clf=GradientBoostingRegressor(learning_rate=a)
    cvpredgb = cross_val_predict(clf, Xs, \
        Ys, cv=10)
    cvpredscores.loc[a]=metrics.r2_score(Ys,cvpredgb)
    print(a)
    print('cv test score: ' '{:.4f}'.format(cvpredscores.loc[a]))

f=plt.figure(figsize=(4,4))
plt.plot(cvpredscores)
plt.xlabel('learning_rate')
plt.ylabel('r2')
f.savefig(folder + 'XXX_YYY_gbr_learningrateoptimization_10fcv.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

bestlr=

# test different n_estimators - with cv

cvpredscores=pd.Series()
for a in [10,50,100,150,200,250,300,500,700]:
    clf=GradientBoostingRegressor(learning_rate=bestlr, n_estimators=a)
    cvpredgb = cross_val_predict(clf, Xs, \
        Ys, cv=10)
    cvpredscores.loc[a]=metrics.r2_score(Ys,cvpredgb)
    print(a)
    print('cv test score: ' '{:.4f}'.format(cvpredscores.loc[a]))

f=plt.figure(figsize=(4,4))
plt.plot(cvpredscores)
plt.xlabel('n_estimators')
plt.ylabel('r2')
f.savefig(folder + 'XXX_YYY_gbr_nestimatorsoptimization_10fcv_learningrate0_35.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

bestnest=500

# test different max_depth - with cv

cvpredscores=pd.Series()
for a in range(2,8):
    clf=GradientBoostingRegressor(learning_rate=bestlr, n_estimators=bestnest, max_depth=a)
    cvpredgb = cross_val_predict(clf, Xs, \
        Ys, cv=10)
    cvpredscores.loc[a]=metrics.r2_score(Ys,cvpredgb)
    print(a)
    print('cv test score: ' '{:.4f}'.format(cvpredscores.loc[a]))

f=plt.figure(figsize=(4,4))
plt.plot(cvpredscores)
plt.xlabel('max_depth')
plt.ylabel('r2')
f.savefig(folder + 'XXX_YYY_gbr_maxdepthoptimization_10fcv_learningrate0_35_nestimator500.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

bestdepth=3

 
### best found: 
clf=GradientBoostingRegressor(learning_rate=bestlr,n_estimators=bestnest, max_depth=bestdepth)
cvscores = cross_val_score(clf, Xs, \
        Ys, cv=10)

cvpredgb = cross_val_predict(clf, Xs, \
        Ys, cv=10)

f=plt.figure(figsize=(6,6))
plt.scatter(Ys,cvpredgb,s=3)
plt.xlabel('splicing ratio - measured')
plt.ylabel('splicing ratio - predicted')
plt.ylim(-15,15)
f.savefig(folder + 'XXX_YYY_gbr_10fcv_nestimators500_learningrate0_35.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

metrics.r2_score(Ys,cvpredgb)

r2=0.71974086880468757


 with cv
[ 0.71633977,  0.74690238,  0.65698592,  0.75243669,  0.75791617,
        0.75600233,  0.76188335,  0.64902551,  0.69886203,  0.67043018]

mean= 0.71667843195964542

best for noise: 0.035 (training set 0.15)


# get cv feature importances

kf=KFold(len(Xs),n_folds=10)

featimp=pd.DataFrame()
for train_index, test_index in kf:
    clf=GradientBoostingRegressor(learning_rate=0.35, n_estimators=500, max_depth=3).fit(Xs.iloc[train_index], Ys.iloc[train_index])
    print(clf.score(Xs.iloc[test_index], Ys.iloc[test_index]))
    featimp=featimp.append(pd.Series(clf.feature_importances_),ignore_index=True)

featimp.columns=Xs.columns
featimp.drop('meanfeatimp', inplace=True)
meanfeatimp=pd.DataFrame(pd.Series(featimp.mean(),index=Xs.columns, name='mean_feature_imp'))
meanfeatimp.sort_values(by='mean_feature_imp',ascending=False, inplace=True)
meanfeatimp['polypyrimidine']=meanfeatimp.index.map(lambda x: bool(('intron' in x) & (x[0:4].find('G')==-1) & (x[0:4].find('A')==-1)))
meanfeatimp['fourmer']=meanfeatimp.index.map(lambda x: x)


f=plt.figure(figsize=(6,4))
ax=sns.barplot(data=meanfeatimp,x='fourmer', y='mean_feature_imp',palette=meanfeatimp['polypyrimidine'].map({False:'Blue',True:'Red'}))
ax.set_xticklabels([])
plt.xlabel('4mers')
plt.ylabel('feature importance')
f.show()

#################
# with random forests
xtrain, xtest, ytrain, ytest = train_test_split(Xs, Ys,\
        test_size=0.2, random_state=0)


clf=RandomForestRegressor().fit(xtrain,ytrain)

clf.score(xtest,ytest)



# test different n_estimators - with cv

cvpredscores=pd.Series()
for a in [10,20,50,70,100]:
for a in [200,300,500,1000]:
    clf=RandomForestRegressor(n_estimators=a)
    cvpredgb = cross_val_predict(clf, Xs, \
        Ys, cv=10)
    cvpredscores.loc[a]=metrics.r2_score(Ys,cvpredgb)
    print(a)
    print('cv test score: ' '{:.4f}'.format(cvpredscores.loc[a]))

f=plt.figure(figsize=(4,4))
plt.plot(cvpredscores)
plt.xlabel('n_estimators')
plt.ylabel('r2')
f.savefig(folder + 'XXX_YYY_randomforest_nestimatorsoptimization_10fcv.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

bestnest=70

# test different max_depth - with cv

cvpredscores=pd.Series()
for a in range(11,20):
    clf=RandomForestRegressor(n_estimators=70, max_depth=a)
    cvpredgb = cross_val_predict(clf, Xs, \
        Ys, cv=10)
    cvpredscores.loc[a]=metrics.r2_score(Ys,cvpredgb)
    print(a)
    print('cv test score: ' '{:.4f}'.format(cvpredscores.loc[a]))

f=plt.figure(figsize=(4,4))
plt.plot(cvpredscores)
plt.xlabel('max_depth')
plt.ylabel('r2')
f.savefig(folder + 'XXX_YYY_randomforest_maxdepthoptimization_10fcv_nestimator70.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

bestdepth=16

 
### best found: learning_rate=0.35 (xxx for noise), n_estimators=500 (xxx for noise)
clf=RandomForestRegressor(n_estimators=bestnest, max_depth=bestdepth)
cvscores = cross_val_score(clf, Xs, \
        Ys, cv=10)

cvpredgb = cross_val_predict(clf, Xs, \
        Ys, cv=10)

f=plt.figure(figsize=(6,6))
plt.scatter(Ys,cvpredgb,s=3)
plt.xlabel('splicing ratio - measured')
plt.ylabel('splicing ratio - predicted')
plt.ylim(-15,15)
f.savefig(folder + 'XXX_YYY_randomforest_10fcv_nestimators70_maxdepth16.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

metrics.r2_score(Ys,cvpredgb)

r2=0.77484653586187013


 with cv

0.784747952011
0.807717424817
0.744824522631
0.793842345372
0.795546996832
0.81626911728
0.798915216216
0.701573237707
0.735225156469
0.766088852353



# get cv feature importances

kf=KFold(len(Xs),n_folds=10)

featimpRF=pd.DataFrame()
for train_index, test_index in kf:
    clf=RandomForestRegressor(n_estimators=bestnest, max_depth=bestdepth).fit(Xs.iloc[train_index], Ys.iloc[train_index])
    print(clf.score(Xs.iloc[test_index], Ys.iloc[test_index]))
    featimpRF=featimpRF.append(pd.Series(clf.feature_importances_),ignore_index=True)

featimpRF.columns=Xs.columns
meanfeatimpRF=pd.DataFrame(pd.Series(featimpRF.mean(),index=Xs.columns, name='mean_feature_imp'))
meanfeatimpRF.sort_values(by='mean_feature_imp',ascending=False, inplace=True)
meanfeatimpRF['polypyrimidine']=meanfeatimpRF.index.map(lambda x: bool(('intron' in x) & (x[0:4].find('G')==-1) & (x[0:4].find('A')==-1)))
meanfeatimpRF['fourmer']=meanfeatimpRF.index.map(lambda x: x)


f=plt.figure(figsize=(6,4))
ax=sns.barplot(data=meanfeatimpRF.iloc[0:50],x='fourmer', y='mean_feature_imp',palette=meanfeatimpRF['polypyrimidine'].map({False:'Blue',True:'Red'}))
ax.set_xticklabels([])
plt.xlabel('4mers')
plt.ylabel('feature importance')
f.savefig(folder + 'featureimp_XXX_YYY_randomforest_10fcv_nestimators70_maxdepth16.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

'''