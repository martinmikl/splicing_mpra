#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 15:22:47 2019

@author: martinm
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
#import shap ### TO USE SHAP REVERT TO REVISION 25 (conda install --revision 25), THEN CHANGE BACK TO REVISION 26
import pickle
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from sklearn.ensemble.partial_dependence import plot_partial_dependence
from sklearn.ensemble.partial_dependence import partial_dependence
from mpl_toolkits.mplot3d import Axes3D
sns.set_context('poster')
from sklearn.cross_validation import cross_val_score

database=pd.read_table('../additional/ATtRACT/ATtRACT_db.txt')
database.drop_duplicates('Matrix_id', inplace=True)
database.set_index('Matrix_id', inplace=True)



#%% Intron retention

Xs=pd.read_pickle('./dataframes/ml/ir/Xs_trainnew_ir.pkl')

cols={'maxent':['maxent5','maxent3'],
      'secondary':['exon1','donor','intron5','introncenter','intron3','acceptor','exon2'],
      'hexamers':[x for x in Xs.columns if (x.split('_')[-1] in ['exon1','intron','exon2'])&(len(x.split('_'))==2)],
      'motifs':[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                'all features':':'}


for feat in ['maxent','secondary','hexamers']:
    mdl=pickle.load(open('./ml/models/irmodel_logratio_trained_'+feat+'.sav', 'rb'))

    explainer = shap.TreeExplainer(mdl)
    shap_values=explainer.shap_values(Xs[cols[feat]])
    f=plt.figure()
    shap.summary_plot(shap_values, Xs[cols[feat]], feature_names=cols[feat])
    f.savefig('./ml/plots/shap/ir_'+feat+'_summaryplot_shapvalues.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/irmodel_logratio_trained_motifs.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols['motifs']])
f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+x.split('_')[-1] for x in cols['motifs']])
f.savefig('./ml/plots/shap/ir_motifs_summaryplot_shapvalues.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/irmodel_logratio_trained_motifs.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols['motifs']])
f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+x for x in cols['motifs']])
f.savefig('./ml/plots/shap/ir_motifs_summaryplot_shapvalues_inclattractID.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+database.loc[x.split('__')[0],'Motif'] + ' ' + x.split('_')[-1] for x in cols['motifs']])
f.savefig('./ml/plots/shap/ir_motifs_summaryplot_shapvalues_inclmotifs.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/irmodel_logratio_trained.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs)
f=plt.figure()
shap.summary_plot(shap_values, Xs, feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+database.loc[x.split('__')[0],'Motif']+' '+x.split('_')[-1] if 'score_motifs' in x else x for x in Xs.columns])
f.savefig('./ml/plots/shap/ir_allfeatures_summaryplot_shapvalues.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
#######################    
'''
feat='maxent'
shap.dependence_plot('maxent5',shap_values,Xs[cols[feat]])
shap.dependence_plot('maxent3',shap_values,Xs[cols[feat]])

shap_interaction_values=explainer.shap_interaction_values(Xs[cols[feat]])
shap.summary_plot(shap_interaction_values, Xs[cols[feat]])

shap.dependence_plot(('maxent5','maxent3'), shap_interaction_values,Xs[cols[feat]])

shap.dependence_plot('exon1',shap_values,Xs[cols[feat]])

shap.dependence_plot('maxent3',shap_values,Xs[cols[feat]])

shap_interaction_values=explainer.shap_interaction_values(Xs[cols[feat]])
shap.summary_plot(shap_interaction_values, Xs[cols[feat]])

shap.dependence_plot(('acceptor','exon2'), shap_interaction_values,Xs[cols[feat]])

shap.force_plot(explainer.expected_value, shap_values[0,:], Xs[cols[feat]].iloc[0,:])


feat='motifs'
mdl=pickle.load(open('./ml/models/irmodel_logratio_trained_'+feat+'.sav', 'rb'))

explainer = shap.TreeExplainer(mdl)
shap_values_motifs=explainer.shap_values(Xs[cols[feat]])
shap_interaction_values_motifs=explainer.shap_interaction_values(Xs[cols[feat]])
'''

#####

mdl=pickle.load(open('./ml/models/irmodel_logratio_trained.sav', 'rb'))

#maxent3second vs alt3/acceptor2
plot_partial_dependence(mdl, Xs, [(0,2),(0,3)], grid_resolution=100)
plot_partial_dependence(mdl, Xs, [0,2,3], grid_resolution=100)

#maxent3first vs intron/acceptor1
plot_partial_dependence(mdl, Xs, [(1,6),(1,7)], grid_resolution=150)
plot_partial_dependence(mdl, Xs, [1,6,7,8], grid_resolution=100)

plot_partial_dependence(mdl, Xs, [0,1,(0,1)], grid_resolution=150)
plot_partial_dependence(mdl, Xs, [3,7,(3,7)], grid_resolution=150)
plot_partial_dependence(mdl, Xs, [(0,1),(3,7)], grid_resolution=200)
plot_partial_dependence(mdl, Xs, [(0,1),(3,7)], grid_resolution=100)

#maxent3first vs acceptor1/maxent3second vs acceptor2

f, axs = plot_partial_dependence(mdl, Xs, [(0,3),(1,7)], grid_resolution=200)
f.savefig('./ml/plots/ir_partialdependence_maxentvssecondary_acceptordonor.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f, axs = plot_partial_dependence(mdl, Xs, [(0,1),(3,7)], grid_resolution=100)
f.savefig('./ml/plots/ir_partialdependence_maxentdonorvsacceptor_secondarydonorvsacceptor.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



f=plt.figure()
target_feature=(0,1)
pdp, axes = partial_dependence(mdl, target_feature, X=Xs, grid_resolution=100)
XX, YY = np.meshgrid(axes[0], axes[1])
Z = pdp[0].reshape(list(map(np.size, axes))).T
ax=Axes3D(fig)
surf=ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=plt.cm.BuPu, edgecolor='k')
f.savefig('./ml/plots/ir_3D_partialdependence_maxentfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure()
target_feature=(3,7)
pdp, axes = partial_dependence(mdl, target_feature, X=Xs, grid_resolution=100)
XX, YY = np.meshgrid(axes[0], axes[1])
Z = pdp[0].reshape(list(map(np.size, axes))).T
ax=Axes3D(fig)
surf=ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=plt.cm.BuPu, edgecolor='k')
f.savefig('./ml/plots/ir_3D_partialdependence_secondaryfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



pdp_1, axes_1 = partial_dependence(mdl, 0, X=Xs, grid_resolution=100)
pdp_2, axes_2 = partial_dependence(mdl, 1, X=Xs, grid_resolution=100)

f=plt.figure(figsize=(4,3))
plt.plot(axes_1[0],pdp_1[0])
plt.plot(axes_2[0],pdp_2[0])
plt.axhline(y=0,linewidth=2, color='gray', alpha=0.5)
f.savefig('./ml/plots/ir_partialdependence_maxentfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

pdp_1, axes_1 = partial_dependence(mdl, 3, X=Xs, grid_resolution=100)
pdp_2, axes_2 = partial_dependence(mdl, 7, X=Xs, grid_resolution=100)

f=plt.figure(figsize=(4,3))
plt.plot(axes_1[0],pdp_1[0])
plt.plot(axes_2[0],pdp_2[0])
plt.axhline(y=0,linewidth=2, color='gray', alpha=0.5)
f.savefig('./ml/plots/ir_partialdependence_secondaryfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



#####

mdl=pickle.load(open('./ml/models/irmodel_logratio_trained_maxent.sav', 'rb'))


#maxent3first vs acceptor1/maxent3second vs acceptor2



f=plt.figure()
target_feature=(0,1)
pdp, axes = partial_dependence(mdl, target_feature, X=Xs, grid_resolution=20)
XX, YY = np.meshgrid(axes[0], axes[1])
Z = pdp[0].reshape(list(map(np.size, axes))).T
ax=Axes3D(f)
surf=ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=plt.cm.BuPu, edgecolor='k')
f.savefig('./ml/plots/ir_3D_partialdependence_maxentfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



pdp_1, axes_1 = partial_dependence(mdl, 0, X=Xs, grid_resolution=100)
pdp_2, axes_2 = partial_dependence(mdl, 1, X=Xs, grid_resolution=100)

f=plt.figure(figsize=(4,3))
plt.plot(axes_1[0],pdp_1[0])
plt.plot(axes_2[0],pdp_2[0])
plt.axhline(y=0,linewidth=2, color='gray', alpha=0.5)
f.savefig('./ml/plots/ir_partialdependence_maxentfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


#%% Cassette exons

Xs=pd.read_pickle('./dataframes/ml/cas/Xs_train_cas.pkl')

cols={'maxent':['maxent5','maxent3'],
      'secondary':['intron1','acceptor','exon5','exoncenter','exon3','donor','intron2'],
      'hexamers':[x for x in Xs.columns if (x.split('_')[-1] in ['intron1','exon','intron2'])&(len(x.split('_'))==2)],
      'motifs':[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                'all features':':'}


for feat in ['maxent','secondary','hexamers']:
    mdl=pickle.load(open('./ml/models/casmodel_logratio_trained_'+feat+'.sav', 'rb'))

    explainer = shap.TreeExplainer(mdl)
    shap_values=explainer.shap_values(Xs[cols[feat]])
    
    print feat 
    f=plt.figure()
    shap.summary_plot(shap_values, Xs[cols[feat]], feature_names=cols[feat])
    f.savefig('./ml/plots/shap/cas_'+feat+'_summaryplot_shapvalues.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/casmodel_logratio_trained_motifs.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols['motifs']])
f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+x.split('_')[-1] for x in cols['motifs']])
f.savefig('./ml/plots/shap/cas_motifs_summaryplot_shapvalues.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/casmodel_logratio_trained_motifs.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols['motifs']])
f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+x for x in cols['motifs']])
f.savefig('./ml/plots/shap/cas_motifs_summaryplot_shapvalues_inclattractID.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/casmodel_logratio_trained_motifs.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols['motifs']])
f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+database.loc[x.split('__')[0],'Motif'] + ' ' + x.split('_')[-1] for x in cols['motifs']])
f.savefig('./ml/plots/shap/cas_motifs_summaryplot_shapvalues_inclmotifs.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/casmodel_logratio_trained.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs)
f=plt.figure()
shap.summary_plot(shap_values, Xs, feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+database.loc[x.split('__')[0],'Motif']+' '+x.split('_')[-1] if 'score_motifs' in x else x for x in Xs.columns])
f.savefig('./ml/plots/shap/cas_allfeatures_summaryplot_shapvalues.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


#####

mdl=pickle.load(open('./ml/models/casmodel_logratio_trained.sav', 'rb'))

#maxent3second vs alt3/acceptor2
plot_partial_dependence(mdl, Xs, [(1,2),(1,3)], grid_resolution=150)

#maxent3first vs intron/acceptor1
plot_partial_dependence(mdl, Xs, [(0,6),(0,7)], grid_resolution=150)

#maxent3first vs acceptor1/maxent3second vs acceptor2

f, axs = plot_partial_dependence(mdl, Xs, [(1,3),(0,7)], grid_resolution=150)
f.savefig('./ml/plots/cas_partialdependence_maxentvssecondary_acceptordonor.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f, axs = plot_partial_dependence(mdl, Xs, [(1,0),(3,7)], grid_resolution=100)
f.savefig('./ml/plots/cas_partialdependence_maxentdonorvsacceptor_secondarydonorvsacceptor.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


fig=plt.figure()
target_feature=(1,0)
pdp, axes = partial_dependence(mdl, target_feature, X=Xs, grid_resolution=100)
XX, YY = np.meshgrid(axes[0], axes[1])
Z = pdp[0].reshape(list(map(np.size, axes))).T
ax=Axes3D(fig)
surf=ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=plt.cm.BuPu, edgecolor='k')

fig=plt.figure()
target_feature=(3,7)
pdp, axes = partial_dependence(mdl, target_feature, X=Xs, grid_resolution=100)
XX, YY = np.meshgrid(axes[0], axes[1])
Z = pdp[0].reshape(list(map(np.size, axes))).T
ax=Axes3D(fig)
surf=ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=plt.cm.BuPu, edgecolor='k')



pdp_1, axes_1 = partial_dependence(mdl, 1, X=Xs, grid_resolution=100)
pdp_2, axes_2 = partial_dependence(mdl, 0, X=Xs, grid_resolution=100)

f=plt.figure(figsize=(4,3))
plt.plot(axes_1[0],pdp_1[0])
plt.plot(axes_2[0],pdp_2[0])
plt.axhline(y=0,linewidth=2, color='gray', alpha=0.5)
f.savefig('./ml/plots/cassette_partialdependence_maxentfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

pdp_1, axes_1 = partial_dependence(mdl, 3, X=Xs, grid_resolution=100)
pdp_2, axes_2 = partial_dependence(mdl, 7, X=Xs, grid_resolution=100)

f=plt.figure(figsize=(4,3))
plt.plot(axes_1[0],pdp_1[0])
plt.plot(axes_2[0],pdp_2[0])
plt.axhline(y=0,linewidth=2, color='gray', alpha=0.5)
f.savefig('./ml/plots/cassette_partialdependence_secondaryfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


fi=pd.DataFrame(index=Xs.columns)
fi['featimp']=mdl.feature_importances_
fi['factorname']=fi.index.map(lambda x: database.loc[x.split('__')[0],'Gene_name']\
+' '+x.split('_')[-1] if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2) else np.nan)

fisort=fi.sort_values('featimp', ascending=False)

fisort[(fisort.factorname=='SRSF1 alt')|(fisort.factorname=='HNRNPA1 alt')]
fisort[(fisort.factorname=='SRSF1 alt')|\
        ((['HNRNPA' in x for x in fisort.factorname.dropna()])&(['alt' in x for x in fisort.factorname.dropna()]))]

list(Xs.columns).index('M106_0.6__score_motifs_alt')
list(Xs.columns).index('1397__score_motifs_alt')

pdp_1, axes_1 = partial_dependence(mdl, list(Xs.columns).index('M106_0.6__score_motifs_alt'), X=Xs, grid_resolution=100)
pdp_2, axes_2 = partial_dependence(mdl, list(Xs.columns).index('1397__score_motifs_alt'), X=Xs, grid_resolution=100)

f=plt.figure(figsize=(4,3))
plt.plot(axes_1[0],pdp_1[0])
plt.plot(axes_2[0],pdp_2[0])
plt.axhline(y=0,linewidth=2, color='gray', alpha=0.5)
f.savefig('./ml/plots/cassette_partialdependence_secondaryfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


fig=plt.figure()
target_feature=(list(Xs.columns).index('1397__score_motifs_alt'),
                list(Xs.columns).index('s48__score_motifs_alt'))
pdp, axes = partial_dependence(mdl, target_feature, X=Xs, grid_resolution=100)
XX, YY = np.meshgrid(axes[0], axes[1])
Z = pdp[0].reshape(list(map(np.size, axes))).T
ax=Axes3D(fig)
surf=ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=plt.cm.BuPu, edgecolor='k')

plot_partial_dependence(mdl, Xs, target_feature,  grid_resolution=100)
plot_partial_dependence(mdl, Xs, (882,278,(882,278)),  grid_resolution=100)

#%% Alternative 5' splice sites

Xs=pd.read_pickle('./dataframes/ml/five/Xs_train_five.pkl')

cols={'maxent':['maxent5first','maxent5second'],
      'secondary':['exon','donor1','alt5','altcenter','alt3','donor2','intron'],
      'hexamers':[x for x in Xs.columns if (x.split('_')[-1] in ['intron','alt','exon'])&(len(x.split('_'))==2)],
      'motifs':[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                'all features':':'}


for feat in ['maxent','secondary','hexamers']:
    mdl=pickle.load(open('./ml/models/fivemodel_logratio_trained_'+feat+'.sav', 'rb'))

    explainer = shap.TreeExplainer(mdl)
    shap_values=explainer.shap_values(Xs[cols[feat]])
    
    print feat 
    f=plt.figure()
    shap.summary_plot(shap_values, Xs[cols[feat]], feature_names=cols[feat])
    f.savefig('./ml/plots/shap/five_'+feat+'_summaryplot_shapvalues.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/fivemodel_logratio_trained_motifs.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols['motifs']])
f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+x.split('_')[-1] for x in cols['motifs']])
f.savefig('./ml/plots/shap/five_motifs_summaryplot_shapvalues.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/threemodel_logratio_trained_motifs.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols['motifs']])
f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+x for x in cols['motifs']])
f.savefig('./ml/plots/shap/three_motifs_summaryplot_shapvalues_inclattractID.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/fivemodel_logratio_trained_motifs.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols['motifs']])
f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+database.loc[x.split('__')[0],'Motif'] + ' ' + x.split('_')[-1] for x in cols['motifs']])
f.savefig('./ml/plots/shap/five_motifs_summaryplot_shapvalues_inclmotifs.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/fivemodel_logratio_trained.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs)
f=plt.figure()
shap.summary_plot(shap_values, Xs, feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+database.loc[x.split('__')[0],'Motif']+' '+x.split('_')[-1] if 'score_motifs' in x else x for x in Xs.columns])
f.savefig('./ml/plots/shap/five_allfeatures_summaryplot_shapvalues.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


#####

mdl=pickle.load(open('./ml/models/fivemodel_logratio_trained.sav', 'rb'))

#maxent3second vs alt3/acceptor2
plot_partial_dependence(mdl, Xs, [(1,6),(1,7)], grid_resolution=100)

#maxent3first vs intron/acceptor1
plot_partial_dependence(mdl, Xs, [(0,2),(0,3)], grid_resolution=100)

#maxent3first vs acceptor1/maxent3second vs acceptor2

f, axs = plot_partial_dependence(mdl, Xs, [(0,3),(1,7)], grid_resolution=100)
f.savefig('./ml/plots/five_partialdependence_maxentvssecondary_firstandsecond.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f, axs = plot_partial_dependence(mdl, Xs, [(0,1),(3,7)], grid_resolution=100)
f.savefig('./ml/plots/five_partialdependence_maxentfirstvssecond_secondaryfirstvssecond.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f, axs = plot_partial_dependence(mdl, Xs, [0,1,(0,1)], grid_resolution=200)
f, axs = plot_partial_dependence(mdl, Xs, [3,7,(3,7)], grid_resolution=200)

###

pdp_1, axes_1 = partial_dependence(mdl, 0, X=Xs, grid_resolution=100)
pdp_2, axes_2 = partial_dependence(mdl, 1, X=Xs, grid_resolution=100)

f=plt.figure(figsize=(4,3))
plt.plot(axes_1[0],pdp_1[0])
plt.plot(axes_2[0],pdp_2[0])
plt.axhline(y=0,linewidth=2, color='gray', alpha=0.5)
f.savefig('./ml/plots/five_partialdependence_maxentfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)




mdl=pickle.load(open('./ml/models/fivemodel_logratio_trained_maxent.sav', 'rb'))


pdp_1, axes_1 = partial_dependence(mdl, 0, X=Xs, grid_resolution=100)
pdp_2, axes_2 = partial_dependence(mdl, 1, X=Xs, grid_resolution=100)

f=plt.figure(figsize=(4,3))
plt.plot(axes_1[0],pdp_1[0])
plt.plot(axes_2[0],pdp_2[0])
plt.axhline(y=0,linewidth=2, color='gray', alpha=0.5)
f.savefig('./ml/plots/five_partialdependence_maxentmodel_maxentfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


#%% Alternative 3' splice sites

Xs=pd.read_pickle('./dataframes/ml/three/Xs_train_three.pkl')

cols={'maxent':['maxent3first','maxent3second'],
      'secondary':['intron','acceptor1','alt5','altcenter','alt3','acceptor2','exon'],
      'hexamers':[x for x in Xs.columns if (x.split('_')[-1] in ['intron','alt','exon'])&(len(x.split('_'))==2)],
      'motifs':[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                'all features':':'}


for feat in ['maxent','secondary','hexamers']:
    mdl=pickle.load(open('./ml/models/threemodel_logratio_trained_'+feat+'.sav', 'rb'))

    explainer = shap.TreeExplainer(mdl)
    shap_values=explainer.shap_values(Xs[cols[feat]])
    
    print feat 
    f=plt.figure()
    shap.summary_plot(shap_values, Xs[cols[feat]], feature_names=cols[feat])
    f.savefig('./ml/plots/shap/three_'+feat+'_summaryplot_shapvalues.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/threemodel_logratio_trained_motifs.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols['motifs']])
f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+x.split('_')[-1] for x in cols['motifs']])
f.savefig('./ml/plots/shap/three_motifs_summaryplot_shapvalues.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/threemodel_logratio_trained_motifs.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols['motifs']])
f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+x for x in cols['motifs']])
f.savefig('./ml/plots/shap/three_motifs_summaryplot_shapvalues_inclattractID.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/threemodel_logratio_trained_motifs.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols['motifs']])
f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+database.loc[x.split('__')[0],'Motif'] + ' ' + x.split('_')[-1] for x in cols['motifs']])
f.savefig('./ml/plots/shap/three_motifs_summaryplot_shapvalues_inclmotifs.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/threemodel_logratio_trained.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs)
f=plt.figure()
shap.summary_plot(shap_values, Xs, feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+database.loc[x.split('__')[0],'Motif']+' '+x.split('_')[-1] if 'score_motifs' in x else x for x in Xs.columns])
f.savefig('./ml/plots/shap/three_allfeatures_summaryplot_shapvalues.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


####
'''
feat='maxent'
mdl=pickle.load(open('./ml/models/threemodel_logratio_trained_'+feat+'.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols[feat]])

shap.dependence_plot('maxent3first',shap_values,Xs[cols[feat]], show=False)
plt.xlim(-5,15)
plt.show()
shap.dependence_plot('maxent3second',shap_values,Xs[cols[feat]], show=False)
plt.xlim(-5,15)
plt.show()

shap_interaction_values=explainer.shap_interaction_values(Xs[cols[feat]])
shap.summary_plot(shap_interaction_values, Xs[cols[feat]])

shap.dependence_plot(('maxent3first','maxent3second'), shap_interaction_values,Xs[cols[feat]], show=False)
plt.xlim(-5,15)
plt.show()

feat='secondary'
mdl=pickle.load(open('./ml/models/threemodel_logratio_trained_'+feat+'.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols[feat]])

shap.dependence_plot('acceptor1',shap_values,Xs[cols[feat]])
shap.dependence_plot('alt3',shap_values,Xs[cols[feat]], alpha=0.2)
shap.dependence_plot('alt3',shap_values,Xs[cols[feat]], display_features=Xs['exon'])

shap_interaction_values=explainer.shap_interaction_values(Xs[cols[feat]])
shap.summary_plot(shap_interaction_values, Xs[cols[feat]])

shap.dependence_plot(('donor2','exon'), shap_interaction_values,Xs[cols[feat]])

shap.force_plot(explainer.expected_value, shap_values[0,:], Xs[cols[feat]].iloc[0,:], matplotlib=True)


from sklearn.ensemble.partial_dependence import plot_partial_dependence
from sklearn.ensemble.partial_dependence import partial_dependence
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d

fig=plt.figure()
target_feature=(0,1)
pdp, axes = partial_dependence(mdl, target_feature, X=Xs[cols[feat]], grid_resolution=50)
XX, YY = np.meshgrid(axes[0], axes[1])
Z = pdp[0].reshape(list(map(np.size, axes))).T
ax=Axes3D(fig)
surf=ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=plt.cm.BuPu, edgecolor='k')

fig, axs = plot_partial_dependence(mdl, Xs[cols[feat]], [3,4,(3,4)])
fig, axs = plot_partial_dependence(mdl, Xs[cols[feat]], [1,5,(1,5)])
fig, axs = plot_partial_dependence(mdl, Xs[cols[feat]], [(1,2),(1,3),(1,4),(1,5),(1,6)], grid_resolution=10)
plot_partial_dependence(mdl, Xs[cols[feat]], [(0,1)], grid_resolution=10)
plot_partial_dependence(mdl, Xs[cols[feat]], [(1,4),(5,4)], grid_resolution=20)
'''

######

mdl=pickle.load(open('./ml/models/threemodel_logratio_trained.sav', 'rb'))

#maxent3second vs alt3/acceptor2
plot_partial_dependence(mdl, Xs, [(1,6),(1,7)], grid_resolution=100)

#maxent3first vs intron/acceptor1
plot_partial_dependence(mdl, Xs, [(0,2),(0,3)], grid_resolution=100)

#maxent3first vs acceptor1/maxent3second vs acceptor2

f, axs = plot_partial_dependence(mdl, Xs, [(0,3),(1,7)], grid_resolution=150)
f.savefig('./ml/plots/three_partialdependence_maxentvssecondary_firstandsecond.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f, axs = plot_partial_dependence(mdl, Xs, [(0,1),(3,7)], grid_resolution=100)
f.savefig('./ml/plots/three_partialdependence_maxentfirstvssecond_secondaryfirstvssecond.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f, axs = plot_partial_dependence(mdl, Xs, [0,1,(0,1)], grid_resolution=100)
f, axs = plot_partial_dependence(mdl, Xs, [3,7,(3,7)], grid_resolution=100)


pdp_1, axes_1 = partial_dependence(mdl, 0, X=Xs, grid_resolution=100)
pdp_2, axes_2 = partial_dependence(mdl, 1, X=Xs, grid_resolution=100)

f=plt.figure(figsize=(4,3))
plt.plot(axes_1[0],pdp_1[0])
plt.plot(axes_2[0],pdp_2[0])
plt.axhline(y=0,linewidth=2, color='gray', alpha=0.5)
f.savefig('./ml/plots/three_partialdependence_maxentfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

pdp_1, axes_1 = partial_dependence(mdl, 3, X=Xs, grid_resolution=100)
pdp_2, axes_2 = partial_dependence(mdl, 7, X=Xs, grid_resolution=100)

f=plt.figure(figsize=(4,3))
plt.plot(axes_1[0],pdp_1[0])
plt.plot(axes_2[0],pdp_2[0])
plt.axhline(y=0,linewidth=2, color='gray', alpha=0.5)
f.savefig('./ml/plots/three_partialdependence_secondaryfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



mdl=pickle.load(open('./ml/models/threemodel_logratio_trained_maxent.sav', 'rb'))


pdp_1, axes_1 = partial_dependence(mdl, 0, X=Xs, grid_resolution=100)
pdp_2, axes_2 = partial_dependence(mdl, 1, X=Xs, grid_resolution=100)

f=plt.figure(figsize=(4,3))
plt.plot(axes_1[0],pdp_1[0])
plt.plot(axes_2[0],pdp_2[0])
plt.axhline(y=0,linewidth=2, color='gray', alpha=0.5)
plt.xlim(0,12)
f.savefig('./ml/plots/three_partialdependence_maxentmodel_maxentfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



#%% Intron retention - noise
Xs=pd.read_pickle('./dataframes/ml/noise/Xs_train_irnoise.pkl')

####

cols={'maxent':['maxent5','maxent3'],
      'secondary':['exon1','donor','intron5','introncenter','intron3','acceptor','exon2'],
      'hexamers':[x for x in Xs.columns if (x.split('_')[-1] in ['exon1','intron','exon2'])&(len(x.split('_'))==2)],
      'motifs':[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                'all features':':'}


for feat in ['maxent','secondary','hexamers']:
    mdl=pickle.load(open('./ml/models/irmodelsel_noiseresgam_trained_'+feat+'.sav', 'rb'))

    explainer = shap.TreeExplainer(mdl)
    shap_values=explainer.shap_values(Xs[cols[feat]])
    f=plt.figure()
    shap.summary_plot(shap_values, Xs[cols[feat]], feature_names=cols[feat])
    f.savefig('./ml/plots/shap/irnoise_'+feat+'_summaryplot_shapvalues.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/irmodelsel_noiseresgam_trained_motifs.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols['motifs']])
f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+x.split('_')[-1] for x in cols['motifs']])
f.savefig('./ml/plots/shap/irnoise_motifs_summaryplot_shapvalues.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/irmodelsel_noiseresgam_trained_motifs.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs[cols['motifs']])
f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+x for x in cols['motifs']])
f.savefig('./ml/plots/shap/irnoise_motifs_summaryplot_shapvalues_inclattractID.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure()
shap.summary_plot(shap_values, Xs[cols['motifs']], feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+database.loc[x.split('__')[0],'Motif'] + ' ' + x.split('_')[-1] for x in cols['motifs']])
f.savefig('./ml/plots/shap/irnoise_motifs_summaryplot_shapvalues_inclmotifs.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

mdl=pickle.load(open('./ml/models/irmodelsel_noiseresgam_trained.sav', 'rb'))
explainer = shap.TreeExplainer(mdl)
shap_values=explainer.shap_values(Xs)
f=plt.figure()
shap.summary_plot(shap_values, Xs, feature_names=[database.loc[x.split('__')[0],'Gene_name']+' '+database.loc[x.split('__')[0],'Motif']+' '+x.split('_')[-1] if 'score_motifs' in x else x for x in Xs.columns])
f.savefig('./ml/plots/shap/irnoise_allfeatures_summaryplot_shapvalues.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
#######################   
''' 
shap.dependence_plot('maxent5',shap_values,Xs)
shap.dependence_plot('maxent3',shap_values,Xs[cols[feat]])

shap_interaction_values=explainer.shap_interaction_values(Xs[cols[feat]])
shap.summary_plot(shap_interaction_values, Xs[cols[feat]])

shap.dependence_plot(('maxent5','maxent3'), shap_interaction_values,Xs)


shap.dependence_plot('exon1',shap_values,Xs[cols[feat]])

shap.dependence_plot('maxent3',shap_values,Xs[cols[feat]])

shap_interaction_values=explainer.shap_interaction_values(Xs[cols[feat]])
shap.summary_plot(shap_interaction_values, Xs[cols[feat]])

shap.dependence_plot(('acceptor','exon2'), shap_interaction_values,Xs[cols[feat]])

shap.force_plot(explainer.expected_value, shap_values[0,:], Xs[cols[feat]].iloc[0,:])


feat='motifs'
mdl=pickle.load(open('./ml/models/irmodel_noiseresgam_trained_'+feat+'.sav', 'rb'))
shap.dependence_plot('exon1',shap_values,Xs[cols[feat]])

explainer = shap.TreeExplainer(mdl)
shap_values_motifs=explainer.shap_values(Xs[cols[feat]])
shap_interaction_values_motifs=explainer.shap_interaction_values(Xs[cols[feat]])

'''
#####

mdl=pickle.load(open('./ml/models/irmodel_noiseres_gam_trained.sav', 'rb'))
mdl=pickle.load(open('./ml/models/irmodelsel_noiseresgam_trained.sav', 'rb'))
'''
#maxent3second vs alt3/acceptor2
plot_partial_dependence(mdl, Xs, [(0,2),(0,3)], grid_resolution=100)
plot_partial_dependence(mdl, Xs, [0,2,3], grid_resolution=100)

#maxent3first vs intron/acceptor1
plot_partial_dependence(mdl, Xs, [(1,6),(1,7)], grid_resolution=150)
plot_partial_dependence(mdl, Xs, [1,6,7,8], grid_resolution=100)

plot_partial_dependence(mdl, Xs, [0,1,(0,1)], grid_resolution=150)
plot_partial_dependence(mdl, Xs, [3,7,(3,7)], grid_resolution=150)
plot_partial_dependence(mdl, Xs, [(0,1),(3,7)], grid_resolution=200)
plot_partial_dependence(mdl, Xs, [(0,1),(3,7)], grid_resolution=100)
'''
#maxent3first vs acceptor1/maxent3second vs acceptor2

f, axs = plot_partial_dependence(mdl, Xs, [(0,3),(1,7)], grid_resolution=200)
f.savefig('./ml/plots/ir_partialdependence_maxentvssecondary_acceptordonor.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f, axs = plot_partial_dependence(mdl, Xs, [(0,1),(3,7)], grid_resolution=100)
f.savefig('./ml/plots/ir_partialdependence_maxentdonorvsacceptor_secondarydonorvsacceptor.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


pdp_1, axes_1 = partial_dependence(mdl, 0, X=Xs, grid_resolution=100)
pdp_2, axes_2 = partial_dependence(mdl, 1, X=Xs, grid_resolution=100)

f=plt.figure(figsize=(4,3))
plt.plot(axes_1[0],pdp_1[0])
plt.plot(axes_2[0],pdp_2[0])
plt.axhline(y=0,linewidth=2, color='gray', alpha=0.5)
f.savefig('./ml/plots/irnoise_partialdependence_maxentfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

pdp_1, axes_1 = partial_dependence(mdl, 3, X=Xs, grid_resolution=100)
pdp_2, axes_2 = partial_dependence(mdl, 7, X=Xs, grid_resolution=100)

f=plt.figure(figsize=(4,3))
plt.plot(axes_1[0],pdp_1[0])
plt.plot(axes_2[0],pdp_2[0])
plt.axhline(y=0,linewidth=2, color='gray', alpha=0.5)
f.savefig('./ml/plots/irnoise_partialdependence_secondaryfirst_and_second.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


#%% compare feature importances for prediction
########

df_features_irtrain=pd.read_pickle('./dataframes/ml/ir/Xs_trainnew_ir.pkl')
Xsnoise=pd.read_pickle('./dataframes/ml/noise/Xs_train_irnoise.pkl')
Xsnoise_test=pd.read_pickle('./dataframes/ml/noise/Xs_test_irnoise.pkl')

sharedcols=list(dict.fromkeys([col if col in Xsnoise.columns else '' for col in df_features_irtrain]))[1:]

forprediction_irnoise.train_gbr_model(Xsnoise[sharedcols], irdf.levelratiospl[Xsnoise.index], 'ir_logratio_sharedfeatures_trained')
forprediction_irnoise.make_prediction(Xsnoise_test[sharedcols], irdf.levelratiospl[Xsnoise_test.index], mode='ir_logratio_sharedfeatures_trained')
#r2=0.787

forprediction_irnoise.train_gbr_model(Xsnoise[sharedcols], irdf.wav_stats[Xsnoise.index], 'ir_wavstats_sharedfeatures_trained')
forprediction_irnoise.make_prediction(Xsnoise_test[sharedcols], irdf.wav_stats[Xsnoise_test.index], mode='ir_wavstats_sharedfeatures_trained')
#r2=0.684

forprediction_irnoise.train_gbr_model(Xsnoise[sharedcols], irdf.noiseresgam[Xsnoise.index], 'ir_noiseresgam_sharedfeatures_trained')
forprediction_irnoise.make_prediction(Xsnoise_test[sharedcols], irdf.noiseresgam[Xsnoise_test.index], mode='ir_noiseresgam_sharedfeatures_trained')
#r2=0.059

mdl_rna=pickle.load(open('./ml/models/ir_logratio_sharedfeatures_trained.sav', 'rb'))
mdl_protein=pickle.load(open('./ml/models/ir_wavstats_sharedfeatures_trained.sav', 'rb'))
mdl_noise=pickle.load(open('./ml/models/ir_noiseresgam_sharedfeatures_trained.sav', 'rb'))

fi=pd.DataFrame([mdl_rna.feature_importances_, mdl_protein.feature_importances_, mdl_noise.feature_importances_]).transpose()
fi.columns=['RNA','protein','noise']
fi.index=sharedcols

cols={'maxent':['maxent5','maxent3'],
      'secondary':['exon1','donor','intron5','introncenter','intron3','acceptor','exon2'],
      'hexamers':[x for x in sharedcols if (x.split('_')[-1] in ['exon1','intron','exon2'])&(len(x.split('_'))==2)],
      'motifs':[x for x in sharedcols if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                'all features':':'}


f=plt.figure(figsize=(3,3))
plt.scatter(fi.RNA, fi.protein, alpha=0.5)
plt.scatter(fi.RNA, fi.noise, alpha=0.5)

plt.scatter(fi.protein, fi.noise, alpha=0.2)
plt.scatter(fi.protein[cols['maxent']], fi.noise[cols['maxent']], color='r')
plt.scatter(fi.protein[cols['secondary']], fi.noise[cols['secondary']], color='magenta')

plt.scatter(fi.RNA, fi.noise, alpha=0.2)
plt.scatter(fi.RNA[cols['maxent']], fi.noise[cols['maxent']], color='r')
plt.scatter(fi.RNA[cols['secondary']], fi.noise[cols['secondary']], color='magenta')

toannot=set(cols['maxent']+list(fi.RNA.sort_values(ascending=False).index)[:5] + list(fi.protein.sort_values(ascending=False).index)[:5])
f=plt.figure(figsize=(4,4))
plt.scatter(fi.RNA, fi.protein, alpha=0.2)
plt.scatter(fi.RNA[cols['maxent']], fi.protein[cols['maxent']], color='r')
plt.scatter(fi.RNA[cols['secondary']], fi.protein[cols['secondary']], color='magenta')
plt.xlim(-0.002, 0.102)
plt.ylim(-0.002, 0.062)
plt.title('Pearson r='+'{:.2f}'.format(pearsonr(fi.RNA, fi.protein)[0]) + \
          '\nSpearman r='+'{:.2f}'.format(spearmanr(fi.RNA, fi.protein)[0]), fontsize=12)
for x in toannot:
    if 'score_motifs' in x:
        plt.annotate(database.loc[x.split('__')[0],'Gene_name']+' '+x.split('_')[-1], xy=(fi.RNA[x], fi.protein[x]), fontsize=11)
    else:
        plt.annotate(x, xy=(fi.RNA[x], fi.protein[x]), fontsize=11)
f.savefig('./ml/plots/ir_featureimportances_sharedfeatures_rnavsprotein.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


toannot=set(cols['maxent']+list(fi.RNA.sort_values(ascending=False).index)[:5] + list(fi.noise.sort_values(ascending=False).index)[:5])
f=plt.figure(figsize=(4,4))
plt.scatter(fi.RNA, fi.noise, alpha=0.2)
plt.scatter(fi.RNA[cols['maxent']], fi.noise[cols['maxent']], color='r')
plt.scatter(fi.RNA[cols['secondary']], fi.noise[cols['secondary']], color='magenta')
plt.xlim(-0.002, 0.102)
plt.ylim(-0.0005, 0.0205)
plt.title('Pearson r='+'{:.2f}'.format(pearsonr(fi.RNA, fi.noise)[0]) + \
          '\nSpearman r='+'{:.2f}'.format(spearmanr(fi.RNA, fi.noise)[0]), fontsize=12)
for x in toannot:
    if 'score_motifs' in x:
        plt.annotate(database.loc[x.split('__')[0],'Gene_name']+' '+x.split('_')[-1], xy=(fi.RNA[x], fi.noise[x]), fontsize=11)
    else:
        plt.annotate(x, xy=(fi.RNA[x], fi.noise[x]), fontsize=11)
f.savefig('./ml/plots/ir_featureimportances_sharedfeatures_rnavsnoise.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


toannot=set(cols['maxent']+list(fi.protein.sort_values(ascending=False).index)[:5] + list(fi.noise.sort_values(ascending=False).index)[:5])
f=plt.figure(figsize=(4,4))
plt.scatter(fi.protein, fi.noise, alpha=0.2)
plt.scatter(fi.protein[cols['maxent']], fi.noise[cols['maxent']], color='r')
plt.scatter(fi.protein[cols['secondary']], fi.noise[cols['secondary']], color='magenta')
plt.xlim(-0.002, 0.062)
plt.ylim(-0.0005, 0.0205)
plt.title('Pearson r='+'{:.2f}'.format(pearsonr(fi.protein, fi.noise)[0]) + \
          '\nSpearman r='+'{:.2f}'.format(spearmanr(fi.protein, fi.noise)[0]), fontsize=12)
for x in toannot:
    if 'score_motifs' in x:
        plt.annotate(database.loc[x.split('__')[0],'Gene_name']+' '+x.split('_')[-1], xy=(fi.protein[x], fi.noise[x]), fontsize=11)
    else:
        plt.annotate(x, xy=(fi.protein[x], fi.noise[x]), fontsize=11)
f.savefig('./ml/plots/ir_featureimportances_sharedfeatures_proteinvsnoise.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



toannot=[]
f=plt.figure(figsize=(4,4))
plt.scatter(fi.RNA, fi.protein, alpha=0.2)
plt.scatter(fi.RNA[cols['maxent']], fi.protein[cols['maxent']], color='r')
plt.scatter(fi.RNA[cols['secondary']], fi.protein[cols['secondary']], color='magenta')
plt.xlim(-0.002, 0.102)
plt.ylim(-0.002, 0.062)
plt.title('Pearson r='+'{:.2f}'.format(pearsonr(fi.RNA, fi.protein)[0]) + \
          '\nSpearman r='+'{:.2f}'.format(spearmanr(fi.RNA, fi.protein)[0]), fontsize=12)
for x in toannot:
    if 'score_motifs' in x:
        plt.annotate(database.loc[x.split('__')[0],'Gene_name']+' '+x.split('_')[-1], xy=(fi.RNA[x], fi.protein[x]), fontsize=11)
    else:
        plt.annotate(x, xy=(fi.RNA[x], fi.protein[x]), fontsize=11)
f.savefig('./ml/plots/ir_featureimportances_sharedfeatures_rnavsprotein_withoutannot.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(4,4))
plt.scatter(fi.RNA, fi.noise, alpha=0.2)
plt.scatter(fi.RNA[cols['maxent']], fi.noise[cols['maxent']], color='r')
plt.scatter(fi.RNA[cols['secondary']], fi.noise[cols['secondary']], color='magenta')
plt.xlim(-0.002, 0.102)
plt.ylim(-0.0005, 0.0205)
plt.title('Pearson r='+'{:.2f}'.format(pearsonr(fi.RNA, fi.noise)[0]) + \
          '\nSpearman r='+'{:.2f}'.format(spearmanr(fi.RNA, fi.noise)[0]), fontsize=12)
for x in toannot:
    if 'score_motifs' in x:
        plt.annotate(database.loc[x.split('__')[0],'Gene_name']+' '+x.split('_')[-1], xy=(fi.RNA[x], fi.noise[x]), fontsize=11)
    else:
        plt.annotate(x, xy=(fi.RNA[x], fi.noise[x]), fontsize=11)
f.savefig('./ml/plots/ir_featureimportances_sharedfeatures_rnavsnoise_withoutannot.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(4,4))
plt.scatter(fi.protein, fi.noise, alpha=0.2)
plt.scatter(fi.protein[cols['maxent']], fi.noise[cols['maxent']], color='r')
plt.scatter(fi.protein[cols['secondary']], fi.noise[cols['secondary']], color='magenta')
plt.xlim(-0.002, 0.062)
plt.ylim(-0.0005, 0.0205)
plt.title('Pearson r='+'{:.2f}'.format(pearsonr(fi.protein, fi.noise)[0]) + \
          '\nSpearman r='+'{:.2f}'.format(spearmanr(fi.protein, fi.noise)[0]), fontsize=12)
for x in toannot:
    if 'score_motifs' in x:
        plt.annotate(database.loc[x.split('__')[0],'Gene_name']+' '+x.split('_')[-1], xy=(fi.protein[x], fi.noise[x]), fontsize=11)
    else:
        plt.annotate(x, xy=(fi.protein[x], fi.noise[x]), fontsize=11)
f.savefig('./ml/plots/ir_featureimportances_sharedfeatures_proteinvsnoise_withoutannot.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



toannot=set(cols['maxent']+cols['secondary']+list(fi.RNA.sort_values(ascending=False).index)[:10] + list(fi.protein.sort_values(ascending=False).index)[:10])
f=plt.figure(figsize=(8,8))
plt.scatter(fi.RNA, fi.protein, alpha=0.2)
plt.scatter(fi.RNA[cols['maxent']], fi.protein[cols['maxent']], color='r')
plt.scatter(fi.RNA[cols['secondary']], fi.protein[cols['secondary']], color='magenta')
plt.xlim(-0.002, 0.102)
plt.ylim(-0.002, 0.062)
plt.title('Pearson r='+'{:.2f}'.format(pearsonr(fi.RNA, fi.protein)[0]) + \
          '\nSpearman r='+'{:.2f}'.format(spearmanr(fi.RNA, fi.protein)[0]), fontsize=12)
for x in toannot:
    if 'score_motifs' in x:
        plt.annotate(database.loc[x.split('__')[0],'Gene_name']+' '+x.split('_')[-1], xy=(fi.RNA[x], fi.protein[x]), fontsize=11)
    else:
        plt.annotate(x, xy=(fi.RNA[x], fi.protein[x]), fontsize=11)
f.savefig('./ml/plots/ir_featureimportances_sharedfeatures_rnavsprotein_large.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


toannot=set(cols['maxent']+cols['secondary']+list(fi.RNA.sort_values(ascending=False).index)[:10] + list(fi.noise.sort_values(ascending=False).index)[:10])
f=plt.figure(figsize=(8,8))
plt.scatter(fi.RNA, fi.noise, alpha=0.2)
plt.scatter(fi.RNA[cols['maxent']], fi.noise[cols['maxent']], color='r')
plt.scatter(fi.RNA[cols['secondary']], fi.noise[cols['secondary']], color='magenta')
plt.xlim(-0.002, 0.102)
plt.ylim(-0.0005, 0.0205)
plt.title('Pearson r='+'{:.2f}'.format(pearsonr(fi.RNA, fi.noise)[0]) + \
          '\nSpearman r='+'{:.2f}'.format(spearmanr(fi.RNA, fi.noise)[0]), fontsize=12)
for x in toannot:
    if 'score_motifs' in x:
        plt.annotate(database.loc[x.split('__')[0],'Gene_name']+' '+x.split('_')[-1], xy=(fi.RNA[x], fi.noise[x]), fontsize=11)
    else:
        plt.annotate(x, xy=(fi.RNA[x], fi.noise[x]), fontsize=11)
f.savefig('./ml/plots/ir_featureimportances_sharedfeatures_rnavsnoise_large.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


toannot=set(cols['maxent']+cols['secondary']+list(fi.protein.sort_values(ascending=False).index)[:10] + list(fi.noise.sort_values(ascending=False).index)[:10])
f=plt.figure(figsize=(8,8))
plt.scatter(fi.protein, fi.noise, alpha=0.2)
plt.scatter(fi.protein[cols['maxent']], fi.noise[cols['maxent']], color='r')
plt.scatter(fi.protein[cols['secondary']], fi.noise[cols['secondary']], color='magenta')
plt.xlim(-0.002, 0.062)
plt.ylim(-0.0005, 0.0205)
plt.title('Pearson r='+'{:.2f}'.format(pearsonr(fi.protein, fi.noise)[0]) + \
          '\nSpearman r='+'{:.2f}'.format(spearmanr(fi.protein, fi.noise)[0]), fontsize=12)
for x in toannot:
    if 'score_motifs' in x:
        plt.annotate(database.loc[x.split('__')[0],'Gene_name']+' '+x.split('_')[-1], xy=(fi.protein[x], fi.noise[x]), fontsize=11)
    else:
        plt.annotate(x, xy=(fi.protein[x], fi.noise[x]), fontsize=11)
f.savefig('./ml/plots/ir_featureimportances_sharedfeatures_proteinvsnoise_large.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


#%%

#############
# COMPARE FEATURE IMPORTANCES BETWEEN SPLICING TYPES
#############

Xs=pd.read_pickle('./dataframes/ml/ir/Xs_trainnew_ir.pkl')

cols={'maxent':['maxent5','maxent3'],
      'secondary':['exon1','donor','intron5','introncenter','intron3','acceptor','exon2'],
      'hexamers':[x for x in Xs.columns if (x.split('_')[-1] in ['exon1','intron','exon2'])&(len(x.split('_'))==2)],
      'motifs':[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                'all features':':'}

mdl=pickle.load(open('./ml/models/irmodel_logratio_trained_hexamers.sav', 'rb'))

featimpir=pd.DataFrame(index=cols['hexamers'])
featimpir['featimp']=mdl.feature_importances_
featimpir['sixmerloc']=featimpir.index.map(lambda idx: idx.replace('exon1','up').replace('intron','alt').replace('exon2','down'))    
featimpir['sixmer']=featimpir.index.map(lambda idx: idx.split('_')[0])


Xs=pd.read_pickle('./dataframes/ml/cas/Xs_train_cas.pkl')

cols={'maxent':['maxent5','maxent3'],
      'secondary':['intron1','acceptor','exon5','exoncenter','exon3','donor','intron2'],
      'hexamers':[x for x in Xs.columns if (x.split('_')[-1] in ['intron1','exon','intron2'])&(len(x.split('_'))==2)],
      'motifs':[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                'all features':':'}

mdl=pickle.load(open('./ml/models/casmodel_logratio_trained_hexamers.sav', 'rb'))

featimpcas=pd.DataFrame(index=cols['hexamers'])
featimpcas['featimp']=mdl.feature_importances_
featimpcas['sixmerloc']=featimpcas.index.map(lambda idx: idx.replace('intron1','up').replace('exon','alt').replace('intron2','down'))    
featimpcas['sixmer']=featimpcas.index.map(lambda idx: idx.split('_')[0])


Xs=pd.read_pickle('./dataframes/ml/five/Xs_train_five.pkl')

cols={'maxent':['maxent5first','maxent5second'],
      'secondary':['exon','donor1','alt5','altcenter','alt3','donor2','intron'],
      'hexamers':[x for x in Xs.columns if (x.split('_')[-1] in ['intron','alt','exon'])&(len(x.split('_'))==2)],
      'motifs':[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                'all features':':'}

mdl=pickle.load(open('./ml/models/fivemodel_logratio_trained_hexamers.sav', 'rb'))

featimpfive=pd.DataFrame(index=cols['hexamers'])
featimpfive['featimp']=mdl.feature_importances_
featimpfive['sixmerloc']=featimpfive.index.map(lambda idx: idx.replace('exon','up').replace('alt','alt').replace('intron','down'))    
featimpfive['sixmer']=featimpfive.index.map(lambda idx: idx.split('_')[0])


Xs=pd.read_pickle('./dataframes/ml/three/Xs_train_three.pkl')

cols={'maxent':['maxent3first','maxent3second'],
      'secondary':['intron','acceptor1','alt5','altcenter','alt3','acceptor2','exon'],
      'hexamers':[x for x in Xs.columns if (x.split('_')[-1] in ['intron','alt','exon'])&(len(x.split('_'))==2)],
      'motifs':[x for x in Xs.columns if (x.split('_')[-1] in ['up','alt','down'])&(len(x.split('_'))>2)],
                'all features':':'}

mdl=pickle.load(open('./ml/models/threemodel_logratio_trained_hexamers.sav', 'rb'))

featimpthree=pd.DataFrame(index=cols['hexamers'])
featimpthree['featimp']=mdl.feature_importances_
featimpthree['sixmerloc']=featimpthree.index.map(lambda idx: idx.replace('intron','up').replace('alt','alt').replace('exon','down'))    
featimpthree['sixmer']=featimpthree.index.map(lambda idx: idx.split('_')[0])

featimp=featimpir.groupby(by='sixmer').sum().join(pd.Series(['intron ret.']*len(featimpir.sixmer.unique()),name='spltype', index=featimpir.sixmer.unique())).reset_index()\
    .append(featimpcas.groupby(by='sixmer').sum().join(pd.Series(['cas. exon']*len(featimpcas.sixmer.unique()),name='spltype', index=featimpcas.sixmer.unique())).reset_index(), ignore_index=True)\
    .append(featimpfive.groupby(by='sixmer').sum().join(pd.Series(["tandem 5'"]*len(featimpfive.sixmer.unique()),name='spltype', index=featimpfive.sixmer.unique())).reset_index(), ignore_index=True)\
    .append(featimpthree.groupby(by='sixmer').sum().join(pd.Series(["tandem 3'"]*len(featimpthree.sixmer.unique()),name='spltype', index=featimpthree.sixmer.unique())).reset_index(), ignore_index=True)

spltypes=['intron ret.','cas. exon',"tandem 5'","tandem 3'"]
featimpcorrpearson=pd.DataFrame('',index=spltypes, columns=spltypes)
for i in spltypes:
    for j in spltypes:
        if i==j:
            featimpcorrpearson.loc[i,j]=1.0
        else:
            featimpcorrpearson.loc[i,j]=featimp[(featimp.spltype==i)|(featimp.spltype==j)].pivot(index='sixmer', columns='spltype', values='featimp').corr().values[0][1]

featimpcorrspearman=pd.DataFrame('',index=spltypes, columns=spltypes)
for i in spltypes:
    for j in spltypes:
        if i==j:
            featimpcorrspearman.loc[i,j]=np.float(1.0)
        else:
            featimpcorrspearman.loc[i,j]=np.float(featimp[(featimp.spltype==i)|(featimp.spltype==j)].pivot(index='sixmer', columns='spltype', values='featimp').corr(method='spearman').values[0][1])


for i in spltypes:
    for j in spltypes:
        if i!=j:
            f=plt.figure(figsize=(3,3))
            plt.scatter(featimp[(featimp.spltype==i)|(featimp.spltype==j)].pivot(index='sixmer', columns='spltype', values='featimp')[i],\
                                featimp[(featimp.spltype==i)|(featimp.spltype==j)].pivot(index='sixmer', columns='spltype', values='featimp')[j],\
                                        s=20,alpha=0.2)
            plt.xlim(-0.0002,0.0152)
            plt.ylim(-0.0002,0.0152)
            plt.xticks([0,0.005,0.01,0.015])
            plt.yticks([0,0.005,0.01,0.015])
            plt.xlabel(i)
            plt.ylabel(j)
            for six in ['TTTCAG','GTAAGT','TTCCTT','CCTTTC']:
                try:
                    plt.annotate(six, xy=(featimp[(featimp.spltype==i)&(featimp.sixmer==six)].featimp.values[0], featimp[(featimp.spltype==j)&(featimp.sixmer==six)].featimp.values[0]), fontsize=12, color='r')
                except:
                    pass
            f.savefig('./figures/Fig4/FigS4X_compare_hexamer_importance_'+i+'_'+j+'.png',
                      dpi = 300, format='png', bbox_inches='tight', frameon=True)



