#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 18:13:06 2017

@author: martinm
"""


import pandas as pd
import numpy as np

#%%

ir_cov300raw=pd.read_table('../rawdata/ir_protein/CoverageIR.tab', header=None)

ir_cov300raw.columns=['binprimer','variant','number_reads']
    
ir_cov300split=ir_cov300raw.variant.apply(lambda x: np.int(x.split('|')[2]))

ir_cov300bin=ir_cov300raw.binprimer.apply(lambda x: np.int(x.split('#')[1].split('_')[0]))
ir_cov300binrep=ir_cov300raw.binprimer.apply(lambda x: np.int(x.split('#')[1].split('_')[1]))

ir_cov300=pd.concat([ir_cov300raw,ir_cov300split,ir_cov300bin,ir_cov300binrep],axis=1)
ir_cov300.columns=['binprimer','variant','number_reads','varid','bin','primernr']
ir_cov300.sort_values(['varid','bin'],inplace=True)
ir_cov300.set_index(['varid','binprimer'],inplace=True)

ir_cov300.to_pickle('../rawdata/ir_protein/ir_cov300.pkl')

#%%

ir_cov300=pd.read_pickle('../rawdata/ir_protein/ir_cov300.pkl')
ir_cov300_totalreadnumber=ir_cov300.apply(lambda x: x.sum(axis=0,level=0))
ir_cov300_totalreadnumber.drop(['variant','bin'],axis=1,inplace=True)


ir_cov300_totalreadnumber.sum() #6,090,074
ir_cov300_totalreadnumber[ir_cov300_totalreadnumber.number_reads==0].count() #310
ir_cov300_totalreadnumber[ir_cov300_totalreadnumber.number_reads<50].count() #952
ir_cov300_totalreadnumber[ir_cov300_totalreadnumber.number_reads<100].count() #1215

                     
                     #%%

### primer effects - correlation between reads from different barcodes for the same bin

#%% Create dataframe for summing over primers 1+2 or over all 3, including total readnumber

# for primers 1 and 2 only

cov300=pd.DataFrame()
for var, group in ir_cov300.groupby(level=0):
    for binnumber in range(1,17):
        summed=group[(group.bin == binnumber) & (group.primernr != 3)].number_reads.sum()
        cov300=cov300.append(pd.Series([var,binnumber,summed]),ignore_index=True)

cov300.columns=['varid','bin','summedreads']
cov300.set_index(['varid','bin'],inplace=True)
    
cov300.to_pickle('../rawdata/ir_protein/cov300ir.pkl')

#%% APPLY FILTERS - only primers 1 and 2


import forexpressionanalysis

cov300f1=forexpressionanalysis.FilterMinVariantCoverage(cov300,50)    # 62346 datapoints==0
   
cov300f2=forexpressionanalysis.FilterMinBinCoverage(cov300f1,5)    # 87294 datapoints==0

cov300f3=forexpressionanalysis.FilterPercentPerBin(cov300f2,2)    # 98791 datapoints==0

#cov300f35=forexpressionanalysis.FilterExtremeBins(cov300f3,0.2)                                                   

cov300f4=forexpressionanalysis.FilterIsolated(cov300f3)    # 103326 datapoints==0

cov300f5=forexpressionanalysis.FilterPercentKept(cov300, cov300f4, 30)    # 103427 datapoints==0

cov300f5.to_pickle('../rawdata/ir_protein/cov300filteredir.pkl')   

#%% check the relative number of reads from each bin and compare with the percentage of cells


binreads300=forexpressionanalysis.BinReadNumbers(cov300)

binpercentages=pd.read_excel('../rawdata/ir_protein/percentofcellssorted_intronretention.xlsx', index_col='bin') 


#%% correct for number of reads per bin and percent sorted

correctionfactor12 = binpercentages.percentcells.divide(100*binreads300/binreads300.sum())

cov300fnorm = forexpressionanalysis.NormalizeByBinPercentages(cov300f5, correctionfactor12)

#sanity check
binreads300norm=forexpressionanalysis.BinReadNumbers(cov300fnorm)


#%% expression estimates
ratiocalcbins=pd.read_excel('../rawdata/ir_protein/calculatedGFPmCherryratiosforbins.xlsx')
xval=np.log2(ratiocalcbins.loc[:,'median']*100)
xvalwidth=[]

for i in range(len(xval)):
    if (i==15):
        xvalwidth.append(0.8)
    else:
        xvalwidth.append(xval.loc[i+2] - xval.loc[i+1])
        

expression300 = forexpressionanalysis.ExpCalcStats(cov300fnorm, xval)


#%%
# normalization, peak detection and smoothing for primers 1 and 2 alone

data300 = forexpressionanalysis.NormalizeToReadNumbers(cov300fnorm)
data300.to_pickle('../rawdata/ir_protein/data300.pkl')

from peakdetect import peakdet
from savitzky_golay import savitzky_golay

data300smoothed = data300.copy()
peaks300=pd.DataFrame()
peaksinfo300=[]
for var,group in data300.groupby(level=0):
    peaks300.loc[int(var), 'rawnumberofpeaks'] = len(peakdet(list(group.summedreads), 0.02, list(xval))[0])
    smoothed = savitzky_golay(list(group.summedreads), 3, 1)
    for binn in np.arange(1,17):
        data300smoothed.loc[int(var)].loc[binn] = smoothed[binn-1]
    peaks300.loc[int(var), 'smoothednumberofpeaks'] = len(peakdet(smoothed, 0.02, list(xval))[0])
    peaksinfo300.append(peakdet(smoothed, 0.02, list(xval)))

peakspositions300=pd.DataFrame(peaksinfo300)
peakspositions300.columns=['maxima','minima']
peakspositions300.index= peaks300.index


for var in peakspositions300.index:
    if (peaks300.loc[var,'smoothednumberofpeaks']>0):
        try:
            peaks300.loc[var, 'xpeak1'] = peakspositions300.ix[var,0][0][0]
            peaks300.loc[var, 'ypeak1'] = peakspositions300.ix[var,0][0][1]
        except:
            peaks300.loc[var, 'xpeak1'] = xval[data300.loc[var].summedreads.tolist().index(data300.loc[var].summedreads.max()) + 1]
            peaks300.loc[var, 'ypeak1'] = data300.loc[var].summedreads.max()            
    if (peaks300.loc[var,'smoothednumberofpeaks']>1):
        peaks300.loc[var, 'xpeak2'] = peakspositions300.ix[var,0][1][0]
        peaks300.loc[var, 'ypeak2'] = peakspositions300.ix[var,0][1][1]
    if (peaks300.loc[var,'smoothednumberofpeaks']>2):
        peaks300.loc[var, 'xpeak3'] = peakspositions300.ix[var,0][2][0]
        peaks300.loc[var, 'ypeak3'] = peakspositions300.ix[var,0][2][1]


#%% 
### Construct one dataframe with all the relevant information

expr300 = expression300.join(ir_cov300_totalreadnumber)
expr300= expr300.join(peaks300)
#expr300=expr300.join(expression300notlog, rsuffix='_notlog')


library=pd.read_pickle('../design/LIBRARY/alllibraries210.pkl')

lib300=library[library.library=='ir']
lib300 = lib300.dropna(axis=1, how='all')

lib300=lib300.join(expr300)

lib300.to_pickle('./dataframes/IR_lib300.pkl')