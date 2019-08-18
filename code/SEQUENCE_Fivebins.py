#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 16:25:06 2017

@author: martinm
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
# from Bio.Seq import Seq

#dropboxprefix=str('/Users/miklm/')
dropboxprefix='C:/Users/martinm/'
fivebins = '/net/mraid08/export/genie/Runs/Martin/fivebins/'


#%%

fivecovraw=pd.read_table(homedir + 'MappingProtein/CoverageFIVE.tab', header=None)


fivecovraw.columns=['binprimer','variant','number_reads']
    
fivecovsplit=fivecovraw.variant.apply(lambda x: np.int(x.split('|')[2]))

fivecovbin=fivecovraw.binprimer.apply(lambda x: np.int(x.split('#')[1].split('_')[0]))
fivecovbinrep=fivecovraw.binprimer.apply(lambda x: np.int(x.split('#')[1].split('_')[1]))

fivecov=pd.concat([fivecovraw,fivecovsplit,fivecovbin,fivecovbinrep],axis=1)
fivecov.columns=['binprimer','variant','number_reads','varid','bin','primernr']
fivecov.sort(['varid','bin'],inplace=True)
fivecov.set_index(['varid','binprimer'],inplace=True)

fivecov.to_pickle(homedir + 'MappingProtein/fivecov.pkl')

#%%

fivecov_totalreadnumber=fivecov.apply(lambda x: x.sum(axis=0,level=0))
fivecov_totalreadnumber.drop(['variant','bin'],axis=1,inplace=True)



fivecov_totalreadnumber.sum() #11984106
fivecov_totalreadnumber[fivecov_totalreadnumber.number_reads==0].count() #130
fivecov_totalreadnumber[fivecov_totalreadnumber.number_reads<50].count() #982
fivecov_totalreadnumber[fivecov_totalreadnumber.number_reads<100].count() #1454




fivecov123=pd.DataFrame()
for var, group in fivecov.groupby(level=0):
    for binnumber in range(1,17):
        summed=group[group.bin==binnumber].number_reads.sum()
        fivecov123=fivecov123.append(pd.Series([var,binnumber,summed]),ignore_index=True)


fivecov123.columns=['varid','bin','summedreads']
fivecov123.set_index(['varid','bin'],inplace=True)

fivecov123.to_pickle(homedir + 'MappingProtein/fivecov123.pkl')

#%% APPLY FILTERS - all 3 primers


import forexpressionanalysis

fivecov123f1=forexpressionanalysis.FilterMinVariantCoverage(fivecov123,60)    # 30255 datapoints==0
   
fivecov123f2=forexpressionanalysis.FilterMinBinCoverage(fivecov123f1,5)    # 46739 datapoints==0

fivecov123f3=forexpressionanalysis.FilterPercentPerBin(fivecov123f2,2)    # 67839 datapoints==0

fivecov123f35=forexpressionanalysis.FilterExtremeBins(fivecov123f3,0.2)    # 68822 datapoints==0                                                   
                                                    
fivecov123f4=forexpressionanalysis.FilterIsolated(fivecov123f35)    # 71972 datapoints==0
    
fivecov123f5=forexpressionanalysis.FilterPercentKept(fivecov123, fivecov123f4, 30)    # 72074 datapoints==0

fivecov123f5.to_pickle(homedir + 'MappingProtein/fivecov123filtered.pkl')   


#%% check how much is left

fivecov123f_totalreadnumber=fivecov123f5.apply(lambda x: x.sum(axis=0,level=0))
fivecov123f_totalreadnumber[fivecov123f_totalreadnumber==0].count()
fivecov123f_totalreadnumber[fivecov123f_totalreadnumber>0].count()



#%% check the relative number of reads from each bin and compare with the percentage of cells


binreads123=forexpressionanalysis.BinReadNumbers(fivecov123)

binpercentages=pd.read_excel(homedir + 'MappingProtein/percentofcellssorted_five.xlsx', index_col='bin') 


#%% correct for number of reads per bin and percent sorted

correctionfactor = binpercentages.percentcells.divide(100*binreads123/binreads123.sum())

fivecov123fnorm = forexpressionanalysis.NormalizeByBinPercentages(fivecov123f5, correctionfactor)

#sanity check

binreads123norm=forexpressionanalysis.BinReadNumbers(fivecov123fnorm)

#%% get ratiocalcbins from FACS csv file (median(mCh)/median(GFP))

ratiocalcbins=pd.DataFrame()
ratiocalcbins.loc[1,'median'] = 171/14123.0
ratiocalcbins.loc[2,'median'] = 215/10859.0
ratiocalcbins.loc[3,'median'] = 236/8077.0
ratiocalcbins.loc[4,'median'] = 250/6119.0
ratiocalcbins.loc[5,'median'] = 269/4889.0
ratiocalcbins.loc[6,'median'] = 286/4082.0
ratiocalcbins.loc[7,'median'] = 314/3571.0
ratiocalcbins.loc[8,'median'] = 390/3365.0
ratiocalcbins.loc[9,'median'] = 746/4238.0
ratiocalcbins.loc[10,'median'] = 1676/5986.0
ratiocalcbins.loc[11,'median'] = 2721/6548.0
ratiocalcbins.loc[12,'median'] = 4161/6904.0
ratiocalcbins.loc[13,'median'] = 5237/6403.0
ratiocalcbins.loc[14,'median'] = 6352/6343.0
ratiocalcbins.loc[15,'median'] = 7002/6153.0
ratiocalcbins.loc[16,'median'] = 7561/4988.0
                 
xvalues = np.log2(ratiocalcbins.loc[:,'median'])
xval=np.log2(ratiocalcbins.loc[:,'median']*100)

#%% expression estimates

expression123 = forexpressionanalysis.ExpCalc(fivecov123fnorm, xval)
expression123.to_pickle(homedir + 'MappingProtein/expression123FIVE.pkl')


#%% Retrieve information from alllibraries

library=pd.read_pickle(homedir + 'dataframes/alllibraries210.pkl')

fivelib=library[library.library=='five']

fivelib['varseq162']=fivelib.varseq.apply(lambda x: x[30:192])


#%% Normalize to number of reads/variant, peak detect on raw and smoothed data

fivedata123 = forexpressionanalysis.NormalizeToReadNumbers(fivecov123fnorm)

fivedata123.to_pickle(homedir + 'MappingProtein/fivedata123.pkl')

from peakdetect import peakdet
from savitzky_golay import savitzky_golay


fivedata123smoothed = fivedata123.copy()
peaks123=pd.DataFrame()
peaksinfo123=[]
for var,group in fivedata123.groupby(level=0):
    peaks123.loc[int(var), 'rawnumberofpeaks'] = len(peakdet(list(group.summedreads), 0.02, list(xval))[0])
    smoothed = savitzky_golay(list(group.summedreads), 3, 1)
    for binn in np.arange(1,17):
        fivedata123smoothed.loc[int(var)].loc[binn] = smoothed[binn-1]
    peaks123.loc[int(var), 'smoothednumberofpeaks'] = len(peakdet(smoothed, 0.02, list(xval))[0])
    peaksinfo123.append(peakdet(smoothed, 0.02, list(xval)))

peakspositions123=pd.DataFrame(peaksinfo123)
peakspositions123.columns=['maxima','minima']
peakspositions123.index = peaks123.index

for var in peakspositions123.index:
    if (peaks123.loc[var,'smoothednumberofpeaks']>0):
        peaks123.loc[var, 'xpeak1'] = peakspositions123.ix[var,0][0][0]
        peaks123.loc[var, 'ypeak1'] = peakspositions123.ix[var,0][0][1]
    if (peaks123.loc[var,'smoothednumberofpeaks']>1):
        peaks123.loc[var, 'xpeak2'] = peakspositions123.ix[var,0][1][0]
        peaks123.loc[var, 'ypeak2'] = peakspositions123.ix[var,0][1][1]
    if (peaks123.loc[var,'smoothednumberofpeaks']>2):
        peaks123.loc[var, 'xpeak3'] = peakspositions123.ix[var,0][2][0]
        peaks123.loc[var, 'ypeak3'] = peakspositions123.ix[var,0][2][1]
        
        

###


xvalwidth=[]

for i in range(len(xval)):
    if (i==15):
        xvalwidth.append(0.8)
    else:
        xvalwidth.append(xval.loc[i+2] - xval.loc[i+1])
        



import random
var = random.choice(fitsp123.index)

f=plt.figure()
plt.bar(xval, fivedata123.loc[var].summedreads, width=xvalwidth)
plt.yticks([0,0.2,0.4,0.6])
plt.title(var)
plt.axvline(x=fits123.loc[var, 'GammaMean'], color = 'r')
plt.axvline(x=fitsp123.loc[var,'GammaMean'], color = 'g')
f.set_size_inches(3,2)
f.show()

for varid in peaks123.index:
    if ((peaks123.loc[varid,'smoothednumberofpeaks'] == 0) & (peaks123.loc[varid,'rawnumberofpeaks'] > 0)):
        peaks123.loc[varid,'smoothednumberofpeaks'] =1
                    



#%% Deal with multiple peaks

### Construct one big dataframe with all the relevant information


fiveexpr=expression123.join(fivecov123f_totalreadnumber)
fiveexpr= fiveexpr.join(peaks123)

fiveexpr.to_pickle(homedir + 'MappingProtein/fiveexpr.pkl')

fivelib['varseq162']=fivelib.varseq.apply(lambda x: x[30:192])

fivelib=fivelib.join(fiveexpr)

fiveexpr.to_pickle(homedir + 'dataframes/fivelib.pkl')


