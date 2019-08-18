# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 16:55:16 2016

@author: martinm
"""
import pandas as pd
import numpy as np

def SumOverBinPrimers(inp):
    group=pd.DataFrame(inp)
    group.columns=['variant','number_reads','bin','primernr']
    binreads=[]
    for binnumber in range(1,17):
        summed=group[group.bin==binnumber].number_reads.sum()
        binreads=binreads.append(summed)
    return binreads


def FilterMinVariantCoverage(data, threshold):
    df=data.copy()
    occur=0
    for var, group in df.groupby(level=0):
        if (group.summedreads.sum()<threshold):
            df.loc[var]=0
            occur +=1
    print(occur)
    return df    


def FilterMinBinCoverage(data, threshold):
    df=data.copy()
    df[df.summedreads < threshold] = 0
    return df


def FilterPercentPerBin(data, threshold):
    df=data.copy()
    occur=0
    for var, group in df.groupby(level=0):
        binn=0
        for row in group.iterrows():
            if (row[1].summedreads/group.summedreads.sum()*100 < threshold):
                df.loc[var].iloc[binn,0]=0
                occur +=1
            binn=binn+1
    print(occur)
    return df       


def FilterExtremeBins(data, threshold):
    df=data.copy()
    occur = 0
    for var, group in df.groupby(level=0):
        if ((group.loc[var].loc[1,'summedreads'] > group.loc[var].loc[2, 'summedreads']) & (group.loc[var].loc[1,'summedreads'] > group.loc[var].loc[3,'summedreads']/threshold)):
            df.loc[var].loc[1,'summedreads']=0
            occur += 1
        if ((group.loc[var].loc[16,'summedreads'] > group.loc[var].loc[15, 'summedreads']) & (group.loc[var].loc[16,'summedreads'] > group.loc[var].loc[14,'summedreads']/threshold)):
            df.loc[var].loc[16,'summedreads']=0
            occur += 1
    print(occur)
    return df         
    

def FilterIsolated(data):
    df=data.copy()
    occur = 0
    for var, group in df.groupby(level=0):
        binn=1
        for row in group.iterrows():
            if (binn==1):
                if (row[1].summedreads>0)&(group.loc[var].loc[2,'summedreads']==0):
                    df.loc[var].loc[binn,'summedreads']=0
                    occur += 1
            elif (binn==16):
                if (row[1].summedreads>0)&(group.loc[var].loc[15,'summedreads']==0):
                    df.loc[var].loc[binn,'summedreads']=0
                    occur += 1
            else:
                if (row[1].summedreads>0)&(group.loc[var].loc[binn+1,'summedreads']==0)&(group.loc[var].loc[binn-1,'summedreads']==0):
                    df.loc[var].loc[binn,'summedreads']=0                    
                    occur += 1
            binn=binn+1
    print(occur)
    return df            




def FilterPercentKept(original, filtered, threshold):
    df=filtered.copy()
    occur = 0
    for var, group in df.groupby(level=0):
        if (group.summedreads.sum()>0):
            if (group.summedreads.sum()/original.loc[var].summedreads.sum()*100 <threshold):
                df.loc[var]=0
                occur +=1
    print(occur)
    return df            


def FilterAll(data, thr_bin_coverage, thr_bin_percent, thr_percent_kept):
    dataf = FilterMinBinCoverage(data, thr_bin_coverage)
    dataf = FilterPercentPerBin(dataf, thr_bin_percent)
    dataf = FilterIsolated(dataf)
    dataf = FilterPercentKept(data, dataf, thr_percent_kept)
    return dataf


def FilterAll2(data, thr_bin_coverage, thr_percent_kept):
    dataf = FilterMinBinCoverage(data, thr_bin_coverage)
    dataf = FilterIsolated(dataf)
    dataf = FilterPercentKept(data, dataf, thr_percent_kept)
    return dataf

def BinReadNumbers(df):
    binreads=pd.Series()
    for binnr in range(1,17):
        binstring='bin == ' + str(binnr)
        binreads.loc[binnr]=df.query(binstring).summedreads.sum()
    return binreads

def NormalizeByBinPercentages(data, correctionfactor):
    df=data.copy()
    for var, group in df.groupby(level=0):
        for binn in range(1,17):
            df.loc[var].loc[binn,'summedreads'] = data.loc[var].loc[binn,'summedreads'] * correctionfactor.loc[binn]
    return df           

def NormalizeToReadNumbers(data):
    df = data.copy()
    for var, group in df.groupby(level=0):
        for binn in range(1,17):
            df.loc[var].loc[binn,'summedreads'] = data.loc[var].loc[binn,'summedreads'] /data.loc[var,'summedreads'].sum()
    return df           
        

def ExpCalc(df, estimate):
    expr=pd.Series()
    for var, group in df.groupby(level=0):
        expr.loc[int(var)]= np.ma.average(estimate, weights=group.summedreads)
    return expr

def ExpCalcStats(df, estimate):
    from statsmodels.stats.weightstats import DescrStatsW
    expr=pd.DataFrame()
    for var, group in df.groupby(level=0):
        st=DescrStatsW(estimate, weights=group.summedreads)
        expr.loc[int(var),'wav_stats']= st.mean
        expr.loc[int(var),'wstd']= st.std
    return expr

def ExpCalcStatsReformatted(df, estimate):
    from statsmodels.stats.weightstats import DescrStatsW
    expr=pd.DataFrame()
    for var in df.index:
        st=DescrStatsW(estimate, weights=df.loc[var].values)
        expr.loc[int(var),'wav_stats']= st.mean
        expr.loc[int(var),'wstd']= st.std
    return expr
   
