#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 11:18:23 2017

@author: martinm
"""
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from scipy.stats.stats import pearsonr


def unbiased_mapping_ir(rnareads):    
    lib300=pd.read_pickle('./dataframes/ir_lib300.pkl')
    rna=pd.DataFrame()
    for var in lib300[lib300.subset!='irmelanomavariants'].index:
        varsq=lib300.varseq162[var]
        ind=1
        intronstarts=[]
        intronends=[]
        reads=[]
        rr=pd.Series(rnareads.loc[var].split(' '))
        readpos=0
        for test in rr:
        #    print(test)
            vspos=0
            test=str(Seq(test).reverse_complement())
            startpos=80
            testsub=test[startpos:90]
            while True:
        #        print(startpos)
                vspos=varsq.find(testsub)
                startpos-=1
                testsub=test[startpos:90]
                if (varsq.find(testsub)==-1)|(startpos<0):
                    if (startpos<0):
                        vspos=0
                    break
            intronends.append(vspos)
            if (vspos==0)|(varsq.find(test[startpos-9:startpos+1])==-1):
                intronstarts.append(0)
            elif (vspos==-1):
                intronstarts.append(-150)
            else:
                intronstarts.append(varsq.find(test[startpos-9:startpos+1])+10)
            reads.append(test)
            readpos+=1
            ind+=1    
        
        introns=pd.DataFrame([intronstarts,intronends, reads]).transpose()
        introns.columns=['intronstart','intronend', 'read']
        
        introns['intlen']=introns.index.map(lambda x: introns.intronend[x]-introns.intronstart[x])
    
        if (len(introns)>0):    
            rna.loc[var,'rnareads']=len(introns)
            rna.loc[var,'unspliced']=len(introns[introns.intlen==0])
            rna.loc[var,'spliced']=len(introns[introns.intlen==float(lib300.intronlength[var])])
            rna.loc[var,'fraction_canonical']=(rna.loc[var,'spliced']+rna.loc[var,'unspliced'])/float(rna.loc[var,'rnareads'])
            rna.loc[var,'fraction_unknowndonor']=len(introns[introns.intlen==-1])/float(rna.loc[var,'rnareads'])
            introns.drop(introns[(introns.intlen<11)|(introns.intlen>140)].index, inplace=True)
            vals=introns.intlen.value_counts()
            if (float(lib300.intronlength[var]) in vals.index):
                vals.drop(float(lib300.intronlength[var]), inplace=True)
            rna.loc[var,'numberofcryptic']=len(vals[vals>rna.loc[var,'rnareads']*0.1])
            crypticreads=[]
            crypticintronlengths=[]
            for i in range(int(rna.loc[var,'numberofcryptic'])):
                crypticreads.append(introns[introns.intlen==vals[vals>rna.loc[var,'rnareads']*0.1].index[i]].read.value_counts().index[0])
                crypticintronlengths.append(str(vals.index[i]))
            rna.loc[var,'cryptic_read']=' '.join(crypticreads)
            rna.loc[var,'cryptic_intronlength']=' '.join(crypticintronlengths)

    rna['fraction_spliced_total']=rna.index.map(lambda x: rna.spliced[x]/float(rna.rnareads[x]))
    rna['fraction_unspliced_total']=rna.index.map(lambda x: rna.unspliced[x]/float(rna.rnareads[x]))
    rna['fraction_spliced_canonical']=rna.index.map(lambda x: rna.spliced[x]/float(rna.unspliced[x] + rna.spliced[x]))
    rna['levelratiospl_unbiased']=rna.index.map(lambda x: np.log2(rna.spliced[x]/float(rna.unspliced[x])))
    rna['levelratiospl_againstallothers']=rna.index.map(lambda x: np.log2(rna.spliced[x]/float(rna.rnareads[x] - rna.spliced[x])))
    return rna


def unbiased_mapping_cas(rnareads):    
    lib300=pd.read_pickle('./dataframes/caslib.pkl')
    rna=pd.DataFrame()
    for var in lib300[lib300.subset!='cassettemelanomavariants'].index:
        varsq='ACTAGTTTACGACGGGTTGGGGATGGCGTGCAGCGCAACCACGAGACGGCCTTCCAAGGTAAGGGGGTTCATTAATCGCCAAGGCCTCACTCCCTTTTTTCCATCTCTCCCCGGACTCACCCGCCAAGGGTGGGTTGGAAACCGAAACGAGTCAGTGTTGAAACGTGTCTCATCCTATTCCTGAAGCCAGAATATTCTGGCCATGAGTCATTGTTTCCGCCCATCTTGATTCTTTTGGAAATGGCAGCTCTTGTTCAAAGACCGGAAAGGGTGGGATGTCAAGACGTC' + \
            lib300.varseq162[var] + 'GGCGCGCCtTGGAGTGGAAGTAGAATGAAGGATTTTTTTTAGAGAGGTGGGGATATCTAAAGGTTTTTATGACGCACGGCTGTTTGCAGGCTCTAACTAAAGGACCATTGTTTATTTGATGTTGATTTAAGTAGTGGATCCTTAGAGATAGTGGTATGGCGGTCTTGAATTGTATCAAAAATCTTGGTTTTCTCTAGGCAATTTTTTGTTCCAATTCAGTTGAATACTCTTCAGTGGATTCAAACCATGAAAAAATAAGTCACCAGGGGAGGATAGCTGAAATAATTCCTAAGGCGGTGCCTGTTTTAATGGAGAAGATATGGGGTGGAGCCTGCGTTTTAAACAAACCCAGATCTGATGCAGGATGTACTTAACTACGTTGAGAAAAACTGATCTGCGCAATTGAGGCGTTACTGAAATATTAGGTGGTGGAGATTTGAGAATAAGGGTTTTCGTCTTTTACCTCATGGGAACTCTGGAAGTCCTTTTGTTAGGATAAATCCTAATAAGACCAAGATAGTACTGTAAAATGAAGTTTAATTATCATGGGTCCCCGCTTAAGAAACTGAAGAACTTATTTTCTTTTTTTGCCCCGGGGTGAATAATAATTGGTTTACTATTGCTTTAGGGGGAAACCTTAGATATTTTAATTTACCTTCTCTCTGGATAGTAGTGTTGTAAGAGAGCAGAAACCCATACTTGAAAATGTGCTTTTCTTTTTTGTTTTCTAGGATGGGTTTGTGGAGTTCTTCCATGTAGAGGACCTAGAAGGTGGCATCAGGAATGTGCTGCTGGCTTTTGCAGGTGTTGCTGGAGTAGGAGCTGGTTTGGCATCTAGA'
        intronstarts=[]
        intronends=[]
        reads=[]
        rr=pd.Series(rnareads.loc[var].split(' '))
        readpos=0
        for test in rr:
            vspos=0
            test=str(Seq(test).reverse_complement())
            startpos=100
            testsub=test[startpos:120]
            while True:
                vspos=varsq.find(testsub)
                startpos-=1
                testsub=test[startpos:120]
                if (varsq.find(testsub)==-1)|(startpos<0):
                    if (startpos<0):
                        vspos=0
                    break
            intronends.append(vspos)
            if (vspos==0)|(varsq.find(test[startpos-9:startpos+1])==-1):
                intronstarts.append(0)
            elif (vspos==-1):
                intronstarts.append(-150)
            else:
                intronstarts.append(varsq.find(test[startpos-9:startpos+1])+10)
            reads.append(test)
            readpos+=1
        
        introns=pd.DataFrame([intronstarts,intronends, reads]).transpose()
        introns.columns=['intronstart','intronend', 'read']
        
        introns['intlen']=introns.index.map(lambda x: introns.intronend[x]-introns.intronstart[x])

        if (len(introns)>0):    
            rna.loc[var,'rnareads']=len(introns)
            rna.loc[var,'unspliced']=len(introns[introns.intlen==1108])
            rna.loc[var,'spliced']=len(introns[introns.intlen==761])
            rna.loc[var,'spliced_downstreamdonor2']=len(introns[introns.intlen==685])
            rna.loc[var,'unspliced_downstreamacceptor2']=len(introns[introns.intlen==1181])
            rna.loc[var,'fraction_canonical']=(rna.loc[var,'spliced']+rna.loc[var,'unspliced'])/float(rna.loc[var,'rnareads'])
            rna.loc[var,'fraction_known']=(rna.loc[var,'spliced']+rna.loc[var,'unspliced'] + \
                   rna.loc[var,'spliced_downstreamdonor2']+rna.loc[var,'unspliced_downstreamacceptor2'])/float(rna.loc[var,'rnareads'])
            vals=introns.intlen.value_counts()
            for knownlength in [1108,761,685,1181]:
                if (knownlength in vals.index):
                    vals.drop(knownlength, inplace=True)
            rna.loc[var,'numberofcryptic']=len(vals[vals>rna.loc[var,'rnareads']*0.1])
            crypticreads=[]
            crypticintronlengths=[]
            for i in range(int(rna.loc[var,'numberofcryptic'])):
                crypticreads.append(introns[introns.intlen==vals[vals>rna.loc[var,'rnareads']*0.1].index[i]].read.value_counts().index[0])
                crypticintronlengths.append(str(vals.index[i]))
            rna.loc[var,'cryptic_read']=' '.join(crypticreads)
            rna.loc[var,'cryptic_intronlength']=' '.join(crypticintronlengths)

    rna['fraction_spliced_total']=rna.index.map(lambda x: rna.spliced[x]/float(rna.rnareads[x]))
    rna['fraction_unspliced_total']=rna.index.map(lambda x: rna.unspliced[x]/float(rna.rnareads[x]))
    rna['fraction_spliced_canonical']=rna.index.map(lambda x: rna.spliced[x]/float(rna.unspliced[x] + rna.spliced[x]))
    rna['levelratiospl_unbiased']=rna.index.map(lambda x: np.log2(rna.spliced[x]/float(rna.unspliced[x])))
    rna['levelratiospl_againstallothers']=rna.index.map(lambda x: np.log2(rna.spliced[x]/float(rna.rnareads[x] - rna.spliced[x])))
    return rna


        
def prepare_rnadf_ir(df_to_process):
    lib300=pd.read_pickle('./dataframes/ir_lib300.pkl')
    df_to_process=df_to_process.join(lib300[lib300.subset!='irmelanomavariants'], rsuffix='_from_lib300')
    df_to_process['rnareads_effective']=df_to_process.index.map(lambda x: df_to_process.spliced[x] + df_to_process.unspliced[x])
    
    for x in df_to_process[(df_to_process.rnareads_effective>100)].index:
        if (df_to_process.spliced[x]==0)|(df_to_process.unspliced[x]==0):
            df_to_process.loc[x,'levelratiospl']=np.log2((df_to_process.spliced[x]+1)/(df_to_process.unspliced[x]+1))
    
        else:
            df_to_process.loc[x,'levelratiospl']=np.log2(df_to_process.spliced[x]/df_to_process.unspliced[x])
    
    df_to_process['wtratio']=df_to_process.name2.apply(lambda x: df_to_process[(df_to_process.subset=='ir_filtered') & \
            (df_to_process.name2 == x)].levelratiospl.mean())
    df_to_process['normratio'] = df_to_process.index.map(lambda x: df_to_process.loc[x, 'levelratiospl'] - df_to_process.loc[x, 'wtratio'])
    return df_to_process

def prepare_rnadf_cas(df_to_process):
    lib300=pd.read_pickle('./dataframes/caslib.pkl')
    lib300.drop('levelratiospl', axis=1, inplace=True)
    df_to_process=df_to_process.join(lib300[lib300.subset!='cassettemelanomavariants'], rsuffix='_fromcaslib')
    df_to_process['rnareads_effective']=df_to_process.index.map(lambda x: df_to_process.spliced[x] + df_to_process.unspliced[x])
    
    for x in df_to_process[(df_to_process.rnareads_effective>100)].index:
        if (df_to_process.spliced[x]==0)|(df_to_process.unspliced[x]==0):
            df_to_process.loc[x,'levelratiospl']=np.log2((df_to_process.spliced[x]+1)/(df_to_process.unspliced[x]+1))
    
        else:
            df_to_process.loc[x,'levelratiospl']=np.log2(df_to_process.spliced[x]/df_to_process.unspliced[x])
    
    df_to_process['wtratio']=df_to_process.commonname.apply(lambda x: df_to_process[(df_to_process.subset=='cassette_filtered') & \
            (df_to_process.commonname == x)].levelratiospl.mean())
    df_to_process['normratio'] = df_to_process.index.map(lambda x: df_to_process.loc[x, 'levelratiospl'] - df_to_process.loc[x, 'wtratio'])
    return df_to_process

def prepare_rnadf_five(df_to_process):
    lib300=pd.read_pickle('/net/mraid08/export/genie/Runs/Martin/fivebins/fivelib.pkl')
    df_to_process=df_to_process.join(lib300[lib300.subset!='fivemelanomavariants'], rsuffix='_fromfivelib')
    df_to_process=df_to_process.drop('levelratiospl',axis=1)
    df_to_process['rnareads_effective']=df_to_process.index.map(lambda x: df_to_process.spliced1[x] + df_to_process.spliced2[x])
    
    for x in df_to_process[(df_to_process.rnareads_effective>100)].index:
        if (df_to_process.spliced1[x]==0)|(df_to_process.spliced2[x]==0):
            df_to_process.loc[x,'levelratiospl']=np.log2((df_to_process.spliced2[x]+1)/(df_to_process.spliced1[x]+1))
    
        else:
            df_to_process.loc[x,'levelratiospl']=np.log2(df_to_process.spliced2[x]/df_to_process.spliced1[x])
    
    df_to_process['wtratio']=df_to_process.commonname.apply(lambda x: df_to_process[(df_to_process.subset=='five_filtered') & \
            (df_to_process.commonname == x)].levelratiospl.mean())
    df_to_process['normratio'] = df_to_process.index.map(lambda x: df_to_process.loc[x, 'levelratiospl'] - df_to_process.loc[x, 'wtratio'])
    return df_to_process

def prepare_rnadf_three(df_to_process):
    lib300=pd.read_pickle('/net/mraid08/export/genie/Runs/Martin/threerna/threelib.pkl')
    df_to_process=df_to_process.join(lib300[lib300.subset!='threemelanomavariants'], rsuffix='_fromthreelib')
    df_to_process['rnareads_effective']=df_to_process.index.map(lambda x: df_to_process.spliced1[x] + df_to_process.spliced2[x])
    
    for x in df_to_process[(df_to_process.rnareads_effective>100)].index:
        if (df_to_process.spliced1[x]==0)|(df_to_process.spliced2[x]==0):
            df_to_process.loc[x,'levelratiospl']=np.log2((df_to_process.spliced2[x]+1)/(df_to_process.spliced1[x]+1))
    
        else:
            df_to_process.loc[x,'levelratiospl']=np.log2(df_to_process.spliced2[x]/df_to_process.spliced1[x])
    
    df_to_process['wtratio']=df_to_process.commonname.apply(lambda x: df_to_process[(df_to_process.subset=='three_filtered') & \
            (df_to_process.commonname == x)].levelratiospl.mean())
    df_to_process['normratio'] = df_to_process.index.map(lambda x: df_to_process.loc[x, 'levelratiospl'] - df_to_process.loc[x, 'wtratio'])
    return df_to_process

def prepare_rnadf_umis_ir(df_to_process):
    df_to_process['unsplicedumis_total']=df_to_process.umis_unspliced.apply(lambda x: len(pd.Series(x.split(' '))))
    df_to_process['unsplicedumis_unique']=df_to_process.umis_unspliced.apply(lambda x: len(pd.Series(x.split(' ')).value_counts()))
    df_to_process['unsplicedumis_unique_per_total']=df_to_process.index.map(lambda x: float(df_to_process.unsplicedumis_unique[x])/df_to_process.unsplicedumis_total[x])
    df_to_process['unsplicedumis_uniquef']=df_to_process.umis_unspliced.apply(lambda x: len(pd.Series(x.split(' ')).value_counts()[pd.Series(x.split(' ')).value_counts()>0.2*np.max(pd.Series(x.split(' ')).value_counts())]))
    df_to_process['unsplicedumis_uniquef_per_total']=df_to_process.index.map(lambda x: float(df_to_process.unsplicedumis_uniquef[x])/df_to_process.unsplicedumis_total[x])
    df_to_process['unsplicedumis_uniquemin']=df_to_process.umis_unspliced.apply(lambda x: len(pd.Series(x.split(' ')).value_counts()[pd.Series(x.split(' ')).value_counts()>1]))
    
    df_to_process['splicedumis_total']=df_to_process.umis_spliced.apply(lambda x: len(pd.Series(x.split(' '))))
    df_to_process['splicedumis_unique']=df_to_process.umis_spliced.apply(lambda x: len(pd.Series(x.split(' ')).value_counts()))
    df_to_process['splicedumis_unique_per_total']=df_to_process.index.map(lambda x: float(df_to_process.splicedumis_unique[x])/df_to_process.splicedumis_total[x])
    df_to_process['splicedumis_uniquef']=df_to_process.umis_spliced.apply(lambda x: len(pd.Series(x.split(' ')).value_counts()[pd.Series(x.split(' ')).value_counts()>0.2*np.max(pd.Series(x.split(' ')).value_counts())]))
    df_to_process['splicedumis_uniquef_per_total']=df_to_process.index.map(lambda x: float(df_to_process.splicedumis_uniquef[x])/df_to_process.splicedumis_total[x])
    df_to_process['splicedumis_uniquemin']=df_to_process.umis_spliced.apply(lambda x: len(pd.Series(x.split(' ')).value_counts()[pd.Series(x.split(' ')).value_counts()>1]))
    
    df_to_process['levelratiospl_unique']=df_to_process.index.map(lambda x: np.log2(df_to_process.splicedumis_unique[x]/float(df_to_process.unsplicedumis_unique[x])))
    df_to_process['levelratiospl_uniquef']=df_to_process.index.map(lambda x: np.log2(df_to_process.splicedumis_uniquef[x]/float(df_to_process.unsplicedumis_uniquef[x])))
    df_to_process['levelratiospl_fractionunique']=df_to_process.index.map(lambda x: np.log2(df_to_process.splicedumis_unique_per_total[x]/float(df_to_process.unsplicedumis_unique_per_total[x])))
    df_to_process['levelratiospl_total']=df_to_process.index.map(lambda x: np.log2(df_to_process.splicedumis_total[x]/float(df_to_process.unsplicedumis_total[x])))
    df_to_process['levelratiospl_min']=df_to_process.index.map(lambda x: np.log2(df_to_process.splicedumis_uniquemin[x]/float(df_to_process.unsplicedumis_uniquemin[x])))
    return df_to_process

def unify_dfs_wide(foldername, lib300):
    files=os.listdir(foldername)
    df=lib300
    for filename in files:
        if ('_from_unbiased_mapping_umis.pkl' in filename):
            conditiondf=pd.read_pickle(foldername + filename)
            if (filename.split('_')[0]=='control'):
                suffix='_from_control'
            else:
                suffix='_from_'+ '_'.join(np.array(filename.split('_'))[[0,2,4]])
            df=df.join(conditiondf[[x for x in conditiondf.columns if ('levelratiospl' in x)]], rsuffix=suffix)
    return df
    

def unify_dfs_long_ir(foldername):
    lib300=pd.read_pickle('/net/mraid08/export/genie/Runs/Martin/IRbins300/lib300.pkl')
    files=os.listdir(foldername)
    df=pd.DataFrame()
    for filename in files:
        if ('_from_unbiased_mapping_umis.pkl' in filename):
            conditiondf=pd.read_pickle(foldername + filename)
            conditiondf['libindex']=conditiondf.index
            conditiondf=conditiondf.join(lib300[lib300.subset!='irmelanomavariants'], rsuffix='_fromlib300')
            if (filename.split('_')[0]=='control'):
                conditiondf.loc[:,'exp']='control'
                conditiondf.loc[:,'rep']='control'
                conditiondf.loc[:,'condition']='control'
            else:
                conditiondf.loc[:,'exp']=filename.split('_')[0]
                conditiondf.loc[:,'rep']=filename.split('_')[2]
                conditiondf.loc[:,'condition']=filename.split('_')[4]    
            df=pd.concat([df,conditiondf], ignore_index=True)
    return df

def unify_dfs_long_cas(foldername):
    lib300=pd.read_pickle('/net/mraid08/export/genie/Runs/Martin/casrna/caslib.pkl')
    files=os.listdir(foldername)
    df=pd.DataFrame()
    for filename in files:
        if ('_from_unbiased_mapping_umis.pkl' in filename):
            conditiondf=pd.read_pickle(foldername + filename)
            conditiondf['libindex']=conditiondf.index
            conditiondf=conditiondf.join(lib300[lib300.subset!='cassettemelanomavariants'], rsuffix='_fromlib300')
            if (filename.split('_')[0]=='control'):
                conditiondf.loc[:,'exp']='control'
                conditiondf.loc[:,'rep']='control'
                conditiondf.loc[:,'condition']='control'
            else:
                conditiondf.loc[:,'exp']=filename.split('_')[0]
                conditiondf.loc[:,'rep']=filename.split('_')[2]
                conditiondf.loc[:,'condition']=filename.split('_')[4]    
            df=pd.concat([df,conditiondf], ignore_index=True)
    return df
    
def add_pearsonr(x,y, **kws):
    r, _ = pearsonr(x,y)
    ax=plt.gca()
    ax.annotate("r={:.4f}".format(r),xy=(-10,8))
    
