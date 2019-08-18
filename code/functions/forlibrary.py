# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 15:47:16 2016

@author: martinm
"""

import pandas as pd



def make_combinatorial(df1,df2):
    '''
    
    '''
    
    comb=pd.DataFrame(columns=list(df1.columns))
    k=0
    for i in df1.index:
        exonlength=int(df1.exonlength[i])
        altintronup=df1.varseq[i][:-30-exonlength]
        altexon=df1.varseq[i][-30-exonlength:-30]
        altintrondown=df1.varseq[i][-30:]

        comb=comb.append(df1.loc[i],ignore_index=True)
        comb.loc[k,'combination']=str('alt' + '\t' + 'alt' + '\t' + 'alt')
        k=k+1

        comb=comb.append(df1.loc[i],ignore_index=True)
        comb.loc[k,'combination']=str('alt' + '\t' + 'alt' + '\t' + 'alt')
        k=k+1

        comb=comb.append(df1.loc[i],ignore_index=True)
        comb.loc[k,'combination']=str('alt' + '\t' + 'alt' + '\t' + 'alt')
        k=k+1

        for j in df2.index:
            
            comb=comb.append(df1.loc[i],ignore_index=True)
            comb.varseq[k]=altintronup + altexon +  df2.introndown[j] 
            comb.loc[k,'combination']=str('alt' + '\t' + 'alt' + '\t' + str(j))
            k=k+1
            
            comb=comb.append(df1.loc[i],ignore_index=True)
            comb.varseq[k]=altintronup + df2.exon[j] +  df2.introndown[j] 
            comb.loc[k,'combination']='alt' + '\t' + str(j) + '\t' + str(j)
            k=k+1
            
            comb=comb.append(df1.loc[i],ignore_index=True)
            comb.varseq[k]=df2.intronup[j] + df2.exon[j] +  df2.introndown[j] 
            comb.loc[k,'combination']=str(j) + '\t' + str(j) + '\t' + str(j)
            k=k+1
            
            comb=comb.append(df1.loc[i],ignore_index=True)
            comb.varseq[k]=altintronup + df2.exon[j] +  altintrondown
            comb.loc[k,'combination']='alt' + '\t' + str(j) + '\t' + 'alt'
            k=k+1
            
            comb=comb.append(df1.loc[i],ignore_index=True)
            comb.varseq[k]=df2.intronup[j] + df2.exon[j] +  altintrondown
            comb.loc[k,'combination']=str(j) + '\t' + str(j) + '\t' + 'alt'
            k=k+1
            
            comb=comb.append(df1.loc[i],ignore_index=True)
            comb.varseq[k]=df2.intronup[j] + altexon +  altintrondown
            comb.loc[k,'combination']=str(j) + '\t' + 'alt' + '\t' + 'alt'
            k=k+1
            
            comb=comb.append(df1.loc[i],ignore_index=True)
            comb.varseq[k]=df2.intronup[j] + altexon +  df2.introndown[j] 
            comb.loc[k,'combination']=str(j) + '\t' + 'alt' + '\t' + str(j)
            k=k+1
            
    
    return comb


    
def make_combinatorialintron(df1,df2):
    '''
    
    '''
    
    comb=pd.DataFrame(columns=list(df1.columns))
    k=0
    for i in df1.index:
        intronlength=int(df1.intronlength[i])
        altintronup=df1.varseq[i][:int(df1.intronstart_varseq[i])]
        altexon=df1.varseq[i][int(df1.intronstart_varseq[i]):int(df1.intronend_varseqnew[i])]
        altintrondown=df1.varseq[i][int(df1.intronend_varseqnew[i]):]

        comb=comb.append(df1.loc[i],ignore_index=True)
        comb.loc[k,'combination']=str('alt' + '\t' + 'alt' + '\t' + 'alt')
        k=k+1

        comb=comb.append(df1.loc[i],ignore_index=True)
        comb.loc[k,'combination']=str('alt' + '\t' + 'alt' + '\t' + 'alt')
        k=k+1

        comb=comb.append(df1.loc[i],ignore_index=True)
        comb.loc[k,'combination']=str('alt' + '\t' + 'alt' + '\t' + 'alt')
        k=k+1
        counter=0
        for j in df2.index:
            if (df2.intronlength[j]==intronlength)&(counter<5):
                spliced=df2.sequence[j][150-intronlength-int(df1.intronstart_varseq[i]):150-intronlength] + df2.sequence[j][150:150+162-int(df1.intronend_varseqnew[i])]
                if (spliced.translate().find("*")==-1):
                    comb=comb.append(df1.loc[i],ignore_index=True)
                    comb.varseq[k]=altintronup + altexon +  df2.sequence[j][150:150+162-int(df1.intronend_varseqnew[i])] 
                    comb.loc[k,'combination']=str('alt' + '\t' + 'alt' + '\t' + str(j))
                    k=k+1
                    
                    comb=comb.append(df1.loc[i],ignore_index=True)
                    comb.varseq[k]=altintronup + df2.sequence[j][150-intronlength:150]  +  df2.sequence[j][150:150+162-int(df1.intronend_varseqnew[i])] 
                    comb.loc[k,'combination']='alt' + '\t' + str(j) + '\t' + str(j)
                    k=k+1
                    
                    comb=comb.append(df1.loc[i],ignore_index=True)
                    comb.varseq[k]=df2.sequence[j][150-intronlength-int(df1.intronstart_varseq[i]):150-intronlength] + df2.sequence[j][150-intronlength:150] +  df2.sequence[j][150:150+162-int(df1.intronend_varseqnew[i])] 
                    comb.loc[k,'combination']=str(j) + '\t' + str(j) + '\t' + str(j)
                    k=k+1
                    
                    comb=comb.append(df1.loc[i],ignore_index=True)
                    comb.varseq[k]=altintronup + df2.sequence[j][150-intronlength:150] +  altintrondown
                    comb.loc[k,'combination']='alt' + '\t' + str(j) + '\t' + 'alt'
                    k=k+1
                    
                    comb=comb.append(df1.loc[i],ignore_index=True)
                    comb.varseq[k]=df2.sequence[j][150-intronlength-int(df1.intronstart_varseq[i]):150-intronlength] + df2.sequence[j][150-intronlength:150] +  altintrondown
                    comb.loc[k,'combination']=str(j) + '\t' + str(j) + '\t' + 'alt'
                    k=k+1
                    
                    comb=comb.append(df1.loc[i],ignore_index=True)
                    comb.varseq[k]=df2.sequence[j][150-intronlength-int(df1.intronstart_varseq[i]):150-intronlength] + altexon +  altintrondown
                    comb.loc[k,'combination']=str(j) + '\t' + 'alt' + '\t' + 'alt'
                    k=k+1
                    
                    comb=comb.append(df1.loc[i],ignore_index=True)
                    comb.varseq[k]=df2.sequence[j][150-intronlength-int(df1.intronstart_varseq[i]):150-intronlength] + altexon +  df2.sequence[j][150:150+162-int(df1.intronend_varseqnew[i])] 
                    comb.loc[k,'combination']=str(j) + '\t' + 'alt' + '\t' + str(j)
                    k=k+1
                    
                    counter=counter+1
    return comb

   
def make_combinatorialintron_alt(df1,df2):
    '''
    
    '''
    
    comb=pd.DataFrame(columns=list(df1.columns))
    k=0
    for i in df1.index:
        altintronup=df1.varseq[i][:int(df1.intronstart_varseq[i])]
        altexon=df1.varseq[i][int(df1.intronstart_varseq[i]):int(df1.intronend_varseqnew[i])]
        altintrondown=df1.varseq[i][int(df1.intronend_varseqnew[i]):]

        for j in df2.index:
            comb=comb.append(df1.loc[i],ignore_index=True)
            comb.varseq[k]=altintronup + altexon +  df2.varseq[i][int(df1.intronend_varseqnew[i]):]
            comb.loc[k,'combination']=str('alt' + '\t' + 'alt' + '\t' + str(j))
            k=k+1
            
            comb=comb.append(df1.loc[i],ignore_index=True)
            comb.varseq[k]=altintronup + df2.varseq[i][int(df1.intronstart_varseq[i]):int(df1.intronend_varseqnew[i])]  +  df2.varseq[i][int(df1.intronend_varseqnew[i]):]
            comb.loc[k,'combination']='alt' + '\t' + str(j) + '\t' + str(j)
            k=k+1
            
            comb=comb.append(df1.loc[i],ignore_index=True)
            comb.varseq[k]=df2.varseq[i][:int(df1.intronstart_varseq[i])] + df2.varseq[i][int(df1.intronstart_varseq[i]):int(df1.intronend_varseqnew[i])] +  df2.varseq[i][int(df1.intronend_varseqnew[i]):]
            comb.loc[k,'combination']=str(j) + '\t' + str(j) + '\t' + str(j)
            k=k+1
            
            comb=comb.append(df1.loc[i],ignore_index=True)
            comb.varseq[k]=altintronup + df2.varseq[i][int(df1.intronstart_varseq[i]):int(df1.intronend_varseqnew[i])] +  altintrondown
            comb.loc[k,'combination']='alt' + '\t' + str(j) + '\t' + 'alt'
            k=k+1
            
            comb=comb.append(df1.loc[i],ignore_index=True)
            comb.varseq[k]=df2.varseq[i][:int(df1.intronstart_varseq[i])] + df2.varseq[i][int(df1.intronstart_varseq[i]):int(df1.intronend_varseqnew[i])] +  altintrondown
            comb.loc[k,'combination']=str(j) + '\t' + str(j) + '\t' + 'alt'
            k=k+1
            
            comb=comb.append(df1.loc[i],ignore_index=True)
            comb.varseq[k]=df2.varseq[i][:int(df1.intronstart_varseq[i])] + altexon +  altintrondown
            comb.loc[k,'combination']=str(j) + '\t' + 'alt' + '\t' + 'alt'
            k=k+1
            
            comb=comb.append(df1.loc[i],ignore_index=True)
            comb.varseq[k]=df2.varseq[i][:int(df1.intronstart_varseq[i])] + altexon +  df2.varseq[i][int(df1.intronend_varseqnew[i]):]
            comb.loc[k,'combination']=str(j) + '\t' + 'alt' + '\t' + str(j)
            k=k+1
            
    return comb