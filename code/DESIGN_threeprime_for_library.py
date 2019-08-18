# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 19:02:47 2016

@author: miklm
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq
#%%


three=pd.read_pickle('./design/jun3prselectedNEW.pkl')



three[['score','score2','score_alt','level','level_alt']]= \
    three[['score','score2','score_alt','level','level_alt']].astype(int)

for i in three.index:
    three.loc[i,'levelratio']=np.log2(float(three.level_alt[i])/float(three.level[i]))

three_filtered=three[(three.diff_nt < 60) & (three.level_alt>20) & (three.level>20) & (three.levelratio>-5) & (three.levelratio<5)]

three_filtered.drop_duplicates('commonname',inplace=True)

for i in three_filtered.index:
    three_filtered.loc[i,'varseq']=three_filtered.sequence[i][71:221]

three_filtered.to_pickle(homedir + 'dataframes/three_filtered.pkl')

three_filtered['acceptor1']=79
three_filtered['acceptor2']=79+three_filtered.diff_nt.astype(int)

three_filtered[['commonname','coordinates','end','end_alt','diff_nt','acceptor1','acceptor2','varseq']].to_csv('./tables/TableS4_threefiltered.csv')


#%%

'''
make 4 groups, within which "exonic - alternative - intronic" sequences are
varied by permutation

The four groups are: number of genes in the group - diff_nt
5x 8
8x 20
7x 26
4x 38
'''

three5all=three_filtered[three_filtered.diff_nt==5]
three26=three_filtered[three_filtered.diff_nt==26]
three38=three_filtered[three_filtered.diff_nt==38]
three44=three_filtered[three_filtered.diff_nt==44]
three20=three_filtered[three_filtered.diff_nt==20]
three20.drop([348443],inplace=True)
three5=three5all.loc[[94383,216874,569557,356837,142122,317803,669635]]
three26.drop([681753,192245],inplace=True)
three32.drop([349526],inplace=True)
three38.drop([79937],inplace=True)


df=three5.copy()
for i in df.index:
    df.loc[i,'intronseq']=df.varseq[i][:79]
    df.loc[i,'altseq']=df.varseq[i][79:79+ int(df.diff_nt[i])]
    df.loc[i,'exonseq']=df.varseq[i][79+ int(df.diff_nt[i]):]
three5=df.copy()

df=three26.copy()
for i in df.index:
    df.loc[i,'intronseq']=df.varseq[i][:79]
    df.loc[i,'altseq']=df.varseq[i][79:79+ int(df.diff_nt[i])]
    df.loc[i,'exonseq']=df.varseq[i][79+ int(df.diff_nt[i]):]
three26=df.copy()

df=three20.copy()
for i in df.index:
    df.loc[i,'intronseq']=df.varseq[i][:79]
    df.loc[i,'altseq']=df.varseq[i][79:79+ int(df.diff_nt[i])]
    df.loc[i,'exonseq']=df.varseq[i][79+ int(df.diff_nt[i]):]
three20=df.copy()

df=three38.copy()
for i in df.index:
    df.loc[i,'intronseq']=df.varseq[i][:79]
    df.loc[i,'altseq']=df.varseq[i][79:79+ int(df.diff_nt[i])]
    df.loc[i,'exonseq']=df.varseq[i][79+ int(df.diff_nt[i]):]
three38=df.copy()

df=three44.copy()
for i in df.index:
    df.loc[i,'intronseq']=df.varseq[i][:79]
    df.loc[i,'altseq']=df.varseq[i][79:79+ int(df.diff_nt[i])]
    df.loc[i,'exonseq']=df.varseq[i][79+ int(df.diff_nt[i]):]
three44=df.copy()


#%%

# create all possible combinations within each group
df=three5.copy()
out=pd.DataFrame()
l=0
for i in df.index:
    for j in df.index:
        for k in df.index:
            out.loc[l,'varseq']=df.loc[i,'intronseq']+df.loc[j,'altseq']+df.loc[k,'exonseq']
            out.loc[l,'geneintronseq']=df.commonname[i]
            out.loc[l,'genealtseq']=df.commonname[j]
            out.loc[l,'geneexonseq']=df.commonname[k]
            out.loc[l,'diff_nt']=df.diff_nt[j]
            l=l+1
three5combvar=out.copy()

df=three26.copy()
out=pd.DataFrame()
l=0
for i in df.index:
    for j in df.index:
        for k in df.index:
            out.loc[l,'varseq']=df.loc[i,'intronseq']+df.loc[j,'altseq']+df.loc[k,'exonseq']
            out.loc[l,'geneintronseq']=df.commonname[i]
            out.loc[l,'genealtseq']=df.commonname[j]
            out.loc[l,'geneexonseq']=df.commonname[k]
            out.loc[l,'diff_nt']=df.diff_nt[j]
            l=l+1
three26combvar=out.copy()

df=three20.copy()
out=pd.DataFrame()
l=0
for i in df.index:
    for j in df.index:
        for k in df.index:
            out.loc[l,'varseq']=df.loc[i,'intronseq']+df.loc[j,'altseq']+df.loc[k,'exonseq']
            out.loc[l,'geneintronseq']=df.commonname[i]
            out.loc[l,'genealtseq']=df.commonname[j]
            out.loc[l,'geneexonseq']=df.commonname[k]
            out.loc[l,'diff_nt']=df.diff_nt[j]
            l=l+1
three20combvar=out.copy()

df=three38.copy()
out=pd.DataFrame()
l=0
for i in df.index:
    for j in df.index:
        for k in df.index:
            out.loc[l,'varseq']=df.loc[i,'intronseq']+df.loc[j,'altseq']+df.loc[k,'exonseq']
            out.loc[l,'geneintronseq']=df.commonname[i]
            out.loc[l,'genealtseq']=df.commonname[j]
            out.loc[l,'geneexonseq']=df.commonname[k]
            out.loc[l,'diff_nt']=df.diff_nt[j]
            l=l+1
three38combvar=out.copy()

df=three44.copy()
out=pd.DataFrame()
l=0
for i in df.index:
    for j in df.index:
        for k in df.index:
            out.loc[l,'varseq']=df.loc[i,'intronseq']+df.loc[j,'altseq']+df.loc[k,'exonseq']
            out.loc[l,'geneintronseq']=df.commonname[i]
            out.loc[l,'genealtseq']=df.commonname[j]
            out.loc[l,'geneexonseq']=df.commonname[k]
            out.loc[l,'diff_nt']=df.diff_nt[j]
            l=l+1
three44combvar=out.copy()

threeprime_combinatorial_variations=pd.concat([three5combvar,three26combvar,three20combvar,three38combvar,three44combvar],ignore_index=True)

for i in threeprime_combinatorial_variations.index:
    threeprime_combinatorial_variations.loc[i,'nostopinspliced']=bool(threeprime_combinatorial_variations.varseq[i][80:].translate().find("*")==-1)

### PICKLE
threeprime_combinatorial_variations.to_pickle('./design/threeprime_combinatorial_variations.pkl')

#%% replace with constitutive splice sites

#select everything larger than 15 diff_nt

threeconst = three_filtered[three_filtered.diff_nt > 20]

#canonical_acceptor = Seq('CTTTCCTCCTTTCCTTTCAGGTC')
canonical_acceptor = Seq('CTCCTTTCCTTTCAGGTC').tomutable()
no_acceptor = Seq('CAGAGAGGA').tomutable()
branchpoint = Seq('CTCAC').tomutable()

threeconstvar=pd.DataFrame(columns=list(threeconst.columns))
j=0
for i in threeconst.index:
    diff=int(threeconst.diff_nt[i])
    threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
    threeconstvar.loc[j,'first_ss']=str('endogenous')
    threeconstvar.loc[j,'second_ss']=str('endogenous')
    j=j+1
    threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
    sequence = threeconst.varseq[i].tomutable()
    sequence[64:82]=canonical_acceptor
    threeconstvar.loc[j,'varseq']=sequence.toseq()
    threeconstvar.loc[j,'first_ss']=str('constitutive')
    threeconstvar.loc[j,'second_ss']=str('endogenous')
    j=j+1
    threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
    sequence = threeconst.varseq[i].tomutable()
    sequence[64+diff:82+diff]=canonical_acceptor
    threeconstvar.loc[j,'varseq']=sequence.toseq()
    threeconstvar.loc[j,'first_ss']=str('endogenous')
    threeconstvar.loc[j,'second_ss']=str('constitutive')
    j=j+1
    threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
    sequence = threeconst.varseq[i].tomutable()
    sequence[64:82]=canonical_acceptor
    sequence[64+diff:82+diff]=canonical_acceptor
    threeconstvar.loc[j,'varseq']=sequence.toseq()
    threeconstvar.loc[j,'first_ss']=str('constitutive')
    threeconstvar.loc[j,'second_ss']=str('constitutive')
    j=j+1
    threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
    sequence = threeconst.varseq[i].tomutable()
    sequence[70:79]=no_acceptor
    threeconstvar.loc[j,'varseq']=sequence.toseq()
    threeconstvar.loc[j,'first_ss']=str('unsplicable')
    threeconstvar.loc[j,'second_ss']=str('endogenous')
    j=j+1
    threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
    sequence = threeconst.varseq[i].tomutable()
    sequence[70+diff:79+diff]=no_acceptor
    threeconstvar.loc[j,'varseq']=sequence.toseq()
    threeconstvar.loc[j,'first_ss']=str('endogenous')
    threeconstvar.loc[j,'second_ss']=str('unsplicable')
    j=j+1


for i in threeconst.index:
    if (threeconst.diff_nt[i]>28):
        diff=int(threeconst.diff_nt[i])
        threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
        sequence = threeconst.varseq[i].tomutable()
        sequence[53:58]=branchpoint
        threeconstvar.loc[j,'varseq']=sequence.toseq()
        threeconstvar.loc[j,'first_ss']=str('constBP')
        threeconstvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
        sequence = threeconst.varseq[i].tomutable()
        sequence[53+diff:58+diff]=branchpoint
        threeconstvar.loc[j,'varseq']=sequence.toseq()
        threeconstvar.loc[j,'first_ss']=str('endogenous')
        threeconstvar.loc[j,'second_ss']=str('constBP')
        j=j+1 
        threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
        sequence = threeconst.varseq[i].tomutable()
        sequence[53:58]=branchpoint
        sequence[53+diff:58+diff]=branchpoint
        threeconstvar.loc[j,'varseq']=sequence.toseq()
        threeconstvar.loc[j,'first_ss']=str('constBP')
        threeconstvar.loc[j,'second_ss']=str('constBP')
        j=j+1 
    
   
threeconst = three_filtered[(three_filtered.diff_nt > 15)&(three_filtered.diff_nt < 21)]

canonical_acceptor = Seq('TTCCTTTCAGGTC').tomutable()


for i in threeconst.index:
        diff=int(threeconst.diff_nt[i])
        threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
        threeconstvar.loc[j,'first_ss']=str('endogenous')
        threeconstvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
        sequence = threeconst.varseq[i].tomutable()
        sequence[69:82]=canonical_acceptor
        threeconstvar.loc[j,'varseq']=sequence.toseq()
        threeconstvar.loc[j,'first_ss']=str('constitutive')
        threeconstvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
        sequence = threeconst.varseq[i].tomutable()
        sequence[69+diff:82+diff]=canonical_acceptor
        threeconstvar.loc[j,'varseq']=sequence.toseq()
        threeconstvar.loc[j,'first_ss']=str('endogenous')
        threeconstvar.loc[j,'second_ss']=str('constitutive')
        j=j+1
        threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
        sequence = threeconst.varseq[i].tomutable()
        sequence[69:82]=canonical_acceptor
        sequence[69+diff:82+diff]=canonical_acceptor
        threeconstvar.loc[j,'varseq']=sequence.toseq()
        threeconstvar.loc[j,'first_ss']=str('constitutive')
        threeconstvar.loc[j,'second_ss']=str('constitutive')
        j=j+1
        threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
        sequence = threeconst.varseq[i].tomutable()
        sequence[70:79]=no_acceptor
        threeconstvar.loc[j,'varseq']=sequence.toseq()
        threeconstvar.loc[j,'first_ss']=str('unsplicable')
        threeconstvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        threeconstvar=threeconstvar.append(threeconst.loc[i], ignore_index=True)
        sequence = threeconst.varseq[i].tomutable()
        sequence[70+diff:79+diff]=no_acceptor
        threeconstvar.loc[j,'varseq']=sequence.toseq()
        threeconstvar.loc[j,'first_ss']=str('endogenous')
        threeconstvar.loc[j,'second_ss']=str('unsplicable')
        j=j+1


### PICKLE

threeconstvar.to_pickle('./design/threeprime_constitutiveandunsplicable_variations.pkl')


#%% SF binding sites


#%%
threeSF=three_filtered.loc[[553047,227218,55191,257190,718514,305051,124697,555659,350853]]


'''
for binding site activity and positioning:				
				potential creation of stop
GROUP1

SRSF1	CACACGA	enh	TCACACGAC	
SRSF2	GGCCTCTG	enh	TGGCCTCTG	(plus2)
SRSF5	TCACAGG	enh	TTCACAGGC	
SRSF6	TGCGTC	enh	CTGCGTCGA	
KeEnh	GACGTC	enh plus 5 to plus 10 in the exon	CGACGTCGA	(plus1)
KeSil	CCAGCA	sil plus 5 to plus 10 in the exon	CCCAGCAGA	
KeNeu	AAAGAG	neu plus 5 to plus 10 in the exon	CAAAGAGGA	(plus1)
hnRNPA1	TAGGGA	sil	TTAGGGAAC	zero (-> PLUS2) stop for sure
hnRNP G	AAGTGTT	sil	CAAGTGTTC	(plus1,zero)
hnRNP U	TGTATTG	sil	TTGTATTGC	(plus2)

GROUP2

RosenbergGENsil	GTGGGG	sil		
RosenbergE5enh	CACCGC	enh		
RosenbergE5sil	GGTGGG	sil		
RosenbergI5enh	TTGTTC	enh		
RosenbergI5sil	CGAACC	sil		
RosenbergE3enh	CGAAGA	enh		
RosenbergE3sil	GGGGGG	sil		
RosenbergI3enh	TCTAAC	enh		plus2 stop for sure
RosenbergI3sil	CCAAGC	sil		


GROUP1:
relative to first splice site
[-30:-21]
[-21:-12]
[-12:-3]
[6:15]
[15:24] --- not for the group with only 26nt diff_nt
[diff-14:diff-5]
[diff+6:diff+15]
[diff+15:diff+24]
[diff+24:diff+33]
'''

splicingfactors1=pd.Series()
splicingfactors1['SRSF1']=Seq('TCACACGAC').tomutable()
splicingfactors1['SRSF2']=Seq('TGGCCTCTG').tomutable()
splicingfactors1['SRSF5']=Seq('TTCACAGGC').tomutable()
splicingfactors1['SRSF6']=Seq('CTGCGTCGA').tomutable()
splicingfactors1['KeEnh1']=Seq('CAGAAGAGT').tomutable()
splicingfactors1['KeSil1']=Seq('CCCAGCAGT').tomutable()
splicingfactors1['KeNeu1']=Seq('CAAAGAGGT').tomutable()
splicingfactors1['KeEnh2']=Seq('CGAAGATGT').tomutable()
splicingfactors1['KeSil2']=Seq('CCTTTTAGT').tomutable()
splicingfactors1['KeNeu2']=Seq('CAAACTTGT').tomutable()
splicingfactors1['KeEnh3']=Seq('CGCAAGAGT').tomutable()
splicingfactors1['KeSil3']=Seq('CCTAGTAGT').tomutable()
splicingfactors1['KeNeu3']=Seq('CAACCTTGT').tomutable()
splicingfactors1['hnRNPA1']=Seq('TTAGGGAAC').tomutable()
splicingfactors1['hnRNPG']=Seq('CAAGTGTTC').tomutable()
splicingfactors1['hnRNPU']=Seq('TTGTATTGC').tomutable()


threeSFvar=pd.DataFrame(columns=list(threeSF.columns))

j=0
for i in threeSF.index:
    diff=int(threeSF.diff_nt[i])
    for sf in splicingfactors1.index:
        threeSFvar=threeSFvar.append(threeSF.loc[i], ignore_index=True)
        sequence = threeSF.varseq[i].tomutable()
        sequence[21:30]=splicingfactors1[sf]
        threeSFvar.loc[j,'varseq']=sequence.toseq()
        threeSFvar.loc[j,'first_ss']=str(sf + ' [-58:-49]')
        threeSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        threeSFvar=threeSFvar.append(threeSF.loc[i], ignore_index=True)
        sequence = threeSF.varseq[i].tomutable()
        sequence[30:39]=splicingfactors1[sf]
        threeSFvar.loc[j,'varseq']=sequence.toseq()
        threeSFvar.loc[j,'first_ss']=str(sf + ' [-49:--40]')
        threeSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        threeSFvar=threeSFvar.append(threeSF.loc[i], ignore_index=True)
        sequence = threeSF.varseq[i].tomutable()
        sequence[39:48]=splicingfactors1[sf]
        threeSFvar.loc[j,'varseq']=sequence.toseq()
        threeSFvar.loc[j,'first_ss']=str(sf + ' [-40:-31]')
        threeSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        threeSFvar=threeSFvar.append(threeSF.loc[i], ignore_index=True)
        sequence = threeSF.varseq[i].tomutable()
        sequence[83:92]=splicingfactors1[sf]
        threeSFvar.loc[j,'varseq']=sequence.toseq()
        threeSFvar.loc[j,'first_ss']=str(sf + ' [4:13]')
        threeSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        if (diff>30):
            threeSFvar=threeSFvar.append(threeSF.loc[i], ignore_index=True)
            sequence = threeSF.varseq[i].tomutable()
            sequence[92:101]=splicingfactors1[sf]
            threeSFvar.loc[j,'varseq']=sequence.toseq()
            threeSFvar.loc[j,'first_ss']=str(sf + ' [13:22]')
            threeSFvar.loc[j,'second_ss']=str('endogenous')
            j=j+1
        threeSFvar=threeSFvar.append(threeSF.loc[i], ignore_index=True)
        sequence = threeSF.varseq[i].tomutable()
        sequence[79+diff+5:79+diff+14]=splicingfactors1[sf]
        threeSFvar.loc[j,'varseq']=sequence.toseq()
        threeSFvar.loc[j,'first_ss']=str('endogenous')
        threeSFvar.loc[j,'second_ss']=str(sf + ' [diff+3:diff+12]')
        j=j+1
        if (diff < 51):
            threeSFvar=threeSFvar.append(threeSF.loc[i], ignore_index=True)
            sequence = threeSF.varseq[i].tomutable()
            sequence[79+diff+14:79+diff+23]=splicingfactors1[sf]
            threeSFvar.loc[j,'varseq']=sequence.toseq()
            threeSFvar.loc[j,'first_ss']=str('endogenous')
            threeSFvar.loc[j,'second_ss']=str(sf + ' [diff+12:diff+21]')
            j=j+1
        if (diff < 40):
            threeSFvar=threeSFvar.append(threeSF.loc[i], ignore_index=True)
            sequence = threeSF.varseq[i].tomutable()
            sequence[79+diff+23:79+diff+32]=splicingfactors1[sf]
            threeSFvar.loc[j,'varseq']=sequence.toseq()
            threeSFvar.loc[j,'first_ss']=str('endogenous')
            threeSFvar.loc[j,'second_ss']=str(sf + ' [diff+21:diff+30]')
            j=j+1

for i in threeSFvar.index:
    if (threeSFvar.varseq[i][80:].translate().find("*")>-1):
        print(str(i))
        print(threeSFvar.varseq[i][80:].translate())
        
from Bio.Restriction import RsrII
from Bio.Restriction import AscI
from Bio.Restriction import SpeI
from Bio.Restriction import AatII


for i in threeSFvar.index:
    if (RsrII.search(threeSFvar.varseq[i])!=[])|(AscI.search(threeSFvar.varseq[i])!=[])|(SpeI.search(threeSFvar.varseq[i])!=[])|(AatII.search(threeSFvar.varseq[i])!=[]):
        threeSFvar.drop(i,inplace=True)
        

threeSFvar.to_pickle('./design/three_splicingfactors_location.pkl')

#%% 

threeSFcomb=threeSF.loc[(threeSF.diff_nt > 35) & (threeSF.diff_nt < 51)]

threeSFcombvar=pd.DataFrame(columns=list(threeSFcomb.columns))


splicingfactorscomb=pd.Series()
splicingfactorscomb['native']=Seq('NNNNNNNNN').tomutable()
splicingfactorscomb['SRSF1']=Seq('TCACACGAC').tomutable()
splicingfactorscomb['SRSF5']=Seq('TTCACAGGC').tomutable()
splicingfactorscomb['hnRNPA1']=Seq('TTAGGGAAC').tomutable()
splicingfactorscomb['hnRNPU']=Seq('TTGTATTGC').tomutable()

k=0


for i in threeSFcomb.index:
    startpositions=[24,42,83,79 + int(threeSFcomb.diff_nt[i]) + 5]
    for sf in splicingfactorscomb.index: 
        sequence = threeSFcomb.varseq[i].tomutable()   
        if (sf!="native"):
            first=sf
            firstpos=startpositions[0]
            sequence[startpositions[0]:startpositions[0]+9] = splicingfactorscomb.loc[sf]
            for sf in splicingfactorscomb.index:
                if (sf!="native"):
                    sequence[startpositions[1]:startpositions[1]+9] = splicingfactorscomb.loc[sf]
                    second=sf
                    secondpos=startpositions[1]
                    threeSFcombvar=threeSFcombvar.append(threeSFcomb.loc[i], ignore_index=True)        
                    threeSFcombvar.loc[k,'varseq']=sequence.toseq()
                    threeSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                    threeSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                    k=k+1
                    
                else:
                    for sf in splicingfactorscomb.index:
                        if (sf!="native"):
                            sequence[startpositions[2]:startpositions[2]+9] = splicingfactorscomb.loc[sf]
                            second=sf
                            secondpos=startpositions[2]
                            threeSFcombvar=threeSFcombvar.append(threeSFcomb.loc[i], ignore_index=True)        
                            threeSFcombvar.loc[k,'varseq']=sequence.toseq()
                            threeSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                            threeSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                            k=k+1
                    else:
                        for sf in splicingfactorscomb.index:
                            if (sf!="native"):
                                sequence[startpositions[3]:startpositions[3]+9] = splicingfactorscomb.loc[sf]
                                second=sf
                                secondpos=startpositions[3]
                                threeSFcombvar=threeSFcombvar.append(threeSFcomb.loc[i], ignore_index=True)        
                                threeSFcombvar.loc[k,'varseq']=sequence.toseq()
                                threeSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                                threeSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                                k=k+1    
        else:       
            for sf in splicingfactorscomb.index:
                if (sf!="native"):
                    first=sf
                    firstpos=startpositions[1]
                    sequence[startpositions[1]:startpositions[1]+9] = splicingfactorscomb.loc[sf]
                    for sf in splicingfactorscomb.index:
                        if (sf!="native"):
                            second=sf
                            secondpos=startpositions[2]
                            sequence[startpositions[2]:startpositions[2]+9] = splicingfactorscomb.loc[sf]
                            threeSFcombvar=threeSFcombvar.append(threeSFcomb.loc[i], ignore_index=True)        
                            threeSFcombvar.loc[k,'varseq']=sequence.toseq()
                            threeSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                            threeSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                            k=k+1
                        else:
                            for sf in splicingfactorscomb.index:
                                if (sf!="native"):
                                    sequence[startpositions[3]:startpositions[3]+9] = splicingfactorscomb.loc[sf]
                                    second=sf
                                    secondpos=startpositions[3]
                                    threeSFcombvar=threeSFcombvar.append(threeSFcomb.loc[i], ignore_index=True)        
                                    threeSFcombvar.loc[k,'varseq']=sequence.toseq()
                                    threeSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                                    threeSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                                    k=k+1    

                else:       
                    for sf in splicingfactorscomb.index:
                        if (sf!="native"):
                            first=sf
                            firstpos=startpositions[2]
                            sequence[startpositions[2]:startpositions[2]+9] = splicingfactorscomb.loc[sf]
                            for sf in splicingfactorscomb.index:
                                if (sf!="native"):
                                    second=sf
                                    secondpos=startpositions[3]
                                    sequence[startpositions[3]:startpositions[3]+9] = splicingfactorscomb.loc[sf]
                                    threeSFcombvar=threeSFcombvar.append(threeSFcomb.loc[i], ignore_index=True)        
                                    threeSFcombvar.loc[k,'varseq']=sequence.toseq()
                                    threeSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                                    threeSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                                    k=k+1
            
for i in threeSFcombvar.index:
    if (threeSFcombvar.varseq[i][80:].translate().find("*")>-1):
        print(str(i))
        print(threeSFcombvar.varseq[i][80:].translate())


            
            
threeSFcombvar.to_pickle('./design/three_splicingfactors_combinatorial_exononly.pkl')
 
#%% Test Rosenberg features

'''
RosenbergGENsil	GTGGGG	sil		
RosenbergE5enh	CACCGC	enh		
RosenbergE5sil	GGTGGG	sil		
RosenbergI5enh	TTGTTC	enh		
RosenbergI5sil	CGAACC	sil		
RosenbergE3enh	CGAAGA	enh		
RosenbergE3sil	GGGGGG	sil		
RosenbergI3enh	TCTAAC	enh		plus2 stop for sure
RosenbergI3sil	CCAAGC	sil		

'''
  
splicingfactors2=pd.Series()
splicingfactors2['RosenbergGENsil']=Seq('GTGGGG').tomutable()
splicingfactors2['RosenbergE5enh']=Seq('CACCGC').tomutable()
splicingfactors2['RosenbergE5sil']=Seq('GGTGGG').tomutable()
splicingfactors2['RosenbergI5enh']=Seq('TTGTTC').tomutable()
splicingfactors2['RosenbergI5sil']=Seq('CGAACC').tomutable()
splicingfactors2['RosenbergE3enh']=Seq('CGAAGA').tomutable()
splicingfactors2['RosenbergE3sil']=Seq('GGGGGG').tomutable()
splicingfactors2['RosenbergI3enh']=Seq('TCTAAC').tomutable()
splicingfactors2['RosenbergI3sil']=Seq('CCAAGC').tomutable()

threeSFRosenberg=pd.DataFrame(columns=list(threeSFcomb.columns))
j=0
for i in threeSFcomb.index:
    diff=int(threeSFcomb.diff_nt[i])
    for sf in splicingfactors2.index:
        threeSFRosenberg=threeSFRosenberg.append(threeSFcomb.loc[i], ignore_index=True)
        sequence = threeSFcomb.varseq[i].tomutable()
        sequence[24:30]=splicingfactors2[sf]
        threeSFRosenberg.loc[j,'varseq']=sequence.toseq()
        threeSFRosenberg.loc[j,'first_ss']=str(sf + ' [-58:-49]')
        threeSFRosenberg.loc[j,'second_ss']=str('endogenous')
        j=j+1
        threeSFRosenberg=threeSFRosenberg.append(threeSFcomb.loc[i], ignore_index=True)
        sequence = threeSFcomb.varseq[i].tomutable()
        sequence[33:39]=splicingfactors2[sf]
        threeSFRosenberg.loc[j,'varseq']=sequence.toseq()
        threeSFRosenberg.loc[j,'first_ss']=str(sf + ' [-49:--40]')
        threeSFRosenberg.loc[j,'second_ss']=str('endogenous')
        j=j+1
        threeSFRosenberg=threeSFRosenberg.append(threeSFcomb.loc[i], ignore_index=True)
        sequence = threeSFcomb.varseq[i].tomutable()
        sequence[42:48]=splicingfactors2[sf]
        threeSFRosenberg.loc[j,'varseq']=sequence.toseq()
        threeSFRosenberg.loc[j,'first_ss']=str(sf + ' [-40:-31]')
        threeSFRosenberg.loc[j,'second_ss']=str('endogenous')
        j=j+1
        threeSFRosenberg=threeSFRosenberg.append(threeSFcomb.loc[i], ignore_index=True)
        sequence = threeSFcomb.varseq[i].tomutable()
        sequence[86:92]=splicingfactors2[sf]
        threeSFRosenberg.loc[j,'varseq']=sequence.toseq()
        threeSFRosenberg.loc[j,'first_ss']=str(sf + ' [4:13]')
        threeSFRosenberg.loc[j,'second_ss']=str('endogenous')
        j=j+1
        if (diff>30):
            threeSFRosenberg=threeSFRosenberg.append(threeSFcomb.loc[i], ignore_index=True)
            sequence = threeSFcomb.varseq[i].tomutable()
            sequence[95:101]=splicingfactors2[sf]
            threeSFRosenberg.loc[j,'varseq']=sequence.toseq()
            threeSFRosenberg.loc[j,'first_ss']=str(sf + ' [13:22]')
            threeSFRosenberg.loc[j,'second_ss']=str('endogenous')
            j=j+1
        threeSFRosenberg=threeSFRosenberg.append(threeSFcomb.loc[i], ignore_index=True)
        sequence = threeSFcomb.varseq[i].tomutable()
        sequence[79+diff+8:79+diff+14]=splicingfactors2[sf]
        threeSFRosenberg.loc[j,'varseq']=sequence.toseq()
        threeSFRosenberg.loc[j,'first_ss']=str('endogenous')
        threeSFRosenberg.loc[j,'second_ss']=str(sf + ' [diff+3:diff+12]')
        j=j+1
        if (diff < 51):
            threeSFRosenberg=threeSFRosenberg.append(threeSFcomb.loc[i], ignore_index=True)
            sequence = threeSFcomb.varseq[i].tomutable()
            sequence[79+diff+17:79+diff+23]=splicingfactors2[sf]
            threeSFRosenberg.loc[j,'varseq']=sequence.toseq()
            threeSFRosenberg.loc[j,'first_ss']=str('endogenous')
            threeSFRosenberg.loc[j,'second_ss']=str(sf + ' [diff+12:diff+21]')
            j=j+1

for i in threeSFRosenberg.index:
    if (threeSFRosenberg.varseq[i][80:].translate().find("*")>-1):
        print(str(i))
        print(threeSFRosenberg.varseq[i][80:].translate())
        
threeSFRosenberg.to_pickle('./design/three_splicingfactors_Rosenberg.pkl')

#%% RNA secondary structure
    
threesec=three_filtered[three_filtered.diff_nt>20]

'''
83:92
83:89
86:92
'''

threesecvar=pd.DataFrame(columns=list(threesec.columns))
k=0
for i in threesec.index:
    diff=int(threesec.diff_nt[i])
    firstssrev = threesec.varseq[i][70:79].reverse_complement()
    firstsscomp = threesec.varseq[i][70:79].complement()
    firstssrev2 = threesec.varseq[i][67:76].reverse_complement()
    firstsscomp2 = threesec.varseq[i][67:76].complement()

    altexonrev=threesec.varseq[i][95:104].reverse_complement()
    altexoncomp=threesec.varseq[i][95:104].complement()
    secondssrev = threesec.varseq[i][70+diff:79+diff].reverse_complement()
    secondsscomp = threesec.varseq[i][70+diff:79+diff].complement()
    secondssrev2 = threesec.varseq[i][67+diff:76+diff].reverse_complement()
    secondsscomp2 = threesec.varseq[i][67+diff:76+diff].complement()
    
    exonrev=threesec.varseq[i][96+diff:105+diff].reverse_complement() #for diff_nt<45
    exoncomp=threesec.varseq[i][96+diff:105+diff].complement()

    threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
    sequence=threesec.varseq[i].tomutable()
    sequence[42:51]=firstssrev
    threesecvar.varseq[k]=sequence.toseq()
    threesecvar.loc[k,'firstss']='rev at [42:51]'
    threesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
    sequence=threesec.varseq[i].tomutable()
    sequence[42:51]=firstsscomp
    threesecvar.varseq[k]=sequence.toseq()
    threesecvar.loc[k,'firstss']='comp at [42:51]'
    threesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
    sequence=threesec.varseq[i].tomutable()
    sequence[42:51]=firstssrev2
    threesecvar.varseq[k]=sequence.toseq()
    threesecvar.loc[k,'firstss']='rev2 at [42:51]'
    threesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
    sequence=threesec.varseq[i].tomutable()
    sequence[42:51]=firstsscomp2
    threesecvar.varseq[k]=sequence.toseq()
    threesecvar.loc[k,'firstss']='comp2 at [42:51]'
    threesecvar.loc[k,'secondss']='endogenous'
    k=k+1

    threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
    sequence=threesec.varseq[i].tomutable()
    sequence[83:92]=firstssrev
    threesecvar.varseq[k]=sequence.toseq()
    threesecvar.loc[k,'firstss']='rev at [83:92]'
    threesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
    sequence=threesec.varseq[i].tomutable()
    sequence[83:92]=firstsscomp
    threesecvar.varseq[k]=sequence.toseq()
    threesecvar.loc[k,'firstss']='comp at [83:92]'
    threesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
    sequence=threesec.varseq[i].tomutable()
    sequence[83:92]=firstssrev2
    threesecvar.varseq[k]=sequence.toseq()
    threesecvar.loc[k,'firstss']='rev2 at [83:92]'
    threesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
    sequence=threesec.varseq[i].tomutable()
    sequence[83:92]=firstsscomp2
    threesecvar.varseq[k]=sequence.toseq()
    threesecvar.loc[k,'firstss']='comp2 at [83:92]'
    threesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
    sequence=threesec.varseq[i].tomutable()
    sequence[83:92]=altexonrev
    threesecvar.varseq[k]=sequence.toseq()
    threesecvar.loc[k,'firstss']='altexonrev at [83:92]'
    threesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
    sequence=threesec.varseq[i].tomutable()
    sequence[83:92]=altexoncomp
    threesecvar.varseq[k]=sequence.toseq()
    threesecvar.loc[k,'firstss']='altexoncomp at [83:92]'
    threesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
    sequence=threesec.varseq[i].tomutable()
    sequence[84+diff:93+diff]=secondssrev
    threesecvar.varseq[k]=sequence.toseq()
    threesecvar.loc[k,'firstss']='endogenous'
    threesecvar.loc[k,'secondss']='rev at [84+diff:93+diff]'
    k=k+1
    threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
    sequence=threesec.varseq[i].tomutable()
    sequence[84+diff:93+diff]=secondsscomp
    threesecvar.varseq[k]=sequence.toseq()
    threesecvar.loc[k,'firstss']='endogenous'
    threesecvar.loc[k,'secondss']='comp at [84+diff:93+diff]'
    k=k+1
    threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
    sequence=threesec.varseq[i].tomutable()
    sequence[84+diff:93+diff]=secondssrev2
    threesecvar.varseq[k]=sequence.toseq()
    threesecvar.loc[k,'firstss']='endogenous'
    threesecvar.loc[k,'secondss']='rev2 at [84+diff:93+diff]'
    k=k+1
    threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
    sequence=threesec.varseq[i].tomutable()
    sequence[84+diff:93+diff]=secondsscomp2
    threesecvar.varseq[k]=sequence.toseq()
    threesecvar.loc[k,'firstss']='endogenous'
    threesecvar.loc[k,'secondss']='comp2 at [84+diff:93+diff]'
    k=k+1
    if (diff<45):
        threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
        sequence=threesec.varseq[i].tomutable()
        sequence[83:92]=exonrev
        threesecvar.varseq[k]=sequence.toseq()
        threesecvar.loc[k,'firstss']='exonrev at [83:92]'
        threesecvar.loc[k,'secondss']='endogenous'
        k=k+1
        threesecvar=threesecvar.append(threesec.loc[i],ignore_index=True)    
        sequence=threesec.varseq[i].tomutable()
        sequence[83:92]=exoncomp
        threesecvar.varseq[k]=sequence.toseq()
        threesecvar.loc[k,'firstss']='exoncomp at [83:92]'
        threesecvar.loc[k,'secondss']='endogenous'
        k=k+1
        
threesecvar.to_pickle('./design/three_secondarystructure_variants.pkl')


#%% nucleosome positioning/recoding

import random
from Bio import SeqUtils

SynonymousCodons = { 
    'CYS': ['TGT', 'TGC'], 
    'ASP': ['GAT', 'GAC'], 
    'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'], 
    'GLN': ['CAA', 'CAG'], 
    'MET': ['ATG'], 
    'ASN': ['AAC', 'AAT'], 
    'PRO': ['CCT', 'CCG', 'CCA', 'CCC'], 
    'LYS': ['AAG', 'AAA'], 
    'STOP': ['TAG', 'TGA', 'TAA'], 
    'THR': ['ACC', 'ACA', 'ACG', 'ACT'], 
    'PHE': ['TTT', 'TTC'], 
    'ALA': ['GCA', 'GCC', 'GCG', 'GCT'], 
    'GLY': ['GGT', 'GGG', 'GGA', 'GGC'], 
    'ILE': ['ATC', 'ATA', 'ATT'], 
    'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'], 
    'HIS': ['CAT', 'CAC'], 
    'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'], 
    'TRP': ['TGG'], 
    'VAL': ['GTA', 'GTC', 'GTG', 'GTT'], 
    'GLU': ['GAG', 'GAA'], 
    'TYR': ['TAT', 'TAC'] 
    } 

gencode={}
for i in SynonymousCodons.keys():
    gencode[SeqUtils.seq1(i)]=SynonymousCodons[i]

gencode['*']=['TAG', 'TGA', 'TAA']

codonmaxgc={
 'A': ['GCG'],
 'C': ['TGC'],
 'D': ['GAC'],
 'E': ['GAG'],
 'F': ['TTC'],
 'G': ['GGC'],
 'H': ['CAC'],
 'I': ['ATC'],
 'K': ['AAG'],
 'L': ['CTG'],
 'M': ['ATG'],
 'N': ['AAC'],
 'P': ['CCG'],
 'Q': ['CAG'],
 'R': ['CGG'],
 'S': ['TCC'],
 'T': ['ACG'],
 'V': ['GTC'],
 'W': ['TGG'],
 '*': ['TGA'],
 'Y': ['TAC']}

codonmingc={
 'A': ['GCA'],
 'C': ['TGT'],
 'D': ['GAT'],
 'E': ['GAA'],
 'F': ['TTT'],
 'G': ['GGA'],
 'H': ['CAT'],
 'I': ['ATA'],
 'K': ['AAA'],
 'L': ['TTA'],
 'M': ['ATG'],
 'N': ['AAT'],
 'P': ['CCT'],
 'Q': ['CAA'],
 'R': ['AGA'],
 'S': ['AGT'],
 'T': ['ACT'],
 'V': ['GTA'],
 'W': ['TGG'],
 '*': ['TAA'],
 'Y': ['TAT']}

for i in three_filtered.index:
    three_filtered.loc[i,'aaseq']=three_filtered.loc[i,'varseq'][80:].translate()

threenuc=pd.DataFrame(columns=list(three_filtered))

k=0
for i in three_filtered.index:
    diff=int(three_filtered.diff_nt[i])
    if (diff<35):
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        threenuc.loc[k,'changes']='endogenous'
        k=k+1
        
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + codonmaxgc[j][0]
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon recoded maximal gc'
        k=k+1
        
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + codonmaxgc[j][0]
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon recoded minimal gc'
        k=k+1
        
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + random.choice(gencode[j])
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon recoded randomly 1'
        k=k+1
    
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + random.choice(gencode[j])
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon recoded randomly 2'
        k=k+1
    
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + random.choice(gencode[j])
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon recoded randomly 3'
        k=k+1
    
    else:
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        threenuc.loc[k,'changes']='endogenous'
        k=k+1
        
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + codonmaxgc[j][0]
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon recoded maximal gc'
        k=k+1
        
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + codonmaxgc[j][0]
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon recoded minimal gc'
        k=k+1
        
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + random.choice(gencode[j])
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon recoded randomly 1'
        k=k+1
    
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + random.choice(gencode[j])
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon recoded randomly 2'
        k=k+1
    
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + random.choice(gencode[j])
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon recoded randomly 3'
        k=k+1
    




        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1:diff//3-6]:
            recod=recod + codonmaxgc[j][0]
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83:83+len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='altexon recoded maximal gc'
        k=k+1
        
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + codonmaxgc[j][0]
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon and altexon recoded maximal gc'
        k=k+1
        
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1:diff//3-6]:
            recod=recod + codonmingc[j][0]
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83:83+len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='altexon recoded minimal gc'
        k=k+1
        
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + codonmaxgc[j][0]
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon and altexon recoded minimal gc'
        k=k+1
        

        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1:diff//3-6]:
            recod=recod + codonmingc[j][0]
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83:83+len(recod)]=recod
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + codonmaxgc[j][0]
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon recoded maximal and altexon recoded minimal gc'
        k=k+1
        
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1:diff//3-6]:
            recod=recod + codonmaxgc[j][0]
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83:83+len(recod)]=recod
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + codonmingc[j][0]
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon recoded minimal and and altexon recoded maximal gc'
        k=k+1
        

        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1:diff//3-6]:
            recod=recod + random.choice(gencode[j])
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83:83+len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='altexon recoded randomly 1'
        k=k+1
    
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + random.choice(gencode[j])
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon and altexon recoded randomly 1'
        k=k+1
    
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1:diff//3-6]:
            recod=recod + random.choice(gencode[j])
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83:83+len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='altexon recoded randomly 2'
        k=k+1
    
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + random.choice(gencode[j])
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon and altexon recoded randomly 2'
        k=k+1
    
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1:diff//3-6]:
            recod=recod + random.choice(gencode[j])
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83:83+len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='altexon recoded randomly 3'
        k=k+1
    
        threenuc=threenuc.append(three_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in three_filtered.aaseq[i][1+diff//3:]:
            recod=recod + random.choice(gencode[j])
        sequence[83+(diff//3)*3:83+(diff//3)*3 + len(recod)]=recod
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='exon and altexon recoded randomly 3'
        k=k+1
    
threedisf=three_filtered[(three_filtered.diff_nt>5)]     
    
imp22dT=Seq('TTTTTTTTTCATTTTTTTTTTT')
perf12dT=Seq('TTTTTTTTTTTT')
perf5dT=Seq('TTTTTC')
imp22ctrl=Seq('TTGGACTGTCATTCACGTCGTT')
perf12ctrl=Seq('TTGGACTGTCAT')
perf5ctrl=Seq('TTGGAC')


for i in threedisf.index:
    threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
    sequence=three_filtered.varseq[i].tomutable()
    sequence[3:25]=imp22dT
    threenuc.varseq[k]=sequence.toseq()
    threenuc.loc[k,'changes']='imperfect 22nt dT tract at position 0'
    k=k+1

    threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
    sequence=three_filtered.varseq[i].tomutable()
    sequence[3:15]=perf12dT
    threenuc.varseq[k]=sequence.toseq()
    threenuc.loc[k,'changes']='perfect 12nt dT tract at position 0'
    k=k+1

    threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
    sequence=three_filtered.varseq[i].tomutable()
    sequence[9:21]=perf12dT
    threenuc.varseq[k]=sequence.toseq()
    threenuc.loc[k,'changes']='perfect 12nt dT tract at position 6'
    k=k+1

    threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
    sequence=three_filtered.varseq[i].tomutable()
    sequence[15:27]=perf12dT
    threenuc.varseq[k]=sequence.toseq()
    threenuc.loc[k,'changes']='perfect 12nt dT tract at position 12'
    k=k+1
    
    threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
    sequence=three_filtered.varseq[i].tomutable()
    sequence[21:33]=perf12dT
    threenuc.varseq[k]=sequence.toseq()
    threenuc.loc[k,'changes']='perfect 12nt dT tract at position 18'
    k=k+1
    
    threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
    sequence=three_filtered.varseq[i].tomutable()
    sequence[27:39]=perf12dT
    threenuc.varseq[k]=sequence.toseq()
    threenuc.loc[k,'changes']='perfect 12nt dT tract at position 24'
    k=k+1
    
    
    threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
    sequence=three_filtered.varseq[i].tomutable()
    sequence[3:25]=imp22ctrl
    threenuc.varseq[k]=sequence.toseq()
    threenuc.loc[k,'changes']='control 22nt at position 0'
    k=k+1

    threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
    sequence=three_filtered.varseq[i].tomutable()
    sequence[3:15]=perf12ctrl
    threenuc.varseq[k]=sequence.toseq()
    threenuc.loc[k,'changes']='control 12nt at position 0'
    k=k+1

    threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
    sequence=three_filtered.varseq[i].tomutable()
    sequence[9:21]=perf12ctrl
    threenuc.varseq[k]=sequence.toseq()
    threenuc.loc[k,'changes']='control 12nt at position 6'
    k=k+1

    threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
    sequence=three_filtered.varseq[i].tomutable()
    sequence[15:27]=perf12ctrl
    threenuc.varseq[k]=sequence.toseq()
    threenuc.loc[k,'changes']='control 12nt at position 12'
    k=k+1

    threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
    sequence=three_filtered.varseq[i].tomutable()
    sequence[21:33]=perf12ctrl
    threenuc.varseq[k]=sequence.toseq()
    threenuc.loc[k,'changes']='control 12nt at position 18'
    k=k+1

    threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
    sequence=three_filtered.varseq[i].tomutable()
    sequence[27:39]=perf12ctrl
    threenuc.varseq[k]=sequence.toseq()
    threenuc.loc[k,'changes']='control 12nt at position 24'
    k=k+1

    diff=int(threedisf.diff_nt[i])
    if (diff>33):
        threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83:89]=perf5dT
        sequence[-30-12:-30-6]=perf5dT
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='perfect 5nt dT tract at altexon position 4'
        k=k+1
    
        threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
        sequence=three_filtered.varseq[i].tomutable()
        sequence[84+diff:84+diff+6]=perf5dT
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='perfect 5nt dT tract at position exon +5'
        k=k+1
    
        threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
        sequence=three_filtered.varseq[i].tomutable()
        sequence[83:89]=perf5dT
        sequence[-30-12:-30-6]=perf5ctrl
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='control 5nt at altexon position 4'
        k=k+1
    
        threenuc=threenuc.append(threedisf.loc[i],ignore_index=True)
        sequence=three_filtered.varseq[i].tomutable()
        sequence[84+diff:84+diff+6]=perf5ctrl
        threenuc.varseq[k]=sequence.toseq()
        threenuc.loc[k,'changes']='control 5nt at position exon +5'
        k=k+1
    
for i in threenuc.index:
    if (threenuc.varseq[i][80:].translate().find("*")>-1):
        print(str(i))
        print(threenuc.changes[i])
        print(threenuc.varseq[i][80:].translate())
#        threesecvar.drop(i,inplace=True)
 
threenuc.drop_duplicates('varseq',inplace=True)

from Bio.Restriction import RsrII
from Bio.Restriction import AscI
from Bio.Restriction import SpeI
from Bio.Restriction import AatII


for i in threenuc.index:
    if (RsrII.search(threenuc.varseq[i])!=[])|(AscI.search(threenuc.varseq[i])!=[])|(SpeI.search(threenuc.varseq[i])!=[])|(AatII.search(threenuc.varseq[i])!=[]):
        threenuc.drop(i,inplace=True)
        

   
threenuc.to_pickle('./design/three_nucleosome_recoding_methylation.pkl')

#%% combine with constitutive exon and intron

threeconstforcomb=three_filtered[(three_filtered.diff_nt>11)&(three_filtered.diff_nt<50)]
const=pd.read_pickle('./design/constexon62.pkl')

threeconstcomb=pd.DataFrame()

k=0
for i in threeconstforcomb.index:
    diff=int(threeconstforcomb.diff_nt[i])
    for j in const.index:
        
        threeconstcomb=threeconstcomb.append(threeconstforcomb.loc[i],ignore_index=True)
        sequence=threeconstforcomb.varseq[i].tomutable()
        sequence[:79]=const.sequence[j][-100-int(const.exonlength[j])-79:-100-int(const.exonlength[j])]
        threeconstcomb.varseq[k]=sequence.toseq()
        threeconstcomb.loc[k,'changes']='intron ' + str(j)
        k=k+1
        
        threeconstcomb=threeconstcomb.append(threeconstforcomb.loc[i],ignore_index=True)
        sequence=threeconstforcomb.varseq[i].tomutable()
        sequence[79+diff:]=const.sequence[j][-100-int(const.exonlength[j]):-100-int(const.exonlength[j])+150-79-diff]
        threeconstcomb.varseq[k]=sequence.toseq()
        threeconstcomb.loc[k,'changes']='exon ' + str(j)
        k=k+1


for i in threeconstcomb.index:
    if (threeconstcomb.varseq[i][80:].translate().find("*")>-1):
        print(str(i))
        print(threeconstcomb.changes[i])
        print(threeconstcomb.varseq[i][80:].translate())
#        threeconstcomb.drop(i,inplace=True)

for i in threeconstcomb.index:
    if (len(threeconstcomb.varseq[i])!=150):
        print(str(i))
        print(threeconstcomb.changes[i])
        print(threeconstcomb.varseq[i][80:].translate())
#        threeconstcomb.drop(i,inplace=True)
  
threeconstcomb.to_pickle('./design/three_combinatorial_with_constitutive.pkl')



#%% Switch elements between the two splice sites

threeswitch=three_filtered[(three_filtered.diff_nt>20)]

#threeswitch.drop(410585,inplace=True)

threeswitchvar=pd.DataFrame(columns=list(threeswitch.columns))

k=0
for i in threeswitch.index:
    diff=int(threeswitch.diff_nt[i])
    
    if (threeswitch.varseq[i][79-16:78].translate().find("*")==-1):
        firstss=threeswitch.varseq[i][64:79]
        secondss=threeswitch.varseq[i][64+diff:79+diff]
        
        threeswitchvar=threeswitchvar.append(threeswitch.loc[i],ignore_index=True)
        threeswitchvar.loc[k,'firstss']='firstss -15nt'
        threeswitchvar.loc[k,'secondss']='secondss -15nt'
        k=k+1
        
        threeswitchvar=threeswitchvar.append(threeswitch.loc[i],ignore_index=True)
        sequence=threeswitch.varseq[i].tomutable()
        sequence[64+diff:79+diff]=firstss
        threeswitchvar.loc[k,'varseq']=sequence.toseq()
        threeswitchvar.loc[k,'firstss']='firstss -15nt'
        threeswitchvar.loc[k,'secondss']='firstss -15nt'
        k=k+1
            
        threeswitchvar=threeswitchvar.append(threeswitch.loc[i],ignore_index=True)
        sequence=threeswitch.varseq[i].tomutable()
        sequence[64:79]=secondss
        threeswitchvar.loc[k,'varseq']=sequence.toseq()
        threeswitchvar.loc[k,'firstss']='secondss -15nt'
        threeswitchvar.loc[k,'secondss']='secondss -15nt'
        k=k+1
        
        threeswitchvar=threeswitchvar.append(threeswitch.loc[i],ignore_index=True)
        sequence=threeswitch.varseq[i].tomutable()
        sequence[64:79]=secondss
        sequence[64+diff:79+diff]=firstss
        threeswitchvar.loc[k,'varseq']=sequence.toseq()
        threeswitchvar.loc[k,'firstss']='secondss -15nt'
        threeswitchvar.loc[k,'secondss']='firstss -15nt'
        k=k+1
    
#three_filtered.drop(410585,inplace=True)

threeswitch=three_filtered[(three_filtered.diff_nt>16)&(three_filtered.diff_nt<55)]
     
for i in threeswitch.index:
    if (threeswitch.varseq[i][78-int(threeswitch.diff_nt[i])+6-1:78].translate().find('*')>-1):
        threeswitch.drop(i,inplace=True)
        

for i in threeswitch.index:
    diff=int(threeswitch.diff_nt[i])
    if (threeswitch.varseq[i][80+diff:80+diff+9].translate().find("*")==-1) & (threeswitch.varseq[i][78:87].translate().find("*")==-1):
        for pos in range(0,3,1):
            threeswitchvar=threeswitchvar.append(threeswitch.loc[i],ignore_index=True)
            sequence=threeswitch.varseq[i].tomutable()
            sequence[79 +diff:82+diff + pos*3]=sequence[79:82 + pos*3]
            threeswitchvar.loc[k,'varseq']=sequence.toseq()
            threeswitchvar.loc[k,'firstss']='firstss ' + str(3+pos*3) + 'nt'
            threeswitchvar.loc[k,'secondss']='firstss ' + str(3+pos*3) + 'nt'
            k=k+1
                
            threeswitchvar=threeswitchvar.append(threeswitch.loc[i],ignore_index=True)
            sequence=threeswitch.varseq[i].tomutable()
            sequence[79:82 + pos*3]=sequence[79 +diff:82+diff + pos*3]
            threeswitchvar.loc[k,'varseq']=sequence.toseq()
            threeswitchvar.loc[k,'firstss']='secondss ' + str(3+pos*3) + 'nt'
            threeswitchvar.loc[k,'secondss']='secondss ' + str(3+pos*3) + 'nt'
            k=k+1
            
            threeswitchvar=threeswitchvar.append(threeswitch.loc[i],ignore_index=True)
            sequence=threeswitch.varseq[i].tomutable()
            sequence[79:82 + pos*3]=sequence[79 +diff:82+diff + pos*3]
            sequence[79 +diff:82+diff + pos*3]= sequence[79:82 + pos*3]
            threeswitchvar.loc[k,'varseq']=sequence.toseq()
            threeswitchvar.loc[k,'firstss']='secondss ' + str(3+pos*3) + 'nt'
            threeswitchvar.loc[k,'secondss']='firstss ' + str(3+pos*3) + 'nt'
            k=k+1
        
    for pos in range(0,(diff-9)/3,1):
        threeswitchvar=threeswitchvar.append(threeswitch.loc[i],ignore_index=True)
        sequence=threeswitch.varseq[i].tomutable()
        sequence[76 - pos*3 +diff:79 + diff]=sequence[76 - pos*3 :79]
        threeswitchvar.loc[k,'varseq']=sequence.toseq()
        threeswitchvar.loc[k,'firstss']='firstss -' + str(3+pos*3) + 'nt'
        threeswitchvar.loc[k,'secondss']='firstss -' + str(3+pos*3) + 'nt'
        k=k+1
            
        threeswitchvar=threeswitchvar.append(threeswitch.loc[i],ignore_index=True)
        sequence=threeswitch.varseq[i].tomutable()
        sequence[76 - pos*3 :79]=sequence[76 - pos*3 +diff:79 + diff]
        threeswitchvar.loc[k,'varseq']=sequence.toseq()
        threeswitchvar.loc[k,'firstss']='secondss -' + str(3+pos*3) + 'nt'
        threeswitchvar.loc[k,'secondss']='secondss -' + str(3+pos*3) + 'nt'
        k=k+1
        
        threeswitchvar=threeswitchvar.append(threeswitch.loc[i],ignore_index=True)
        sequence=threeswitch.varseq[i].tomutable()
        sequence[76 - pos*3 :79]=sequence[76 - pos*3 +diff:79 + diff]
        sequence[76 - pos*3 +diff:79 + diff]=sequence[76 - pos*3 :79]
        threeswitchvar.loc[k,'varseq']=sequence.toseq()
        threeswitchvar.loc[k,'firstss']='secondss -' + str(3+pos*3) + 'nt'
        threeswitchvar.loc[k,'secondss']='firstss -' + str(3+pos*3) + 'nt'
        k=k+1
    
for i in threeswitchvar.index:
    if (threeswitchvar.varseq[i][80:].translate().find("*")>-1):
        print(str(i))
        print(threeswitchvar.firstss[i])
        print(threeswitchvar.secondss[i])
        print(threeswitchvar.varseq[i][80:].translate())
        threeswitchvar.drop(i,inplace=True)

for i in threeswitchvar.index:
    if (len(threeswitchvar.varseq[i])!=150):
        print(str(i))
        print(threeswitchvar.changes[i])
        print(threeswitchvar.varseq[i][:51+int(threeswitchvar.diff_nt[i])].translate())
#        threeswitchvar.drop(i,inplace=True)
  
threeswitchvar.to_pickle('./design/three_switch_splice_site_sequences.pkl')


#%% length

threelength= three_filtered[three_filtered.diff_nt>33]
threelengthvar=pd.DataFrame()
linker=Seq('gaagctgctgcaagagaagctgcagctagggaggctgcagctagggaggctgctgcaaga')

k=0

for i in threelength.index:
    diff=int(threelength.diff_nt[i])
    seqext=threelength.varseq[i]
        
    
    count = (diff-4-8)//6
    print(count)
    for ex in range(1,count,1):
        
        threelengthvar=threelengthvar.append(threelength.loc[i],ignore_index=True)
        threelengthvar.varseq[k]=seqext[:83] + seqext[83+ex*6:148] + linker[:2+ex*6]
        threelengthvar.loc[k,'diff_nt']=str(diff - ex*6)
        threelengthvar.loc[k,'changes']=str(ex*6) + 'nt deletion at position 4'
        k=k+1

        threelengthvar=threelengthvar.append(threelength.loc[i],ignore_index=True)
        middle = 83 + (((diff-4-8)//3)//2)*3
        threelengthvar.varseq[k]=seqext[:middle-ex*3] + seqext[middle+ex*3:148] + linker[:2+ex*6]
        threelengthvar.loc[k,'diff_nt']=str(diff - ex*6)
        threelengthvar.loc[k,'changes']=str(ex*6) + 'nt deletion at center'
        k=k+1
        
        threelengthvar=threelengthvar.append(threelength.loc[i],ignore_index=True)
        threelengthvar.varseq[k]=seqext[:79+diff-8 - ex*6] + seqext[79+diff-8:148] + linker[:2+ex*6]
        threelengthvar.loc[k,'diff_nt']=str(diff - ex*6)
        threelengthvar.loc[k,'changes']=str(ex*6) + 'nt deletion at position -8'
        k=k+1

for i in threelengthvar.index:
    if (threelengthvar.varseq[i][80:].translate().find("*")>-1):
        print(str(i))
        print(threelengthvar.varseq[i][80:].translate())
        threelengthvar.drop(i,inplace=True)

for i in threelengthvar.index:
    if (len(threelengthvar.varseq[i])!=150):
        print(str(i))
        print(threelengthvar.varseq[i][:51+int(threeswitchvar.diff_nt[i])].translate())
#        threeswitchvar.drop(i,inplace=True)
    
threelengthvar.to_pickle('./design/three_exonlength.pkl')
    
