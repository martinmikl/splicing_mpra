# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 21:43:04 2016

@author: miklm
"""


import pandas as pd
import numpy as np
from Bio.Seq import Seq

#%%

cassette=pd.read_pickle('./design/juncasselected.pkl')

cassette[['level','level_alt']]= \
    cassette[['level','level_alt']].astype(int)


for i in cassette.index:
    cassette.loc[i,'levelratio']=np.log2(float(cassette.level_alt[i])/float(cassette.level[i]))

for i in cassette.index:
    if (cassette.strand[i]=='+'):
        cassette.loc[i,'intronup_length']=int(cassette.end[i]-cassette.start[i])
        cassette.loc[i,'introndown_length']=int(cassette.end2[i]-cassette.start2[i])
    else:
        cassette.loc[i,'intronup_length']=int(cassette.end2[i]-cassette.start2[i])
        cassette.loc[i,'introndown_length']=int(cassette.end[i]-cassette.start[i])


cassette_filtered=cassette[(cassette.level_alt>20) & (cassette.level>20) & (cassette.levelratio>-4) & (cassette.levelratio<4)]
#
cassette_filtered.drop_duplicates('commonname',inplace=True)


for i in cassette_filtered.index:
    cassette_filtered.loc[i,'varseq']=str(cassette_filtered.sequence[i])\
        [100+int(cassette_filtered.intronup_length[i]) \
        + int(cassette_filtered.exonlength[i]) +30 \
        - 150:100+int(cassette_filtered.intronup_length[i]) \
        + int(cassette_filtered.exonlength[i]) +30].upper()

cassette_filtered.to_pickle('./design/cassette_filtered.pkl')

cassette_filtered['exonstart_varseq']=120-cassette_filtered.exonlength.astype(int)
cassette_filtered['exonend_varseq']=120

cassette_filtered[['commonname','coordinates','end','start2','exonlength','exonstart_varseq','exonend_varseq','varseq']].to_csv('./tables/TableS2_casfiltered.csv')

#%%

'''
make 4 groups, within which "exonic - alternative - intronic" sequences are
varied by permutation

The four groups are: number of genes in the group - diff_nt
5x 35
6x 47
8x 53
8x 68
'''


df=cassette_filtered.copy()
for i in df.index:
    df.loc[i,'intronup']=df.sequence[i][100+int(df.intronup_length[i]) + int(df.exonlength[i]) +30 - 150 :100+int(df.intronup_length[i])]
    df.loc[i,'altexon']=df.sequence[i][100+int(df.intronup_length[i]):100+ int(df.intronup_length[i]) + int(df.exonlength[i])]
    df.loc[i,'introndown']=df.sequence[i][100+int(df.intronup_length[i]) + int(df.exonlength[i]):100+int(df.intronup_length[i]) + int(df.exonlength[i])+30]
cassette_filtered=df.copy()

cassette35=cassette_filtered[cassette_filtered.exonlength==35]
cassette47=cassette_filtered[cassette_filtered.exonlength==47]
cassette50=cassette_filtered[cassette_filtered.exonlength==50]
#cassette68=cassette_filtered[cassette_filtered.exonlength==68]
cassette71=cassette_filtered[cassette_filtered.exonlength==71]
cassette50.drop(492702,inplace=True)
cassette71.drop(561749,inplace=True)


#%%

# create all possible combinations within each group
df=cassette35.copy()
out=pd.DataFrame()
l=0
for i in df.index:
    for j in df.index:
        for k in df.index:
            out.loc[l,'varseq']=df.loc[i,'intronup']+df.loc[j,'altexon']+df.loc[k,'introndown']
            out.loc[l,'geneintronup']=df.commonname[i]
            out.loc[l,'genealtexon']=df.commonname[j]
            out.loc[l,'geneintrondown']=df.commonname[k]
            out.loc[l,'exonlength']=df.exonlength[j]
            l=l+1
cassette35combvar=out.copy()

df=cassette47.copy()
out=pd.DataFrame()
l=0
for i in df.index:
    for j in df.index:
        for k in df.index:
            out.loc[l,'varseq']=df.loc[i,'intronup']+df.loc[j,'altexon']+df.loc[k,'introndown']
            out.loc[l,'geneintronup']=df.commonname[i]
            out.loc[l,'genealtexon']=df.commonname[j]
            out.loc[l,'geneintrondown']=df.commonname[k]
            out.loc[l,'exonlength']=df.exonlength[j]
            l=l+1
cassette47combvar=out.copy()

df=cassette50.copy()
out=pd.DataFrame()
l=0
for i in df.index:
    for j in df.index:
        for k in df.index:
            out.loc[l,'varseq']=df.loc[i,'intronup']+df.loc[j,'altexon']+df.loc[k,'introndown']
            out.loc[l,'geneintronup']=df.commonname[i]
            out.loc[l,'genealtexon']=df.commonname[j]
            out.loc[l,'geneintrondown']=df.commonname[k]
            out.loc[l,'exonlength']=df.exonlength[j]
            l=l+1
cassette50combvar=out.copy()

df=cassette71.copy()
out=pd.DataFrame()
l=0
for i in df.index:
    for j in df.index:
        for k in df.index:
            out.loc[l,'varseq']=df.loc[i,'intronup']+df.loc[j,'altexon']+df.loc[k,'introndown']
            out.loc[l,'geneintronup']=df.commonname[i]
            out.loc[l,'genealtexon']=df.commonname[j]
            out.loc[l,'geneintrondown']=df.commonname[k]
            out.loc[l,'exonlength']=df.exonlength[j]
            l=l+1
cassette71combvar=out.copy()



cassetteprime_combinatorial_variations=pd.concat([cassette35combvar,cassette47combvar,cassette50combvar,cassette71combvar],ignore_index=True)

### PICKLE
cassetteprime_combinatorial_variations.to_pickle('./design/cassette_combinatorial_variations.pkl')

#%% replace with constitutive splice sites

#select everything larger than 15 diff_nt and smaller than 50 diff_nt and change to constitutive splice sites and splice sites mutants

cassetteconst = cassette_filtered[(cassette_filtered.exonlength > 34) & (cassette_filtered.exonlength < 72)]
cassetteconst.drop([492702, 561749],inplace=True)

canonical_donor = Seq('CAGGTAAGT').tomutable()
no_donor = Seq('CTGCTC').tomutable()
canonical_acceptor = Seq('CTCCTTTCCTTTCAGGTC').tomutable()
no_acceptor = Seq('CAGAGAGGA').tomutable()
GC_donor = Seq('CAGGCAAGT').tomutable()
branchpoint = Seq('CTCAC').tomutable()


cassetteconstvar=pd.DataFrame(columns=list(cassetteconst.columns))
j=0
for i in cassetteconst.index:
    exonlength=int(cassetteconst.exonlength[i])
    cassetteconstvar=cassetteconstvar.append(cassetteconst.loc[i], ignore_index=True)
    cassetteconstvar.loc[j,'varseq'] = cassetteconst.sequence[i]\
    [100+int(cassetteconst.intronup_length[i]) + int(cassetteconst.exonlength[i]) +30 - 150:100+int(cassetteconst.intronup_length[i]) + int(cassetteconst.exonlength[i]) +30]
    cassetteconstvar.loc[j,'first_ss']=str('endogenous')
    cassetteconstvar.loc[j,'second_ss']=str('endogenous')
    j=j+1
    cassetteconstvar=cassetteconstvar.append(cassetteconst.loc[i], ignore_index=True)
    sequence = cassetteconst.sequence[i]\
        [100+int(cassetteconst.intronup_length[i]) + int(cassetteconst.exonlength[i]) +30 - 150:100+int(cassetteconst.intronup_length[i]) + int(cassetteconst.exonlength[i]) +30]\
        .tomutable()
    sequence[-33:-24]=canonical_donor
    sequence[-(30+exonlength+15):-(30+exonlength-3)]=canonical_acceptor
    cassetteconstvar.loc[j,'varseq']=sequence.toseq()
    cassetteconstvar.loc[j,'first_ss']=str('constitutive')
    cassetteconstvar.loc[j,'second_ss']=str('constitutive')
    j=j+1
    cassetteconstvar=cassetteconstvar.append(cassetteconst.loc[i], ignore_index=True)
    sequence = cassetteconst.sequence[i]\
        [100+int(cassetteconst.intronup_length[i]) + int(cassetteconst.exonlength[i]) +30 - 150:100+int(cassetteconst.intronup_length[i]) + int(cassetteconst.exonlength[i]) +30]\
        .tomutable()
    sequence[-33:-24]=GC_donor
    sequence[-(30+exonlength+15):-(30+exonlength-3)]=canonical_acceptor
    cassetteconstvar.loc[j,'varseq']=sequence.toseq()
    cassetteconstvar.loc[j,'first_ss']=str('constitutive')
    cassetteconstvar.loc[j,'second_ss']=str('GC')
    j=j+1
    cassetteconstvar=cassetteconstvar.append(cassetteconst.loc[i], ignore_index=True)
    sequence = cassetteconst.sequence[i]\
        [100+int(cassetteconst.intronup_length[i]) + int(cassetteconst.exonlength[i]) +30 - 150:100+int(cassetteconst.intronup_length[i]) + int(cassetteconst.exonlength[i]) +30]\
        .tomutable()
    sequence[-30:-24]=no_donor
    sequence[-(30+exonlength+9):-(30+exonlength)]=no_acceptor
    cassetteconstvar.loc[j,'varseq']=sequence.toseq()
    cassetteconstvar.loc[j,'first_ss']=str('unsplicable')
    cassetteconstvar.loc[j,'second_ss']=str('unsplicable')
    j=j+1
    cassetteconstvar=cassetteconstvar.append(cassetteconst.loc[i], ignore_index=True)
    sequence = cassetteconst.sequence[i]\
        [100+int(cassetteconst.intronup_length[i]) + int(cassetteconst.exonlength[i]) +30 - 150:100+int(cassetteconst.intronup_length[i]) + int(cassetteconst.exonlength[i]) +30]\
        .tomutable()
    sequence[-(30+exonlength+26):-(30+exonlength+21)]=branchpoint
    cassetteconstvar.loc[j,'varseq']=sequence.toseq()
    cassetteconstvar.loc[j,'first_ss']=str('constBP')
    cassetteconstvar.loc[j,'second_ss']=str('endogenous')
    j=j+1
    
cassetteconstvar.to_pickle('./design/cassetteprime_constitutiveandunsplicable_variations.pkl')

#%%

cassetteSF = cassette_filtered.loc[[293188,174196,559501,636056,11570,238913,493093,649564,661524,666274,80120,212181,566550]]

cassetteSFvar=pd.DataFrame(columns=list(cassetteSF.columns))

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
[-58:-49] --- not for exonlength >51
[-49:-40]
[-40:-31]
[5:14]
[14:23] 
[diff-21:diff-12] --- not for exonlength <40
[diff-12:diff-3]
[diff+6:diff+15]
[diff+15:diff+24]
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


j=0
for i in cassetteSF.index:
    diff=int(cassetteSF.exonlength[i])
    for sf in splicingfactors1.index:
        if (diff < 51):
            cassetteSFvar=cassetteSFvar.append(cassetteSF.loc[i], ignore_index=True)
            sequence = cassette_filtered.sequence[i]\
                [100+int(cassette_filtered.intronup_length[i]) + diff +30 - 150: \
                100+int(cassette_filtered.intronup_length[i]) + diff +30]\
                .tomutable()
            sequence[120-diff-58:120-diff-49]=splicingfactors1[sf]
            cassetteSFvar.loc[j,'varseq']=sequence.toseq()
            cassetteSFvar.loc[j,'first_ss']=str(sf + ' [-58:-49]')
            cassetteSFvar.loc[j,'second_ss']=str('endogenous')
            j=j+1
    
        cassetteSFvar=cassetteSFvar.append(cassetteSF.loc[i], ignore_index=True)
        sequence = cassette_filtered.sequence[i]\
            [100+int(cassette_filtered.intronup_length[i]) + diff +30 - 150: \
            100+int(cassette_filtered.intronup_length[i]) + diff +30]\
            .tomutable()
        sequence[120-diff-49:120-diff-40]=splicingfactors1[sf]
        cassetteSFvar.loc[j,'varseq']=sequence.toseq()
        cassetteSFvar.loc[j,'first_ss']=str(sf + ' [-49:-40]')
        cassetteSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1

        cassetteSFvar=cassetteSFvar.append(cassetteSF.loc[i], ignore_index=True)
        sequence = cassette_filtered.sequence[i]\
            [100+int(cassette_filtered.intronup_length[i]) + diff +30 - 150: \
            100+int(cassette_filtered.intronup_length[i]) + diff +30]\
            .tomutable()
        sequence[120-diff-40:120-diff-31]=splicingfactors1[sf]
        cassetteSFvar.loc[j,'varseq']=sequence.toseq()
        cassetteSFvar.loc[j,'first_ss']=str(sf + ' [-40:-31]')
        cassetteSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1

        cassetteSFvar=cassetteSFvar.append(cassetteSF.loc[i], ignore_index=True)
        sequence = cassette_filtered.sequence[i]\
            [100+int(cassette_filtered.intronup_length[i]) + diff +30 - 150: \
            100+int(cassette_filtered.intronup_length[i]) + diff +30]\
            .tomutable()
        sequence[120-diff+5:120-diff+14]=splicingfactors1[sf]
        cassetteSFvar.loc[j,'varseq']=sequence.toseq()
        cassetteSFvar.loc[j,'first_ss']=str(sf + ' [5:14]')
        cassetteSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1

        cassetteSFvar=cassetteSFvar.append(cassetteSF.loc[i], ignore_index=True)
        sequence = cassette_filtered.sequence[i]\
            [100+int(cassette_filtered.intronup_length[i]) + diff +30 - 150: \
            100+int(cassette_filtered.intronup_length[i]) + diff +30]\
            .tomutable()
        sequence[120-diff+14:120-diff+23]=splicingfactors1[sf]
        cassetteSFvar.loc[j,'varseq']=sequence.toseq()
        cassetteSFvar.loc[j,'first_ss']=str(sf + ' [+14:+23]')
        cassetteSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        if (diff>40):
            cassetteSFvar=cassetteSFvar.append(cassetteSF.loc[i], ignore_index=True)
            sequence = cassette_filtered.sequence[i]\
                [100+int(cassette_filtered.intronup_length[i]) + diff +30 - 150: \
                100+int(cassette_filtered.intronup_length[i]) + diff +30]\
                .tomutable()
            sequence[120-21:120-12]=splicingfactors1[sf]
            cassetteSFvar.loc[j,'varseq']=sequence.toseq()
            cassetteSFvar.loc[j,'first_ss']=str('endogenous')
            cassetteSFvar.loc[j,'second_ss']=str(sf + ' [-21:-12]')
            j=j+1

        cassetteSFvar=cassetteSFvar.append(cassetteSF.loc[i], ignore_index=True)
        sequence = cassette_filtered.sequence[i]\
            [100+int(cassette_filtered.intronup_length[i]) + diff +30 - 150: \
            100+int(cassette_filtered.intronup_length[i]) + diff +30]\
            .tomutable()
        sequence[120-12:120-3]=splicingfactors1[sf]
        cassetteSFvar.loc[j,'varseq']=sequence.toseq()
        cassetteSFvar.loc[j,'first_ss']=str('endogenous')
        cassetteSFvar.loc[j,'second_ss']=str(sf + ' [-12:-3]')
        j=j+1

        cassetteSFvar=cassetteSFvar.append(cassetteSF.loc[i], ignore_index=True)
        sequence = cassette_filtered.sequence[i]\
            [100+int(cassette_filtered.intronup_length[i]) + diff +30 - 150: \
            100+int(cassette_filtered.intronup_length[i]) + diff +30]\
            .tomutable()
        sequence[120 + 6:120 + 15]=splicingfactors1[sf]
        cassetteSFvar.loc[j,'varseq']=sequence.toseq()
        cassetteSFvar.loc[j,'first_ss']=str('endogenous')
        cassetteSFvar.loc[j,'second_ss']=str(sf + ' [+6:+15]')
        j=j+1

        cassetteSFvar=cassetteSFvar.append(cassetteSF.loc[i], ignore_index=True)
        sequence = cassette_filtered.sequence[i]\
            [100+int(cassette_filtered.intronup_length[i]) + diff +30 - 150: \
            100+int(cassette_filtered.intronup_length[i]) + diff +30]\
            .tomutable()
        sequence[120 + 15:120 + 24]=splicingfactors1[sf]
        cassetteSFvar.loc[j,'varseq']=sequence.toseq()
        cassetteSFvar.loc[j,'first_ss']=str('endogenous')
        cassetteSFvar.loc[j,'second_ss']=str(sf + ' [+15:+24]')
        j=j+1



for i in cassetteSFvar.index:
    if (cassetteSFvar.varseq[i][122-int(cassetteSFvar.exonlength[i]):120].translate().find("*")>-1):
        print(i)
        print(cassetteSFvar.varseq[i][122-int(cassetteSFvar.exonlength[i]):120].translate())


from Bio.Restriction import RsrII
from Bio.Restriction import AscI
from Bio.Restriction import SpeI
from Bio.Restriction import AatII


for i in cassetteSFvar.index:
    if (RsrII.search(cassetteSFvar.varseq[i])!=[])|(AscI.search(cassetteSFvar.varseq[i])!=[])|(SpeI.search(cassetteSFvar.varseq[i])!=[])|(AatII.search(cassetteSFvar.varseq[i])!=[]):
        cassetteSFvar.drop(i,inplace=True)
        


cassetteSFvar.to_pickle('./design/cassette_splicingfactors_location.pkl')


#%% 

cassetteSFcomb=cassetteSF.loc[(cassetteSF.exonlength==50)]

cassetteSFcombvar=pd.DataFrame(columns=list(cassetteSFcomb.columns))


splicingfactorscomb=pd.Series()
splicingfactorscomb['native']=Seq('NNNNNNNNN').tomutable()
splicingfactorscomb['SRSF1']=Seq('TCACACGAC').tomutable()
splicingfactorscomb['SRSF5']=Seq('TTCACAGGC').tomutable()
splicingfactorscomb['hnRNPA1']=Seq('TTAGGGAAC').tomutable()
splicingfactorscomb['hnRNPU']=Seq('TTGTATTGC').tomutable()

k=0


for i in cassetteSFcomb.index:
    startpositions=[12,30,75,99]
    for sf in splicingfactorscomb.index: 
        sequence = cassetteconst.sequence[i]\
            [100+int(cassetteconst.intronup_length[i]) + 50 +30 - 150: \
            100+int(cassetteconst.intronup_length[i]) + 50 +30]\
            .tomutable()
        if (sf!="native"):
            first=sf
            firstpos=startpositions[0]
            sequence[startpositions[0]:startpositions[0]+9] = splicingfactorscomb.loc[sf]
            for sf in splicingfactorscomb.index:
                if (sf!="native"):
                    sequence[startpositions[1]:startpositions[1]+9] = splicingfactorscomb.loc[sf]
                    second=sf
                    secondpos=startpositions[1]
                    cassetteSFcombvar=cassetteSFcombvar.append(cassetteSFcomb.loc[i], ignore_index=True)        
                    cassetteSFcombvar.loc[k,'varseq']=sequence.toseq()
                    cassetteSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                    cassetteSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                    k=k+1
                    
                else:
                    for sf in splicingfactorscomb.index:
                        if (sf!="native"):
                            sequence[startpositions[2]:startpositions[2]+9] = splicingfactorscomb.loc[sf]
                            second=sf
                            secondpos=startpositions[2]
                            cassetteSFcombvar=cassetteSFcombvar.append(cassetteSFcomb.loc[i], ignore_index=True)        
                            cassetteSFcombvar.loc[k,'varseq']=sequence.toseq()
                            cassetteSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                            cassetteSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                            k=k+1
                    else:
                        for sf in splicingfactorscomb.index:
                            if (sf!="native"):
                                sequence[startpositions[3]:startpositions[3]+9] = splicingfactorscomb.loc[sf]
                                second=sf
                                secondpos=startpositions[3]
                                cassetteSFcombvar=cassetteSFcombvar.append(cassetteSFcomb.loc[i], ignore_index=True)        
                                cassetteSFcombvar.loc[k,'varseq']=sequence.toseq()
                                cassetteSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                                cassetteSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
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
                            cassetteSFcombvar=cassetteSFcombvar.append(cassetteSFcomb.loc[i], ignore_index=True)        
                            cassetteSFcombvar.loc[k,'varseq']=sequence.toseq()
                            cassetteSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                            cassetteSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                            k=k+1
                        else:
                            for sf in splicingfactorscomb.index:
                                if (sf!="native"):
                                    sequence[startpositions[3]:startpositions[3]+9] = splicingfactorscomb.loc[sf]
                                    second=sf
                                    secondpos=startpositions[3]
                                    cassetteSFcombvar=cassetteSFcombvar.append(cassetteSFcomb.loc[i], ignore_index=True)        
                                    cassetteSFcombvar.loc[k,'varseq']=sequence.toseq()
                                    cassetteSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                                    cassetteSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
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
                                    cassetteSFcombvar=cassetteSFcombvar.append(cassetteSFcomb.loc[i], ignore_index=True)        
                                    cassetteSFcombvar.loc[k,'varseq']=sequence.toseq()
                                    cassetteSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                                    cassetteSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                                    k=k+1
            
for i in cassetteSFcombvar.index:
    if (cassetteSFcombvar.varseq[i][72:120].translate().find("*")>-1):
        print(str(i))
        print(cassetteSFcombvar.varseq[i][72:120].translate())


            
            
cassetteSFcombvar.to_pickle('./design/cassette_splicingfactors_combinatorial_exononly.pkl')


 
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

relative to first splice site
[-46:-40]
[-37:-31]
[8:14]
[20:23] 
[diff-18:diff-12] 
[diff-9:diff-3]
[diff+9:diff+15]
[diff+18:diff+24]	

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

cassetteSFRosenberg=pd.DataFrame(columns=list(cassetteSFcomb.columns))
j=0
for i in cassetteSFcomb.index:
    for sf in splicingfactors2.index:
        cassetteSFRosenberg=cassetteSFRosenberg.append(cassetteSFcomb.loc[i], ignore_index=True)
        sequence = cassetteconst.sequence[i]\
            [100+int(cassetteconst.intronup_length[i]) + 50 +30 - 150: \
            100+int(cassetteconst.intronup_length[i]) + 50 +30]\
            .tomutable()
        sequence[24:30]=splicingfactors2[sf]
        cassetteSFRosenberg.loc[j,'varseq']=sequence.toseq()
        cassetteSFRosenberg.loc[j,'first_ss']=str(sf + ' [-46:-40]')
        cassetteSFRosenberg.loc[j,'second_ss']=str('endogenous')
        j=j+1
        cassetteSFRosenberg=cassetteSFRosenberg.append(cassetteSFcomb.loc[i], ignore_index=True)
        sequence = cassetteconst.sequence[i]\
            [100+int(cassetteconst.intronup_length[i]) + 50 +30 - 150: \
            100+int(cassetteconst.intronup_length[i]) + 50 +30]\
            .tomutable()
        sequence[33:39]=splicingfactors2[sf]
        cassetteSFRosenberg.loc[j,'varseq']=sequence.toseq()
        cassetteSFRosenberg.loc[j,'first_ss']=str(sf + ' [-37:-31]')
        cassetteSFRosenberg.loc[j,'second_ss']=str('endogenous')
        j=j+1
        cassetteSFRosenberg=cassetteSFRosenberg.append(cassetteSFcomb.loc[i], ignore_index=True)
        sequence = cassetteconst.sequence[i]\
            [100+int(cassetteconst.intronup_length[i]) + 50 +30 - 150: \
            100+int(cassetteconst.intronup_length[i]) + 50 +30]\
            .tomutable()
        sequence[78:84]=splicingfactors2[sf]
        cassetteSFRosenberg.loc[j,'varseq']=sequence.toseq()
        cassetteSFRosenberg.loc[j,'first_ss']=str(sf + ' [8:14]')
        cassetteSFRosenberg.loc[j,'second_ss']=str('endogenous')
        j=j+1
        cassetteSFRosenberg=cassetteSFRosenberg.append(cassetteSFcomb.loc[i], ignore_index=True)
        sequence = cassetteconst.sequence[i]\
            [100+int(cassetteconst.intronup_length[i]) + 50 +30 - 150: \
            100+int(cassetteconst.intronup_length[i]) + 50 +30]\
            .tomutable()
        sequence[87:93]=splicingfactors2[sf]
        cassetteSFRosenberg.loc[j,'varseq']=sequence.toseq()
        cassetteSFRosenberg.loc[j,'first_ss']=str(sf + ' [20:23]')
        cassetteSFRosenberg.loc[j,'second_ss']=str('endogenous')
        j=j+1
        cassetteSFRosenberg=cassetteSFRosenberg.append(cassetteSFcomb.loc[i], ignore_index=True)
        sequence = cassetteconst.sequence[i]\
            [100+int(cassetteconst.intronup_length[i]) + 50 +30 - 150: \
            100+int(cassetteconst.intronup_length[i]) + 50 +30]\
            .tomutable()
        sequence[102:108]=splicingfactors2[sf]
        cassetteSFRosenberg.loc[j,'varseq']=sequence.toseq()
        cassetteSFRosenberg.loc[j,'first_ss']=str('endogenous')
        cassetteSFRosenberg.loc[j,'second_ss']=str(sf + ' [diff-18:diff-12]')
        j=j+1
        cassetteSFRosenberg=cassetteSFRosenberg.append(cassetteSFcomb.loc[i], ignore_index=True)
        sequence = cassetteconst.sequence[i]\
            [100+int(cassetteconst.intronup_length[i]) + 50 +30 - 150: \
            100+int(cassetteconst.intronup_length[i]) + 50 +30]\
            .tomutable()
        sequence[111:117]=splicingfactors2[sf]
        cassetteSFRosenberg.loc[j,'varseq']=sequence.toseq()
        cassetteSFRosenberg.loc[j,'first_ss']=str('endogenous')
        cassetteSFRosenberg.loc[j,'second_ss']=str(sf + ' [diff-9:diff-3]')
        j=j+1
        cassetteSFRosenberg=cassetteSFRosenberg.append(cassetteSFcomb.loc[i], ignore_index=True)
        sequence = cassetteconst.sequence[i]\
            [100+int(cassetteconst.intronup_length[i]) + 50 +30 - 150: \
            100+int(cassetteconst.intronup_length[i]) + 50 +30]\
            .tomutable()
        sequence[129:135]=splicingfactors2[sf]
        cassetteSFRosenberg.loc[j,'varseq']=sequence.toseq()
        cassetteSFRosenberg.loc[j,'first_ss']=str('endogenous')
        cassetteSFRosenberg.loc[j,'second_ss']=str(sf + ' [diff+9:diff+15]')
        j=j+1
        cassetteSFRosenberg=cassetteSFRosenberg.append(cassetteSFcomb.loc[i], ignore_index=True)
        sequence = cassetteconst.sequence[i]\
            [100+int(cassetteconst.intronup_length[i]) + 50 +30 - 150: \
            100+int(cassetteconst.intronup_length[i]) + 50 +30]\
            .tomutable()
        sequence[138:144]=splicingfactors2[sf]
        cassetteSFRosenberg.loc[j,'varseq']=sequence.toseq()
        cassetteSFRosenberg.loc[j,'first_ss']=str('endogenous')
        cassetteSFRosenberg.loc[j,'second_ss']=str(sf + ' [diff+18:diff+24]')
        j=j+1

for i in cassetteSFRosenberg.index:
    if (cassetteSFRosenberg.varseq[i][72:120].translate().find("*")>-1):
        print(str(i))
        print(cassetteSFRosenberg.varseq[i][72:120].translate())
        
cassetteSFRosenberg.to_pickle('./design/cassette_splicingfactors_Rosenberg.pkl')

#%% RNA secondary structure

addforsec=cassette_filtered.loc[[28911,693744,397830,647738,416417,247730,172080,364621,219263]]

cassettesec=pd.concat([cassette35,cassette47,cassette50,cassette71,addforsec])


'''
-30-37:-30-28
-30-diff+5:-30-diff+14
-30-14:-30-5
-15:-6
'''

cassettesecvar=pd.DataFrame(columns=list(cassettesec.columns))
k=0
for i in cassettesec.index:
    diff=int(cassettesec.exonlength[i])
    firstssrev = cassettesec.varseq[i][-30-diff-9:-30-diff].reverse_complement()
    firstsscomp = cassettesec.varseq[i][-30-diff-9:-30-diff].complement()
    firstssrev2 = cassettesec.varseq[i][-30-diff-12:-30-diff-3].reverse_complement()
    firstsscomp2 = cassettesec.varseq[i][-30-diff-12:-30-diff-3].complement()

    altexonrev=cassettesec.varseq[i][-30-diff + 17:-30-diff + 26].reverse_complement()
    altexoncomp=cassettesec.varseq[i][-30-diff + 17:-30-diff + 26].complement()
    exonrev=cassettesec.varseq[i][-30-24:-30-15].reverse_complement() 
    exoncomp=cassettesec.varseq[i][-30-24:-30-15].complement()
    
    secondssrev = cassettesec.varseq[i][-30:-21].reverse_complement()
    secondsscomp = cassettesec.varseq[i][-30:-21].complement()
    secondssrev2 = cassettesec.varseq[i][-27:-18].reverse_complement()
    secondsscomp2 = cassettesec.varseq[i][-27:-18].complement()
    

    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-diff-37:-30-diff-28]=firstssrev
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='rev at [-30-37:-30-28]'
    cassettesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-diff-37:-30-diff-28]=firstsscomp
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='comp at [-30-37:-30-28]'
    cassettesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-diff-37:-30-diff-28]=firstssrev2
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='rev2 at [-30-37:-30-28]'
    cassettesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-diff-37:-30-diff-28]=firstsscomp2
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='comp2 at [-30-37:-30-28]'
    cassettesecvar.loc[k,'secondss']='endogenous'
    k=k+1

    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-diff+5:-30-diff+14]=firstssrev
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='rev at [-30-diff+5:-30-diff+14]'
    cassettesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-diff+5:-30-diff+14]=firstsscomp
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='comp at [-30-diff+5:-30-diff+14]'
    cassettesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-diff+5:-30-diff+14]=firstssrev2
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='rev2 at [-30-diff+5:-30-diff+14]'
    cassettesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-diff+5:-30-diff+14]=firstsscomp2
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='comp2 at [-30-diff+5:-30-diff+14]'
    cassettesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-diff+5:-30-diff+14]=altexonrev
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='altexonrev at [-30-diff+5:-30-diff+14]'
    cassettesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-diff+5:-30-diff+14]=altexoncomp
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='altexoncomp at [-30-diff+5:-30-diff+14]'
    cassettesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    
    
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-12:-30-3]=secondssrev
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='endogenous'
    cassettesecvar.loc[k,'secondss']='rev at [-30-14:-30-5]'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-12:-30-3]=secondsscomp
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='endogenous'
    cassettesecvar.loc[k,'secondss']='comp at [-30-14:-30-5]'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-12:-30-3]=secondssrev2
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='endogenous'
    cassettesecvar.loc[k,'secondss']='rev2 at [-30-14:-30-5]'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-12:-30-3]=secondsscomp2
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='endogenous'
    cassettesecvar.loc[k,'secondss']='comp2 at [-30-14:-30-5]'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-12:-30-3]=exonrev
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='endogenous'
    cassettesecvar.loc[k,'secondss']='exonrev at [-30-14:-30-5]'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-30-12:-30-3]=exoncomp
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='endogenous'
    cassettesecvar.loc[k,'secondss']='exoncomp at [-30-14:-30-5]'
    k=k+1

    
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-15:-6]=secondssrev
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='endogenous'
    cassettesecvar.loc[k,'secondss']='rev at [-15:-6]'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-15:-6]=secondsscomp
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='endogenous'
    cassettesecvar.loc[k,'secondss']='comp at [-15:-6]'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-15:-6]=secondssrev2
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='endogenous'
    cassettesecvar.loc[k,'secondss']='rev2 at [-15:-6]'
    k=k+1
    cassettesecvar=cassettesecvar.append(cassettesec.loc[i],ignore_index=True)    
    sequence=cassettesec.varseq[i].tomutable()
    sequence[-15:-6]=secondsscomp2
    cassettesecvar.varseq[k]=sequence.toseq()
    cassettesecvar.loc[k,'firstss']='endogenous'
    cassettesecvar.loc[k,'secondss']='comp2 at [-15:-6]'
    k=k+1
    
for i in cassettesecvar.index:
    if (cassettesecvar.varseq[i][122-int(cassettesecvar.exonlength[i]):120].translate().find("*")>-1):
        print(str(i))
        print(cassettesecvar.varseq[i][122-int(cassettesecvar.exonlength[i]):120].translate())
        cassettesecvar.drop(i,inplace=True)
        
        
       
cassettesecvar.to_pickle('./design/cassette_secondarystructure_variants.pkl')

#%% constitutive exons

junctions=pd.read_pickle('./design/junctions.pkl')

constexons=pd.DataFrame()


for i in range(len(junctions)-1):
    if (junctions.loc[junctions.index[i],'end']<junctions.loc[junctions.index[i+1],'start']):
        constexons= constexons.append(junctions.iloc[i,:])
        constexons.loc[junctions.index[i],'exonstart']= \
            junctions.loc[junctions.index[i],'end']
        constexons.loc[junctions.index[i],'exonend']= \
            junctions.loc[junctions.index[i+1],'start']
            
        constexons.loc[junctions.index[i],'score1']= \
            junctions.loc[junctions.index[i],'score']
        constexons.loc[junctions.index[i],'score2']= \
            junctions.loc[junctions.index[i+1],'score']
        constexons.loc[junctions.index[i],'level1']= \
            junctions.loc[junctions.index[i],'level']
        constexons.loc[junctions.index[i],'level2']= \
            junctions.loc[junctions.index[i+1],'level']

        constexons.loc[junctions.index[i],'sig1']= \
            junctions.loc[junctions.index[i],'sig']
        constexons.loc[junctions.index[i],'sig2']= \
            junctions.loc[junctions.index[i+1],'sig']

        constexons.loc[junctions.index[i],'exonlength']= \
            constexons.loc[junctions.index[i],'exonend']- \
            constexons.loc[junctions.index[i],'exonstart']
            
constexons.to_pickle('./design/constexons.pkl')
constexons=pd.read_pickle('./design/constexons.pkl')

constexonssmall=constexons[constexons.exonlength<90]

constexon35=constexons[constexons.exonlength==35]
constexon47=constexons[constexons.exonlength==47]
constexon50=constexons[constexons.exonlength==50]
constexon71=constexons[constexons.exonlength==71]

cassette35=cassette_filtered[cassette_filtered.exonlength==35]
cassette47=cassette_filtered[cassette_filtered.exonlength==47]
cassette50=cassette_filtered[cassette_filtered.exonlength==50]
#cassette68=cassette_filtered[cassette_filtered.exonlength==68]
cassette71=cassette_filtered[cassette_filtered.exonlength==71]
cassette50.drop(492702,inplace=True)
cassette71.drop(561749,inplace=True)

cassette74=cassette_filtered[cassette_filtered.exonlength==74]
constexon74=constexons[constexons.exonlength==74]

cassette62=cassette_filtered[cassette_filtered.exonlength==62]
constexon62=constexons[constexons.exonlength==62]


f=open('./design/constcas.bed', 'w')
for i in constexonssmall.index:
    if (constexonssmall.loc[i,'strand']=='+'):
        f.write(constexonssmall.loc[i,'chr'] + "\t" + str(int(constexonssmall.loc[i,'exonend'])-200) + \
            "\t" + str(int(constexonssmall.loc[i,'exonend'])+100) + "\t" + constexonssmall.loc[i,'name'] + \
            "\t" + str(constexonssmall.loc[i,'exonlength']) + "\t" + constexonssmall.loc[i,'strand'] + "\n")
    else:
        f.write(constexonssmall.loc[i,'chr'] + "\t" + str(int(constexonssmall.loc[i,'exonstart'])-100) + \
            "\t" + str(int(constexonssmall.loc[i,'exonstart'])+200) + "\t" + constexonssmall.loc[i,'name'] + \
            "\t" + str(constexonssmall.loc[i,'exonlength']) + "\t" + constexonssmall.loc[i,'strand'] + "\n")

f.close()


'''
bedtools getfasta -fi hg19.fa -bed constcas.bed -fo constcasseqaroundsplicesite.fa -s

'''     

#parse fasta file

from Bio import SeqIO

f=open('./design/constcasseqaroundsplicesite.fa','r')
records=list(SeqIO.parse(f, 'fasta'))

c=0
for i in constexonssmall.index:
    constexonssmall.loc[i,'coordinates']=records[c].id
    constexonssmall.loc[i,'sequence']=records[c].seq
    c=c+1


cassette35=cassette_filtered[cassette_filtered.exonlength==35]
cassette47=cassette_filtered[cassette_filtered.exonlength==47]
cassette50=cassette_filtered[cassette_filtered.exonlength==50]
#cassette68=cassette_filtered[cassette_filtered.exonlength==68]
cassette71=cassette_filtered[cassette_filtered.exonlength==71]
cassette50.drop(492702,inplace=True)
cassette71.drop(561749,inplace=True)

cassette74=cassette_filtered[cassette_filtered.exonlength==74]
cassette62=cassette_filtered[cassette_filtered.exonlength==62]


for i in constexonssmall.index:
    if (constexonssmall.exonlength[i]%3==2):
        constexonssmall.loc[i,'nostopinplus2']=bool(constexonssmall.sequence[i][-100-int(constexonssmall.exonlength[i])+2:-100].translate().find("*")==-1)

for i in constexonssmall.index:
    if (constexonssmall.nostopinplus2[i]==True):
        exonlength=int(constexonssmall.exonlength[i])
        constexonssmall.loc[i,'varseq']=constexonssmall.sequence[i][-220:-70]
        constexonssmall.loc[i,'intronup']=constexonssmall.sequence[i][-220:-100-exonlength]
        constexonssmall.loc[i,'exon']=constexonssmall.sequence[i][-100-exonlength:-100]
        constexonssmall.loc[i,'introndown']=constexonssmall.sequence[i][-100:-70]
    else:
        constexonssmall.drop(i,inplace=True)

constexon62=constexonssmall[(constexonssmall.exonlength==62)&(constexonssmall.nostopinplus2==True)]

constexon74=constexonssmall[(constexonssmall.exonlength==74)&(constexonssmall.nostopinplus2==True)]

constexon71=constexonssmall[(constexonssmall.exonlength==71)&(constexonssmall.nostopinplus2==True)]
constexon50=constexonssmall[(constexonssmall.exonlength==50)&(constexonssmall.nostopinplus2==True)]
constexon47=constexonssmall[(constexonssmall.exonlength==47)&(constexonssmall.nostopinplus2==True)]
constexon50.drop(8995,inplace=True)
constexon47.drop(11096,inplace=True)
constexon71.drop([3971,5279,7641,28555],inplace=True)
constexon35=constexonssmall[(constexonssmall.exonlength==35)&(constexonssmall.nostopinplus2==True)]


import forlibrary

comb62=forlibrary.make_combinatorial(cassette62,constexon62)

comb50=forlibrary.make_combinatorial(cassette50,constexon50)

comb47=forlibrary.make_combinatorial(cassette47,constexon47)

comb71=forlibrary.make_combinatorial(cassette71,constexon71)

comb35=forlibrary.make_combinatorial(cassette35,constexon35)

cassette_comb_with_constitutive=pd.concat([comb35,comb47,comb50,comb62,comb71],ignore_index=True)

for i in cassette_comb_with_constitutive.index:
    if (cassette_comb_with_constitutive.varseq[i][122-int(cassette_comb_with_constitutive.exonlength[i]):120].translate().find("*")>-1):
        print(str(i))
        print(cassette_comb_with_constitutive.varseq[i][122-int(cassette_comb_with_constitutive.exonlength[i]):120].translate())
#        cassette_comb_with_constitutive.drop(i,inplace=True)
    
cassette_comb_with_constitutive.to_pickle('./design/cassette_comb_with_constitutive.pkl')

constidxs=[]
for i in cassette_comb_with_constitutive.combination.dropna().drop_duplicates():
    constidxs+=i.split('\t')
    
constexonssmall[[str(x) in pd.Series(constidxs).unique()[1:] for x in constexonssmall.index]][['chr','start','end','name','strand','exonlength']].to_csv('./tables/TableS10.csv')



#%% nucleosome positioning

# get codon table from script forrecoding.py

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

for i in cassette_filtered.index:
    cassette_filtered.loc[i,'aaseq']=cassette_filtered.loc[i,'varseq'][-30-int(cassette_filtered.exonlength[i])+2:-30].translate()

cassettenuc=pd.DataFrame(columns=list(cassette_filtered))

k=0
for i in cassette_filtered.index:
    cassettenuc=cassettenuc.append(cassette_filtered.loc[i],ignore_index=True)
    cassettenuc.loc[k,'changes']='endogenous'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassette_filtered.loc[i],ignore_index=True)
    recod=str()
    for j in cassette_filtered.aaseq[i][1:-2]:
        recod=recod + codonmaxgc[j][0]
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[-30- int(cassette_filtered.exonlength[i])+2+3:-30-6]=recod
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded maximal gc'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassette_filtered.loc[i],ignore_index=True)
    recod=str()
    for j in cassette_filtered.aaseq[i][1:-2]:
        recod=recod + codonmingc[j][0]
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[-30-int(cassette_filtered.exonlength[i])+2+3:-30-6]=recod
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded minimal gc'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassette_filtered.loc[i],ignore_index=True)
    recod=str()
    for j in cassette_filtered.aaseq[i][1:-2]:
        recod=recod + random.choice(gencode[j])
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[-30- int(cassette_filtered.exonlength[i])+2+3:-30-6]=recod
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded randomly 1'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassette_filtered.loc[i],ignore_index=True)
    recod=str()
    for j in cassette_filtered.aaseq[i][1:-2]:
        recod=recod + random.choice(gencode[j])
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[-30- int(cassette_filtered.exonlength[i])+2+3:-30-6]=recod
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded randomly 2'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassette_filtered.loc[i],ignore_index=True)
    recod=str()
    for j in cassette_filtered.aaseq[i][1:-2]:
        recod=recod + random.choice(gencode[j])
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[-30- int(cassette_filtered.exonlength[i])+2+3:-30-6]=recod
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded randomly 3'
    k=k+1
    

addfornuc=cassette_filtered.loc[[28911,693744,397830,647738,416417,247730,172080,364621,219263]]

cassettedisf=pd.concat([cassette35,cassette47,cassette50,cassette71,addfornuc])

imp22dT=Seq('TTTTTTTTTCATTTTTTTTTTT')
perf12dT=Seq('TTTTTTTTTTTT')
perf5dT=Seq('TTTTTC')
imp22ctrl=Seq('TTGGACTGTCATTCACGTCGTT')
perf12ctrl=Seq('TTGGACTGTCAT')
perf5ctrl=Seq('TTGGAC')


for i in cassettedisf.index:
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[3:25]=imp22dT
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='imperfect 22nt dT tract at position 0'
    k=k+1

    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[3:15]=perf12dT
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='perfect 12nt dT tract at position 0'
    k=k+1

    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[9:21]=perf12dT
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='perfect 12nt dT tract at position 6'
    k=k+1

    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[15:27]=perf12dT
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='perfect 12nt dT tract at position 12'
    k=k+1
    
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[3:25]=imp22ctrl
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='control 22nt at position 0'
    k=k+1

    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[3:15]=perf12ctrl
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='control 12nt at position 0'
    k=k+1

    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[9:21]=perf12ctrl
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='control 12nt at position 6'
    k=k+1

    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[15:27]=perf12ctrl
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='control 12nt at position 12'
    k=k+1


    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+8:-30-int(cassettedisf.exonlength[i])+14]=perf5dT
    sequence[-30-12:-30-6]=perf5dT
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='perfect 5nt dT tract at position exonstart+8 and exondend-14'
    k=k+1

    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+11:-30-int(cassettedisf.exonlength[i])+17]=perf5dT
    sequence[-30-15:-30-9]=perf5dT
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='perfect 5nt dT tract at position exonstart+11 and exondend-17'
    k=k+1

    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+8:-30-int(cassettedisf.exonlength[i])+14]=perf5dT
    sequence[-30-12:-30-6]=perf5ctrl
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='control 5nt at position exonstart+8 and exondend-14'
    k=k+1

    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassette_filtered.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+11:-30-int(cassettedisf.exonlength[i])+17]=perf5dT
    sequence[-30-15:-30-9]=perf5ctrl
    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='control 5nt at position exonstart+11 and exondend-17'
    k=k+1


perT=Seq('TTC')
perC=Seq('CCA')

for i in cassettedisf.index:
    sequenceT=cassette_filtered.varseq[i].tomutable()
    sequenceC=cassette_filtered.varseq[i].tomutable()
    count= int((cassettedisf.exonlength[i]-2)//9)
    for j in range(0,count,1):
        sequenceT[-30-int(cassettedisf.exonlength[i])+5 + 9*j:-30-int(cassettedisf.exonlength[i])+5 + 9*j+3]=perT
        sequenceC[-30-int(cassettedisf.exonlength[i])+5 + 9*j:-30-int(cassettedisf.exonlength[i])+5 + 9*j+3]=perC
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    cassettenuc.varseq[k]=sequenceT.toseq()
    cassettenuc.loc[k,'changes']='9nt periodicity of TT'
    k=k+1
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    cassettenuc.varseq[k]=sequenceC.toseq()
    cassettenuc.loc[k,'changes']='9nt periodicity of CC'
    k=k+1

    sequenceT=cassette_filtered.varseq[i].tomutable()
    sequenceC=cassette_filtered.varseq[i].tomutable()
    count= int((cassettedisf.exonlength[i]-2)//10)
    for j in range(0,count,1):
        sequenceT[-30-int(cassettedisf.exonlength[i])+5 + 10*j:-30-int(cassettedisf.exonlength[i])+5 + 10*j+3]=perT
        sequenceC[-30-int(cassettedisf.exonlength[i])+5 + 10*j:-30-int(cassettedisf.exonlength[i])+5 + 10*j+3]=perC
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    cassettenuc.varseq[k]=sequenceT.toseq()
    cassettenuc.loc[k,'changes']='10nt periodicity of TT'
    k=k+1
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    cassettenuc.varseq[k]=sequenceC.toseq()
    cassettenuc.loc[k,'changes']='10nt periodicity of CC'
    k=k+1

    sequenceT=cassette_filtered.varseq[i].tomutable()
    sequenceC=cassette_filtered.varseq[i].tomutable()
    count= int((cassettedisf.exonlength[i]-2)//11)
    for j in range(0,count,1):
        sequenceT[-30-int(cassettedisf.exonlength[i])+5 + 11*j:-30-int(cassettedisf.exonlength[i])+5 + 11*j+3]=perT
        sequenceC[-30-int(cassettedisf.exonlength[i])+5 + 11*j:-30-int(cassettedisf.exonlength[i])+5 + 11*j+3]=perC
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    cassettenuc.varseq[k]=sequenceT.toseq()
    cassettenuc.loc[k,'changes']='11nt periodicity of TT'
    k=k+1
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    cassettenuc.varseq[k]=sequenceC.toseq()
    cassettenuc.loc[k,'changes']='11nt periodicity of CC'
    k=k+1


# if something needs to be removed, get rid of the "second half of exon" one, and potentially make the list of contexts smaller
for i in cassettedisf.index:
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    count= int((cassettedisf.exonlength[i]-2)//9)
    for j in range(0,count,1):
        sequence[-30-int(cassettedisf.exonlength[i])+6 + 9*j:-30-int(cassettedisf.exonlength[i])+6 + 9*j+2]=Seq('CG')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='CG in exon every 9nt'
    k=k+1
    
   
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    recod=str()
    for j in cassettedisf.aaseq[i][1:-2]:
        recod=recod + codonmingc[j][0]
    sequence=cassettedisf.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+2+3:-30-6]=recod
    count= int((cassettedisf.exonlength[i]-2)//9)
    for j in range(0,count,1):
        sequence[-30-int(cassettedisf.exonlength[i])+6 + 9*j:-30-int(cassettedisf.exonlength[i])+6 + 9*j+2]=Seq('CG')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded minimal gc, CG in exon every 9nt'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    count= int((cassettedisf.exonlength[i]-2)//6)
    for j in range(0,count,1):
        sequence[-30-int(cassettedisf.exonlength[i])+6 + 6*j:-30-int(cassettedisf.exonlength[i])+6 + 6*j+2]=Seq('CG')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='CG in exon every 6nt'
    k=k+1
    
   
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    recod=str()
    for j in cassettedisf.aaseq[i][1:-2]:
        recod=recod + codonmingc[j][0]
    sequence=cassettedisf.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+2+3:-30-6]=recod
    count= int((cassettedisf.exonlength[i]-2)//6)
    for j in range(0,count,1):
        sequence[-30-int(cassettedisf.exonlength[i])+6 + 6*j:-30-int(cassettedisf.exonlength[i])+6 + 6*j+2]=Seq('CG')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded minimal gc, CG in exon every 6nt'
    k=k+1

    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    for j in range(0,4,1):
        sequence[-30-5-18 + 6*j:-30-5-18 + 6*j+2]=Seq('CG')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='CG in last 20 nt of the exon every 6nt'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+2+3:-30-6]=recod
    for j in range(0,4,1):
        sequence[-30-5-18 + 6*j:-30-5-18 + 6*j+2]=Seq('CG')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded minimal GC, CG in last 20 nt of the exon every 6nt'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    for j in range(0,7,1):
        sequence[-30-5-18 + 3*j:-30-5-18 + 3*j+2]=Seq('CG')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='CG in last 20 nt of the exon every 3nt'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+2+3:-30-6]=recod
    for j in range(0,4,1):
        sequence[-30-5-18 + 3*j:-30-5-18 + 3*j+2]=Seq('CG')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded minimal GC, CG in last 20 nt of the exon every 3nt'
    k=k+1

    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    count= int((150-30-cassettedisf.exonlength[i]-28)//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('CG')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='CG in upstream intron every 6nt'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+2+3:-30-6]=recod
    count= int((150-30-cassettedisf.exonlength[i]-28)//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('CG')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded minimal GC, CG in upstream intron every 6nt'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    count= int((150-30-cassettedisf.exonlength[i]-28)//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('CG')
    count= int((cassettedisf.exonlength[i]-2)//6)
    for j in range(0,count,1):
        sequence[-30-int(cassettedisf.exonlength[i])+6 + 6*j:-30-int(cassettedisf.exonlength[i])+6 + 6*j+2]=Seq('CG')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='CG in upstream intron and in exon every 6nt'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+2+3:-30-6]=recod
    count= int((150-30-cassettedisf.exonlength[i]-28)//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('CG')
    count= int((cassettedisf.exonlength[i]-2)//6)
    for j in range(0,count,1):
        sequence[-30-int(cassettedisf.exonlength[i])+6 + 6*j:-30-int(cassettedisf.exonlength[i])+6 + 6*j+2]=Seq('CG')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded minimal GC, CG in upstream intron and in exon every 6nt'
    k=k+1
    
# same for GC

    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    count= int((cassettedisf.exonlength[i]-2)//9)
    for j in range(0,count,1):
        sequence[-30-int(cassettedisf.exonlength[i])+6 + 9*j:-30-int(cassettedisf.exonlength[i])+6 + 9*j+2]=Seq('GC')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='GC in exon every 9nt'
    k=k+1
    
   
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    recod=str()
    for j in cassettedisf.aaseq[i][1:-2]:
        recod=recod + codonmingc[j][0]
    sequence=cassettedisf.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+2+3:-30-6]=recod
    count= int((cassettedisf.exonlength[i]-2)//9)
    for j in range(0,count,1):
        sequence[-30-int(cassettedisf.exonlength[i])+6 + 9*j:-30-int(cassettedisf.exonlength[i])+6 + 9*j+2]=Seq('GC')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded minimal gc, GC in exon every 9nt'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    count= int((cassettedisf.exonlength[i]-2)//6)
    for j in range(0,count,1):
        sequence[-30-int(cassettedisf.exonlength[i])+6 + 6*j:-30-int(cassettedisf.exonlength[i])+6 + 6*j+2]=Seq('GC')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='GC in exon every 6nt'
    k=k+1
    
   
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    recod=str()
    for j in cassettedisf.aaseq[i][1:-2]:
        recod=recod + codonmingc[j][0]
    sequence=cassettedisf.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+2+3:-30-6]=recod
    count= int((cassettedisf.exonlength[i]-2)//6)
    for j in range(0,count,1):
        sequence[-30-int(cassettedisf.exonlength[i])+6 + 6*j:-30-int(cassettedisf.exonlength[i])+6 + 6*j+2]=Seq('GC')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded minimal gc, GC in exon every 6nt'
    k=k+1

    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    for j in range(0,4,1):
        sequence[-30-5-18 + 6*j:-30-5-18 + 6*j+2]=Seq('GC')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='GC in last 20 nt of the exon every 6nt'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+2+3:-30-6]=recod
    for j in range(0,4,1):
        sequence[-30-5-18 + 6*j:-30-5-18 + 6*j+2]=Seq('GC')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded minimal GC, GC in last 20 nt of the exon every 6nt'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    for j in range(0,7,1):
        sequence[-30-5-18 + 3*j:-30-5-18 + 3*j+2]=Seq('GC')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='GC in last 20 nt of the exon every 3nt'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+2+3:-30-6]=recod
    for j in range(0,4,1):
        sequence[-30-5-18 + 3*j:-30-5-18 + 3*j+2]=Seq('GC')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded minimal GC, GC in last 20 nt of the exon every 3nt'
    k=k+1

    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    count= int((150-30-cassettedisf.exonlength[i]-28)//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('GC')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='GC in upstream intron every 6nt'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+2+3:-30-6]=recod
    count= int((150-30-cassettedisf.exonlength[i]-28)//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('GC')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded minimal GC, GC in upstream intron every 6nt'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    count= int((150-30-cassettedisf.exonlength[i]-28)//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('GC')
    count= int((cassettedisf.exonlength[i]-2)//6)
    for j in range(0,count,1):
        sequence[-30-int(cassettedisf.exonlength[i])+6 + 6*j:-30-int(cassettedisf.exonlength[i])+6 + 6*j+2]=Seq('GC')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='GC in upstream intron and in exon every 6nt'
    k=k+1
    
    cassettenuc=cassettenuc.append(cassettedisf.loc[i],ignore_index=True)
    sequence=cassettedisf.varseq[i].tomutable()
    sequence[-30-int(cassettedisf.exonlength[i])+2+3:-30-6]=recod
    count= int((150-30-cassettedisf.exonlength[i]-28)//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('GC')
    count= int((cassettedisf.exonlength[i]-2)//6)
    for j in range(0,count,1):
        sequence[-30-int(cassettedisf.exonlength[i])+6 + 6*j:-30-int(cassettedisf.exonlength[i])+6 + 6*j+2]=Seq('GC')

    cassettenuc.varseq[k]=sequence.toseq()
    cassettenuc.loc[k,'changes']='exon recoded minimal GC, GC in upstream intron and in exon every 6nt'
    k=k+1
    
for i in cassettenuc.index:
    if (cassettenuc.varseq[i][122-int(cassettenuc.exonlength[i]):120].translate().find("*")>-1):
        print(str(i))
        print(cassettenuc.changes[i])
        print(cassettenuc.varseq[i][122-int(cassettenuc.exonlength[i]):120].translate())
#        cassettesecvar.drop(i,inplace=True)
    
cassettenuc.drop_duplicates('varseq',inplace=True)

from Bio.Restriction import RsrII
from Bio.Restriction import AscI
from Bio.Restriction import SpeI
from Bio.Restriction import AatII


for i in cassettenuc.index:
    if (RsrII.search(cassettenuc.varseq[i])!=[])|(AscI.search(cassettenuc.varseq[i])!=[])|(SpeI.search(cassettenuc.varseq[i])!=[])|(AatII.search(cassettenuc.varseq[i])!=[]):
        cassettenuc.drop(i,inplace=True)
        

cassettenuc.to_pickle('./design/cassette_nucleosome_recoding_methylation.pkl')

#%% test different splice site sequences

cassetteswitch=pd.concat([cassette35,cassette47,cassette50,cassette71])

cassetteswitch.drop(325472,inplace=True)

cassetteswitchvar=pd.DataFrame()

k=0

for i in cassetteswitch.index:
    exonlength=int(cassetteswitch.exonlength[i])
    
    for j in cassetteswitch.index:
        cassetteswitchvar=cassetteswitchvar.append(cassetteswitch.loc[i],ignore_index=True)
        sequence=cassetteswitch.varseq[i].tomutable()
        sequence[-30-exonlength-15:-30-exonlength+3]= \
            cassetteswitch.varseq[j][-30-int(cassetteswitch.exonlength[j])-15:\
            -30-int(cassetteswitch.exonlength[j]) +3]
        cassetteswitchvar.varseq[k]=sequence.toseq()
        cassetteswitchvar.loc[k,'changes']='acceptor from ' + str(j)
        k=k+1
        
        cassetteswitchvar=cassetteswitchvar.append(cassetteswitch.loc[i],ignore_index=True)
        sequence=cassetteswitch.varseq[i].tomutable()
        sequence[-33:-24]=cassetteswitch.varseq[j][-33:-24]
        cassetteswitchvar.varseq[k]=sequence.toseq()
        cassetteswitchvar.loc[k,'changes']='donor from ' + str(j)
        k=k+1
        
        cassetteswitchvar=cassetteswitchvar.append(cassetteswitch.loc[i],ignore_index=True)
        sequence=cassetteswitch.varseq[i].tomutable()
        sequence[-30-exonlength-15:-30-exonlength+3]= \
            cassetteswitch.varseq[j][-30-int(cassetteswitch.exonlength[j])-15:\
            -30-int(cassetteswitch.exonlength[j]) +3]
        sequence[-33:-24]=cassetteswitch.varseq[j][-33:-24]
        cassetteswitchvar.varseq[k]=sequence.toseq()
        cassetteswitchvar.loc[k,'changes']='acceptor and donor from ' + str(j)
        k=k+1
        
for i in cassetteswitchvar.index:
    if (cassetteswitchvar.varseq[i][122-int(cassetteswitchvar.exonlength[i]):120].translate().find("*")>-1):
        print(str(i))
        print(cassetteswitchvar.changes[i])
        print(cassetteswitchvar.varseq[i][122-int(cassetteswitchvar.exonlength[i]):120].translate())
        cassetteswitchvar.drop(i,inplace=True)
    
cassetteswitchvar.to_pickle('./design/cassette_switch_splice_site_sequences.pkl')

#%% length

cassettelength= pd.concat([cassette35,cassette47,cassette50,cassette71])

cassettelengthvar=pd.DataFrame()

k=0

for i in cassettelength.index:
    exonlength=int(cassettelength.exonlength[i])
    seqext = cassettelength.sequence[i]\
        [100+int(cassettelength.intronup_length[i]) + int(cassettelength.exonlength[i]) +30 - 200: \
        100+int(cassettelength.intronup_length[i]) + int(cassettelength.exonlength[i]) +30]
        
    
    count = (exonlength-5-6)//6
    print(count)
    for ex in range(1,count,1):
        
        cassettelengthvar=cassettelengthvar.append(cassettelength.loc[i],ignore_index=True)
        cassettelengthvar.varseq[k]=seqext[50-ex*6:-30-exonlength+5] + seqext[-30-exonlength+5 +ex*6:]
        cassettelengthvar.loc[k,'exonlength']=str(exonlength - ex*6)
        cassettelengthvar.loc[k,'changes']=str(ex*6) + 'nt deletion at position 5'
        k=k+1

        cassettelengthvar=cassettelengthvar.append(cassettelength.loc[i],ignore_index=True)
        middle = -30 - ((exonlength//3)//2)*3
        cassettelengthvar.varseq[k]=seqext[50-ex*6:middle-ex*3] + seqext[middle+ex*3:]
        cassettelengthvar.loc[k,'exonlength']=str(exonlength - ex*6)
        cassettelengthvar.loc[k,'changes']=str(ex*6) + 'nt deletion at center'
        k=k+1
        
        cassettelengthvar=cassettelengthvar.append(cassettelength.loc[i],ignore_index=True)
        cassettelengthvar.varseq[k]=seqext[50-ex*6:-30-6 -ex*6] + seqext[-30-6:]
        cassettelengthvar.loc[k,'exonlength']=str(exonlength - ex*6)
        cassettelengthvar.loc[k,'changes']=str(ex*6) + 'nt deletion at position -6'
        k=k+1

for i in cassettelengthvar.index:
    if (cassettelengthvar.varseq[i][122-int(cassettelengthvar.exonlength[i]):120].translate().find("*")>-1):
        print(str(i))
        print(cassettelengthvar.changes[i])
        print(cassettelengthvar.varseq[i][122-int(cassettelengthvar.exonlength[i]):120].translate())
#        cassettelengthvar.drop(i,inplace=True)

for i in cassettelengthvar.index:
    if (len(cassettelengthvar.varseq[i])!=150):
        print(i)
        print(cassettelengthvar.subset[i])        
        cassettelengthvar.drop(i,inplace=True)
        
    
cassettelengthvar.to_pickle('./design/cassette_exonlength.pkl')
