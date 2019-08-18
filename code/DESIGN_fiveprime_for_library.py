# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 12:48:10 2016

@author: miklm
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq


#%%

five=pd.read_pickle('./design/jun5prselected.pkl')

five[['score','score2','score_alt','level','level_alt']]= \
    five[['score','score2','score_alt','level','level_alt']].astype(int)

for i in five.index:
    five.loc[i,'levelratio']=np.log2(float(five.level_alt[i])/float(five.level[i]))

five_filtered=five[(five.level_alt>20) & (five.level>20) & (five.levelratio>-5) & (five.levelratio<5)]

five_filtered.drop_duplicates('commonname',inplace=True)

five_filtered.drop(476039,inplace=True) #no evidence for splicing

for i in five_filtered.index:
    five_filtered.loc[i,'varseq']=str(five_filtered.sequence[i])[99:261].upper()

for i in five_filtered.index:
    if (five_filtered.loc[i,'varseq'][:51+int(five_filtered.diff_nt[i])].translate().find("*") > -1):
        print(str(i))
        print(five_filtered.loc[i,'varseq'][:51+int(five_filtered.diff_nt[i])].translate())

five_filtered.to_pickle('./design/five_filtered.pkl')

five_filtered['donor1']=51
five_filtered['donor2']=51+five_filtered.diff_nt.astype(int)

five_filtered[['commonname','coordinates','start','start_alt','diff_nt','donor1','donor2','varseq']].to_csv('./tables/TableS3_fivefiltered.csv')


#%%

'''
make 4 groups, within which "exonic - alternative - intronic" sequences are
varied by permutation

The four groups are: number of genes in the group - diff_nt
7x 8
7x 17
5x 26
4x 41
'''

five8=five_filtered[five_filtered.diff_nt==8]
five8.drop(476039,inplace=True) #no evidence for splicing
five17=five_filtered[five_filtered.diff_nt==17]
five26=five_filtered[five_filtered.diff_nt==26]
five41=five_filtered[five_filtered.diff_nt==41]

df=five8.copy()
for i in df.index:
    df.loc[i,'exonseq']=df.sequence[i][99:150]
    df.loc[i,'altseq']=df.sequence[i][150:150+ int(df.diff_nt[i])]
    df.loc[i,'intronseq']=df.sequence[i][150+ int(df.diff_nt[i]):261]
five8=df.copy()

df=five17.copy()
for i in df.index:
    df.loc[i,'exonseq']=df.sequence[i][99:150]
    df.loc[i,'altseq']=df.sequence[i][150:150+ int(df.diff_nt[i])]
    df.loc[i,'intronseq']=df.sequence[i][150+ int(df.diff_nt[i]):261]

five17=df.copy()

df=five26.copy()
for i in df.index:
    df.loc[i,'exonseq']=df.sequence[i][99:150]
    df.loc[i,'altseq']=df.sequence[i][150:150+ int(df.diff_nt[i])]
    df.loc[i,'intronseq']=df.sequence[i][150+ int(df.diff_nt[i]):261]
five26=df.copy()

df=five41.copy()
for i in df.index:
    df.loc[i,'exonseq']=df.sequence[i][99:150]
    df.loc[i,'altseq']=df.sequence[i][150:150+ int(df.diff_nt[i])]
    df.loc[i,'intronseq']=df.sequence[i][150+ int(df.diff_nt[i]):261]
five41=df.copy()


# mutate downstream 3'ss from CAG to GAC
five8.intronseq[351211]=Seq('GTTGCGTGGTGGAGTGAGGGCAGGGCCTGGGTCTGGGGGCGGGGCCTGGGTCTGACCCCAGGGACTGGTTTTTGTGTTGACGTGTCTCTGTGGTGACGAGCAG')

# drop variant with downstream 3'ss
five17.drop(114892,inplace=True)
#%%
# create all possible combinations within each group
df=five8.copy()
out=pd.DataFrame()
l=0
for i in df.index:
    for j in df.index:
        for k in df.index:
            out.loc[l,'varseq']=df.loc[i,'exonseq']+df.loc[j,'altseq']+df.loc[k,'intronseq']
            out.loc[l,'geneexonseq']=df.commonname[i]
            out.loc[l,'genealtseq']=df.commonname[j]
            out.loc[l,'geneintronseq']=df.commonname[k]
            out.loc[l,'diff_nt']=df.diff_nt[j]
            l=l+1
five8combvar=out.copy()

df=five17.copy()
out=pd.DataFrame()
l=0
for i in df.index:
    for j in df.index:
        for k in df.index:
            out.loc[l,'varseq']=df.loc[i,'exonseq']+df.loc[j,'altseq']+df.loc[k,'intronseq']
            out.loc[l,'geneexonseq']=df.commonname[i]
            out.loc[l,'genealtseq']=df.commonname[j]
            out.loc[l,'geneintronseq']=df.commonname[k]
            out.loc[l,'diff_nt']=df.diff_nt[j]
            l=l+1
five17combvar=out.copy()

df=five26.copy()
out=pd.DataFrame()
l=0
for i in df.index:
    for j in df.index:
        for k in df.index:
            out.loc[l,'varseq']=df.loc[i,'exonseq']+df.loc[j,'altseq']+df.loc[k,'intronseq']
            out.loc[l,'geneexonseq']=df.commonname[i]
            out.loc[l,'genealtseq']=df.commonname[j]
            out.loc[l,'geneintronseq']=df.commonname[k]
            out.loc[l,'diff_nt']=df.diff_nt[j]

            l=l+1
five26combvar=out.copy()

df=five41.copy()
out=pd.DataFrame()
l=0
for i in df.index:
    for j in df.index:
        for k in df.index:
            out.loc[l,'varseq']=df.loc[i,'exonseq']+df.loc[j,'altseq']+df.loc[k,'intronseq']
            out.loc[l,'geneexonseq']=df.commonname[i]
            out.loc[l,'genealtseq']=df.commonname[j]
            out.loc[l,'geneintronseq']=df.commonname[k]
            out.loc[l,'diff_nt']=df.diff_nt[j]
            l=l+1
five41combvar=out.copy()

fiveprime_combinatorial_variations=pd.concat([five8combvar,five17combvar,five26combvar,five41combvar],ignore_index=True)

for i in fiveprime_combinatorial_variations.index:
    fiveprime_combinatorial_variations.loc[i,'nostopbeforesplicesite']=bool(fiveprime_combinatorial_variations.varseq[i][:50+int(fiveprime_combinatorial_variations.loc[i,'diff_nt'])].translate().find("*")==-1)

### PICKLE
fiveprime_combinatorial_variations.to_pickle('./design/fiveprime_combinatorial_variations.pkl')

#%%
def parse_five(df):
    '''parse a sequence from teh 5prime collection and return the dataframe 
    including sequences for the exonic (exonseq), the alternative part (altseq)
    and the intronic part (intronseq)
    '''
    
    for i in df.index:
        df.loc[i,'exonseq']=df.sequence[i][100:150]
        df.loc[i,'altseq']=df.sequence[i][150:150+ int(df.diff_nt[i])]
        df.loc[i,'intronseq']=df.sequence[i][150+ int(df.diff_nt[i]):262]
    return df[:]

#%% replace with constitutive splice sites

#select everything larger than 15 diff_nt and smaller than 50 diff_nt and change to constitutive splice sites and splice sites mutants

fiveconst = five_filtered[(five_filtered.diff_nt > 15) & (five_filtered.diff_nt < 50)]
fiveconst.drop(54910,inplace=True)

canonical_donor = Seq('CAGGTAAGT').tomutable()
no_donor = Seq('CTGCTC').tomutable()
GC_donor = Seq('CAGGCAAGT').tomutable()


fiveconstvar=pd.DataFrame(columns=list(fiveconst.columns))
j=0
for i in fiveconst.index:
    diff=int(fiveconst.diff_nt[i])
    fiveconstvar=fiveconstvar.append(fiveconst.loc[i], ignore_index=True)
    fiveconstvar.loc[j,'varseq'] = fiveconst.sequence[i][99:261]
    fiveconstvar.loc[j,'first_ss']=str('endogenous')
    fiveconstvar.loc[j,'second_ss']=str('endogenous')
    j=j+1
    fiveconstvar=fiveconstvar.append(fiveconst.loc[i], ignore_index=True)
    sequence = fiveconst.sequence[i][99:261].tomutable()
    sequence[48:57]=canonical_donor
    fiveconstvar.loc[j,'varseq']=sequence.toseq()
    fiveconstvar.loc[j,'first_ss']=str('constitutive')
    fiveconstvar.loc[j,'second_ss']=str('endogenous')
    j=j+1
    fiveconstvar=fiveconstvar.append(fiveconst.loc[i], ignore_index=True)
    sequence = fiveconst.sequence[i][99:261].tomutable()
    sequence[48+diff:57+diff]=canonical_donor
    fiveconstvar.loc[j,'varseq']=sequence.toseq()
    fiveconstvar.loc[j,'first_ss']=str('endogenous')
    fiveconstvar.loc[j,'second_ss']=str('constitutive')
    j=j+1
    fiveconstvar=fiveconstvar.append(fiveconst.loc[i], ignore_index=True)
    sequence = fiveconst.sequence[i][99:261].tomutable()
    sequence[48:57]=canonical_donor
    sequence[48+diff:57+diff]=canonical_donor
    fiveconstvar.loc[j,'varseq']=sequence.toseq()
    fiveconstvar.loc[j,'first_ss']=str('constitutive')
    fiveconstvar.loc[j,'second_ss']=str('constitutive')
    j=j+1
    fiveconstvar=fiveconstvar.append(fiveconst.loc[i], ignore_index=True)
    sequence = fiveconst.sequence[i][99:261].tomutable()
    sequence[51:57]=no_donor
    fiveconstvar.loc[j,'varseq']=sequence.toseq()
    fiveconstvar.loc[j,'first_ss']=str('unsplicable')
    fiveconstvar.loc[j,'second_ss']=str('endogenous')
    j=j+1
    fiveconstvar=fiveconstvar.append(fiveconst.loc[i], ignore_index=True)
    sequence = fiveconst.sequence[i][99:261].tomutable()
    sequence[51+diff:57+diff]=no_donor
    fiveconstvar.loc[j,'varseq']=sequence.toseq()
    fiveconstvar.loc[j,'first_ss']=str('endogenous')
    fiveconstvar.loc[j,'second_ss']=str('unsplicable')
    j=j+1
    fiveconstvar=fiveconstvar.append(fiveconst.loc[i], ignore_index=True)
    sequence = fiveconst.sequence[i][99:261].tomutable()
    sequence[48:57]=GC_donor
    fiveconstvar.loc[j,'varseq']=sequence.toseq()
    fiveconstvar.loc[j,'first_ss']=str('GC')
    fiveconstvar.loc[j,'second_ss']=str('endogenous')
    j=j+1
    fiveconstvar=fiveconstvar.append(fiveconst.loc[i], ignore_index=True)
    sequence = fiveconst.sequence[i][99:261].tomutable()
    sequence[48+diff:57+diff]=GC_donor
    fiveconstvar.loc[j,'varseq']=sequence.toseq()
    fiveconstvar.loc[j,'first_ss']=str('endogenous')
    fiveconstvar.loc[j,'second_ss']=str('GC')
    j=j+1


fiveconstvar.to_pickle('./design/fiveprime_constitutiveandunsplicable_variations.pkl')

#%% SF binding sites

fiveSF=five_filtered.loc[[242748, 309946,541250, 344238,305982,481167,596032,388363,70023,195774,320413]]

'''
YuI1	D: GGGCCACTTGGA	sil	INTRONIC (+11 to +22)
YuI2	B: CGCTGGTCATTC	sil	
YuI3	C: GAGGTATCAGCTT	sil	
YuI4	A: CGTTAGAGTAGC	sil	
YuE1	F: CTTAATTTTAGT	sil	EXONIC (-18 to -7)
YuE2	E: TAGTTTAGTTAG	sil	
'''

YuI1 = Seq('GGGCCACTTGGA').tomutable()
YuI2 = Seq('CGCTGGTCATTC').tomutable()
#YuI3 = Seq('GAGGTATCAGCTT').tomutable()
#YuI4 = Seq('CGTTAGAGTAGC').tomutable()
YuE1 = Seq('CTTAATTTTAGT').tomutable()
#YuE2 = Seq('TAGTTTAGTTAG').tomutable()

# I excluded the ones that (can) lead to a stop codon when incorporated at the indicated position


fiveSFvar=pd.DataFrame(columns=list(fiveSF.columns))

j=0
for i in fiveSF.index:
    diff=int(fiveSF.diff_nt[i])
    fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
    fiveSFvar.loc[j,'varseq'] = fiveSF.sequence[i][99:261]
    fiveSFvar.loc[j,'first_ss']=str('endogenous')
    fiveSFvar.loc[j,'second_ss']=str('endogenous')
    j=j+1
    fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
    sequence = fiveSF.sequence[i][99:261].tomutable()
    sequence[61:73]=YuI1
    fiveSFvar.loc[j,'varseq']=sequence.toseq()
    fiveSFvar.loc[j,'first_ss']=str('YuI1')
    fiveSFvar.loc[j,'second_ss']=str('endogenous')
    j=j+1
    fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
    sequence = fiveSF.sequence[i][99:261].tomutable()
    sequence[61:73]=YuI2
    fiveSFvar.loc[j,'varseq']=sequence.toseq()
    fiveSFvar.loc[j,'first_ss']=str('YuI2')
    fiveSFvar.loc[j,'second_ss']=str('endogenous')
    j=j+1
    fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
    sequence = fiveSF.sequence[i][99:261].tomutable()
    sequence[61:73]=YuI1
    sequence[61+diff:73+diff]=YuI1
    fiveSFvar.loc[j,'varseq']=sequence.toseq()
    fiveSFvar.loc[j,'first_ss']=str('YuI1')
    fiveSFvar.loc[j,'second_ss']=str('YuI1')
    j=j+1
    fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
    sequence = fiveSF.sequence[i][99:261].tomutable()
    sequence[61:73]=YuI2
    sequence[61+diff:73+diff]=YuI1
    fiveSFvar.loc[j,'varseq']=sequence.toseq()
    fiveSFvar.loc[j,'first_ss']=str('YuI2')
    fiveSFvar.loc[j,'second_ss']=str('YuI1')
    j=j+1
    fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
    sequence = fiveSF.sequence[i][99:261].tomutable()
    sequence[61:73]=YuI1
    sequence[61+diff:73+diff]=YuI2
    fiveSFvar.loc[j,'varseq']=sequence.toseq()
    fiveSFvar.loc[j,'first_ss']=str('YuI1')
    fiveSFvar.loc[j,'second_ss']=str('YuI2')
    j=j+1
    fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
    sequence = fiveSF.sequence[i][99:261].tomutable()
    sequence[61:73]=YuI2
    sequence[61+diff:73+diff]=YuI2
    fiveSFvar.loc[j,'varseq']=sequence.toseq()
    fiveSFvar.loc[j,'first_ss']=str('YuI2')
    fiveSFvar.loc[j,'second_ss']=str('YuI2')
    j=j+1
    fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
    sequence = fiveSF.sequence[i][99:261].tomutable()
    sequence[33:45]=YuE1
    fiveSFvar.loc[j,'varseq']=sequence.toseq()
    fiveSFvar.loc[j,'first_ss']=str('YuE1')
    fiveSFvar.loc[j,'second_ss']=str('endogenous')
    j=j+1
    fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
    sequence = fiveSF.sequence[i][99:261].tomutable()
    sequence[33:45]=YuE1
    sequence[33+diff:45+diff]=YuE1
    fiveSFvar.loc[j,'varseq']=sequence.toseq()
    fiveSFvar.loc[j,'first_ss']=str('YuE1')
    fiveSFvar.loc[j,'second_ss']=str('YuE1')
    j=j+1
    fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
    sequence = fiveSF.sequence[i][99:261].tomutable()
    sequence[33:45]=YuE1
    sequence[61+diff:73+diff]=YuI1
    fiveSFvar.loc[j,'varseq']=sequence.toseq()
    fiveSFvar.loc[j,'first_ss']=str('YuE1')
    fiveSFvar.loc[j,'second_ss']=str('YuI1')
    j=j+1
    fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
    sequence = fiveSF.sequence[i][99:261].tomutable()
    sequence[33:45]=YuE1
    sequence[61+diff:73+diff]=YuI2
    fiveSFvar.loc[j,'varseq']=sequence.toseq()
    fiveSFvar.loc[j,'first_ss']=str('YuE1')
    fiveSFvar.loc[j,'second_ss']=str('YuI2')
    j=j+1
    
    fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
    sequence = fiveSF.sequence[i][99:261].tomutable()
    sequence[61+diff:73+diff]=YuI1
    fiveSFvar.loc[j,'varseq']=sequence.toseq()
    fiveSFvar.loc[j,'first_ss']=str('endogenous')
    fiveSFvar.loc[j,'second_ss']=str('YuI1')
    j=j+1
    fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
    sequence = fiveSF.sequence[i][99:261].tomutable()
    sequence[61+diff:73+diff]=YuI2
    fiveSFvar.loc[j,'varseq']=sequence.toseq()
    fiveSFvar.loc[j,'first_ss']=str('endogenous')
    fiveSFvar.loc[j,'second_ss']=str('YuI2')
    j=j+1
    fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
    sequence = fiveSF.sequence[i][99:261].tomutable()
    sequence[33+diff:45+diff]=YuE1
    fiveSFvar.loc[j,'varseq']=sequence.toseq()
    fiveSFvar.loc[j,'first_ss']=str('endogenous')
    fiveSFvar.loc[j,'second_ss']=str('YuE1')
    j=j+1

#for i in fiveSFvar.index:
#    print(fiveSFvar.varseq[i][:51+int(fiveSFvar.diff_nt[i])].translate())

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
splicingfactors1['KeEnh']=Seq('CGACGTCGA').tomutable()
splicingfactors1['KeSil']=Seq('CCCAGCAGA').tomutable()
splicingfactors1['KeNeu']=Seq('CAAAGAGGA').tomutable()
splicingfactors1['hnRNPA1']=Seq('TTAGGGAAC').tomutable()
splicingfactors1['hnRNPG']=Seq('CAAGTGTTC').tomutable()
splicingfactors1['hnRNPU']=Seq('TTGTATTGC').tomutable()



for i in fiveSF.index:
    diff=int(fiveSF.diff_nt[i])
    for sf in splicingfactors1.index:
        fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
        sequence = fiveSF.sequence[i][99:261].tomutable()
        sequence[21:30]=splicingfactors1[sf]
        fiveSFvar.loc[j,'varseq']=sequence.toseq()
        fiveSFvar.loc[j,'first_ss']=str(sf + ' [-30:-21]')
        fiveSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
        sequence = fiveSF.sequence[i][99:261].tomutable()
        sequence[30:39]=splicingfactors1[sf]
        fiveSFvar.loc[j,'varseq']=sequence.toseq()
        fiveSFvar.loc[j,'first_ss']=str(sf + ' [-21:-12]')
        fiveSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
        sequence = fiveSF.sequence[i][99:261].tomutable()
        sequence[39:48]=splicingfactors1[sf]
        fiveSFvar.loc[j,'varseq']=sequence.toseq()
        fiveSFvar.loc[j,'first_ss']=str(sf + ' [-12:-3]')
        fiveSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
        sequence = fiveSF.sequence[i][99:261].tomutable()
        sequence[57:66]=splicingfactors1[sf]
        fiveSFvar.loc[j,'varseq']=sequence.toseq()
        fiveSFvar.loc[j,'first_ss']=str(sf + ' [6:15]')
        fiveSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        if (diff>30):
            fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
            sequence = fiveSF.sequence[i][99:261].tomutable()
            sequence[66:75]=splicingfactors1[sf]
            fiveSFvar.loc[j,'varseq']=sequence.toseq()
            fiveSFvar.loc[j,'first_ss']=str(sf + ' [15:24]')
            fiveSFvar.loc[j,'second_ss']=str('endogenous')
            j=j+1
        fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
        sequence = fiveSF.sequence[i][99:261].tomutable()
        sequence[51+diff-14:51+diff-5]=splicingfactors1[sf]
        fiveSFvar.loc[j,'varseq']=sequence.toseq()
        fiveSFvar.loc[j,'first_ss']=str('endogenous')
        fiveSFvar.loc[j,'second_ss']=str(sf + ' [diff-14:diff-5]')
        j=j+1
        fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
        sequence = fiveSF.sequence[i][99:261].tomutable()
        sequence[51+diff+6:51+diff+15]=splicingfactors1[sf]
        fiveSFvar.loc[j,'varseq']=sequence.toseq()
        fiveSFvar.loc[j,'first_ss']=str('endogenous')
        fiveSFvar.loc[j,'second_ss']=str(sf + ' [diff+6:diff+15]')
        j=j+1
        fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
        sequence = fiveSF.sequence[i][99:261].tomutable()
        sequence[51+diff+15:51+diff+24]=splicingfactors1[sf]
        fiveSFvar.loc[j,'varseq']=sequence.toseq()
        fiveSFvar.loc[j,'first_ss']=str('endogenous')
        fiveSFvar.loc[j,'second_ss']=str(sf + ' [diff+15:diff+24]')
        j=j+1
        fiveSFvar=fiveSFvar.append(fiveSF.loc[i], ignore_index=True)
        sequence = fiveSF.sequence[i][99:261].tomutable()
        sequence[51+diff+24:51+diff+33]=splicingfactors1[sf]
        fiveSFvar.loc[j,'varseq']=sequence.toseq()
        fiveSFvar.loc[j,'first_ss']=str('endogenous')
        fiveSFvar.loc[j,'second_ss']=str(sf + ' [diff+24:diff+33]')
        j=j+1


fiveSFvar.to_pickle('./design/five_splicingfactors_location.pkl')

#%%

fiveSFcomb=fiveSF.loc[[70023,195774,388363,481167,596032]]

fiveSFcombvar=pd.DataFrame(columns=list(fiveSFcomb.columns))

startpositions=[21,30,39,57,66,75]

splicingfactorscomb=pd.Series()
splicingfactorscomb['native']=Seq('NNNNNNNNN').tomutable()
splicingfactorscomb['SRSF1']=Seq('TCACACGAC').tomutable()
splicingfactorscomb['SRSF5']=Seq('TTCACAGGC').tomutable()
splicingfactorscomb['hnRNPA1']=Seq('TTAGGGAAC').tomutable()
splicingfactorscomb['hnRNPU']=Seq('TTGTATTGC').tomutable()

k=0
 
for i in fiveSFcomb.index:
    for sf in splicingfactorscomb.index: 
        sequence = fiveSFcomb.sequence[i][99:261].tomutable()   
        if (sf!="native"):
            first=sf
            firstpos=startpositions[0]
            sequence[startpositions[0]:startpositions[0]+9] = splicingfactorscomb.loc[sf]
            for sf in splicingfactorscomb.index:
                if (sf!="native"):
                    sequence[startpositions[1]:startpositions[1]+9] = splicingfactorscomb.loc[sf]
                    second=sf
                    secondpos=startpositions[1]
                    fiveSFcombvar=fiveSFcombvar.append(fiveSFcomb.loc[i], ignore_index=True)        
                    fiveSFcombvar.loc[k,'varseq']=sequence.toseq()
                    fiveSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                    fiveSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                    k=k+1
                    
                else:
                    for sf in splicingfactorscomb.index:
                        if (sf!="native"):
                            sequence[startpositions[2]:startpositions[2]+9] = splicingfactorscomb.loc[sf]
                            second=sf
                            secondpos=startpositions[2]
                            fiveSFcombvar=fiveSFcombvar.append(fiveSFcomb.loc[i], ignore_index=True)        
                            fiveSFcombvar.loc[k,'varseq']=sequence.toseq()
                            fiveSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                            fiveSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
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
                            fiveSFcombvar=fiveSFcombvar.append(fiveSFcomb.loc[i], ignore_index=True)        
                            fiveSFcombvar.loc[k,'varseq']=sequence.toseq()
                            fiveSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                            fiveSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                            k=k+1


for i in fiveSFcomb.index:
    startpositions=[57,66,51 + int(fiveSFcomb.diff_nt[i]) - 14]
    for sf in splicingfactorscomb.index: 
        sequence = fiveSFcomb.sequence[i][99:261].tomutable()   
        if (sf!="native"):
            first=sf
            firstpos=startpositions[0]
            sequence[startpositions[0]:startpositions[0]+9] = splicingfactorscomb.loc[sf]
            for sf in splicingfactorscomb.index:
                if (sf!="native"):
                    sequence[startpositions[1]:startpositions[1]+9] = splicingfactorscomb.loc[sf]
                    second=sf
                    secondpos=startpositions[1]
                    fiveSFcombvar=fiveSFcombvar.append(fiveSFcomb.loc[i], ignore_index=True)        
                    fiveSFcombvar.loc[k,'varseq']=sequence.toseq()
                    fiveSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                    fiveSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                    k=k+1
                    
                else:
                    for sf in splicingfactorscomb.index:
                        if (sf!="native"):
                            sequence[startpositions[2]:startpositions[2]+9] = splicingfactorscomb.loc[sf]
                            second=sf
                            secondpos=startpositions[2]
                            fiveSFcombvar=fiveSFcombvar.append(fiveSFcomb.loc[i], ignore_index=True)        
                            fiveSFcombvar.loc[k,'varseq']=sequence.toseq()
                            fiveSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                            fiveSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
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
                            fiveSFcombvar=fiveSFcombvar.append(fiveSFcomb.loc[i], ignore_index=True)        
                            fiveSFcombvar.loc[k,'varseq']=sequence.toseq()
                            fiveSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                            fiveSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                            k=k+1
            


for i in fiveSFcomb.index:
    startpositions=[57 + int(fiveSFcomb.diff_nt[i]),66 + int(fiveSFcomb.diff_nt[i]),75 + int(fiveSFcomb.diff_nt[i])]
    for sf in splicingfactorscomb.index: 
        sequence = fiveSFcomb.sequence[i][99:261].tomutable()   
        if (sf!="native"):
            first=sf
            firstpos=startpositions[0]
            sequence[startpositions[0]:startpositions[0]+9] = splicingfactorscomb.loc[sf]
            for sf in splicingfactorscomb.index:
                if (sf!="native"):
                    sequence[startpositions[1]:startpositions[1]+9] = splicingfactorscomb.loc[sf]
                    second=sf
                    secondpos=startpositions[1]
                    fiveSFcombvar=fiveSFcombvar.append(fiveSFcomb.loc[i], ignore_index=True)        
                    fiveSFcombvar.loc[k,'varseq']=sequence.toseq()
                    fiveSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                    fiveSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                    k=k+1
                    
                else:
                    for sf in splicingfactorscomb.index:
                        if (sf!="native"):
                            sequence[startpositions[2]:startpositions[2]+9] = splicingfactorscomb.loc[sf]
                            second=sf
                            secondpos=startpositions[2]
                            fiveSFcombvar=fiveSFcombvar.append(fiveSFcomb.loc[i], ignore_index=True)        
                            fiveSFcombvar.loc[k,'varseq']=sequence.toseq()
                            fiveSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                            fiveSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
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
                            fiveSFcombvar=fiveSFcombvar.append(fiveSFcomb.loc[i], ignore_index=True)        
                            fiveSFcombvar.loc[k,'varseq']=sequence.toseq()
                            fiveSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                            fiveSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                            k=k+1
            
            
fiveSFcombvar.to_pickle('./design/five_splicingfactors_combinatorial.pkl')
   
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

fiveSFRosenberg=pd.DataFrame(columns=list(fiveSFcomb.columns))
j=0
for i in fiveSF.index:
    diff=int(fiveSF.diff_nt[i])
    for sf in splicingfactors2.index:
        fiveSFRosenberg=fiveSFRosenberg.append(fiveSF.loc[i], ignore_index=True)
        sequence = fiveSF.sequence[i][99:261].tomutable()
        sequence[30:36]=splicingfactors2[sf]
        fiveSFRosenberg.loc[j,'varseq']=sequence.toseq()
        fiveSFRosenberg.loc[j,'first_ss']=str(sf + ' [-21:-15]')
        fiveSFRosenberg.loc[j,'second_ss']=str('endogenous')
        j=j+1
        fiveSFRosenberg=fiveSFRosenberg.append(fiveSF.loc[i], ignore_index=True)
        sequence = fiveSF.sequence[i][99:261].tomutable()
        sequence[39:45]=splicingfactors2[sf]
        fiveSFRosenberg.loc[j,'varseq']=sequence.toseq()
        fiveSFRosenberg.loc[j,'first_ss']=str(sf + ' [-12:-6]')
        fiveSFRosenberg.loc[j,'second_ss']=str('endogenous')
        j=j+1
        fiveSFRosenberg=fiveSFRosenberg.append(fiveSF.loc[i], ignore_index=True)
        sequence = fiveSF.sequence[i][99:261].tomutable()
        sequence[60:66]=splicingfactors2[sf]
        fiveSFRosenberg.loc[j,'varseq']=sequence.toseq()
        fiveSFRosenberg.loc[j,'first_ss']=str(sf + ' [9:15]')
        fiveSFRosenberg.loc[j,'second_ss']=str('endogenous')
        j=j+1
        fiveSFRosenberg=fiveSFRosenberg.append(fiveSF.loc[i], ignore_index=True)
        sequence = fiveSF.sequence[i][99:261].tomutable()
        sequence[51+diff-14:51+diff-8]=splicingfactors2[sf]
        fiveSFRosenberg.loc[j,'varseq']=sequence.toseq()
        fiveSFRosenberg.loc[j,'first_ss']=str('endogenous')
        fiveSFRosenberg.loc[j,'second_ss']=str(sf + ' [diff-14:diff-8]')
        j=j+1
        fiveSFRosenberg=fiveSFRosenberg.append(fiveSF.loc[i], ignore_index=True)
        sequence = fiveSF.sequence[i][99:261].tomutable()
        sequence[51+diff+9:51+diff+15]=splicingfactors2[sf]
        fiveSFRosenberg.loc[j,'varseq']=sequence.toseq()
        fiveSFRosenberg.loc[j,'first_ss']=str('endogenous')
        fiveSFRosenberg.loc[j,'second_ss']=str(sf + ' [diff+9:diff+15]')
        j=j+1
        fiveSFRosenberg=fiveSFRosenberg.append(fiveSF.loc[i], ignore_index=True)
        sequence = fiveSF.sequence[i][99:261].tomutable()
        sequence[51+diff+18:51+diff+24]=splicingfactors2[sf]
        fiveSFRosenberg.loc[j,'varseq']=sequence.toseq()
        fiveSFRosenberg.loc[j,'first_ss']=str('endogenous')
        fiveSFRosenberg.loc[j,'second_ss']=str(sf + ' [diff+18:diff+24]')
        j=j+1

fiveSFRosenberg.to_pickle('./design/five_splicingfactors_Rosenberg.pkl')

#%% secondary structure

fivesec = five_filtered[(five_filtered.diff_nt > 30)]


'''
51-12:51-3
51+15:51:51+24 
51+diff-14:51+diff-5
51+diff+15:51+diff+24
'''

fivesecvar=pd.DataFrame(columns=list(fivesec.columns))
k=0
for i in fivesec.index:
    diff=int(fivesec.diff_nt[i])
    exonrev=fivesec.varseq[i][51-24:51-15].reverse_complement() 
    exoncomp=fivesec.varseq[i][51-24:51-15].complement()

    firstssrev = fivesec.varseq[i][51:60].reverse_complement()
    firstsscomp = fivesec.varseq[i][51:60].complement()
    firstssrev2 = fivesec.varseq[i][54:63].reverse_complement()
    firstsscomp2 = fivesec.varseq[i][54:63].complement()

    altexonrev=fivesec.varseq[i][51+diff-26:51+diff-17].reverse_complement() 
    altexoncomp=fivesec.varseq[i][51+diff-26:51+diff-17].complement()
    
    secondssrev = fivesec.varseq[i][51+diff:51+diff+9].reverse_complement()
    secondsscomp = fivesec.varseq[i][51+diff:51+diff+9].complement()
    secondssrev2 = fivesec.varseq[i][51+diff+3:51+diff+12].reverse_complement()
    secondsscomp2 = fivesec.varseq[i][51+diff+3:51+diff+12].complement()
    

    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51-12:51-3]=firstssrev
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='rev at [51-12:51-3]'
    fivesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51-12:51-3]=firstsscomp
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='comp at [51-12:51-3]'
    fivesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51-12:51-3]=firstssrev2
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='rev2 at [51-12:51-3]'
    fivesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51-12:51-3]=firstsscomp2
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='comp2 at [51-12:51-3]'
    fivesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51-12:51-3]=exonrev
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='exonrev at [51-12:51-3]'
    fivesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51-12:51-3]=exoncomp
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='exonrev at [51-12:51-3]'
    fivesecvar.loc[k,'secondss']='endogenous'
    k=k+1


    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51+15:51+24]=firstssrev
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='rev at [51+15:51:51+24]'
    fivesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51+15:51+24]=firstsscomp
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='comp at [51+15:51:51+24]'
    fivesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51+15:51+24]=firstssrev2
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='rev2 at [51+15:51:51+24]'
    fivesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51+15:51+24]=firstsscomp2
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='comp2 at [51+15:51:51+24]'
    fivesecvar.loc[k,'secondss']='endogenous'
    k=k+1
    
    
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51+diff-14:51+diff-5]=altexonrev
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='endogenous'
    fivesecvar.loc[k,'secondss']='altexoncomp at [51+diff-14:51+diff-5]'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51+diff-14:51+diff-5]=altexoncomp
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='endogenous'
    fivesecvar.loc[k,'secondss']='altexoncomp at [51+diff-14:51+diff-5]'
    k=k+1  
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51+diff-14:51+diff-5]=secondssrev
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='endogenous'
    fivesecvar.loc[k,'secondss']='rev at [51+diff-14:51+diff-5]'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51+diff-14:51+diff-5]=secondsscomp
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='endogenous'
    fivesecvar.loc[k,'secondss']='comp at [51+diff-14:51+diff-5]'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51+diff-14:51+diff-5]=secondssrev2
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='endogenous'
    fivesecvar.loc[k,'secondss']='rev2 at [51+diff-14:51+diff-5]'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51+diff-14:51+diff-5]=secondsscomp2
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='endogenous'
    fivesecvar.loc[k,'secondss']='comp2 at [51+diff-14:51+diff-5]'
    k=k+1

    
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51+diff+15:51+diff+24]=secondssrev
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='endogenous'
    fivesecvar.loc[k,'secondss']='rev at [51+diff+15:51+diff+24]'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51+diff+15:51+diff+24]=secondsscomp
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='endogenous'
    fivesecvar.loc[k,'secondss']='comp at [51+diff+15:51+diff+24]'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51+diff+15:51+diff+24]=secondssrev2
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='endogenous'
    fivesecvar.loc[k,'secondss']='rev2 at [51+diff+15:51+diff+24]'
    k=k+1
    fivesecvar=fivesecvar.append(fivesec.loc[i],ignore_index=True)    
    sequence=fivesec.varseq[i].tomutable()
    sequence[51+diff+15:51+diff+24]=secondsscomp2
    fivesecvar.varseq[k]=sequence.toseq()
    fivesecvar.loc[k,'firstss']='endogenous'
    fivesecvar.loc[k,'secondss']='comp2 at [51+diff+15:51+diff+24]'
    k=k+1
    
for i in fivesecvar.index:
    if (fivesecvar.loc[i,'varseq'][:51+int(fivesecvar.diff_nt[i])].translate().find("*") > -1):
        print(str(i))
        print(fivesecvar.loc[i,'varseq'][:51+int(fivesecvar.diff_nt[i])].translate())
        fivesecvar.drop(i,inplace=True)

        
       
fivesecvar.to_pickle('./design/five_secondarystructure_variants.pkl')


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

for i in five_filtered.index:
    five_filtered.loc[i,'aaseq']=five_filtered.loc[i,'varseq'][:51+int(five_filtered.diff_nt[i])].translate()

fivenuc=pd.DataFrame(columns=list(five_filtered))

k=0
for i in five_filtered.index:
    diff=int(five_filtered.diff_nt[i])
    if (diff<25):
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        fivenuc.loc[k,'changes']='endogenous'
        k=k+1
        
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + codonmaxgc[j][0]
        sequence=five_filtered.varseq[i].tomutable()
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon recoded maximal gc'
        k=k+1
        
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + codonmingc[j][0]
        sequence=five_filtered.varseq[i].tomutable()
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon recoded minimal gc'
        k=k+1
        
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + random.choice(gencode[j])
        sequence=five_filtered.varseq[i].tomutable()
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon recoded randomly 1'
        k=k+1
    
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + random.choice(gencode[j])
        sequence=five_filtered.varseq[i].tomutable()
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon recoded randomly 2'
        k=k+1
    
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + random.choice(gencode[j])
        sequence=five_filtered.varseq[i].tomutable()
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon recoded randomly 3'
        k=k+1
    
    else:
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        fivenuc.loc[k,'changes']='endogenous'
        k=k+1
        
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + codonmaxgc[j][0]
        sequence=five_filtered.varseq[i].tomutable()
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon recoded maximal gc'
        k=k+1
        
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + codonmaxgc[j][0]
        sequence=five_filtered.varseq[i].tomutable()
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon recoded minimal gc'
        k=k+1
        
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + random.choice(gencode[j])
        sequence=five_filtered.varseq[i].tomutable()
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon recoded randomly 1'
        k=k+1
    
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + random.choice(gencode[j])
        sequence=five_filtered.varseq[i].tomutable()
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon recoded randomly 2'
        k=k+1
    
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + random.choice(gencode[j])
        sequence=five_filtered.varseq[i].tomutable()
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon recoded randomly 3'
        k=k+1
    




        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][19:19+diff//3-3]:
            recod=recod + codonmaxgc[j][0]
        sequence=five_filtered.varseq[i].tomutable()
        sequence[57:57+len(recod)]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='altexon recoded maximal gc'
        k=k+1
        
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + codonmaxgc[j][0]
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon and altexon recoded maximal gc'
        k=k+1
        
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][19:19+diff//3-3]:
            recod=recod + codonmingc[j][0]
        sequence=five_filtered.varseq[i].tomutable()
        sequence[57:57+len(recod)]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='altexon recoded minimal gc'
        k=k+1
        
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + codonmingc[j][0]
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon and altexon recoded minimal gc'
        k=k+1
        

        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][19:19+diff//3-3]:
            recod=recod + codonmingc[j][0]
        sequence=five_filtered.varseq[i].tomutable()
        sequence[57:57+len(recod)]=recod
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + codonmaxgc[j][0]
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon recoded maximal and altexon recoded minimal gc'
        k=k+1
        
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][19:19+diff//3-3]:
            recod=recod + codonmaxgc[j][0]
        sequence=five_filtered.varseq[i].tomutable()
        sequence[57:57+len(recod)]=recod
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + codonmaxgc[j][0]
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon recoded minimal and and altexon recoded maximal gc'
        k=k+1
        

        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][19:19+diff//3-3]:
            recod=recod + random.choice(gencode[j])
        sequence=five_filtered.varseq[i].tomutable()
        sequence[57:57+len(recod)]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='altexon recoded randomly 1'
        k=k+1
    
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + random.choice(gencode[j])
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon and altexon recoded randomly 1'
        k=k+1
    
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][19:19+diff//3-3]:
            recod=recod + random.choice(gencode[j])
        sequence=five_filtered.varseq[i].tomutable()
        sequence[57:57+len(recod)]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='altexon recoded randomly 2'
        k=k+1
    
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + random.choice(gencode[j])
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon and altexon recoded randomly 2'
        k=k+1
    
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][19:19+diff//3-3]:
            recod=recod + random.choice(gencode[j])
        sequence=five_filtered.varseq[i].tomutable()
        sequence[57:57+len(recod)]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='altexon recoded randomly 3'
        k=k+1
    
        fivenuc=fivenuc.append(five_filtered.loc[i],ignore_index=True)
        recod=str()
        for j in five_filtered.aaseq[i][0:16]:
            recod=recod + random.choice(gencode[j])
        sequence[0:48]=recod
        fivenuc.varseq[k]=sequence.toseq()
        fivenuc.loc[k,'changes']='exon and altexon recoded randomly 3'
        k=k+1
 
for i in fivenuc.index:
    if (fivenuc.varseq[i][:51+int(fivenuc.diff_nt[i])].translate().find("*")>-1):
        print(str(i))
        print(fivenuc.changes[i])
        print(fivenuc.varseq[i][:51+int(fivenuc.diff_nt[i])].translate())
#        threesecvar.drop(i,inplace=True)
 
fivenuc.drop_duplicates('varseq',inplace=True)
   
fivenuc.to_pickle('./design/five_nucleosome_recoding_methylation.pkl')
 

#%% combine with constitutive exon and intron

fiveconstforcomb=five_filtered[(five_filtered.diff_nt>11)&(five_filtered.diff_nt<60)]
const=pd.read_pickle('./design/constexon62.pkl')

const[['chr','exonend','name','strand']].to_csv('./tables/TableS11.csv')
const[['chr','exonstart','name','strand']].to_csv('./tables/TableS12.csv')


fiveconstcomb=pd.DataFrame()

k=0
for i in fiveconstforcomb.index:
    diff=int(fiveconstforcomb.diff_nt[i])
    for j in const.index:
        
        fiveconstcomb=fiveconstcomb.append(fiveconstforcomb.loc[i],ignore_index=True)
        sequence=fiveconstforcomb.varseq[i].tomutable()
        sequence[:51]=const.sequence[j][-100-51:-100]
        fiveconstcomb.varseq[k]=sequence.toseq()
        fiveconstcomb.loc[k,'changes']='exon ' + str(j)
        k=k+1
        
        fiveconstcomb=fiveconstcomb.append(fiveconstforcomb.loc[i],ignore_index=True)
        sequence=fiveconstforcomb.varseq[i].tomutable()
        sequence[51+diff:]=const.sequence[j][-100:-100+162-51-diff]
        fiveconstcomb.varseq[k]=sequence.toseq()
        fiveconstcomb.loc[k,'changes']='intron ' + str(j)
        k=k+1


for i in fiveconstcomb.index:
    if (fiveconstcomb.varseq[i][:51+int(fiveconstcomb.diff_nt[i])].translate().find("*")>-1):
        print(str(i))
        print(fiveconstcomb.changes[i])
        print(fiveconstcomb.varseq[i][:51+int(fiveconstcomb.diff_nt[i])].translate())
#        fiveconstcomb.drop(i,inplace=True)

for i in fiveconstcomb.index:
    if (len(fiveconstcomb.varseq[i])!=162):
        print(str(i))
        print(fiveconstcomb.changes[i])
        print(fiveconstcomb.varseq[i][:51+int(fiveconstcomb.diff_nt[i])].translate())
#        fiveconstcomb.drop(i,inplace=True)
  
fiveconstcomb.to_pickle('./design/five_combinatorial_with_constitutive.pkl')




#%% Switch elements between the two splice sites

fiveswitch=five_filtered[(five_filtered.diff_nt>11)]

#fiveswitch.drop(410585,inplace=True)

fiveswitchvar=pd.DataFrame(columns=list(fiveswitch.columns))

k=0
for i in fiveswitch.index:
    diff=int(fiveswitch.diff_nt[i])
    
    firstss=fiveswitch.varseq[i][51:57]
    secondss=fiveswitch.varseq[i][51+diff:57+diff]
    
    fiveswitchvar=fiveswitchvar.append(fiveswitch.loc[i],ignore_index=True)
    fiveswitchvar.loc[k,'firstss']='firstss 6nt'
    fiveswitchvar.loc[k,'secondss']='secondss 6nt'
    k=k+1
    
    fiveswitchvar=fiveswitchvar.append(fiveswitch.loc[i],ignore_index=True)
    sequence=fiveswitch.varseq[i].tomutable()
    sequence[51+diff:57+diff]=firstss
    fiveswitchvar.loc[k,'varseq']=sequence.toseq()
    fiveswitchvar.loc[k,'firstss']='firstss 6nt'
    fiveswitchvar.loc[k,'secondss']='firstss 6nt'
    k=k+1
        
    fiveswitchvar=fiveswitchvar.append(fiveswitch.loc[i],ignore_index=True)
    sequence=fiveswitch.varseq[i].tomutable()
    sequence[51:57]=secondss
    fiveswitchvar.loc[k,'varseq']=sequence.toseq()
    fiveswitchvar.loc[k,'firstss']='secondss 6nt'
    fiveswitchvar.loc[k,'secondss']='secondss 6nt'
    k=k+1
    
    fiveswitchvar=fiveswitchvar.append(fiveswitch.loc[i],ignore_index=True)
    sequence=fiveswitch.varseq[i].tomutable()
    sequence[51:57]=secondss
    sequence[51+diff:57+diff]=firstss
    fiveswitchvar.loc[k,'varseq']=sequence.toseq()
    fiveswitchvar.loc[k,'firstss']='secondss 6nt'
    fiveswitchvar.loc[k,'secondss']='firstss 6nt'
    k=k+1
    
#five_filtered.drop(410585,inplace=True)

fiveswitch=five_filtered[(five_filtered.diff_nt>16)&(five_filtered.diff_nt<55)]
     
for i in fiveswitch.index:
    if (fiveswitch.varseq[i][51 + int(fiveswitch.diff_nt[i]):51 + int(fiveswitch.diff_nt[i]) + int(fiveswitch.diff_nt[i])-6].translate().find('*')>-1):
        fiveswitch.drop(i,inplace=True)
        

for i in fiveswitch.index:
    diff=int(fiveswitch.diff_nt[i])
    if (fiveswitch.varseq[i][51+diff-9:51+diff].translate().find("*")==-1) & (fiveswitch.varseq[i][51-11:51+1].translate().find("*")==-1):
        for pos in range(0,3,1):
            fiveswitchvar=fiveswitchvar.append(fiveswitch.loc[i],ignore_index=True)
            sequence=fiveswitch.varseq[i].tomutable()
            sequence[48 - pos*3 +diff:51+diff]=sequence[48 - pos*3:51]
            fiveswitchvar.loc[k,'varseq']=sequence.toseq()
            fiveswitchvar.loc[k,'firstss']='firstss -' + str(3+pos*3) + 'nt'
            fiveswitchvar.loc[k,'secondss']='firstss -' + str(3+pos*3) + 'nt'
            k=k+1
                
            fiveswitchvar=fiveswitchvar.append(fiveswitch.loc[i],ignore_index=True)
            sequence=fiveswitch.varseq[i].tomutable()
            sequence[48 - pos*3:51]=sequence[48 - pos*3+diff:51+diff]
            fiveswitchvar.loc[k,'varseq']=sequence.toseq()
            fiveswitchvar.loc[k,'firstss']='secondss -' + str(3+pos*3) + 'nt'
            fiveswitchvar.loc[k,'secondss']='secondss -' + str(3+pos*3) + 'nt'
            k=k+1
            
            fiveswitchvar=fiveswitchvar.append(fiveswitch.loc[i],ignore_index=True)
            sequence=fiveswitch.varseq[i].tomutable()
            sequence[48 - pos*3:51]=sequence[48 - pos*3+diff:51+diff]
            sequence[48 - pos*3+diff:51+diff]=sequence[48 - pos*3:51]
            fiveswitchvar.loc[k,'varseq']=sequence.toseq()
            fiveswitchvar.loc[k,'firstss']='secondss -' + str(3+pos*3) + 'nt'
            fiveswitchvar.loc[k,'secondss']='firstss -' + str(3+pos*3) + 'nt'
            k=k+1
        
    for pos in range(0,(diff-12)/3,1):
        fiveswitchvar=fiveswitchvar.append(fiveswitch.loc[i],ignore_index=True)
        sequence=fiveswitch.varseq[i].tomutable()
        sequence[51+diff:60+pos*3 + diff]=sequence[51:60 + pos*3]
        fiveswitchvar.loc[k,'varseq']=sequence.toseq()
        fiveswitchvar.loc[k,'firstss']='firstss ' + str(9+pos*3) + 'nt'
        fiveswitchvar.loc[k,'secondss']='firstss ' + str(9+pos*3) + 'nt'
        k=k+1
            
        fiveswitchvar=fiveswitchvar.append(fiveswitch.loc[i],ignore_index=True)
        sequence=fiveswitch.varseq[i].tomutable()
        sequence[51:60 + pos*3]=sequence[51+diff:60+pos*3 + diff]
        fiveswitchvar.loc[k,'varseq']=sequence.toseq()
        fiveswitchvar.loc[k,'firstss']='secondss ' + str(9+pos*3) + 'nt'
        fiveswitchvar.loc[k,'secondss']='secondss ' + str(9+pos*3) + 'nt'
        k=k+1
        
        fiveswitchvar=fiveswitchvar.append(fiveswitch.loc[i],ignore_index=True)
        sequence=fiveswitch.varseq[i].tomutable()
        sequence[51:60 + pos*3]=sequence[51+diff:60+pos*3 + diff]
        sequence[51+diff:60+pos*3 + diff]=sequence[51:60 + pos*3]
        fiveswitchvar.loc[k,'varseq']=sequence.toseq()
        fiveswitchvar.loc[k,'firstss']='secondss ' + str(9+pos*3) + 'nt'
        fiveswitchvar.loc[k,'secondss']='firstss ' + str(9+pos*3) + 'nt'
        k=k+1
    
for i in fiveswitchvar.index:
    if (fiveswitchvar.varseq[i][:51+int(fiveswitchvar.diff_nt[i])].translate().find("*")>-1):
        print(str(i))
        print(fiveswitchvar.firstss[i])
        print(fiveswitchvar.secondss[i])
        print(fiveswitchvar.varseq[i][:51+int(fiveswitchvar.diff_nt[i])].translate())
        fiveswitchvar.drop(i,inplace=True)

for i in fiveswitchvar.index:
    if (len(fiveswitchvar.varseq[i])!=162):
        print(str(i))
        print(fiveswitchvar.changes[i])
        print(fiveswitchvar.varseq[i][:51+int(fiveswitchvar.diff_nt[i])].translate())
#        fiveswitchvar.drop(i,inplace=True)
  
fiveswitchvar.to_pickle('./design/five_switch_splice_site_sequences.pkl')



#%% length

fivelength= five_filtered[five_filtered.diff_nt>30]
fivelengthvar=pd.DataFrame()

k=0

for i in fivelength.index:
    diff=int(fivelength.diff_nt[i])
    seqext = fivelength.sequence[i][99:301]
        
    
    count = (diff-5-6)//6
    print(count)
    for ex in range(1,count,1):
        if (count<8):
            fivelengthvar=fivelengthvar.append(fivelength.loc[i],ignore_index=True)
            fivelengthvar.varseq[k]=seqext[:57] + seqext[57+ex*6:162+ex*6]
            fivelengthvar.loc[k,'diff_nt']=str(diff - ex*6)
            fivelengthvar.loc[k,'changes']=str(ex*6) + 'nt deletion at position 6'
            k=k+1
    
            fivelengthvar=fivelengthvar.append(fivelength.loc[i],ignore_index=True)
            middle = 57 + (((diff-5-6)//3)//2)*3
            fivelengthvar.varseq[k]=seqext[:middle-ex*3] + seqext[middle+ex*3:162+ex*6]
            fivelengthvar.loc[k,'diff_nt']=str(diff - ex*6)
            fivelengthvar.loc[k,'changes']=str(ex*6) + 'nt deletion at center'
            k=k+1
            
            fivelengthvar=fivelengthvar.append(fivelength.loc[i],ignore_index=True)
            fivelengthvar.varseq[k]=seqext[:51+diff-5 - ex*6] + seqext[51+diff-5:162+ex*6]
            fivelengthvar.loc[k,'diff_nt']=str(diff - ex*6)
            fivelengthvar.loc[k,'changes']=str(ex*6) + 'nt deletion at position -5'
            k=k+1

for i in fivelengthvar.index:
    if (len(fivelengthvar.varseq[i])!=162):
        print(str(i))
        print(fivelengthvar.changes[i])
        print(fivelengthvar.varseq[i][:51+int(fivelengthvar.diff_nt[i])].translate())
#        fiveswitchvar.drop(i,inplace=True)

for i in fivelengthvar.index:
    if (fivelengthvar.varseq[i][:51+int(fivelengthvar.diff_nt[i])].translate().find("*")>-1):
        print(str(i))
        print(fivelengthvar.varseq[i][:51+int(fivelengthvar.diff_nt[i])].translate())
#        fiveswitchvar.drop(i,inplace=True)
    
fivelengthvar.to_pickle('./design/five_exonlength.pkl')
    
