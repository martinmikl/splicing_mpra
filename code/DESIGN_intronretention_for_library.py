# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 14:45:30 2016

@author: martinm
"""

import pandas as pd

from Bio.Seq import Seq
import forlibrary

#%% Intron retention

ir=pd.read_pickle('./design/intronretentionvariants.pkl')

irfilt=ir[((ir.firstexonlength>80) & (ir.intronlength>58))]

irfilt.drop_duplicates(['name2'],inplace=True)

for i in irfilt.index:
    irfilt.intronlength[i]=int(irfilt.intronend_varseqnew[i])-int(irfilt.intronstart_varseq[i])
    
#for i in irfilt.index:
#    print(irfilt.name2[i])
#    print(irfilt.varseq[i][int(ir.intronstart_varseq[i])-4 : int(ir.intronstart_varseq[i])+4])
    
for i in irfilt.index:
    if ((irfilt.varseq[i][int(ir.intronstart_varseq[i]): int(ir.intronstart_varseq[i])+2]!='GT') \
            & (irfilt.varseq[i][int(ir.intronstart_varseq[i]) : int(ir.intronstart_varseq[i])+2]!='GC') & \
            (irfilt.varseq[i][int(ir.intronstart_varseq[i]) : int(ir.intronstart_varseq[i])+2]!='gc')):
        irfilt.drop(i, inplace=True)
        
irmain=irfilt[irfilt.intronend_varseqnew>135]

irmain.to_pickle('./design/irmain.pkl')

irmain['varseq']=irmain.varseq.apply(lambda x: str(x).upper())
irmain[['name2','coordinates','intronstart','intronend','intronlength','intronstart_varseq','intronend_varseqnew','varseq']].to_csv('./tables/TableS1_irmain.csv')


#%%

ir85=irmain[(irmain.intronlength==85)&(irmain.intronend_varseqnew==138)]
ir82=irmain[(irmain.intronlength==82)&(irmain.intronend_varseqnew==137)]
ir97=irmain[((irmain.intronlength==97)|(irmain.intronlength==100))&(irmain.intronend_varseqnew==137)]
for i in ir97.index:
    if (ir97.intronlength[i]==100):
        varseqvar=str('CAT') + str(ir97.varseq[i][:70]) + str(ir97.varseq[i][73:])
        ir97.varseq[i]=Seq(varseqvar)
        ir97.intronlength[i]==int(97)
        ir97.intronstart_varseq[i]=int(40)
        
irmain.drop(irmain[irmain.intronlength==85].index, inplace=True)
irmain.drop(irmain[irmain.intronlength==82].index, inplace=True)
irmain.drop(irmain[irmain.intronlength==97].index, inplace=True)

irmain_permutations=pd.DataFrame(columns=list(irmain.columns))
k=0

for i in irmain.index:
    intronsequence = irmain.varseq[i][int(irmain.intronstart_varseq[i]):int(irmain.intronend_varseqnew[i])]
    intronlength=len(intronsequence)
    for j in irmain.index:
        if(irmain.intronlength[i]%3==irmain.intronlength[j]%3):
            varseqvar=irmain.sequence[j][int(irmain.firstexonlength[j])-int(irmain.intronend_varseqnew[j]) + int(intronlength) :\
                int(irmain.firstexonlength[j])] \
                + intronsequence + irmain.varseq[j][int(irmain.intronend_varseqnew[j]):]
            irmain_permutations=irmain_permutations.append(irmain.loc[i],ignore_index=True)
            irmain_permutations.varseq[k]=varseqvar
            irmain_permutations.loc[k,'intron_gene']=irmain.loc[i,'name2']
            irmain_permutations.loc[k,'exons_gene']=irmain.loc[j,'name2']
            irmain_permutations.loc[k,'exons_gene_intronend']=irmain.loc[j,'intronend_varseqnew']
            k=k+1


for i in irmain_permutations.index:
    intronsequence = irmain_permutations.varseq[i][int(irmain_permutations.intronstart_varseq[i]):int(irmain_permutations.intronend_varseqnew[i])]
    intronlength=len(intronsequence)    
    irmain_permutations.loc[i,'splicedvarseqperm']=irmain_permutations.loc[i,'varseq'] \
        [:int(irmain_permutations.exons_gene_intronend[i])-intronlength] + irmain_permutations.loc[i,'varseq'] \
        [int(irmain_permutations.exons_gene_intronend[i]):]


for i in irmain_permutations.index:
    irmain_permutations.loc[i,'stopinsplicedperm0']=bool(irmain_permutations.splicedvarseqperm[i].translate().find("*")==-1)
    irmain_permutations.loc[i,'stopinunsplicedperm0']=bool(irmain_permutations.loc[i,'varseq'].translate().find("*")==-1)

irmain_permutations.drop(irmain_permutations[irmain_permutations.stopinsplicedperm0==False].index, inplace=True)

irmain_permutations.to_pickle('./design/irmain_permutations.pkl')
irmain_permutations=pd.read_pickle('./design/irmain_permutations.pkl')


#%% replace with constitutive splice sites


irconst = irmain

canonical_donor = Seq('CAGGTAAGT').tomutable()
no_donor = Seq('CTGCTC').tomutable()
canonical_acceptor = Seq('CTCCTTTCCTTTCAGGTC').tomutable()
no_acceptor = Seq('CAGAGAGGA').tomutable()
GC_donor = Seq('CAGGCAAGT').tomutable()
U12_donor_AT = Seq('ATATCCTTT')
U12_donor_GT = Seq('GTATCCTTT')
U12_branch_acceptor = Seq('TTCCTTAACTTCCTTTCAGATC')
branchpoint = Seq('CTCAC').tomutable()

irconstvar=pd.DataFrame(columns=list(irconst.columns))
j=0
for i in irconst.index:
    intronstart=int(irconst.intronend_varseqnew[i]) - int(irconst.intronlength[i])
    intronend=int(irconst.intronend_varseqnew[i])
    irconstvar=irconstvar.append(irconst.loc[i], ignore_index=True)
    irconstvar.loc[j,'first_ss']=str('endogenous')
    irconstvar.loc[j,'second_ss']=str('endogenous')
    j=j+1
    irconstvar=irconstvar.append(irconst.loc[i], ignore_index=True)
    sequence = irconst.varseq[i].tomutable()
    sequence[intronstart-3:intronstart+6]=canonical_donor
    sequence[intronend-15:intronend+3]=canonical_acceptor
    irconstvar.loc[j,'varseq']=sequence.toseq()
    irconstvar.loc[j,'first_ss']=str('constitutive')
    irconstvar.loc[j,'second_ss']=str('constitutive')
    j=j+1 
    irconstvar=irconstvar.append(irconst.loc[i], ignore_index=True)
    sequence = irconst.varseq[i].tomutable()
    sequence[intronstart-3:intronstart+6]=GC_donor
    sequence[intronend-15:intronend+3]=canonical_acceptor
    irconstvar.loc[j,'varseq']=sequence.toseq()
    irconstvar.loc[j,'first_ss']=str('GC')
    irconstvar.loc[j,'second_ss']=str('constitutive')
    j=j+1 
    irconstvar=irconstvar.append(irconst.loc[i], ignore_index=True)
    sequence = irconst.varseq[i].tomutable()
    sequence[intronstart:intronstart+9]=U12_donor_GT
    sequence[intronend-19:intronend+3]=U12_branch_acceptor
    irconstvar.loc[j,'varseq']=sequence.toseq()
    irconstvar.loc[j,'first_ss']=str('U12_GT')
    irconstvar.loc[j,'second_ss']=str('U12')
    j=j+1 
    irconstvar=irconstvar.append(irconst.loc[i], ignore_index=True)
    sequence = irconst.varseq[i].tomutable()
    sequence[intronstart:intronstart+9]=U12_donor_AT
    sequence[intronend-19:intronend+3]=U12_branch_acceptor
    irconstvar.loc[j,'varseq']=sequence.toseq()
    irconstvar.loc[j,'first_ss']=str('U12_AT')
    irconstvar.loc[j,'second_ss']=str('U12')
    j=j+1 
    irconstvar=irconstvar.append(irconst.loc[i], ignore_index=True)
    sequence = irconst.varseq[i].tomutable()
    sequence[intronstart:intronstart+6]=no_donor
    sequence[intronend-9:intronend]=no_acceptor
    irconstvar.loc[j,'varseq']=sequence.toseq()
    irconstvar.loc[j,'first_ss']=str('nonsplicable')
    irconstvar.loc[j,'second_ss']=str('nonsplicable')
    j=j+1 
    irconstvar=irconstvar.append(irconst.loc[i], ignore_index=True)
    sequence = irconst.varseq[i].tomutable()
    sequence[intronend-26:intronend-21]=branchpoint
    irconstvar.loc[j,'varseq']=sequence.toseq()
    irconstvar.loc[j,'first_ss']=str('endogenous')
    irconstvar.loc[j,'second_ss']=str('constBP')
    j=j+1 
    
irconstvar.to_pickle('./design/irprime_constitutiveandunsplicable_variations.pkl')

#%%

irSF = irmain.loc[[6,59,77,90,104,0,9,46,83,100,124,114]]

irSFvar=pd.DataFrame(columns=list(irSF.columns))

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
RosenbergI3enh	TCTAAC	enh		plus2 stop 
RosenbergI3sil	CCAAGC	sil		


GROUP1:
[-30:-21]
[-21:-12]
[-12:-3]
[5:14]
[14:23]
[-49:-40]
[-40:-31]
[5:14]
[14:23]
'''

### Create sequences containing the motifs and not giving rise to stop codons in frame

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


j=0
for i in irSF.index:
    intstart=int(irSF.intronstart_varseq[i])-int(irSF.intronstart_varseq[i])%3
    intend=int(irSF.intronend_varseqnew[i])-int(irSF.intronstart_varseq[i])%3
    for sf in splicingfactors1.index:
        if (intstart > 30):
            irSFvar=irSFvar.append(irSF.loc[i], ignore_index=True)
            sequence = irSF.varseq[i].tomutable()
            sequence[intstart-30:intstart-21]=splicingfactors1[sf]
            irSFvar.loc[j,'varseq']=sequence.toseq()
            irSFvar.loc[j,'first_ss']=str(sf + ' [-30:-21]')
            irSFvar.loc[j,'second_ss']=str('endogenous')
            j=j+1
        if (intstart > 20):
            irSFvar=irSFvar.append(irSF.loc[i], ignore_index=True)
            sequence = irSF.varseq[i].tomutable()
            sequence[intstart-21:intstart-12]=splicingfactors1[sf]
            irSFvar.loc[j,'varseq']=sequence.toseq()
            irSFvar.loc[j,'first_ss']=str(sf + ' [-21:-12]')
            irSFvar.loc[j,'second_ss']=str('endogenous')
            j=j+1
        irSFvar=irSFvar.append(irSF.loc[i], ignore_index=True)
        sequence = irSF.varseq[i].tomutable()
        sequence[intstart-12:intstart-3]=splicingfactors1[sf]
        irSFvar.loc[j,'varseq']=sequence.toseq()
        irSFvar.loc[j,'first_ss']=str(sf + ' [-12:-3]')
        irSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        irSFvar=irSFvar.append(irSF.loc[i], ignore_index=True)
        sequence = irSF.varseq[i].tomutable()
        sequence[intstart+8:intstart+17]=splicingfactors1[sf]
        irSFvar.loc[j,'varseq']=sequence.toseq()
        irSFvar.loc[j,'first_ss']=str(sf + ' [+6:+15]')
        irSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        irSFvar=irSFvar.append(irSF.loc[i], ignore_index=True)
        sequence = irSF.varseq[i].tomutable()
        sequence[intstart+17:intstart+26]=splicingfactors1[sf]
        irSFvar.loc[j,'varseq']=sequence.toseq()
        irSFvar.loc[j,'first_ss']=str(sf + ' [+15:+24]')
        irSFvar.loc[j,'second_ss']=str('endogenous')
        j=j+1
        irSFvar=irSFvar.append(irSF.loc[i], ignore_index=True)
        sequence = irSF.varseq[i].tomutable()
        sequence[intend-49:intend-40]=splicingfactors1[sf]
        irSFvar.loc[j,'varseq']=sequence.toseq()
        irSFvar.loc[j,'first_ss']=str('endogenous')
        irSFvar.loc[j,'second_ss']=str(sf + ' [-49:-40]')
        j=j+1
        irSFvar=irSFvar.append(irSF.loc[i], ignore_index=True)
        sequence = irSF.varseq[i].tomutable()
        sequence[intend-40:intend-31]=splicingfactors1[sf]
        irSFvar.loc[j,'varseq']=sequence.toseq()
        irSFvar.loc[j,'first_ss']=str('endogenous')
        irSFvar.loc[j,'second_ss']=str(sf + ' [-40:-31]')
        j=j+1

        irSFvar=irSFvar.append(irSF.loc[i], ignore_index=True)
        sequence = irSF.varseq[i].tomutable()
        sequence[intend+9:intend+18]=splicingfactors1[sf]
        irSFvar.loc[j,'varseq']=sequence.toseq()
        irSFvar.loc[j,'first_ss']=str('endogenous')
        irSFvar.loc[j,'second_ss']=str(sf + ' [+9:+18]')
        j=j+1

  


for i in irSFvar.index:
    irSFvar.loc[i,'varseqSF_spliced']= Seq(str(irSFvar.varseq[i][:int(irSFvar.intronstart_varseq[i])]) 
        + str(irSFvar.varseq[i][int(irSFvar.intronend_varseqnew[i]):]))
    if (irSFvar.loc[i,'varseqSF_spliced'].translate().find("*") > -1):
        print(str(i))
        print(irSFvar.loc[i,'varseqSF_spliced'].translate())

irSFvar.drop(list(irSFvar.columns[10:12]),axis=1,inplace=True)

irSFvar.drop(list(irSFvar.columns[11:23]),axis=1,inplace=True)

irSFvar.drop(list(irSFvar.columns[17:21]),axis=1,inplace=True)

irSFvar.drop('splicedvarseq',axis=1,inplace=True)

irSFvar.to_pickle('./design/ir_splicingfactors_location.pkl')


#%% 
irSFcomb = irSF.loc[[59,77,0,83,100,124]]

irSFcombvar=pd.DataFrame(columns=list(irSFcomb.columns))


splicingfactorscomb=pd.Series()
splicingfactorscomb['native']=Seq('NNNNNNNNN').tomutable()
splicingfactorscomb['SRSF1']=Seq('TCACACGAC').tomutable()
splicingfactorscomb['SRSF5']=Seq('TTCACAGGC').tomutable()
splicingfactorscomb['hnRNPA1']=Seq('TTAGGGAAC').tomutable()
splicingfactorscomb['hnRNPU']=Seq('TTGTATTGC').tomutable()

k=0


for i in irSFcomb.index:
    startpositions=[int(irSF.intronstart_varseq[i])-int(irSF.intronstart_varseq[i])%3 -12,\
        int(irSF.intronstart_varseq[i])-int(irSF.intronstart_varseq[i])%3 + 8,\
        int(irSF.intronend_varseqnew[i])-int(irSF.intronstart_varseq[i])%3 - 40,\
        int(irSF.intronend_varseqnew[i])-int(irSF.intronstart_varseq[i])%3 + 9]
    for sf in splicingfactorscomb.index: 
        sequence = irSFcomb.varseq[i].tomutable()
        if (sf!="native"):
            first=sf
            firstpos=startpositions[0]
            sequence[startpositions[0]:startpositions[0]+9] = splicingfactorscomb.loc[sf]
            for sf in splicingfactorscomb.index:
                if (sf!="native"):
                    sequence[startpositions[1]:startpositions[1]+9] = splicingfactorscomb.loc[sf]
                    second=sf
                    secondpos=startpositions[1]
                    irSFcombvar=irSFcombvar.append(irSFcomb.loc[i], ignore_index=True)        
                    irSFcombvar.loc[k,'varseq']=sequence.toseq()
                    irSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                    irSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                    k=k+1
                    
                else:
                    for sf in splicingfactorscomb.index:
                        if (sf!="native"):
                            sequence[startpositions[2]:startpositions[2]+9] = splicingfactorscomb.loc[sf]
                            second=sf
                            secondpos=startpositions[2]
                            irSFcombvar=irSFcombvar.append(irSFcomb.loc[i], ignore_index=True)        
                            irSFcombvar.loc[k,'varseq']=sequence.toseq()
                            irSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                            irSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                            k=k+1
                    else:
                        for sf in splicingfactorscomb.index:
                            if (sf!="native"):
                                sequence[startpositions[3]:startpositions[3]+9] = splicingfactorscomb.loc[sf]
                                second=sf
                                secondpos=startpositions[3]
                                irSFcombvar=irSFcombvar.append(irSFcomb.loc[i], ignore_index=True)        
                                irSFcombvar.loc[k,'varseq']=sequence.toseq()
                                irSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                                irSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
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
                            irSFcombvar=irSFcombvar.append(irSFcomb.loc[i], ignore_index=True)        
                            irSFcombvar.loc[k,'varseq']=sequence.toseq()
                            irSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                            irSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                            k=k+1
                        else:
                            for sf in splicingfactorscomb.index:
                                if (sf!="native"):
                                    sequence[startpositions[3]:startpositions[3]+9] = splicingfactorscomb.loc[sf]
                                    second=sf
                                    secondpos=startpositions[3]
                                    irSFcombvar=irSFcombvar.append(irSFcomb.loc[i], ignore_index=True)        
                                    irSFcombvar.loc[k,'varseq']=sequence.toseq()
                                    irSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                                    irSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
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
                                    irSFcombvar=irSFcombvar.append(irSFcomb.loc[i], ignore_index=True)        
                                    irSFcombvar.loc[k,'varseq']=sequence.toseq()
                                    irSFcombvar.loc[k,'first_SF']=str(str(first) + ' at position ' + str(firstpos))
                                    irSFcombvar.loc[k,'second_SF']=str(str(second) + ' at position ' + str(secondpos))
                                    k=k+1
            

for i in irSFcombvar.index:
    irSFcombvar.loc[i,'varseqSF_spliced']= Seq(str(irSFcombvar.varseq[i][:int(irSFcombvar.intronstart_varseq[i])]) 
        + str(irSFcombvar.varseq[i][int(irSFcombvar.intronend_varseqnew[i]):]))
    if (irSFcombvar.loc[i,'varseqSF_spliced'].translate().find("*") > -1):
        print(str(i))
        print(irSFcombvar.loc[i,'varseqSF_spliced'].translate())

irSFcombvar.drop(list(irSFcombvar.columns[10:12]),axis=1,inplace=True)

irSFcombvar.drop(list(irSFcombvar.columns[11:23]),axis=1,inplace=True)

irSFcombvar.drop(list(irSFcombvar.columns[17:21]),axis=1,inplace=True)

irSFcombvar.drop('splicedvarseq',axis=1,inplace=True)

            
irSFcombvar.to_pickle('./design/ir_splicingfactors_combinatorial.pkl')


 
#%% Test Rosenberg features (From Rosenberg et al., 2015, Cell)

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

irSFRosenberg=pd.DataFrame(columns=list(irSFcomb.columns))

j=0
for i in irSFcomb.index:
    intstart=int(irSFcomb.intronstart_varseq[i])-int(irSFcomb.intronstart_varseq[i])%3
    intend=int(irSFcomb.intronend_varseqnew[i])-int(irSFcomb.intronstart_varseq[i])%3
    for sf in splicingfactors2.index:
        if (intstart > 30):
            irSFRosenberg=irSFRosenberg.append(irSFcomb.loc[i], ignore_index=True)
            sequence = irSFcomb.varseq[i].tomutable()
            sequence[intstart-30:intstart-24]=splicingfactors2[sf]
            irSFRosenberg.loc[j,'varseq']=sequence.toseq()
            irSFRosenberg.loc[j,'first_ss']=str(sf + ' [-27:-21]')
            irSFRosenberg.loc[j,'second_ss']=str('endogenous')
            j=j+1
        if (intstart > 20):
            irSFRosenberg=irSFRosenberg.append(irSFcomb.loc[i], ignore_index=True)
            sequence = irSFcomb.varseq[i].tomutable()
            sequence[intstart-21:intstart-15]=splicingfactors2[sf]
            irSFRosenberg.loc[j,'varseq']=sequence.toseq()
            irSFRosenberg.loc[j,'first_ss']=str(sf + ' [-18:-12]')
            irSFRosenberg.loc[j,'second_ss']=str('endogenous')
            j=j+1
        irSFRosenberg=irSFRosenberg.append(irSFcomb.loc[i], ignore_index=True)
        sequence = irSFcomb.varseq[i].tomutable()
        sequence[intstart-12:intstart-6]=splicingfactors2[sf]
        irSFRosenberg.loc[j,'varseq']=sequence.toseq()
        irSFRosenberg.loc[j,'first_ss']=str(sf + ' [-9:-3]')
        irSFRosenberg.loc[j,'second_ss']=str('endogenous')
        j=j+1
        irSFRosenberg=irSFRosenberg.append(irSFcomb.loc[i], ignore_index=True)
        sequence = irSFcomb.varseq[i].tomutable()
        sequence[intstart+11:intstart+17]=splicingfactors2[sf]
        irSFRosenberg.loc[j,'varseq']=sequence.toseq()
        irSFRosenberg.loc[j,'first_ss']=str(sf + ' [+9:+15]')
        irSFRosenberg.loc[j,'second_ss']=str('endogenous')
        j=j+1
        irSFRosenberg=irSFRosenberg.append(irSFcomb.loc[i], ignore_index=True)
        sequence = irSFcomb.varseq[i].tomutable()
        sequence[intstart+20:intstart+26]=splicingfactors2[sf]
        irSFRosenberg.loc[j,'varseq']=sequence.toseq()
        irSFRosenberg.loc[j,'first_ss']=str(sf + ' [+18:+24]')
        irSFRosenberg.loc[j,'second_ss']=str('endogenous')
        j=j+1
        irSFRosenberg=irSFRosenberg.append(irSFcomb.loc[i], ignore_index=True)
        sequence = irSFcomb.varseq[i].tomutable()
        sequence[intend-49:intend-43]=splicingfactors2[sf]
        irSFRosenberg.loc[j,'varseq']=sequence.toseq()
        irSFRosenberg.loc[j,'first_ss']=str('endogenous')
        irSFRosenberg.loc[j,'second_ss']=str(sf + ' [-52:-40]')
        j=j+1
        irSFRosenberg=irSFRosenberg.append(irSFcomb.loc[i], ignore_index=True)
        sequence = irSFcomb.varseq[i].tomutable()
        sequence[intend-40:intend-34]=splicingfactors2[sf]
        irSFRosenberg.loc[j,'varseq']=sequence.toseq()
        irSFRosenberg.loc[j,'first_ss']=str('endogenous')
        irSFRosenberg.loc[j,'second_ss']=str(sf + ' [-37:-31]')
        j=j+1

        irSFRosenberg=irSFRosenberg.append(irSFcomb.loc[i], ignore_index=True)
        sequence = irSFcomb.varseq[i].tomutable()
        sequence[intend+12:intend+18]=splicingfactors2[sf]
        irSFRosenberg.loc[j,'varseq']=sequence.toseq()
        irSFRosenberg.loc[j,'first_ss']=str('endogenous')
        irSFRosenberg.loc[j,'second_ss']=str(sf + ' [+12:+18]')
        j=j+1

  


for i in irSFRosenberg.index:
    irSFRosenberg.loc[i,'varseqSF_spliced']= Seq(str(irSFRosenberg.varseq[i][:int(irSFRosenberg.intronstart_varseq[i])]) 
        + str(irSFRosenberg.varseq[i][int(irSFRosenberg.intronend_varseqnew[i]):]))
    if (irSFRosenberg.loc[i,'varseqSF_spliced'].translate().find("*") > -1):
        print(str(i))
        print(irSFRosenberg.loc[i,'varseqSF_spliced'].translate())

irSFRosenberg.drop(list(irSFRosenberg.columns[10:12]),axis=1,inplace=True)

irSFRosenberg.drop(list(irSFRosenberg.columns[11:23]),axis=1,inplace=True)

irSFRosenberg.drop(list(irSFRosenberg.columns[17:21]),axis=1,inplace=True)

irSFRosenberg.drop('splicedvarseq',axis=1,inplace=True)

       
irSFRosenberg.to_pickle('./design/ir_splicingfactors_Rosenberg.pkl')

#%% secondary structure

irsec=irmain

'''
start-12:start-3
start+15:start+24
end-37:end-28
end+3:end+12
'''

irsecvar=pd.DataFrame(columns=list(irsec.columns))
k=0
for i in irsec.index:
    start=int(irsec.intronstart_varseq[i])
    end=int(irsec.intronend_varseqnew[i])
    startinframe=int(irsec.intronstart_varseq[i]) - int(irsec.intronstart_varseq[i])%3
    if (int(irsec.intronstart_varseq[i])%3==0):
        endinframe=int(irsec.intronend_varseqnew[i])
    else:
        endinframe=int(irsec.intronend_varseqnew[i])+3 - int(irsec.intronstart_varseq[i])%3
    
    exonrev=irsec.varseq[i][startinframe-24:startinframe-15].reverse_complement() 
    exoncomp=irsec.varseq[i][startinframe-24:startinframe-15].complement()
    firstssrev = irsec.varseq[i][start:start+9].reverse_complement()
    firstsscomp = irsec.varseq[i][start:start+9].complement()
    firstssrev2 = irsec.varseq[i][start+3:start+12].reverse_complement()
    firstsscomp2 = irsec.varseq[i][start+3:start+12].complement()
    
    secondssrev = irsec.varseq[i][end-9:end].reverse_complement()
    secondsscomp = irsec.varseq[i][end-9:end].complement()
    secondssrev2 = irsec.varseq[i][end-12:end-3].reverse_complement()
    secondsscomp2 = irsec.varseq[i][end-12:end-3].complement()
    altexonrev=irsec.varseq[i][endinframe+15:endinframe+24].reverse_complement() #if end<139
    altexoncomp=irsec.varseq[i][endinframe+15:endinframe+24].complement()
    

    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[startinframe-12:startinframe-3]=firstssrev
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='rev at [start-12:start-3]'
    irsecvar.loc[k,'secondss']='endogenous'
    k=k+1
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[startinframe-12:startinframe-3]=firstsscomp
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='comp at [start-12:start-3]'
    irsecvar.loc[k,'secondss']='endogenous'
    k=k+1
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[startinframe-12:startinframe-3]=firstssrev2
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='rev2 at [start-12:start-3]'
    irsecvar.loc[k,'secondss']='endogenous'
    k=k+1
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[startinframe-12:startinframe-3]=firstsscomp2
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='comp2 at [start-12:start-3]'
    irsecvar.loc[k,'secondss']='endogenous'
    k=k+1
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[startinframe-12:startinframe-3]=exonrev
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='exonrev at [start-12:start-3]'
    irsecvar.loc[k,'secondss']='endogenous'
    k=k+1
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[startinframe-12:startinframe-3]=exoncomp
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='exoncomp at [start-12:start-3]'
    irsecvar.loc[k,'secondss']='endogenous'
    k=k+1
    

    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[start+15:start+24]=firstssrev
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='rev at [start+15:start+24]'
    irsecvar.loc[k,'secondss']='endogenous'
    k=k+1
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[start+15:start+24]=firstsscomp
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='comp at [start+15:start+24]'
    irsecvar.loc[k,'secondss']='endogenous'
    k=k+1
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[start+15:start+24]=firstssrev2
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='rev2 at [start+15:start+24]'
    irsecvar.loc[k,'secondss']='endogenous'
    k=k+1
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[start+15:start+24]=firstsscomp2
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='comp2 at [start+15:start+24]'
    irsecvar.loc[k,'secondss']='endogenous'
    k=k+1

    
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[end-37:end-28]=secondssrev
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='endogenous'
    irsecvar.loc[k,'secondss']='rev at [end-37:end-28]'
    k=k+1
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[end-37:end-28]=secondsscomp
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='endogenous'
    irsecvar.loc[k,'secondss']='comp at [end-37:end-28]'
    k=k+1
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[end-37:end-28]=secondssrev2
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='endogenous'
    irsecvar.loc[k,'secondss']='rev2 at [end-37:end-28]'
    k=k+1
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[end-37:end-28]=secondsscomp2
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='endogenous'
    irsecvar.loc[k,'secondss']='comp2 at [end-37:end-28]'
    k=k+1
    
    
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[endinframe+3:endinframe+12]=secondssrev
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='endogenous'
    irsecvar.loc[k,'secondss']='rev at [end+3:end+12]'
    k=k+1
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[endinframe+3:endinframe+12]=secondsscomp
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='endogenous'
    irsecvar.loc[k,'secondss']='comp at [end+3:end+12]'
    k=k+1
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[endinframe+3:endinframe+12]=secondssrev2
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='endogenous'
    irsecvar.loc[k,'secondss']='rev2 at [end+3:end+12]'
    k=k+1
    irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
    sequence=irsec.varseq[i].tomutable()
    sequence[endinframe+3:endinframe+12]=secondsscomp2
    irsecvar.varseq[k]=sequence.toseq()
    irsecvar.loc[k,'firstss']='endogenous'
    irsecvar.loc[k,'secondss']='comp2 at [end+3:end+12]'
    k=k+1
    if (endinframe<139):
        irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
        sequence=irsec.varseq[i].tomutable()
        sequence[endinframe+3:endinframe+12]=exonrev
        irsecvar.varseq[k]=sequence.toseq()
        irsecvar.loc[k,'firstss']='endogenous'
        irsecvar.loc[k,'secondss']='exonrev at [end+3:end+12]'
        k=k+1
        irsecvar=irsecvar.append(irsec.loc[i],ignore_index=True)    
        sequence=irsec.varseq[i].tomutable()
        sequence[endinframe+3:endinframe+12]=exoncomp
        irsecvar.varseq[k]=sequence.toseq()
        irsecvar.loc[k,'firstss']='endogenous'
        irsecvar.loc[k,'secondss']='exoncomp at [end+3:end+12]'
        k=k+1

    

for i in irsecvar.index:
    start=int(irsecvar.intronstart_varseq[i])
    end=int(irsecvar.intronend_varseqnew[i])
    spliced=irsecvar.varseq[i][:start] + irsecvar.varseq[i][end:]
    if (spliced.translate().find("*")>-1):
        print(str(i))
        print(spliced.translate())
        irsecvar.drop(i,inplace=True)
        
        
       
irsecvar.to_pickle('./design/ir_secondarystructure_variants.pkl')

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

irnuc=pd.DataFrame(columns=list(irmain))

k=0
for i in irmain.index:
    start=int(irmain.intronstart_varseq[i])
    end=int(irmain.intronend_varseqnew[i])
    startinframe=int(irmain.intronstart_varseq[i]) - int(irmain.intronstart_varseq[i])%3
    endinframe=int(irmain.intronend_varseqnew[i])+3 - int(irmain.intronstart_varseq[i])%3
    aaseq1=irmain.varseq[i][:startinframe].translate()
    aaseq2=irmain.varseq[i][endinframe:].translate()
    aaintron=irmain.varseq[i][startinframe+6:endinframe-32].translate()
    exon1max=str()
    for j in aaseq1[:-1]:
        exon1max=exon1max + codonmaxgc[j][0]

    exon1min=str()
    for j in aaseq1[:-1]:
        exon1min=exon1min + codonmingc[j][0]
    exon1rand1=str()
    for j in aaseq1[:-1]:
        exon1rand1=exon1rand1 + random.choice(gencode[j])
    exon1rand2=str()
    for j in aaseq1[:-1]:
        exon1rand2=exon1rand2 + random.choice(gencode[j])
    exon1rand3=str()
    for j in aaseq1[:-1]:
        exon1rand3=exon1rand3 + random.choice(gencode[j])
    exon2max=str()
    for j in aaseq2[1:]:
        exon2max=exon2max + codonmaxgc[j][0]
    exon2min=str()
    for j in aaseq2[1:]:
        exon2min=exon2min + codonmingc[j][0]
    exon2rand1=str()
    for j in aaseq2[1:]:
        exon2rand1=exon2rand1 + random.choice(gencode[j])
    exon2rand2=str()
    for j in aaseq2[1:]:
        exon2rand2=exon2rand2 + random.choice(gencode[j])
    exon2rand3=str()
    for j in aaseq2[1:]:
        exon2rand3=exon2rand3 + random.choice(gencode[j])
    intronmax=str()
    for j in aaintron:
        intronmax=intronmax + codonmaxgc[j][0]
    intronmin=str()
    for j in aaintron:
        intronmin=intronmin + codonmingc[j][0]
    intronrand1=str()
    for j in aaintron:
        intronrand1=intronrand1 + random.choice(gencode[j])
    intronrand2=str()
    for j in aaintron:
        intronrand2=intronrand2 + random.choice(gencode[j])
    intronrand3=str()
    for j in aaintron:
        intronrand3=intronrand3 + random.choice(gencode[j])
    
    

    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    irnuc.loc[k,'changes']='endogenous'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1max
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2max
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded maximal gc'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1min
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2min
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded minimal gc'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1rand1
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2rand1
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded rand1'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1rand2
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2rand2
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded rand2'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1rand3
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2rand3
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded rand3'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronmax
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded maximal gc'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronmin
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded minimal gc'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronrand1
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded rand1'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronrand2
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded rand2'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronrand3
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded rand3'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1rand1
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2rand1
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronrand1
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon and intron recoded rand1'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1rand2
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2rand2
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronrand2
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon and intron recoded rand2'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1rand3
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2rand3
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronrand3
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon and intron recoded rand3'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1max
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2max
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronmin
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded maximal gc and intron recoded minimal gc'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1min
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2min
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronmax
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded minimal gc and intron recoded maximal gc'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1max
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2max
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronrand1
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded maximal gc and intron recoded rand1'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1min
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2min
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronrand1
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded minimal gc and intron recoded rand1'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1max
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2max
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronrand2
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded maximal gc and intron recoded rand2'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1min
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2min
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronrand2
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded minimal gc and intron recoded rand2'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1max
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2max
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronrand3
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded maximal gc and intron recoded rand3'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1min
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2min
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronrand3
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded minimal gc and intron recoded rand3'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1rand1
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2rand1
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronmax
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded maximal gc and exon recoded rand1'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1rand1
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2rand1
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronmin
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded minimal gc and exon recoded rand1'
    k=k+1
 
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1rand2
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2rand2
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronmax
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded maximal gc and exon recoded rand2'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1rand2
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2rand2
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronmin
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded minimal gc and exon recoded rand2'
    k=k+1
 
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1rand3
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2rand3
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronmax
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded maximal gc and exon recoded rand3'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1rand3
    sequence[endinframe+3:endinframe+3+len(exon2max)]=exon2rand3
    sequence[startinframe+6:startinframe+6+len(intronmax)]=intronmin
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded minimal gc and exon recoded rand3'
    k=k+1
 

imp22dT=Seq('TTTTTTTTTCATTTTTTTTTTT')
perf12dT=Seq('TTTTTTTTTTTT')
perf5dT=Seq('TTTTTC')
imp22ctrl=Seq('TTGGACTGTCATTCACGTCGTT')
perf12ctrl=Seq('TTGGACTGTCAT')
perf5ctrl=Seq('TTGGAC')


for i in irmain.index:
    start=int(irmain.intronstart_varseq[i])
    if (start>27):
        end=int(irmain.intronend_varseqnew[i])
        startinframe=int(irmain.intronstart_varseq[i]) - int(irmain.intronstart_varseq[i])%3
        endinframe=int(irmain.intronend_varseqnew[i])+3 - int(irmain.intronstart_varseq[i])%3
        
        startpositionsexon=[startinframe-21,startinframe-12,endinframe+3,endinframe+12]
        startpositionsintron=[startinframe+9,startinframe+18,endinframe-48,endinframe-39]
    
        irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
        sequence=irmain.varseq[i].tomutable()
        for s in startpositionsexon:
            sequence[s:s+6]=perf5dT
        irnuc.varseq[k]=sequence.toseq()
        irnuc.loc[k,'changes']='perf5dT in exon'
        k=k+1
        
        irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
        sequence=irmain.varseq[i].tomutable()
        for s in startpositionsexon:
            sequence[s:s+6]=perf5ctrl
        irnuc.varseq[k]=sequence.toseq()
        irnuc.loc[k,'changes']='perf5ctrl in exon'
        k=k+1
        
        irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
        sequence=irmain.varseq[i].tomutable()
        for s in startpositionsintron:
            sequence[s:s+6]=perf5dT
        irnuc.varseq[k]=sequence.toseq()
        irnuc.loc[k,'changes']='perf5dT in intron'
        k=k+1
        
        irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
        sequence=irmain.varseq[i].tomutable()
        for s in startpositionsintron:
            sequence[s:s+6]=perf5ctrl
        irnuc.varseq[k]=sequence.toseq()
        irnuc.loc[k,'changes']='perf5ctrl in intron'
        k=k+1


# methylation

for i in irmain.index:
    start=int(irmain.intronstart_varseq[i])
    end=int(irmain.intronend_varseqnew[i])
    startinframe=int(irmain.intronstart_varseq[i]) - int(irmain.intronstart_varseq[i])%3
    endinframe=int(irmain.intronend_varseqnew[i])+3 - int(irmain.intronstart_varseq[i])%3
    aaseq1=irmain.varseq[i][:startinframe].translate()
    aaseq2=irmain.varseq[i][endinframe:].translate()
    aaintron=irmain.varseq[i][startinframe+6:endinframe-32].translate()
    exon1min=str()
    for j in aaseq1[:-1]:
        exon1min=exon1min + codonmingc[j][0]
    exon1rand1=str()
    for j in aaseq1[:-1]:
        exon1rand1=exon1rand1 + random.choice(gencode[j])
    exon2min=str()
    for j in aaseq2[1:]:
        exon2min=exon2min + codonmingc[j][0]
    exon2rand1=str()
    for j in aaseq2[1:]:
        exon2rand1=exon2rand1 + random.choice(gencode[j])
    intronmin=str()
    for j in aaintron:
        intronmin=intronmin + codonmingc[j][0]
    intronrand1=str()
    for j in aaintron:
        intronrand1=intronrand1 + random.choice(gencode[j])
    
    
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    count= int(startinframe//9)
    for j in range(0,count,1):
        sequence[1 + 9*j:1 + 9*j+2]=Seq('CG')
    count= int((162-endinframe)//9)
    for j in range(0,count,1):
        sequence[endinframe + 1 + 9*j:endinframe + 1 + 9*j+2]=Seq('CG')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='CG in exon every 9nt'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1min
    sequence[endinframe+3:endinframe+3+len(exon2min)]=exon2min
    count= int(startinframe//9)
    for j in range(0,count,1):
        sequence[1 + 9*j:1 + 9*j+2]=Seq('CG')
    count= int((162-endinframe)//9)
    for j in range(0,count,1):
        sequence[endinframe + 1 + 9*j:endinframe + 1 + 9*j+2]=Seq('CG')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded minimal gc, CG in exon every 9nt'
    k=k+1

    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    count= int(startinframe//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('CG')
    count= int((162-endinframe)//6)
    for j in range(0,count,1):
        sequence[endinframe + 1 + 6*j:endinframe + 1 + 6*j+2]=Seq('CG')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='CG in exon every 6nt'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1min
    sequence[endinframe+3:endinframe+3+len(exon2min)]=exon2min
    count= int(startinframe//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('CG')
    count= int((162-endinframe)//6)
    for j in range(0,count,1):
        sequence[endinframe + 1 + 6*j:endinframe + 1 + 6*j+2]=Seq('CG')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded minimal gc, CG in exon every 6nt'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    count= int((endinframe-33-(startinframe+6))//9)
    for j in range(0,count,1):
        sequence[startinframe+6 +1 + 9*j:startinframe+6+1 + 9*j+2]=Seq('CG')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='CG in intron every 9nt'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[startinframe+6:startinframe+6+len(intronmin)]=intronmin
    count= int((endinframe-33-(startinframe+6))//9)
    for j in range(0,count,1):
        sequence[startinframe+6 +1 + 9*j:startinframe+6+1 + 9*j+2]=Seq('CG')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded minimal gc, CG in intron every 9nt'
    k=k+1

    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    count= int((endinframe-33-(startinframe+6))//6)
    for j in range(0,count,1):
        sequence[startinframe+6 +1 + 6*j:startinframe+6+1 + 6*j+2]=Seq('CG')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='CG in intron every 6nt'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[startinframe+6:startinframe+6+len(intronmin)]=intronmin
    count= int((endinframe-33-(startinframe+6))//6)
    for j in range(0,count,1):
        sequence[startinframe+6 +1 + 6*j:startinframe+6+1 + 6*j+2]=Seq('CG')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded minimal gc, CG in intron every 6nt'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    count= int(startinframe//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('CG')
    count= int((162-endinframe)//6)
    for j in range(0,count,1):
        sequence[endinframe + 1 + 6*j:endinframe + 1 + 6*j+2]=Seq('CG')
    count= int((endinframe-33-(startinframe+6))//6)
    for j in range(0,count,1):
        sequence[startinframe+6 +1 + 6*j:startinframe+6+1 + 6*j+2]=Seq('CG')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='CG in exon and intron every 6nt'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1min
    sequence[endinframe+3:endinframe+3+len(exon2min)]=exon2min
    sequence[startinframe+6:startinframe+6+len(intronmin)]=intronmin
    count= int(startinframe//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('CG')
    count= int((162-endinframe)//6)
    for j in range(0,count,1):
        sequence[endinframe + 1 + 6*j:endinframe + 1 + 6*j+2]=Seq('CG')
    count= int((endinframe-33-(startinframe+6))//6)
    for j in range(0,count,1):
        sequence[startinframe+6 +1 + 6*j:startinframe+6+1 + 6*j+2]=Seq('CG')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon and intron recoded minimal gc, CG in exon and intron every 6nt'
    k=k+1
    
   
          
# same for GC

    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    count= int(startinframe//9)
    for j in range(0,count,1):
        sequence[1 + 9*j:1 + 9*j+2]=Seq('GC')
    count= int((162-endinframe)//9)
    for j in range(0,count,1):
        sequence[endinframe + 1 + 9*j:endinframe + 1 + 9*j+2]=Seq('GC')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='GC in exon every 9nt'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1min
    sequence[endinframe+3:endinframe+3+len(exon2min)]=exon2min
    count= int(startinframe//9)
    for j in range(0,count,1):
        sequence[1 + 9*j:1 + 9*j+2]=Seq('GC')
    count= int((162-endinframe)//9)
    for j in range(0,count,1):
        sequence[endinframe + 1 + 9*j:endinframe + 1 + 9*j+2]=Seq('GC')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded minimal gc, GC in exon every 9nt'
    k=k+1

    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    count= int(startinframe//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('GC')
    count= int((162-endinframe)//6)
    for j in range(0,count,1):
        sequence[endinframe + 1 + 6*j:endinframe + 1 + 6*j+2]=Seq('GC')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='GC in exon every 6nt'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1min
    sequence[endinframe+3:endinframe+3+len(exon2min)]=exon2min
    count= int(startinframe//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('GC')
    count= int((162-endinframe)//6)
    for j in range(0,count,1):
        sequence[endinframe + 1 + 6*j:endinframe + 1 + 6*j+2]=Seq('GC')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon recoded minimal gc, GC in exon every 6nt'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    count= int((endinframe-33-(startinframe+6))//9)
    for j in range(0,count,1):
        sequence[startinframe+6 +1 + 9*j:startinframe+6+1 + 9*j+2]=Seq('GC')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='GC in intron every 9nt'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[startinframe+6:startinframe+6+len(intronmin)]=intronmin
    count= int((endinframe-33-(startinframe+6))//9)
    for j in range(0,count,1):
        sequence[startinframe+6 +1 + 9*j:startinframe+6+1 + 9*j+2]=Seq('GC')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded minimal gc, GC in intron every 9nt'
    k=k+1

    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    count= int((endinframe-33-(startinframe+6))//6)
    for j in range(0,count,1):
        sequence[startinframe+6 +1 + 6*j:startinframe+6+1 + 6*j+2]=Seq('GC')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='GC in intron every 6nt'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[startinframe+6:startinframe+6+len(intronmin)]=intronmin
    count= int((endinframe-33-(startinframe+6))//6)
    for j in range(0,count,1):
        sequence[startinframe+6 +1 + 6*j:startinframe+6+1 + 6*j+2]=Seq('GC')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='intron recoded minimal gc, GC in intron every 6nt'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    count= int(startinframe//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('GC')
    count= int((162-endinframe)//6)
    for j in range(0,count,1):
        sequence[endinframe + 1 + 6*j:endinframe + 1 + 6*j+2]=Seq('GC')
    count= int((endinframe-33-(startinframe+6))//6)
    for j in range(0,count,1):
        sequence[startinframe+6 +1 + 6*j:startinframe+6+1 + 6*j+2]=Seq('GC')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='GC in exon and intron every 6nt'
    k=k+1
    
    irnuc=irnuc.append(irmain.loc[i],ignore_index=True)
    sequence=irmain.varseq[i].tomutable()
    sequence[:startinframe-3]=exon1min
    sequence[endinframe+3:endinframe+3+len(exon2min)]=exon2min
    sequence[startinframe+6:startinframe+6+len(intronmin)]=intronmin
    count= int(startinframe//6)
    for j in range(0,count,1):
        sequence[1 + 6*j:1 + 6*j+2]=Seq('GC')
    count= int((162-endinframe)//6)
    for j in range(0,count,1):
        sequence[endinframe + 1 + 6*j:endinframe + 1 + 6*j+2]=Seq('GC')
    count= int((endinframe-33-(startinframe+6))//6)
    for j in range(0,count,1):
        sequence[startinframe+6 +1 + 6*j:startinframe+6+1 + 6*j+2]=Seq('GC')
    irnuc.varseq[k]=sequence.toseq()
    irnuc.loc[k,'changes']='exon and intron recoded minimal gc, GC in exon and intron every 6nt'
    k=k+1
    
  
    
for i in irnuc.index:
    start=int(irnuc.intronstart_varseq[i])
    end=int(irnuc.intronend_varseqnew[i])
    spliced=irnuc.varseq[i][:start] + irnuc.varseq[i][end:]
    if (spliced.translate().find("*")>-1):
        print(str(i))
        print(spliced.translate())
#        irnuc.drop(i,inplace=True)
 
for i in irnuc.index:
    if (len(irnuc.varseq[i])!=162):
        print(str(i))
        print(len(irnuc.varseq[i]))
#        irnuc.drop(i,inplace=True)

       
irnuc.drop_duplicates('varseq',inplace=True)        
       
irnuc.to_pickle('./design/ir__nucleosome_recoding_methylation.pkl')
       
       
#%% constitutive introns
       
       
junctions=pd.read_pickle('./design/junctions.pkl')

constintrons=pd.DataFrame()

for i in range(len(junctions)):
    if (int(junctions.loc[junctions.index[i],'end'])-int(junctions.loc[junctions.index[i],'start'])<120):
        constintrons= constintrons.append(junctions.iloc[i,:])

        constintrons.loc[junctions.index[i],'intronlength']= \
            int(constintrons.loc[junctions.index[i],'end'])- \
            int(constintrons.loc[junctions.index[i],'start'])

constintronssmall=constintrons[(constintrons.intronlength>70)&(constintrons.level>100)&(constintrons.intronlength%3==1)]
   
   


f=open('./design/constintron.bed', 'w')
for i in constintronssmall.index:
    if (constintronssmall.loc[i,'strand']=='+'):
        f.write(constintronssmall.loc[i,'chr'] + "\t" + str(int(constintronssmall.loc[i,'end'])-150) + \
            "\t" + str(int(constintronssmall.loc[i,'end'])+50) + "\t" + constintronssmall.loc[i,'name'] + \
            "\t" + str(constintronssmall.loc[i,'intronlength']) + "\t" + constintronssmall.loc[i,'strand'] + "\n")
    else:
        f.write(constintronssmall.loc[i,'chr'] + "\t" + str(int(constintronssmall.loc[i,'start'])-50) + \
            "\t" + str(int(constintronssmall.loc[i,'start'])+150) + "\t" + constintronssmall.loc[i,'name'] + \
            "\t" + str(constintronssmall.loc[i,'intronlength']) + "\t" + constintronssmall.loc[i,'strand'] + "\n")

f.close()


'''
bedtools getfasta -fi hg19.fa -bed constintron.bed -fo constintronseqaroundsplicesite.fa -s

'''     

#parse fasta file

from Bio import SeqIO

f=open('./design/constintronseqaroundsplicesite.fa','r')
records=list(SeqIO.parse(f, 'fasta'))

c=0
for i in constintronssmall.index:
    constintronssmall.loc[i,'coordinates']=records[c].id
    constintronssmall.loc[i,'sequence']=records[c].seq
    c=c+1
    
    
ir82comb=forlibrary.make_combinatorialintron_alt(ir82,ir82)
ir85comb=forlibrary.make_combinatorialintron_alt(ir85,ir85)
ir97comb=forlibrary.make_combinatorialintron_alt(ir97,ir97)

intronretention_combinatorial_variations_threeway=pd.concat([ir82comb,ir85comb,ir97comb],ignore_index=True)

for i in intronretention_combinatorial_variations_threeway.index:
    start=int(intronretention_combinatorial_variations_threeway.intronstart_varseq[i])
    end=int(intronretention_combinatorial_variations_threeway.intronend_varseqnew[i])
    spliced=intronretention_combinatorial_variations_threeway.varseq[i][:start] + intronretention_combinatorial_variations_threeway.varseq[i][end:]
    if (spliced.translate().find("*")>-1):
        print(str(i))
        print(spliced.translate())
#        intronretention_combinatorial_variations_threeway.drop(i,inplace=True)

### PICKLE
intronretention_combinatorial_variations_threeway.to_pickle('./library/intronretention_combinatorial_variations_threeway.pkl')

    
introncomb=forlibrary.make_combinatorialintron(irmain,constintronssmall)

for i in introncomb.index:
    start=int(introncomb.intronstart_varseq[i])
    end=int(introncomb.intronend_varseqnew[i])
    spliced=introncomb.varseq[i][:start] + introncomb.varseq[i][end:]
    if (spliced.translate().find("*")>-1):
        print(str(i))
        print(spliced.translate())
        introncomb.drop(i,inplace=True)

introncomb.to_pickle('./design/intronretention_combinatorial_variations_threeway_withconstitutive.pkl')

constidxs=[]
for i in introncomb.combination.dropna().drop_duplicates():
    constidxs+=i.split('\t')
    
constintronssmall[[str(x) in pd.Series(constidxs).unique()[1:] for x in constintronssmall.index]][['chr','start','end','name','strand','intronlength']].to_csv('./tables/TableS9.csv')


#%% test different splice site sequences

irswitch=irmain[irmain.intronend_varseqnew==137]

irswitchvar=pd.DataFrame()

k=0

for i in irswitch.index:
    starti=int(irswitch.intronstart_varseq[i])
    
    for j in irswitch.index:
        startj=int(irswitch.intronstart_varseq[j])

        irswitchvar=irswitchvar.append(irswitch.loc[i],ignore_index=True)
        sequence=irswitch.varseq[i].tomutable()
        sequence[137-15:140]=irswitch.varseq[j][137-15:140]
        irswitchvar.varseq[k]=sequence.toseq()
        irswitchvar.loc[k,'changes']='acceptor from ' + str(j)
        k=k+1
        
        irswitchvar=irswitchvar.append(irswitch.loc[i],ignore_index=True)
        sequence=irswitch.varseq[i].tomutable()
        sequence[starti-3:starti+6]=irswitch.varseq[j][startj-3:startj+6]
        irswitchvar.varseq[k]=sequence.toseq()
        irswitchvar.loc[k,'changes']='donor from ' + str(j)
        k=k+1
        
        irswitchvar=irswitchvar.append(irswitch.loc[i],ignore_index=True)
        sequence=irswitch.varseq[i].tomutable()
        sequence[137-15:140]=irswitch.varseq[j][137-15:140]
        sequence[starti-3:starti+6]=irswitch.varseq[j][startj-3:startj+6]
        irswitchvar.varseq[k]=sequence.toseq()
        irswitchvar.loc[k,'changes']='acceptor and donor from ' + str(j)
        k=k+1
        
for i in irswitchvar.index:
    start=int(irswitchvar.intronstart_varseq[i])
    end=int(irswitchvar.intronend_varseqnew[i])
    spliced=irswitchvar.varseq[i][:start] + irswitchvar.varseq[i][end:]
    if (spliced.translate().find("*")>-1):
        print(str(i))
        print(spliced)
        print(spliced.translate())
        irswitchvar.drop(i,inplace=True)

    
irswitchvar.to_pickle('./design/ir_switch_splice_site_sequences.pkl')

#%% length

irlength= irmain

linker=Seq('gaagctgctgcaagagaagctgcagctagggaggctgcagctagggaggctgctgcaaga')
irlengthvar=pd.DataFrame()

k=0

for i in irlength.index:
    intronlength=int(irlength.intronlength[i])
    start=int(irlength.intronstart_varseq[i])+6
    end=int(irlength.intronend_varseqnew[i])-31
    count = (end-start)//6
    print(count)
    for ex in range(0,count-1,1):
        if (ex<9):
            irlengthvar=irlengthvar.append(irlength.loc[i],ignore_index=True)
            irlengthvar.varseq[k]=linker[-6-ex*6:] + irlength.varseq[i][:start] + irlength.varseq[i][start+6 +ex*6:]
            irlengthvar.loc[k,'intronstart_varseq']=str(start + ex*6)
            irlengthvar.loc[k,'intronlength']=str(intronlength - 6 - ex*6)
            irlengthvar.loc[k,'changes']=str(6 + ex*6) + 'nt deletion at position 6'
            k=k+1
    
            irlengthvar=irlengthvar.append(irlength.loc[i],ignore_index=True)
            middle = start + (((end-start)//3)//2)*3
            irlengthvar.varseq[k]=linker[-6-ex*6:] + irlength.varseq[i][:middle-3-ex*3] + irlength.varseq[i][middle+3+ex*3:]
            irlengthvar.loc[k,'intronstart_varseq']=str(start + ex*6)
            irlengthvar.loc[k,'intronlength']=str(intronlength - 6 - ex*6)
            irlengthvar.loc[k,'changes']=str(6 + ex*6) + 'nt deletion at center'
            k=k+1
            
            irlengthvar=irlengthvar.append(irlength.loc[i],ignore_index=True)
            irlengthvar.varseq[k]=linker[-6-ex*6:] + irlength.varseq[i][:end-6-ex*6] + irlength.varseq[i][end:]
            irlengthvar.loc[k,'intronstart_varseq']=str(start + ex*6)
            irlengthvar.loc[k,'intronlength']=str(intronlength - 6 - ex*6)
            irlengthvar.loc[k,'changes']=str(6 + ex*6) + 'nt deletion at position -31'
            k=k+1

for i in irlengthvar.index:
    start=int(irlengthvar.intronstart_varseq[i])
    end=int(irlengthvar.intronend_varseqnew[i])
    spliced=irlengthvar.varseq[i][:start] + irlengthvar.varseq[i][end:]
    if (spliced.translate().find("*")>-1):
        print(str(i))
        print(spliced)
        print(spliced.translate())
#        irlengthvar.drop(i,inplace=True)

    
    
irlengthvar.to_pickle('./design/ir_lengthvars.pkl')
    
