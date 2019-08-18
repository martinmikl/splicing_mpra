# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:36:32 2016

@author: martinm
"""

import pandas as pd
from Bio.Seq import Seq


#%% import splice junctions from ENCODE (Thomas Gingeras Lab, K562 RNA seq, 2011, BED file from ENCODE website, ENCFF000HFA). 


splicejunctionsk562=pd.read_table('../additional/ENCFF000HFAsplicejunctions.bed')

splicejunctionsk562.columns=['chr','start','end','name','score','strand','level','sig','score2']

junctions=splicejunctionsk562[splicejunctionsk562.sig<0.05]


#%% extract alternative splicing events

junctions_alt3prime=pd.DataFrame()

for i in range(len(junctions)):
    if (i<len(junctions)-1):
        if (junctions.loc[junctions.index[i],'start']==junctions.loc[junctions.index[i+1],'start']):
            junctions_alt3prime= junctions_alt3prime.append(junctions.iloc[i,:])
            junctions_alt3prime.loc[junctions.index[i],'end_alt']= \
                junctions.loc[junctions.index[i+1],'end']
            junctions_alt3prime.loc[junctions.index[i],'score_alt']= \
                junctions.loc[junctions.index[i+1],'score']
            junctions_alt3prime.loc[junctions.index[i],'level_alt']= \
                junctions.loc[junctions.index[i+1],'level']
            junctions_alt3prime.loc[junctions.index[i],'sig_alt']= \
                junctions.loc[junctions.index[i+1],'sig']
            junctions_alt3prime.loc[junctions.index[i],'diff_nt']= \
                junctions_alt3prime.loc[junctions.index[i],'end_alt']- \
                junctions_alt3prime.loc[junctions.index[i],'end']
            
junctions_alt5prime=pd.DataFrame()

for i in range(len(junctions)):
    if (i<len(junctions)-1):
        if (junctions.loc[junctions.index[i],'end']==junctions.loc[junctions.index[i+1],'end']):
            junctions_alt5prime= junctions_alt5prime.append(junctions.iloc[i,:])
            junctions_alt5prime.loc[junctions.index[i],'start_alt']= \
                junctions.loc[junctions.index[i+1],'start']
            junctions_alt5prime.loc[junctions.index[i],'score_alt']= \
                junctions.loc[junctions.index[i+1],'score']
            junctions_alt5prime.loc[junctions.index[i],'level_alt']= \
                junctions.loc[junctions.index[i+1],'level']
            junctions_alt5prime.loc[junctions.index[i],'sig_alt']= \
                junctions.loc[junctions.index[i+1],'sig']
            junctions_alt5prime.loc[junctions.index[i],'diff_nt']= \
                junctions_alt5prime.loc[junctions.index[i],'start_alt']- \
                junctions_alt5prime.loc[junctions.index[i],'start']
            
jun5pr_plus=junctions_alt5prime[(junctions_alt5prime.strand=='+')& \
(junctions_alt5prime.diff_nt<80)&(junctions_alt5prime.diff_nt%3>0)]    

jun5pr_minus=junctions_alt3prime[(junctions_alt3prime.strand=='-')& \
(junctions_alt3prime.diff_nt<80)&(junctions_alt3prime.diff_nt%3>0)] 

jun5pr_minus=jun5pr_minus.rename(columns = {'end':'start_alt','start':'end','end_alt':'start'})

jun5pr=pd.DataFrame().append([jun5pr_plus, jun5pr_minus])

### PICKLE
jun5pr.to_pickle('./design/jun5pr.pkl')


#%% define window to get sequence containing the future varseq for 5pr splice sites: 150 nt upstream and downstream of first 5pr splice site


f=open('./design/jun5pr.bed', 'w')
for i in jun5pr.index:
    f.write(jun5pr.loc[i,'chr'] + "\t" + str(int(jun5pr.loc[i,'start'])-150) + \
        "\t" + str(int(jun5pr.loc[i,'start'])+150) + "\t" + jun5pr.loc[i,'name'] + \
        "\t" + str(jun5pr.loc[i,'diff_nt']) + "\t" + jun5pr.loc[i,'strand'] + "\n")
        
f.close()


#%% get sequences in fasta format using bedtools

'''
 bedtools getfasta -fi hg19.fa 
     -bed ./design/jun5pr.bed 
     -fo ./design/jun5prseqaroundsplicesite.fa 
     -s

'''

#%% parse fasta file

from Bio import SeqIO

f=open('./design/jun5prseqaroundsplicesite.fa','r')
records=list(SeqIO.parse(f, 'fasta'))

c=0
for i in jun5pr.index:
    jun5pr.loc[i,'coordinates']=records[c].id
    jun5pr.loc[i,'sequence']=records[c].seq
    c=c+1

for i in jun5pr.index:
    jun5pr.loc[i,'stopinseqbefore_0']=bool(jun5pr.loc[i,'sequence'][90:150].translate().find('*')==-1)
    jun5pr.loc[i,'stopinseqbefore_1']=bool(jun5pr.loc[i,'sequence'][88:150].translate().find('*')==-1)
    jun5pr.loc[i,'stopinseqbefore_2']=bool(jun5pr.loc[i,'sequence'][89:150].translate().find('*')==-1)

    jun5pr.loc[i,'stopinseqbeforealt_0']=bool(jun5pr.loc[i,'sequence'][90:150+int(jun5pr.loc[i,'diff_nt'])].translate().find('*')==-1)
    jun5pr.loc[i,'stopinseqbeforealt_1']=bool(jun5pr.loc[i,'sequence'][88:150+int(jun5pr.loc[i,'diff_nt'])].translate().find('*')==-1)
    jun5pr.loc[i,'stopinseqbeforealt_2']=bool(jun5pr.loc[i,'sequence'][89:150+int(jun5pr.loc[i,'diff_nt'])].translate().find('*')==-1)

jun5prselected=pd.DataFrame(columns=list(jun5pr.columns))

for i in jun5pr.index:
    if (int(jun5pr.loc[i,'diff_nt'])%3==2)&(jun5pr.loc[i,'stopinseqbefore_0']==True)&(jun5pr.loc[i,'stopinseqbeforealt_0']==True):
        jun5prselected.loc[i,:]=jun5pr.loc[i,:]

# Write to bed file


f=open('./design/jun5prselected.bed', 'w')
for i in jun5prselected.index:
    f.write(jun5prselected.loc[i,'chr'] + "\t" + str(int(jun5prselected.loc[i,'start'])-150) + \
        "\t" + str(int(jun5prselected.loc[i,'start'])+150) + "\t" + jun5prselected.loc[i,'name'] + \
        "\t" + str(i) + "\t" + jun5prselected.loc[i,'strand'] + "\n")
        
f.close()

'''
intersectBed -wao -s -a jun5prselected.bed -b hg19refseqgenes.bed > jun5prselected_intersectwithrefseqgenes.bed
'''


jun5prselectedrefseq=pd.read_table('./design/jun5prselected_intersectwithrefseqgenes.bed', header=None)

refseqgenes=pd.read_table('../additional/hg19refseqgenes.txt', header=None)
refseqgenes=refseqgenes.drop_duplicates('name')
refseqgenes.set_index('name', inplace=True)

for i in jun5prselectedrefseq.index:
    if (jun5prselectedrefseq.loc[i,'accession']!='.'):
        jun5prselectedrefseq.loc[i,'commonname']=refseqgenes.loc[str(jun5prselectedrefseq.loc[i,'accession']),'name2']
        
for i in jun5prselected.index:
    jun5prselected.loc[i,'commonname']=jun5prselectedrefseq.loc[jun5prselected.loc[i,'name'],'commonname']
    jun5prselected.loc[i,'accession']=jun5prselectedrefseq.loc[jun5prselected.loc[i,'name'],'accession']

for i in jun5prselected.index:
    if (jun5prselected.loc[i,'accession']=='.'):
        jun5prselected.drop(i,inplace=True)
        
### PICKLE
jun5prselected.to_pickle('./design/jun5prselected.pkl')
jun5prselected=pd.read_pickle('./design/jun5prselected.pkl')
        
#%% for 3prime
        
jun3pr_plus=junctions_alt3prime[(junctions_alt3prime.strand=='+')& \
(junctions_alt3prime.diff_nt<80)&(junctions_alt3prime.diff_nt%3>0)]    

jun3pr_minus=junctions_alt5prime[(junctions_alt5prime.strand=='-')& \
(junctions_alt5prime.diff_nt<80)&(junctions_alt5prime.diff_nt%3>0)] 

jun3pr_minus=jun3pr_minus.rename(columns = {'end':'start','start':'end_alt','start_alt':'end'})

jun3pr=pd.DataFrame().append([jun3pr_plus, jun3pr_minus])
        
        
# define window to get sequence containing the future varseq for 3pr splice sites: 150 nt upstream and downstream of first 3pr splice site


f=open('./design/jun3pr.bed', 'w')
for i in jun3pr.index:
    f.write(jun3pr.loc[i,'chr'] + "\t" + str(int(jun3pr.loc[i,'end'])-150) + \
        "\t" + str(int(jun3pr.loc[i,'end'])+150) + "\t" + jun3pr.loc[i,'name'] + \
        "\t" + str(jun3pr.loc[i,'diff_nt']) + "\t" + jun3pr.loc[i,'strand'] + "\n")
        
f.close()

# get sequences in fasta format using bedtools

'''
bedtools getfasta -fi hg19.fa -bed jun3pr.bed -fo jun3prseqaroundsplicesite.fa -s

intersectBed -wao -s -a jun3pr.bed -b hg19refseqgenes.bed > jun3pr_intersectwithrefseqgenes.bed

'''

f=open('./design/jun3prseqaroundsplicesite.fa','r')
records=list(SeqIO.parse(f, 'fasta'))

c=0
for i in jun3pr.index:
    jun3pr.loc[i,'coordinates']=records[c].id
    jun3pr.loc[i,'sequence']=records[c].seq
    c=c+1


jun3prrefseq=pd.read_table('./design/jun3pr_intersectwithrefseqgenes.bed',header=None)

refseqgenes=pd.read_table('../additional/hg19refseqgenes.txt')
refseqgenes=refseqgenes.drop_duplicates('name')
refseqgenes.set_index('name', inplace=True)


jun3prrefseq=jun3prrefseq.rename(columns = {9:'accession'})
jun3prrefseq=jun3prrefseq.rename(columns = {3:'name'})

for i in jun3prrefseq.index:
    if (jun3prrefseq.loc[i,'accession']!='.'):
        jun3prrefseq.loc[i,'commonname']=refseqgenes.loc[str(jun3prrefseq.loc[i,'accession']),'name2']

jun3prrefseq=jun3prrefseq.drop_duplicates('name')    

jun3prrefseq.set_index('name',inplace=True)    
    
for i in jun3pr.index:
    jun3pr.loc[i,'commonname']=jun3prrefseq.loc[jun3pr.loc[i,'name'],'commonname']
    jun3pr.loc[i,'accession']=jun3prrefseq.loc[jun3pr.loc[i,'name'],'accession']

for i in jun3pr.index:
    if (jun3pr.loc[i,'accession']=='.'):
        jun3pr.drop(i,inplace=True)


#%%

for i in jun3pr.index:

    jun3pr.loc[i,'stopinseqafteralt_0']=bool(jun3pr.loc[i,'sequence'][150:210].translate().find('*')==-1)
    jun3pr.loc[i,'stopinseqafteralt_1']=bool(jun3pr.loc[i,'sequence'][151:210].translate().find('*')==-1)
    jun3pr.loc[i,'stopinseqafteralt_2']=bool(jun3pr.loc[i,'sequence'][152:210].translate().find('*')==-1)

jun3prselected=pd.DataFrame(columns=list(jun3pr.columns))

for i in jun3pr.index:
    if (int(jun3pr.loc[i,'diff_nt'])%3==2)&(jun3pr.loc[i,'stopinseqafteralt_1']==True):
        jun3prselected.loc[i,:]=jun3pr.loc[i,:]


for i in jun3prselected.index:
    jun3prselected.loc[i,'varseq'] = jun3prselected.sequence[i][71:221]

for i in jun3prselected.index:
    jun3prselected.loc[i,'nostopinspliced'] = bool(jun3prselected.varseq[i][80:].translate().find("*")==-1)

for i in jun3prselected.index:
    if (jun3prselected.nostopinspliced[i]==False):
        jun3prselected.drop(i,inplace=True)

### PICKLE
jun3prselected.to_pickle('./design/jun3prselectedNEW.pkl')



#%% cassette exons

junctions_cassette=pd.DataFrame()

for i in range(len(junctions)-10):
    if (i<len(junctions)-1):
        if (junctions.loc[junctions.index[i],'start']==junctions.loc[junctions.index[i+1],'start']):
            cassette=0
            for scan in range(i+2,i+10):
                if (junctions.loc[junctions.index[scan],'end']==junctions.loc[junctions.index[i+1],'end']):
                    cassette=scan
                    break
            
            if (cassette>0) & ((junctions.index[i] in junctions_cassette.index)==False):
                junctions_cassette= junctions_cassette.append(junctions.iloc[i,:])
                junctions_cassette.loc[junctions.index[i],'end2']= \
                    junctions.loc[junctions.index[i+1],'end']
                junctions_cassette.loc[junctions.index[i],'score_skipped']= \
                    junctions.loc[junctions.index[i+1],'score']
                junctions_cassette.loc[junctions.index[i],'level_alt']= \
                    junctions.loc[junctions.index[i+1],'level']
                junctions_cassette.loc[junctions.index[i],'sig_alt']= \
                    junctions.loc[junctions.index[i+1],'sig']
                junctions_cassette.loc[junctions.index[i],'start2']= \
                    junctions.loc[junctions.index[cassette],'start']
                junctions_cassette.loc[junctions.index[i],'score2']= \
                    junctions.loc[junctions.index[cassette],'score']
                junctions_cassette.loc[junctions.index[i],'exonlength']= \
                    junctions_cassette.loc[junctions.index[i],'start2']- \
                    junctions_cassette.loc[junctions.index[i],'end']
            
            
            

            
juncas=junctions_cassette[(junctions_cassette.exonlength>0)&(junctions_cassette.exonlength<90)&(junctions_cassette.exonlength %3>0)]


f=open('./design/juncas.bed', 'w')
for i in juncas.index:
    f.write(juncas.loc[i,'chr'] + "\t" + str(int(juncas.loc[i,'start'])-100) + \
        "\t" + str(int(juncas.loc[i,'end2'])+100) + "\t" + juncas.loc[i,'name'] + \
        "\t" + str(juncas.loc[i,'exonlength']) + "\t" + juncas.loc[i,'strand'] + "\n")
        
f.close()


'''
 bedtools getfasta -fi hg19.fa -bed juncas.bed -fo juncasseqaroundsplicesite.fa -s

'''     

#parse fasta file

from Bio import SeqIO

f=open('./design/juncasseqaroundsplicesite.fa','r')
records=list(SeqIO.parse(f, 'fasta'))

c=0
for i in juncas.index:
    juncas.loc[i,'coordinates']=records[c].id
    juncas.loc[i,'sequence']=records[c].seq
    c=c+1

intron1=juncas.end[i]-juncas.start[i]
intron2=juncas.end2[i]-juncas.start2[i]
skippedlength=juncas.end2[i]-juncas.start[i]
sumincl=intron1 + juncas.exonlength[i] + intron2

print(juncas.loc[i,'sequence'][int(juncas.loc[i,'end']-juncas.loc[i,'start']+100 -10) : int(juncas.loc[i,'end']-juncas.loc[i,'start']+100 + juncas.loc[i,'exonlength']+10)])


for i in juncas.index:
    intron1=juncas.end[i]-juncas.start[i]
    intron2=juncas.end2[i]-juncas.start2[i]
    if (juncas.strand[i]=='+'):
        juncas.loc[i,'stopinexon_0']=bool(juncas.sequence[i][int(intron1 + 100) : int(- intron2-100)].translate().find('*')==-1)
        juncas.loc[i,'stopinexon_1']=bool(juncas.sequence[i][int(intron1 + 101) : int(- intron2-100)].translate().find('*')==-1)
        juncas.loc[i,'stopinexon_2']=bool(juncas.sequence[i][int(intron1 + 102) : int(- intron2-100)].translate().find('*')==-1)
    else:
        juncas.loc[i,'stopinexon_0']=bool(juncas.sequence[i][int(intron2 + 100) : int(- intron1-100)].translate().find('*')==-1)
        juncas.loc[i,'stopinexon_1']=bool(juncas.sequence[i][int(intron2 + 101) : int(- intron1-100)].translate().find('*')==-1)
        juncas.loc[i,'stopinexon_2']=bool(juncas.sequence[i][int(intron2 + 102) : int(- intron1-100)].translate().find('*')==-1)


### PICKLE
juncas.to_pickle('./design/juncas.pkl')



juncasselected=pd.DataFrame(columns=list(juncas.columns))

for i in juncas.index:
    if (int(juncas.loc[i,'exonlength'])%3==2)&(juncas.loc[i,'stopinexon_2']==True):
        juncasselected.loc[i,:]=juncas.loc[i,:]



# Write to bed file


f=open('./design/juncasselected.bed', 'w')
for i in juncasselected.index:
    f.write(juncasselected.loc[i,'chr'] + "\t" + str(int(juncasselected.loc[i,'start'])-150) + \
        "\t" + str(int(juncasselected.loc[i,'start'])+150) + "\t" + juncasselected.loc[i,'name'] + \
        "\t" + str(i) + "\t" + juncasselected.loc[i,'strand'] + "\n")
        
f.close()

'''
intersectBed -wao -s -a juncasselected.bed -b hg19refseqgenes.bed > juncasselected_intersectwithrefseqgenes.bed

'''

juncasselectedrefseq=pd.read_table('./design/juncasselecteddoublefirstline.bed', header=None)


for i in juncasselectedrefseq.index:
    if (juncasselectedrefseq.loc[i,'accession']!='.'):
        juncasselectedrefseq.loc[i,'commonname']=refseqgenes.loc[str(juncasselectedrefseq.loc[i,'accession']),'name2']
    
    
for i in juncasselected.index:
    juncasselected.loc[i,'commonname']=juncasselectedrefseq.loc[juncasselected.loc[i,'name'],'commonname']
    juncasselected.loc[i,'accession']=juncasselectedrefseq.loc[juncasselected.loc[i,'name'],'accession']

for i in juncasselected.index:
    if (juncasselected.loc[i,'accession']=='.'):
        juncasselected.drop(i,inplace=True)

juncasselected.drop_duplicates('sequence',inplace=True)

### PICKLE
juncasselected.to_pickle('./design/juncasselected.pkl')
juncasselected=pd.read_pickle('./design/juncasselected.pkl')

#%% INTRON RETENTION

introns=pd.read_table('../additional/Braunschweig2014_TableS6_IntronsHuman.txt')

import re

for i in range(136918,len(introns)-1):
    introns.loc[i,'strand']=re.split("[0123456789]", introns.Coordinates[i])[1]
    coordlist=re.split("[-_:+]", introns.Coordinates[i])
    introns.loc[i,'chromosome']=coordlist[0]
    introns.loc[i,'start']=int(coordlist[2])
    introns.loc[i,'end']=int(coordlist[3])
    introns.loc[i,'intronlength']=int(coordlist[3])-int(coordlist[2])
    introns.loc[i,'upstreamexon']=str(coordlist[0] + ':' + coordlist[1] + '-' + coordlist[2])
    introns.loc[i,'downstreamexon']=str(coordlist[0] + ':' + coordlist[3] + '-' + coordlist[4])


intronsshort=pd.DataFrame(columns=list(introns.columns))

for i in introns.index:
    if (introns.loc[i,'intronlength']<100):
        intronsshort.loc[i,:]=introns.loc[i,:]

# annotation file downloaded from UCSC table browser

hg19annotation=pd.read_table('../additional/hg19refgeneannotation.txt')

hg19exons=pd.DataFrame()

for i in hg19annotation.index:
    exonstarts=hg19annotation.exonStarts[i].split(',')
    exonends=hg19annotation.exonEnds[i].split(',')
    for j in range(0,len(exonstarts)-1):
        newrow=pd.Series([hg19annotation.chrom[i],exonstarts[j],exonends[j],hg19annotation.name[i],hg19annotation.name2[i],hg19annotation.strand[i]])
        hg19exons=hg19exons.append(newrow, ignore_index=True)


hg19exons.columns=list(['chr','start','end','name','name2','strand'])


### PICKLE

hg19exons.to_pickle('./design/hg19exons.pkl')

hg19exonsalt = hg19exons[((hg19exons.duplicated('start') == True) & \
    (hg19exons.duplicated('end') == False)) | \
    ((hg19exons.duplicated('start') == False) & (hg19exons.duplicated('end') == True))]


hg19exonsaltplus=pd.DataFrame(columns=list(['chr','start','end','name','name2','strand']))

for i in hg19exonsalt.index:
    indexlist = hg19exons[((hg19exons.chr==hg19exonsalt.chr[i]) & \
        ((hg19exons.start==hg19exonsalt.start[i]) | \
        (hg19exons.end==hg19exonsalt.end[i])))].index.tolist()
    hg19exonsaltplus=hg19exonsaltplus.append(hg19exons.loc[indexlist])

hg19exonsaltnondup=hg19exonsaltplus.drop_duplicates(['start','end'])


#test=hg19exonsaltnondup.head(300)

genedf=pd.DataFrame(columns=list(['chr','start','end','name','name2','strand']))
genedf.loc[0,'name2']=str(" ")
count=1

for i in hg19exonsaltnondup.index:
    if (hg19exonsaltnondup.name2[i]!=genedf.iloc[0,4]):
        genedf=hg19exonsaltnondup[hg19exonsaltnondup.name2==hg19exonsaltnondup.name2[i]]
    if ((len(genedf[(genedf.start==hg19exonsaltnondup.start[i])].index.tolist())>1) & \
            (len(genedf[(genedf.end==hg19exonsaltnondup.end[i])].index.tolist())>1)):
        hg19exonsaltnondup.loc[i,'exonwithintronretention']=True
    else:
        hg19exonsaltnondup.loc[i,'exonwithintronretention']=False
    count=count+1


### PICKLE
hg19exonsaltnondup.to_pickle('./design/hg19exonsaltnondup.pkl')

hg19exonswithintronretention=hg19exonsaltnondup[hg19exonsaltnondup.exonwithintronretention==True]
hg19exonsIR=pd.DataFrame()

for i in hg19exonswithintronretention.index:
    indexlist = hg19exonsaltnondup[((hg19exonsaltnondup.start == \
        hg19exonswithintronretention.start[i]) | (hg19exonsaltnondup.end == \
        hg19exonswithintronretention.end[i]))].index.tolist()
    if (len(indexlist)==3):
        hg19exonsIR = hg19exonsIR.append(hg19exonswithintronretention.loc[i])
        starts=pd.Series(hg19exonsaltnondup.loc[indexlist,'start']).drop_duplicates()
        ends=pd.Series(hg19exonsaltnondup.loc[indexlist,'end']).drop_duplicates()
        boundaries=pd.concat([starts,ends],ignore_index=True).order()
        hg19exonsIR.loc[i,'intronstart'] = int(boundaries.values[1])
        hg19exonsIR.loc[i,'intronend'] = int(boundaries.values[2])
        hg19exonsIR.loc[i,'intronlength'] = int(boundaries.values[2])-int(boundaries.values[1])
        hg19exonsIR.loc[i,'indexlist']=str(indexlist)

intronretention=hg19exonsIR[((hg19exonsIR.intronlength<121) & (hg19exonsIR.intronlength>20))]

for i in intronretention.index:
    if ((intronretention.duplicated('name2')[i]) \
            & ((intronretention.duplicated('start')[i]) | (intronretention.duplicated('end')[i]))):
                intronretention.loc[i,'checked']=False
    else:
        intronretention.loc[i,'checked']=True
        
for i in intronretention.index:
    if (len(intronretention.chr[i])>6):
        intronretention.drop(i, inplace=True)
       

f=open('./design/intronretention.bed', 'w')
for i in intronretention.index:
    f.write(intronretention.loc[i,'chr'] + "\t" + str(intronretention.loc[i,'start']) + \
        "\t" + str(intronretention.loc[i,'end']) + "\t" + intronretention.loc[i,'name2'] + \
        "\t" + str(intronretention.loc[i,'intronlength']) + "\t" + intronretention.loc[i,'strand'] + "\n")
        
f.close()


'''
 bedtools getfasta -fi hg19.fa -bed intronretention.bed -fo intronretentionsequences.fa -s

'''     

#parse fasta file

from Bio import SeqIO

f=open('./design/intronretentionsequences.fa','r')
records=list(SeqIO.parse(f, 'fasta'))

c=0
for i in intronretention.index:
    intronretention.loc[i,'coordinates']=records[c].id
    intronretention.loc[i,'sequence']=records[c].seq
    c=c+1

for i in intronretention.index:
    if ((int(intronretention.start[i])<int(intronretention.intronstart[i])) & ((int(intronretention.intronend[i])<int(intronretention.end[i])))):
        print(i)
    else:
        intronretention.drop(i, inplace=True)

for i in intronretention.index:
    intronretention.loc[i,'sequencelength']=len(intronretention.sequence[i])

for i in intronretention.index:
    if (intronretention.sequencelength[i] < 162):
        intronretention.drop(i,inplace=True)

for i in intronretention.index:
    if (intronretention.strand[i]=='+'):
        intronretention.loc[i,'firstexonlength']=int(intronretention.intronstart[i]) - int(intronretention.start[i])
        intronretention.loc[i,'secondexonlength']=int(intronretention.end[i]) - int(intronretention.intronend[i])
    elif (intronretention.strand[i]=='-'):
        intronretention.loc[i,'firstexonlength']=int(intronretention.end[i]) - int(intronretention.intronend[i])
        intronretention.loc[i,'secondexonlength']=int(intronretention.intronstart[i]) - int(intronretention.start[i])
    

for i in intronretention.index:
    intronretention.loc[i,'splicedexonseq']=Seq(str(intronretention.sequence[i][0:int(intronretention.firstexonlength[i])])+str(intronretention.sequence[i][-int(intronretention.secondexonlength[i]):]))
    intronretention.loc[i,'intronseq']=intronretention.sequence[i][int(intronretention.firstexonlength[i]):-int(intronretention.secondexonlength[i])]
        


for i in intronretention.index:
    intronretention.loc[i,'varseq']=intronretention.sequence[i][int(intronretention.firstexonlength[i]) + int(intronretention.intronlength[i]) - 145: int(intronretention.firstexonlength[i]) + int(intronretention.intronlength[i]) + 25]
    intronretention.loc[i,'intronend_varseq']=170-25;
    intronretention.loc[i,'intronstart_varseq']=170-25-intronretention.intronlength[i]


for i in intronretention.index:
    if (len(intronretention.varseq[i])<162):
        intronretention.loc[i,'varseq']=intronretention.sequence[i][:170]
        intronretention.loc[i,'intronend_varseq']=intronretention.firstexonlength[i] + intronretention.intronlength[i]
        intronretention.loc[i,'intronstart_varseq']=intronretention.firstexonlength[i]

for i in intronretention.index:
    intronretention.loc[i,'splicedvarseq']=intronretention.loc[i,'varseq'][:int(intronretention.intronstart_varseq[i])] + intronretention.loc[i,'varseq'][int(intronretention.intronend_varseq[i]):]
    intronretention.loc[i,'stopinspliced0']=bool(intronretention.loc[i,'splicedvarseq'][0:168].translate().find("*")==-1)
    intronretention.loc[i,'stopinspliced1']=bool(intronretention.loc[i,'splicedvarseq'][1:169].translate().find("*")==-1)
    intronretention.loc[i,'stopinspliced2']=bool(intronretention.loc[i,'splicedvarseq'][2:170].translate().find("*")==-1)

for i in intronretention.index:
    intronretention.loc[i,'stopinunspliced0']=bool(intronretention.loc[i,'varseq'].translate().find("*")==-1)
    intronretention.loc[i,'stopinunspliced1']=bool(intronretention.loc[i,'varseq'][1:].translate().find("*")==-1)
    intronretention.loc[i,'stopinunspliced2']=bool(intronretention.loc[i,'varseq'][2:].translate().find("*")==-1)

###PICKLE
intronretention.to_pickle('./design/intronretention.pkl')

intronretention=pd.read_pickle('./design/intronretention.pkl')

for i in intronretention.index:
    if ((intronretention.stopinspliced0[i]==False) & (intronretention.stopinspliced1[i]==False) & (intronretention.stopinspliced2[i]==False)):
        intronretention.drop(i, inplace=True)

intronretentionvariants=pd.DataFrame(columns=list(intronretention.columns))
j=0

for i in intronretention.index:
    if (intronretention.stopinspliced0[i]==True):
        intronretentionvariants=intronretentionvariants.append(intronretention.loc[i], ignore_index=True)
        intronretentionvariants.varseq[j]=intronretention.varseq[i][6:168]
        intronretentionvariants.intronend_varseq[j]=intronretention.intronend_varseq[i]-6
        intronretentionvariants.intronstart_varseq[j]=intronretention.intronstart_varseq[i]-6
        j=j+1
    if (intronretention.stopinspliced1[i]==True):
        intronretentionvariants=intronretentionvariants.append(intronretention.loc[i], ignore_index=True)
        intronretentionvariants.varseq[j]=intronretention.varseq[i][7:169]
        intronretentionvariants.intronend_varseq[j]=intronretention.intronend_varseq[i]-7
        intronretentionvariants.intronstart_varseq[j]=intronretention.intronstart_varseq[i]-7
        j=j+1
    if (intronretention.stopinspliced2[i]==True):
        intronretentionvariants=intronretentionvariants.append(intronretention.loc[i], ignore_index=True)
        intronretentionvariants.varseq[j]=intronretention.varseq[i][8:170]
        intronretentionvariants.intronend_varseq[j]=intronretention.intronend_varseq[i]-8
        intronretentionvariants.intronstart_varseq[j]=intronretention.intronstart_varseq[i]-8
        j=j+1

for i in intronretentionvariants.index:
    intronretentionvariants.loc[i,'splicedvarseq']=intronretentionvariants.loc[i,'varseq'][:int(intronretentionvariants.intronstart_varseq[i])] + intronretentionvariants.loc[i,'varseq'][int(intronretentionvariants.intronend_varseq[i]):]

for i in intronretentionvariants.index:
    if (len(intronretentionvariants.splicedvarseq[i])%3==0):
        stopposition=int(intronretentionvariants.intronstart_varseq[i])+ int(intronretentionvariants.intronlength[i]//2)
        varseqvar=intronretentionvariants.varseq[i][:stopposition]+intronretentionvariants.varseq[i][stopposition+2:] + Seq('CA')

        intronretentionvariants.varseq[i]=varseqvar
        intronretentionvariants.loc[i,'intronend_varseqnew']=intronretentionvariants.intronend_varseq[i]-2

    elif (len(intronretentionvariants.splicedvarseq[i])%3==1):
        stopposition=int(intronretentionvariants.intronstart_varseq[i])+ int(intronretentionvariants.intronlength[i]//2)
        varseqvar=intronretentionvariants.varseq[i][:stopposition]+intronretentionvariants.varseq[i][stopposition+1:] + Seq('C')
        intronretentionvariants.varseq[i]=varseqvar
        intronretentionvariants.loc[i,'intronend_varseqnew']=intronretentionvariants.intronend_varseq[i]-1
    elif (len(intronretentionvariants.splicedvarseq[i])%3==2):
        intronretentionvariants.loc[i,'intronend_varseqnew']=intronretentionvariants.intronend_varseq[i]

for i in intronretentionvariants.index:
    intronretentionvariants.loc[i,'splicedvarseq2']=intronretentionvariants.loc[i,'varseq'][:int(intronretentionvariants.intronstart_varseq[i])] + intronretentionvariants.loc[i,'varseq'][int(intronretentionvariants.intronend_varseqnew[i]):]
    intronretentionvariants.loc[i,'stopinsplicedvarseq']=bool(intronretentionvariants.loc[i,'splicedvarseq2'].translate().find("*")==-1)
    intronretentionvariants.loc[i,'stopinunplicedvarseq2']=bool(intronretentionvariants.varseq[i].translate().find("*")==-1)
    
for i in intronretentionvariants.index:
    intronretentionvariants.loc[i,'splicedexonlength']=len(intronretentionvariants.splicedvarseq[i])%3
    intronretentionvariants.loc[i,'splicedexonlength2']=len(intronretentionvariants.splicedvarseq2[i])%3
###PICKLE
intronretentionvariants.to_pickle('./design/intronretentionvariants.pkl')
      