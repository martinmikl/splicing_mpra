# -*- coding: utf-8 -*-
"""
Created on Wed Apr 06 16:33:31 2016

@author: martinm
"""


import pandas as pd

from Bio.Seq import Seq

import random


#%% cassette

cassettecomb=pd.read_pickle('./design/cassette_combinatorial_variations.pkl')
cassetteconstvar=pd.read_pickle('./design/cassetteprime_constitutiveandunsplicable_variations.pkl')
cassetteSFvar=pd.read_pickle('./design/cassette_splicingfactors_location.pkl')
cassetteSFcombvar=pd.read_pickle('./design/cassette_splicingfactors_combinatorial_exononly.pkl')
cassetteSFRosenberg=pd.read_pickle('./design/cassette_splicingfactors_Rosenberg.pkl')
cassettesecvar=pd.read_pickle('./design/cassette_secondarystructure_variants.pkl')
cassette_comb_with_constitutive=pd.read_pickle('./design/cassette_comb_with_constitutive.pkl')
cassettenuc=pd.read_pickle('./design/cassette_nucleosome_recoding_methylation.pkl')
cassetteswitchvar=pd.read_pickle('./design/cassette_switch_splice_site_sequences.pkl')
cassettelengthvar=pd.read_pickle('./design/cassette_exonlength.pkl')
cassettemelanomavariants=pd.read_pickle('./design/cassettemelanomavariants.pkl')
cassette_filtered=pd.read_pickle('./design/cassette_filtered.pkl')

cassette_filtered.loc[:,'library']='cassette'
cassette_filtered.loc[:,'subset']='cassette_filtered'

cassettecomb.loc[:,'library']='cassette'
cassettecomb.loc[:,'subset']='cassettecomb'

cassetteconstvar.loc[:,'library']='cassette'
cassetteconstvar.loc[:,'subset']='cassetteconstvar'

cassetteSFvar.loc[:,'library']='cassette'
cassetteSFvar.loc[:,'subset']='cassetteSFvar'

cassetteSFcombvar.loc[:,'library']='cassette'
cassetteSFcombvar.loc[:,'subset']='cassetteSFcombvar'

cassetteSFRosenberg.loc[:,'library']='cassette'
cassetteSFRosenberg.loc[:,'subset']='cassetteSFRosenberg'

cassettesecvar.loc[:,'library']='cassette'
cassettesecvar.loc[:,'subset']='cassettesecvar'

cassette_comb_with_constitutive.loc[:,'library']='cassette'
cassette_comb_with_constitutive.loc[:,'subset']='cassette_comb_with_constitutive'

cassettenuc.loc[:,'library']='cassette'
cassettenuc.loc[:,'subset']='cassettenuc'

cassetteswitchvar.loc[:,'library']='cassette'
cassetteswitchvar.loc[:,'subset']='cassetteswitchvar'

cassettelengthvar.loc[:,'library']='cassette'
cassettelengthvar.loc[:,'subset']='cassettelengthvar'

cassettemelanomavariants.loc[:,'library']='cassette'
cassettemelanomavariants.loc[:,'subset']='cassettemelanomavariants'

cassettenuc.drop(cassettenuc[cassettenuc.commonname=='LOC643201'].index,inplace=True)
cassettenuc.drop(cassettenuc[cassettenuc.commonname=='RNASE1'].index,inplace=True)


cassettelibrary=pd.concat([cassettecomb,cassetteconstvar,cassetteSFvar,cassetteSFcombvar,cassetteSFRosenberg,cassettesecvar,cassette_comb_with_constitutive,cassettenuc,cassetteswitchvar,cassettelengthvar,cassettemelanomavariants,cassettemelanomavariants,cassettemelanomavariants,cassette_filtered,cassette_filtered],ignore_index=True)

from Bio.Restriction import RsrII
from Bio.Restriction import AscI
from Bio.Restriction import SpeI
from Bio.Restriction import AatII


for i in cassettelibrary.index:
    if (RsrII.search(cassettelibrary.varseq[i])!=[])|(AscI.search(cassettelibrary.varseq[i])!=[])|(SpeI.search(cassettelibrary.varseq[i])!=[])|(AatII.search(cassettelibrary.varseq[i])!=[]):
        print(i)
        print(cassettelibrary.subset[i])        
        cassettelibrary.drop(i,inplace=True)
        
for i in cassettelibrary.index:
    if (len(cassettelibrary.varseq[i])!=150):
        print(i)
        print(cassettelibrary.subset[i])        
        

for i in cassettelibrary.index:
    if (cassettelibrary.varseq[i][122-int(cassettelibrary.exonlength[i]):120].translate().find("*")>-1):
        print(str(i))
        print(cassettelibrary.subset[i])        
        print(cassettelibrary.varseq[i][122-int(cassettelibrary.exonlength[i]):120].translate())
#        cassettelibrary.drop(i,inplace=True)

cassettelibrary.to_pickle('./design/LIBRARY/cassettelibrary.pkl')

#%% three

threecomb=pd.read_pickle('./design/threeprime_combinatorial_variations.pkl')
threeconstvar=pd.read_pickle('./design/threeprime_constitutiveandunsplicable_variations.pkl')
threeSFvar=pd.read_pickle('./design/three_splicingfactors_location.pkl')
threeSFcombvar=pd.read_pickle('./design/three_splicingfactors_combinatorial_exononly.pkl')
threeSFRosenberg=pd.read_pickle('./design/three_splicingfactors_Rosenberg.pkl')
threesecvar=pd.read_pickle('./design/three_secondarystructure_variants.pkl')
three_comb_with_constitutive=pd.read_pickle('./design/three_combinatorial_with_constitutive.pkl')
threenuc=pd.read_pickle('./design/three_nucleosome_recoding_methylation.pkl')
threeswitchvar=pd.read_pickle('./design/three_switch_splice_site_sequences.pkl')
threelengthvar=pd.read_pickle('./design/three_exonlength.pkl')
threemelanomavariants=pd.read_pickle('./design/threemelanomavariants.pkl')
three_filtered=pd.read_pickle('./design/three_filtered.pkl')

three_filtered.loc[:,'library']='three'
three_filtered.loc[:,'subset']='three_filtered'


threecomb.loc[:,'library']='three'
threecomb.loc[:,'subset']='threecomb'

threeconstvar.loc[:,'library']='three'
threeconstvar.loc[:,'subset']='threeconstvar'

threeSFvar.loc[:,'library']='three'
threeSFvar.loc[:,'subset']='threeSFvar'

threeSFcombvar.loc[:,'library']='three'
threeSFcombvar.loc[:,'subset']='threeSFcombvar'

threeSFRosenberg.loc[:,'library']='three'
threeSFRosenberg.loc[:,'subset']='threeSFRosenberg'

threesecvar.loc[:,'library']='three'
threesecvar.loc[:,'subset']='threesecvar'

three_comb_with_constitutive.loc[:,'library']='three'
three_comb_with_constitutive.loc[:,'subset']='three_comb_with_constitutive'

threenuc.loc[:,'library']='three'
threenuc.loc[:,'subset']='threenuc'

threeswitchvar.loc[:,'library']='three'
threeswitchvar.loc[:,'subset']='threeswitchvar'

threelengthvar.loc[:,'library']='three'
threelengthvar.loc[:,'subset']='threelengthvar'

threemelanomavariants.loc[:,'library']='three'
threemelanomavariants.loc[:,'subset']='threemelanomavariants'


threelibrary=pd.concat([threecomb,threeconstvar,threeSFvar,threeSFcombvar,threeSFRosenberg,threesecvar,three_comb_with_constitutive,threenuc,threeswitchvar,threelengthvar,threemelanomavariants,threemelanomavariants,threemelanomavariants,three_filtered,three_filtered,three_filtered],ignore_index=True)

from Bio.Restriction import RsrII
from Bio.Restriction import AscI
from Bio.Restriction import SpeI
from Bio.Restriction import AatII


for i in threelibrary.index:
    if (RsrII.search(threelibrary.varseq[i])!=[])|(AscI.search(threelibrary.varseq[i])!=[])|(SpeI.search(threelibrary.varseq[i])!=[])|(AatII.search(threelibrary.varseq[i])!=[]):
        print(i)
        print(threelibrary.subset[i])        
        threelibrary.drop(i,inplace=True)
        
for i in threelibrary.index:
    if (len(threelibrary.varseq[i])!=150):
        print(i)
        print(threelibrary.subset[i])        
        threelibrary.drop(i,inplace=True)
        

for i in threelibrary.index:
    if (threelibrary.varseq[i][80:].translate().find("*")>-1):
        print(str(i))
        print(threelibrary.subset[i])        
        print(threelibrary.varseq[i][80:].translate())
        threelibrary.drop(i,inplace=True)

threelibrary.to_pickle('./design/LIBRARY/threelibrary.pkl')

#%% five

fivecomb=pd.read_pickle('./design/fiveprime_combinatorial_variations.pkl')
fiveconstvar=pd.read_pickle('./design/fiveprime_constitutiveandunsplicable_variations.pkl')
fiveSFvar=pd.read_pickle('./design/five_splicingfactors_location.pkl')
fiveSFcombvar=pd.read_pickle('./design/five_splicingfactors_combinatorial_exononly.pkl')
fiveSFRosenberg=pd.read_pickle('./design/five_splicingfactors_Rosenberg.pkl')
fivesecvar=pd.read_pickle('./design/five_secondarystructure_variants.pkl')
five_comb_with_constitutive=pd.read_pickle('./design/five_combinatorial_with_constitutive.pkl')
fivenuc=pd.read_pickle('./design/five_nucleosome_recoding_methylation.pkl')
fiveswitchvar=pd.read_pickle('./design/five_switch_splice_site_sequences.pkl')
fivelengthvar=pd.read_pickle('./design/five_exonlength.pkl')
fivemelanomavariants=pd.read_pickle('./design/fivemelanomavariants.pkl')
five_filtered=pd.read_pickle('./design/five_filtered.pkl')

five_filtered.loc[:,'library']='five'
five_filtered.loc[:,'subset']='five_filtered'


fivecomb.loc[:,'library']='five'
fivecomb.loc[:,'subset']='fivecomb'

fiveconstvar.loc[:,'library']='five'
fiveconstvar.loc[:,'subset']='fiveconstvar'

fiveSFvar.loc[:,'library']='five'
fiveSFvar.loc[:,'subset']='fiveSFvar'

fiveSFcombvar.loc[:,'library']='five'
fiveSFcombvar.loc[:,'subset']='fiveSFcombvar'

fiveSFRosenberg.loc[:,'library']='five'
fiveSFRosenberg.loc[:,'subset']='fiveSFRosenberg'

fivesecvar.loc[:,'library']='five'
fivesecvar.loc[:,'subset']='fivesecvar'

five_comb_with_constitutive.loc[:,'library']='five'
five_comb_with_constitutive.loc[:,'subset']='five_comb_with_constitutive'

fivenuc.loc[:,'library']='five'
fivenuc.loc[:,'subset']='fivenuc'

fiveswitchvar.loc[:,'library']='five'
fiveswitchvar.loc[:,'subset']='fiveswitchvar'

fivelengthvar.loc[:,'library']='five'
fivelengthvar.loc[:,'subset']='fivelengthvar'

fivemelanomavariants.loc[:,'library']='five'
fivemelanomavariants.loc[:,'subset']='fivemelanomavariants'


fivelibrary=pd.concat([fivecomb,fiveconstvar,fiveSFvar,fiveSFcombvar,fiveSFRosenberg,fivesecvar,five_comb_with_constitutive,fivenuc,fiveswitchvar,fivelengthvar,five_filtered,five_filtered,five_filtered,fivemelanomavariants,fivemelanomavariants,fivemelanomavariants],ignore_index=True)

from Bio.Restriction import RsrII
from Bio.Restriction import AscI

for i in fivelibrary.index:
    if (RsrII.search(fivelibrary.varseq[i])!=[])|(AscI.search(fivelibrary.varseq[i])!=[]):
        print(i)
        print(fivelibrary.subset[i])        
        fivelibrary.drop(i,inplace=True)
        
for i in fivelibrary.index:
    if (len(fivelibrary.varseq[i])!=162):
        print(i)
        print(fivelibrary.subset[i])        
        fivelibrary.drop(i,inplace=True)
        

for i in fivelibrary.index:
    if (fivelibrary.varseq[i][:51+int(fivelibrary.diff_nt[i])].translate().find("*")>-1):
        print(str(i))
        print(fivelibrary.subset[i])        
        print(fivelibrary.varseq[i][:51+int(fivelibrary.diff_nt[i])].translate())
        fivelibrary.drop(i,inplace=True)

fivelibrary.to_pickle('./design/LIBRARY/fivelibrary.pkl')


#%% ir

ircomb=pd.read_pickle('./design/irmain_permutations.pkl')
ircombthreeway=pd.read_pickle('./design/intronretention_combinatorial_variations_threeway.pkl')
irconstvar=pd.read_pickle('./design/irprime_constitutiveandunsplicable_variations.pkl')
irSFvar=pd.read_pickle('./design/ir_splicingfactors_location.pkl')
irSFcombvar=pd.read_pickle('./design/ir_splicingfactors_combinatorial.pkl')
irSFRosenberg=pd.read_pickle('./design/ir_splicingfactors_Rosenberg.pkl')
irsecvar=pd.read_pickle('./design/ir_secondarystructure_variants.pkl')
ir_comb_with_constitutive=pd.read_pickle('./design/intronretention_combinatorial_variations_threeway_withconstitutive.pkl')
irnuc=pd.read_pickle('./design/ir__nucleosome_recoding_methylation.pkl')
irswitchvar=pd.read_pickle('./design/ir_switch_splice_site_sequences.pkl')
irlengthvar=pd.read_pickle('./design/ir_lengthvars.pkl')
irmelanomavariants=pd.read_pickle('./design/irmelanomavariants.pkl')
ir_filtered=pd.read_pickle('./design/irmain.pkl')

ir_filtered.loc[:,'library']='ir'
ir_filtered.loc[:,'subset']='ir_filtered'


ircomb.loc[:,'library']='ir'
ircomb.loc[:,'subset']='ircomb'

ircombthreeway.loc[:,'library']='ir'
ircombthreeway.loc[:,'subset']='ircombthreeway'

irconstvar.loc[:,'library']='ir'
irconstvar.loc[:,'subset']='irconstvar'

irSFvar.loc[:,'library']='ir'
irSFvar.loc[:,'subset']='irSFvar'

irSFcombvar.loc[:,'library']='ir'
irSFcombvar.loc[:,'subset']='irSFcombvar'

irSFRosenberg.loc[:,'library']='ir'
irSFRosenberg.loc[:,'subset']='irSFRosenberg'

irsecvar.loc[:,'library']='ir'
irsecvar.loc[:,'subset']='irsecvar'

ir_comb_with_constitutive.loc[:,'library']='ir'
ir_comb_with_constitutive.loc[:,'subset']='ir_comb_with_constitutive'

irnuc.loc[:,'library']='ir'
irnuc.loc[:,'subset']='irnuc'

irswitchvar.loc[:,'library']='ir'
irswitchvar.loc[:,'subset']='irswitchvar'

irlengthvar.loc[:,'library']='ir'
irlengthvar.loc[:,'subset']='irlengthvar'

irmelanomavariants.loc[:,'library']='ir'
irmelanomavariants.loc[:,'subset']='irmelanomavariants'


irlibrary=pd.concat([ircomb,ircombthreeway,irconstvar,irSFvar,irSFcombvar,irSFRosenberg,irsecvar,ir_comb_with_constitutive,irnuc,irswitchvar,irlengthvar,ir_filtered,ir_filtered,ir_filtered,irmelanomavariants,irmelanomavariants,irmelanomavariants],ignore_index=True)

from Bio.Restriction import RsrII
from Bio.Restriction import AscI

for i in irlibrary.index:
    if (RsrII.search(irlibrary.varseq[i])!=[])|(AscI.search(irlibrary.varseq[i])!=[]):
        print(i)
        print(irlibrary.subset[i])        
        irlibrary.drop(i,inplace=True)
        
for i in irlibrary.index:
    if (len(irlibrary.varseq[i])!=162):
        print(i)
        print(irlibrary.subset[i])     
        print(irlibrary.name2[i])        
        irlibrary.drop(i,inplace=True)
        
for i in irlibrary.index:
    start=int(irlibrary.intronstart_varseq[i])
    end=int(irlibrary.intronend_varseqnew[i])
    spliced=irlibrary.varseq[i][:start] + irlibrary.varseq[i][end:]
    if (spliced.translate().find("*")>-1):
        print(str(i))
        print(irlibrary.subset[i])        
        print(spliced.translate())
        irlibrary.drop(i,inplace=True)

irlibrary.to_pickle('./design/LIBRARY/irlibrary.pkl')


#%%


cassettelibrary=pd.read_pickle('./design/LIBRARY/cassettelibrary.pkl')
threelibrary=pd.read_pickle('./design/LIBRARY/threelibrary.pkl')
fivelibrary=pd.read_pickle('./design/LIBRARY/fivelibrary.pkl')
irlibrary=pd.read_pickle('./design/LIBRARY/irlibrary.pkl')

# ADD RESTRICTION SITE LINKER FOR CASSETTE AND THREE

'''
ACTAGTCTTGACGTC 
'''

linker='ACTAGTCTTGACGTC'

for i in cassettelibrary.index:
    cassettelibrary.varseq[i]=linker + str(cassettelibrary.varseq[i][3:])



for i in threelibrary.index:
    threelibrary.varseq[i]=linker + str(threelibrary.varseq[i][3:])


#%% add linker and barcodes

cassettelibrary=pd.read_pickle('./design/LIBRARY/cassettelibrary.pkl')
threelibrary=pd.read_pickle('./design/LIBRARY/threelibrary.pkl')
fivelibrary=pd.read_pickle('./design/LIBRARY/fivelibrary.pkl')
irlibrary=pd.read_pickle('./design/LIBRARY/irlibrary.pkl')



f=open("./design/barcodes_pool_12bps_3SNPs_filtered.tab",'r')

bcspl=f.read().splitlines()    

f.close()


linker='ACTAGTCTTGACGTC' # SpeI and AatII sites


primers=['GCCCCACGGAGGTGCCAC','CCTCCTCACGGCGACGCG','CTCCCGGGCATGCGAATT','TCAACCAGTCGCGGTCCA',
         'TGCGAGTTAGGGGACGGT','ACGGACGCGGGTATAGCA','CGAAATGGGCCGCATTGC','CACTGCGGCTGATGACGA',
            'GACAGATGCGCCGTGGAT','AGCCACCCGATCCAATGC','ATGGGGTTCGGTATGCGC','AAGGCTCCCCGAGACGAT']


'''

ir CGAAATGGGCCGCATTGC CACTGCGGCTGATGACGA
cas ATGGGGTTCGGTATGCGC AAGGCTCCCCGAGACGAT
five GACAGATGCGCCGTGGAT AGCCACCCGATCCAATGC
three TGCGAGTTAGGGGACGGT ACGGACGCGGGTATAGCA

'''


def addbcthree(vs):
    while True:
        bc=random.choice(bcspl)
        vsnew=primers[4] + bc + linker +str(vs[3:]) + primers[5]
        if (RsrII.search(Seq(vsnew))==[])&(AscI.search(Seq(vsnew))==[]):
            bcspl.remove(bc)
            break
    return vsnew;
    
threelibrary.varseq=threelibrary.varseq.apply(lambda x: addbcthree(x))

def addbccas(vs):
    while True:
        bc=random.choice(bcspl)
        vsnew=primers[10] + bc+ linker + str(vs[3:]) + primers[11]
        if (RsrII.search(Seq(vsnew))==[])&(AscI.search(Seq(vsnew))==[]):
            bcspl.remove(bc)
            break
    return vsnew;
    
cassettelibrary.varseq=cassettelibrary.varseq.apply(lambda x: addbccas(x))

def addbcfive(vs):
    while True:
        bc=random.choice(bcspl)
        vsnew=primers[8] + bc+str(vs) + primers[9]
        if (RsrII.search(Seq(vsnew))==[])&(AscI.search(Seq(vsnew))==[]):
            bcspl.remove(bc)
            break
    return vsnew;
    
fivelibrary.varseq=fivelibrary.varseq.apply(lambda x: addbcfive(x))

def addbcir(vs):
    while True:
        bc=random.choice(bcspl)
        vsnew=primers[6] + bc+str(vs) + primers[7]
        if (RsrII.search(Seq(vsnew))==[])&(AscI.search(Seq(vsnew))==[]):
            bcspl.remove(bc)
            break
    return vsnew;
    
irlibrary.varseq=irlibrary.varseq.apply(lambda x: addbcir(x))


threelibrary.to_pickle('./design/LIBRARY/threelibrary210.pkl')
threelibrary=pd.read_pickle('./design/LIBRARY/threelibrary210.pkl')
threelibrary['varseq162']=threelibrary.varseq.apply(lambda x: x[45:192])

cassettelibrary.to_pickle('./design/LIBRARY/cassettelibrary210.pkl')
cassettelibrary=pd.read_pickle('./design/LIBRARY/cassettelibrary210.pkl')
cassettelibrary['varseq162']=cassettelibrary.varseq.apply(lambda x: x[45:192])

fivelibrary.to_pickle('./design/LIBRARY/fivelibrary210.pkl')
fivelibrary=pd.read_pickle('./design/LIBRARY/fivelibrary210.pkl')
fivelibrary['varseq162']=fivelibrary.varseq.apply(lambda x: x[30:192])

irlibrary.to_pickle('./design/LIBRARY/irlibrary210.pkl')
irlibrary=pd.read_pickle('./design/LIBRARY/irlibrary210.pkl')
irlibrary['varseq162']=irlibrary.varseq.apply(lambda x: x[30:192])

#%%

alllibraries=pd.concat([fivelibrary,threelibrary,cassettelibrary,irlibrary],ignore_index=True)
alllibraries.to_pickle('./design/LIBRARY/alllibraries210.pkl')




