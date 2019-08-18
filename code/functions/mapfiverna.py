#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from Bio import SeqIO
import sys

mapfolder=''

#%%
filenumber = sys.argv[1]
fqread1 = SeqIO.to_dict(SeqIO.parse(mapfolder + 'Split1-' + str(filenumber),'fastq'))
fqread2 = SeqIO.to_dict(SeqIO.parse(mapfolder + 'Split2-' + str(filenumber),'fastq'))


sublib=pd.read_pickle('../../code/design/LIBRARY/fivelibrary210.pkl')
sublib['barcode']=sublib.varseq.apply(lambda x: x[18:30])
readsmap=pd.Series([''],index=sublib.index)

liftfor='GACAGATGCGCCGTGGAT'
#upstream='ACTAGTTTACGACGGGTT'

for read in fqread1.keys():
    if (fqread1[read].seq.find(liftfor)>-1):
        testbc=fqread1[read].seq[fqread1[read].seq.find(liftfor)+18:fqread1[read].seq.find(liftfor)+30]
        readsmap.loc[sublib[sublib.barcode==testbc].index]=readsmap.loc[sublib[sublib.barcode==testbc].index] +' '+ [str(fqread2[read].seq)]           
    elif (fqread2[read].seq.find(liftfor)>-1):
        testbc=fqread2[read].seq[fqread2[read].seq.find(liftfor)+18:fqread2[read].seq.find(liftfor)+30]
        readsmap.loc[sublib[sublib.barcode==testbc].index]=readsmap.loc[sublib[sublib.barcode==testbc].index] +' '+ [str(fqread1[read].seq)]           


readsmap.to_pickle(mapfolder + 'coveragePYTHON-' + str(filenumber) + '.pkl')
