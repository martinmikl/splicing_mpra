#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
first argument: 
    ir      retained introns
    cas     cassette exons
    five    tandem 5' splice sites
    three   tandem 3' splice sites
second argument:
    name of the csv file containing four columns: an identifier, the DNA sequence of the region containing a splice site, the first coordinate (first splice site, according to splicing type) and the second coordinate (second splice site, according to splicing type)
third argument:
    name of output file (optional)
'''

import pandas as pd
import sys
import forprediction_ir
import forprediction_cas
import forprediction_five
import forprediction_three

try:
    outputfilename=sys.argv[3]
except:
    outputfilename='./predictions'


if sys.argv[1]=='ir':
    df_features=forprediction_ir.make_features_ir(sys.argv[2])
    forprediction_ir.predict_splicingratio_all_models(df_features, \
                                                  filename=outputfilename)
elif sys.argv[1]=='cas':
    df_features=forprediction_cas.make_features_cas(sys.argv[2])
    forprediction_cas.predict_splicingratio_all_models(df_features, \
                                                  filename=outputfilename)
elif sys.argv[1]=='five':
    df_features=forprediction_five.make_features_five(sys.argv[2])
    forprediction_five.predict_splicingratio_all_models(df_features, \
                                                  filename=outputfilename)
elif sys.argv[1]=='three':
    df_features=forprediction_three.make_features_three(sys.argv[2])
    forprediction_three.predict_splicingratio_all_models(df_features, \
                                                  filename=outputfilename)
else:
    print 'Input type not specified; use ir, cas, five and three only'
    sys.exit()