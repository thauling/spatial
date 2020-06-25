# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 14:29:29 2020

@author: thomas.hauling
"""

#################################################check BARCODES####################################
import os
import pandas as pd
import numpy as np

os.chdir(r"C:\Users\thomas.hauling\WORK\Python_Scripts\padlock filtering\filter barcodes for GC content")

dfb = pd.read_csv("barcodelib.csv")
###remove nonsense
dfb = dfb.drop_duplicates('barcode', keep = 'first')
dfb = dfb['barcode'].dropna()
#######################################
#remove those with terminal Gs?? possibly important if oligos conjugated to dyes - otherwise IGNORE and outcomment this paragrapgh
#termGs = pd.DataFrame()
#termGs['5'] = dfb.str[0]
#termGs[ '3'] = dfb.str[-1]
#termGs = termGs.to_numpy()
#termGs = termGs.astype('U')
#freqGt = np.char.count(termGs, 'G')
#dfb = dfb[(freqGt[:,0] == 0) & (freqGt[:,1] == 0)] ##creates SERIES

#######################################convert back to numpy
dfb = dfb.to_numpy()
dfb = dfb.astype('U')
####################
###remove polyXs
dfpORIGINAL = dfb
#remove = np.sum([(np.char.count(dfb, 'AAAAAA')), (np.char.count(dfb, 'AAAAA')), (np.char.count(dfb, 'GGGG')), (np.char.count(dfb, 'CCCC')),(np.char.count(dfb, 'TTTT'))], axis=0)
remove = np.sum([(np.char.count(dfb, 'AAAAA')), (np.char.count(dfb, 'TTTTT')), (np.char.count(dfb, 'GGGG')), (np.char.count(dfb, 'CCCC'))], axis=0)
dfb = dfb[remove == 0]
#######################################

#######################################
freqAb = np.char.count(dfb, 'A')
freqTb = np.char.count(dfb, 'T')
freqCb = np.char.count(dfb, 'C')
freqGb = np.char.count(dfb, 'G')
atsumb = np.sum([freqAb, freqTb], axis=0)
gcsumb = np.sum([freqCb, freqGb], axis=0)
atgcsumb = np.sum([atsumb, gcsumb], axis=0)
atgcdivb = np.divide(atsumb, gcsumb)
bc50 = dfb[atgcdivb == 1]
bc4060 = dfb[(atgcdivb >= 0.67) & (atgcdivb <= 1.5)]

###now save this stuff

np.savetxt("filtered_barcodes_50.csv", bc50, delimiter=",", fmt='%s') ##fmt='%s' is important for correct dtype format whatever
np.savetxt("filtered_barcodes_4060.csv", bc4060, delimiter=",", fmt='%s') ##fmt='%s' is important for correct dtype format whatever

#np.savetxt("bc50inc_tGs.csv", bc50, delimiter=",", fmt='%s') ##fmt='%s' is important for correct dtype format whatever
#np.savetxt("bc4060inc_tGs.csv", bc4060, delimiter=",", fmt='%s') ##fmt='%s' is important for correct dtype format whatever