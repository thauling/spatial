#this is an attempt to write a script to generate bridge probe libraries, WORKS!
import os
#from sys import maxint
# Load the Pandas libraries with alias 'pd'
import pandas as pd
#import numpy as np
import random

print(os.getcwd())
# load probedesign input file and the n extract barcodes
pdinput = pd.read_csv("WholeBrain_99-5pc.csv")

# Preview the first 5 lines of the loaded data
print(pdinput.head())
#CALCULATE encoding scheme and save as csv file, eg. 'encoding.csv'#

#eg RANDOM, generate random encoder for one gene
print (random.sample(set('0123456'), 7)) # set = channels, ,7 specifies 7 random picks = one/ round
#rndencoding = (random.sample(set('0123456'), 7))

# generate random list for every gene
genenumber = len(pdinput)
genelist = pdinput['gene']
encodingarray = pd.DataFrame()


for x in genelist:
    print (random.sample(set('0123456'), 7))   

##create random encoder seq for n genes (as specified in pdinput file), for nrounds and create encoding array
#rndround = (random.sample(set('0123456'), len(pdinput)))
#encodingarray = pd.concat([encodingarray, genelist])
#nrounds = 7

#for y in range (0, nrounds):
#    encodingarray[y] = (random.sample(set('0123456'), len(pdinput)))
### end of random encoding

### use this for Reed-Solomon encoding, requires csv file with decoding scheme
# load Reed-Solomon (or other encoding scheme), containing gene names (plus barcodes?) plus single digit encoding strings, e.g. gene1, barcode1, 1,3,5,0,2,4,6
#encoding = pd.read_csv("encoding_file_test.csv")
encoding = pd.read_csv("encoding_99-5pc.csv")
print(encoding.head())
encoding_numerical = encoding
# load dye probe library containing dye probe reverse complementary sequences plus corresponding digits,eg
# 0, GCTTGCAA...
# 1, CACGTAAC...
# 2,...
dyelib = pd.read_csv("dyelib.csv")
print(dyelib.head())

#convert integers to strings in encoding
encoding = encoding.astype(str)
print(encoding)

#replace channel numbers with reverse complementary sequences in encoding variable.
#first specifcy corresponding reverse complemntary sequences
ch0 = dyelib.at[0, 'reverse_complementary_seq']
ch1 = dyelib.at[1, 'reverse_complementary_seq']
ch2 = dyelib.at[2, 'reverse_complementary_seq']
ch3 = dyelib.at[3, 'reverse_complementary_seq']
ch4 = dyelib.at[4, 'reverse_complementary_seq']
ch5 = dyelib.at[5, 'reverse_complementary_seq']
ch6 = dyelib.at[6, 'reverse_complementary_seq']

encoding.replace({"0": ch0}, inplace = True)
encoding.replace({"1": ch1}, inplace = True)
encoding.replace({"2": ch2}, inplace = True)
encoding.replace({"3": ch3}, inplace = True)
encoding.replace({"4": ch4}, inplace = True)
encoding.replace({"5": ch5}, inplace = True)
encoding.replace({"6": ch6}, inplace = True)

#now, for every gene, combine (i.e. concatenate barcode and round-specifc dye-sequence
bridge_r0 = encoding.barcode + encoding.r0
bridge_r1 = encoding.barcode + encoding.r1
bridge_r2 = encoding.barcode + encoding.r2
bridge_r3 = encoding.barcode + encoding.r3
bridge_r4 = encoding.barcode + encoding.r4
bridge_r5 = encoding.barcode + encoding.r5
bridge_r6 = encoding.barcode + encoding.r6

bridge_r0 = bridge_r0.to_frame()
bridge_r1 = bridge_r1.to_frame()
bridge_r2 = bridge_r2.to_frame()
bridge_r3 = bridge_r3.to_frame()
bridge_r4 = bridge_r4.to_frame()
bridge_r5 = bridge_r5.to_frame()
bridge_r6 = bridge_r6.to_frame()
#genelist = genelist.tolist()
genelist = 'r0_' + genelist.astype(str)
bridge_r0['gene'] = genelist
bridge_r1['gene'] = genelist
bridge_r2['gene'] = genelist
bridge_r3['gene'] = genelist
bridge_r4['gene'] = genelist
bridge_r5['gene'] = genelist
bridge_r6['gene'] = genelist

#write to csv
bridge_r0.to_csv('bridge_r0.csv', index=False, sep=",")
bridge_r1.to_csv('bridge_r1.csv', index=False, sep=",")
bridge_r2.to_csv('bridge_r2.csv', index=False, sep=",")
bridge_r3.to_csv('bridge_r3.csv', index=False, sep=",")
bridge_r4.to_csv('bridge_r4.csv', index=False, sep=",")
bridge_r5.to_csv('bridge_r5.csv', index=False, sep=",")
bridge_r6.to_csv('bridge_r6.csv', index=False, sep=",")

##or alternatively, for in-script generated random array
encodingarray_num = encodingarray

encodingarray.replace({"0": ch0}, inplace = True)
encodingarray.replace({"1": ch1}, inplace = True)
encodingarray.replace({"2": ch2}, inplace = True)
encodingarray.replace({"3": ch3}, inplace = True)
encodingarray.replace({"4": ch4}, inplace = True)
encodingarray.replace({"5": ch5}, inplace = True)
encodingarray.replace({"6": ch6}, inplace = True)

#make bridge probes by adding barcodes


#print (rndencoding)
print (genenumber)
print (encodingarray)
print (genelist)


   