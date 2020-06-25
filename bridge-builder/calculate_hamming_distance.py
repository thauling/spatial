# -*- coding: utf-8 -*- calculate hamming distance, not working yet
"""
Created on Wed Jan  8 17:02:42 2020

@author: thomas.hauling
"""

import distance
import pandas as pd
distance.hamming("hamming", "hamning")

files = pd.read_csv('7dyes.csv')

seq3 = ['4012','0212','3112','0123', '4321', '3241'] #str has to be used for digit sequences STARTING with 0!

def convert(string):
        li = list(string.split(' '))
        return li
a = seq3[0]
b = seq3[0+1]

a = convert(a)

help(distance)

hammingdist = distance.hamming('01234', '01212', normalized=True)

seq1 = ['0123', '4321', '3241']

seq2 = ['1233', '1234', '2413']



hammingdist2 = distance.hamming(seq1, seq2, normalized=True)


for i in range(len(seq3)):
    print(i)
    print(len(seq3))
    c = []
    d = []
   # a = a.append('blub')
    #c = c.append[seq3[i-1]]
    #d = d.append[seq3[i]]
    #hammingdist = distance.hamming(c,d, normalized=True)
    