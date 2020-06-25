# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 14:32:42 2020

@author: thomas.hauling
"""

###filter padlocks and their primers for GC content

import os
import pandas as pd
import numpy as np

##specify directory that contains padlock sequences 
os.chdir(r"C:\Users\thomas.hauling\WORK\Python_Scripts\padlock filtering")
##specify csv file that contains padlock sequences 
df1 = pd.read_csv("demo_padlock_csv_file.csv")

df2 = df1.drop_duplicates('padlock', keep = 'first')
df3= df2['padlock'].dropna()
df3 = df3[df3.str[0].str.isupper()]
good = df2.dropna() ###remove non-padlock lines from data frame

#############################################################################NEWNEWNEW################
#convert dataframe to numpy arrays for easier processing

goodab = good.to_numpy()
goodab2 = goodab.astype('U')

###remove all the crap### polyA tails, >= quadGs (since difficult to synthesize). How to specify >= AAAAA etc?
remove = np.sum([(np.char.count(goodab2, 'AAAAAAA')), (np.char.count(goodab2, 'AAAAAA')), (np.char.count(goodab2, 'AAAAA')), (np.char.count(goodab2, 'GGGGG')),(np.char.count(goodab2, 'GGGG'))], axis=0)

goodab2 = goodab2[remove[:,1] == 0]
#######done####

### code to filter for GC:AT ratio
##gives 6column arrays (col1 = targets, col5 = padlocks)
freqA = np.char.count(goodab2[:,1], 'A')
freqT = np.char.count(goodab2[:,1], 'T')
freqC = np.char.count(goodab2[:,1], 'C')
freqG = np.char.count(goodab2[:,1], 'G')
###########################################3
atsum = np.sum([freqA, freqT], axis=0)
gcsum = np.sum([freqC, freqG], axis=0)
atgcsum = np.sum([atsum, gcsum], axis=0)
errors = atgcsum[atgcsum<30] ##extract those that have fewer than 30nt = length of seq = sum of As, Ts, Gs and Cs
errors2 = goodab2[atgcsum<30]
fine = goodab2[atgcsum == 30]

goodab2 = goodab2[atgcsum == 30]
freqA = np.char.count(goodab2[:,1], 'A')
freqT = np.char.count(goodab2[:,1], 'T')
freqC = np.char.count(goodab2[:,1], 'C')
freqG = np.char.count(goodab2[:,1], 'G')
atsum = np.sum([freqA, freqT], axis=0)
gcsum = np.sum([freqC, freqG], axis=0)
atgcsum = np.sum([atsum, gcsum], axis=0)
errors3 = atgcsum[atgcsum<30] ##approved fro further processing IF ZERO (0) array 
goodab2 = goodab2[atgcsum == 30]
#### now analysis
atgcdiv = np.divide(atsum, gcsum) ### 1s are 50/50 ratio ###remove any "non-30s before this operation!!
#################
comba = np.column_stack((goodab2,freqA, freqT, freqC, freqG, atgcdiv))
#combadf = pd.DataFrame(comba)
comba2 = np.column_stack((freqA, freqT, freqC, freqG, atgcdiv))

seq50 = comba[comba2[:,4] == 1] ####always 1 (0.5*max length)/(0.5*max length)  
seq4060 = comba[(comba2[:,4] >= 0.67) & (comba2[:,4] <= 1.5)] ##dont forget '=' character

#seq40plus = comba[comba2[:,4] > 0.67] ####max length (30) *0.4 (12) / max length - max length (30) *0.4 (18)
#seq60less = comba[comba2[:,4] < 1.5] ####max length (30) *0.6 (18)/ max length - max length (30) *0.6 (12)

## WORKS!!!
## now EXPORT!
## save padlocks that fullfill the 50:50 and 40:60 or 60:40 GC-to-AT ratio  
np.savetxt("seq50.csv", seq50, delimiter=",", fmt='%s') ##fmt='%s' is important for correct dtype format 
np.savetxt("seq4060.csv", seq4060, delimiter=",", fmt='%s') ##fmt='%s' is important for correct dtype format 

###############now primer based selection###### export those padlocks that have primers with 40 - 60% GC content  

#convert back to numpy arrays
goodab2df = pd.DataFrame(goodab2) 

##specify primer template

primrc = goodab2df[5].str[5:15] #specify length and position of primer based on padlock column (starting from 5prim +n (here 5))

##find primers that contain GC clamp
gcclamp = primrc.str[0:5]#GC clamp = one or two G and/ or C within 5 bases from 3'end 
primrc = primrc.to_numpy()
primrc = primrc.astype('U')
gcclamp = gcclamp.to_numpy()
gcclamp = gcclamp.astype('U')
freqA2 = np.char.count(primrc, 'A')
freqT2 = np.char.count(primrc, 'T')
freqC2 = np.char.count(primrc, 'C')
freqG2 = np.char.count(primrc, 'G')
freqC3 = np.char.count(gcclamp, 'C')
freqG3 = np.char.count(gcclamp, 'G')
atsum2 = np.sum([freqA2, freqT2], axis=0)
gcsum2 = np.sum([freqC2, freqG2], axis=0)
atgcsum2 = np.sum([atsum2, gcsum2], axis=0)
errors4 = atgcsum2[atgcsum2<10] ##adjust ot seq length - 10nt , check ln 83
errors5 = primrc[atsum2 == 0]
errors6 = primrc[gcsum2 == 0]
goodabOLD = goodab2
goodab2 = goodab2[(gcsum2 > 0) & (atsum2 > 0)]
primrcOLD = primrc
primrc = primrc[(gcsum2 > 0) & (atsum2 > 0)] 
gcclampOLD = gcclamp
gcclamp = gcclamp[(gcsum2 > 0) & (atsum2 > 0)] 
freqA2 = freqA2[(gcsum2 > 0) & (atsum2 > 0)]
freqT2 = freqT2[(gcsum2 > 0) & (atsum2 > 0)]
freqC2 = freqC2[(gcsum2 > 0) & (atsum2 > 0)]
freqG2 = freqG2[(gcsum2 > 0) & (atsum2 > 0)]
atsum2 = np.sum([freqA2, freqT2], axis=0)
gcsum2 = np.sum([freqC2, freqG2], axis=0)
atgcsum2 = np.sum([atsum2, gcsum2], axis=0)

#### now analysis
atgcdiv2 = np.divide(atsum2, gcsum2) 
#combadf = pd.DataFrame(comba)
comba3 = np.column_stack((primrc,goodab2, atgcdiv2))
comba4 = np.column_stack((freqA2, freqT2, freqC2, freqG2, atgcdiv2))
seq50p = comba3[comba4[:,4] == 1]
seq4060p = comba3[(comba4[:,4] >= 0.67) & (comba4[:,4] <= 1.5)]


## export probes with 50 and 40to60% GC ratio in target seq
np.savetxt("seq50p.csv", seq50p, delimiter=",", fmt='%s') ##fmt='%s' is important for correct dtype format whatever
np.savetxt("seq4060p.csv", seq4060p, delimiter=",", fmt='%s') ##fmt='%s' is important for correct dtype format whatever

### G/C clamp
#gcclampOLD = gcclamp
#gcclamp = gcclamp[(gcsum2 > 0) & (atsum2 > 0)] 
freqC3 = np.char.count(gcclamp, 'C')
freqG3 = np.char.count(gcclamp, 'G')

gcsum3 = np.sum([freqC3, freqG3], axis=0)

###select all GC clamp primers
gcclampALL = comba3[(gcsum3 >= 1) & (gcsum3 <= 2)]

###from these devide into bins: GC50% and GC40-60%
gcclamp50 = comba3[(gcsum3 >= 1) & (gcsum3 <= 2) & (comba4[:,4] == 1)]

#gcclamp50 = gcclampALL[(gcsum3 >= 1) & (gcsum3 <= 2)]
gcclamp4060 = comba3[(gcsum3 >= 1) & (gcsum3 <= 2) & (comba4[:,4] >= 0.67) & (comba4[:,4] <= 1.5)]
#### works!!
# export padlocks with primers that have a GC content between 40 and 60% and include a GC clamp
np.savetxt("gcclampALL.csv", gcclampALL, delimiter=",", fmt='%s') ##fmt='%s' is important for correct dtype format whatever
np.savetxt("gcclamp50.csv", gcclamp50, delimiter=",", fmt='%s') 
np.savetxt("gcclamp4060.csv", gcclamp4060, delimiter=",", fmt='%s')
