# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 11:54:39 2020

@author: thomas.hauling
"""

#### padlock extraction script. accepts padlock design output file and list with probe numbers per gene -> creates padlock purchase list

import os
import pandas as pd
#import re
import math
import collections
from Bio.Seq import Seq


os.chdir(r"C:\Users\thomas.hauling\WORK\Python_Scripts\padlock filtering\extract numbers of padlocks based on gene name")
##specfify padlock input file

padlocks = pd.read_csv(r"C:\Users\thomas.hauling\WORK\Python_Scripts\padlock filtering\extract numbers of padlocks based on gene name\padlock_input.csv")

##specify probes to extract per gene

probenumbers = pd.read_csv("probestoorder.csv")

##round probenumbers to integer values
##create probenumberdict
df = padlocks
##remove duplicates
df = df.drop_duplicates('padlock', keep = 'first')

#to_drop = ['GGGG', 'CCCC', 'AAAAA', 'AAAA', 'TTTT']

#df1 = df.iloc[:,0:5] ##dont use
index = df.index
columns = df.columns
values = df.values

##################################################################
plseq = df['padlock'].dropna()
target = df['target']
acronym = df['acronym']
combined = df[['acronym', 'padlock']] ##use this
mygenes = probenumbers['Gene']
myprobes = probenumbers['nProbes to buy']

todrop = []
for pls in plseq:
    if pls.find("AAAA") > -1:
        todrop.append(pls)

zdf = df
zdf[~zdf['padlock'].isin(todrop)]
##make gene:number of probes dict## works - keep
probeints = []

for a in myprobes:
    intnumber = math.ceil(a)
    probeints.append(intnumber)
    
probedict = dict(zip(probenumbers[probenumbers.columns[0]], probeints))
#floatdict = dict(zip(probenumbers[probenumbers.columns[0]], probenumbers[probenumbers.columns[1]]))
### this works!!

unique = df.acronym.unique()
unique = unique.tolist()
#### seems fine
fasta = []
genenames = []
for a2 in df['acronym']:
    if a2.startswith('>'):
       fasta.append(a2)
    else:
        genenames.append(a2)

found = []
notfound = []

probecounter=collections.Counter(genenames)
confirmed = []
notwanted = []

for name in genenames:
    if name in probedict:
        confirmed.append(name)
    else:
        notwanted.append(name)

uniconfirmed = list(set(confirmed))

#candidate = 'Gnas'
enoughprobes = []
notenough = []

actualVSwanted = []


for candidate in uniconfirmed:
    
    v1 = probecounter.get(candidate)
    v2 = probedict.get(candidate)
    data = pd.DataFrame([[candidate, v1, v2]], columns=list('GHW'))
    actualVSwanted.append(data)
    
    if v1 >= v2:
        enoughprobes.append(candidate)
    else:
        notenough.append(candidate)
###this works
### EXTRACT n RANDOM probes 
enoughext = []
enoughexta = []
enoughdf = []

for a in enoughprobes:
    b = probedict.get(a)
    if a in unique:
        found.append(a)
        c = df.loc[df['acronym']==a]
        cran = c['padlock'].sample(n=b, random_state=1, replace=False) #True - duplicates can occur but won t crash, False - no duplicates but requires probes >= probe selection 
        acran = c['acronym'].sample(n=b, random_state=1, replace=False)
        enoughext.append(cran)
        enoughexta.append(acran)
        cran = cran.tolist()
        acran = set(acran)
        data2 = pd.DataFrame([[acran, cran]], columns=list('AB'))
        enoughdf.append(data2)
            

notenoughext = []
notenoughexta = []
notenoughdf = []

for a in notenough:
    b = probedict.get(a)
    if a in unique:
        found.append(a)
        c = df.loc[df['acronym']==a]
        cran = c['padlock'].sample(n=b, random_state=1, replace=True) #True - duplicates can occur but won t crash, False - no duplicates but requires probes >= probe selection 
        acran = c['acronym'].sample(n=b, random_state=1, replace=True)
        notenoughext.append(cran)
        notenoughexta.append(acran)   
        cran = cran.tolist()
        acran = set(acran)
        data3 = pd.DataFrame([[acran, cran]], columns=list('AB'))
        notenoughdf.append(data3)
           
        
    else:
        notfound.append(a)           


###clean up data befroe saving
##enoughdf = enoughdf['A'].astype(str) NOT WORKING
        
#######!!!!!!!!!!!!!!!!!! THIS WORKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!      
###link index of randomext to gene name via df(index)

###SAVE all probes to csv file
        
    ##SPECIFY names#
outputpre = "demo"
outputsuf = "output.csv"  

######################################3
#newname = outputpre + "_PROBESenough_" + outputsuf   
        
pd.concat(enoughext).to_csv(outputpre + '_PROBESenough_' + outputsuf, index = False, header = True)
pd.concat(enoughexta).to_csv(outputpre + '_ACRONYMSenough_' + outputsuf, index = False, header = True)
pd.concat(notenoughext).to_csv(outputpre + '_PROBESnotenough_' + outputsuf, index = False, header = True)
pd.concat(notenoughexta).to_csv(outputpre + '_ACRONYMSnotenough_' + outputsuf, index = False, header = True)
pd.concat(actualVSwanted).to_csv(outputpre + '_actualVSwanted_' + outputsuf, index = False, header = True)
#
#
pd.concat(enoughdf).to_csv(outputpre + '_enoughDF_' + outputsuf, index = False, header = True)
#pd.concat(notenoughdf).to_csv(outputpre + '_notenoughDF_' + outputsuf, index = False, header = True)

####MAKE PRIMERS####
##enough padlcoks
expadlocks = pd.read_csv(outputpre + '_PROBESenough_' + outputsuf)
exacr = pd.read_csv(outputpre + '_ACRONYMSenough_' + outputsuf)
primrc = expadlocks['padlock'].str[5:15] #specify length and position of primer
#prim = primrc.reverse_complement()

primlist = []

for dna in primrc:
    myseq = Seq(dna)
    prim = myseq.reverse_complement()
    primlist.append(prim)


primlist2 = pd.DataFrame(map(str, primlist))
primlist2.to_csv(outputpre + '_primers_' + outputsuf, index = False, header = True)
primrc.to_csv(outputpre + '_primers-rc_' + outputsuf, index = False, header = True)
####
#############################################################################################
##NOT enough padlcoks
expadlocks = pd.read_csv(outputpre + '_PROBESnotenough_' + outputsuf)
exacr = pd.read_csv(outputpre + '_ACRONYMSnotenough_' + outputsuf)
primrc = expadlocks['padlock'].str[5:15] #specify length and position of primer
#prim = primrc.reverse_complement()

primlist = []

for dna in primrc:
    myseq = Seq(dna)
    prim = myseq.reverse_complement()
    primlist.append(prim)


primlist2 = pd.DataFrame(map(str, primlist))
primlist2.to_csv(outputpre + '_primers_notEnough_' + outputsuf, index = False, header = True)
primrc.to_csv(outputpre + '_primers-rc_notEnough_' + outputsuf, index = False, header = True)







