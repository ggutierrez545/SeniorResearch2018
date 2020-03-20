# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 13:45:49 2018

@author: gjgut
"""

import numpy as np
import pandas as pd

A = pd.read_csv('110_Senate.csv')
B = A.as_matrix()

reps, bills = np.shape(B)

#find number of primary sponsorships directly from data
primeships = np.zeros(reps)
for ii in range(reps):
    primeships[ii] = np.count_nonzero(B[ii] == 1)

#find each time a bill primarily sponsored by a rep was cosponsored. This list will
#contain duplicates that reflect when a specific rep cosponosored multiple bills
#from that primary sponsor    
all_cospons =  [[] for i in range(reps)]
for ii in range(bills):
    templist = []
    for jj in range(reps):
        #make a list of cosponsors for a certain bill
        if B[jj][ii] == 2:
            templist.append(jj)
        #find the primary sponsor of that bill and keep track of them    
        elif B[jj][ii] == 1:
            prime = jj
    #update the primary sponsors list of cosponsors        
    all_cospons[prime] = all_cospons[prime] + templist

reps_admat = np.zeros([reps,reps])
for ii in range(len(all_cospons)):
    for jj in range(len(all_cospons[ii])):
        reps_admat[ii][all_cospons[ii][jj]] = reps_admat[ii][all_cospons[ii][jj]] + 1
        
reps_weights = reps_admat/np.max(reps_admat)

RepOutList = [[] for i in range(reps)]
RepOutWeights = [[] for i in range(reps)]
for ii in range(reps):
    templist = []
    tempweights = []
    for jj in range(reps):
        if reps_weights[ii][jj] > 0.0:
            templist.append(jj)
            tempweights.append(reps_weights[ii][jj])
    RepOutList[ii] = templist
    RepOutWeights[ii] = tempweights
    
            
        
       



        
    
            
                        
            
        