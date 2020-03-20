# -*- coding: utf-8 -*-
"""
calculates in-degree PageRank (where nodes "score points" for having incoming connections) for a given network
defined by input adjacency matrix. This implementation takes the Markov chain approach (see p. 282 of Lay's second
edition Linear Algebra text), with each column of Atilde giving the probability of transition to other states.
These probabilities are based on either following existing connections with uniform probability (weighted by alpha) 
or jumping to a random node (weighted by 1-alpha). 
I have checked this routine on a few simple 3-node networks (such as this: http://pr.efactory.de/e-pagerank-algorithm.shtml)
and it worked perfectly. Written by Christian G. Fink
"""

def pagerank(A,alpha=0.85,tol=0.001):
    
    import numpy as np
    import time
    
    N=len(A) #number of nodes
    kout=np.sum(A,axis=0) #vector of out-degree values for each node
    #if out-degree is equal to 0, make it equal to 1 so that you don't divide by zero when creating Atilde
    for i in range(len(kout)):
        if kout[i]==0:
            kout[i]=1.0
    
    start = time.time()         
    #create matrix that we can apply power method to in order to compute PR
    Atilde = np.zeros((N,N))
    for j in range(len(kout)):
        Atilde[:,j:j+1]=alpha*A[:,j:j+1]/kout[j] + ((1-alpha)/N) * np.ones((N,1))
        
    err=1.0 #set to nonsense value to ensure the while loop is entered
    PR_prev=(1/N)*np.ones((N,1)) #nonsense values to compute error on first iteration
    while(err>tol):
        PR = np.dot(Atilde,PR_prev)
        PR = PR/np.sum(PR) #it is crucial that all the elements of PR sum to 1, since we are interpreting this as steady-state probability distribution
        err = np.sum(abs(PR-PR_prev))
        PR_prev = np.copy(PR)
        
    PageRank = []
    for ii in range(len(PR)):
        for jj in range(len(PR[ii])):
            PageRank.append(PR[ii][jj])
    
    end = time.time()
    runtime = end - start
    
    return PageRank, runtime