# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 14:31:45 2018

@author: gjgut
"""
import numpy as np
from scipy.stats import stats

A_rank = sorted(range(len(AvgInfs)), key=lambda i: AvgInfs[i])[-453:]
Pagerank_rank = sorted(range(len(PageRank)), key=lambda i: PageRank[i])[-453:]
Algo_rank = sorted(range(len(AlgoInfs)), key=lambda i: AlgoInfs[i])[-453:]
OutDeg_rank = sorted(range(len(OutDegree)), key=lambda i: OutDegree[i])[-453:]

A_rank.reverse()
Pagerank_rank.reverse()
Algo_rank.reverse()
OutDeg_rank.reverse()

positions = list(range(453))

dictionary = dict(zip(A_rank,positions))

A_rankings = [dictionary[i] for i in A_rank]
B_rankings = [dictionary[i] for i in Pagerank_rank]
Algo_rankings = [dictionary[i] for i in Algo_rank]
OutDeg_rankings = [dictionary[i] for i in OutDeg_rank]

concord = np.zeros(453)
discord = np.zeros(453)

for ii in B_rankings:
    for jj in B_rankings:
        if jj > ii:
            if B_rankings[ii] < B_rankings[jj]:
                concord[ii] = concord[ii] + 1
            elif B_rankings[ii] > B_rankings[jj]:
                discord[ii] = discord[ii] + 1
                
concord_total = np.sum(concord)
discord_total = np.sum(discord)

kendall_tau = (concord_total-discord_total)/(concord_total+discord_total)

print(kendall_tau)           
    
tau, p_value = stats.kendalltau(A_rankings,Algo_rankings)

print(tau)

tau2, p_value2 = stats.kendalltau(A_rankings,OutDeg_rankings)

print(tau2)
    