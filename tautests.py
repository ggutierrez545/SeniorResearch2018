# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:46:02 2018

@author: gjgut
"""

import numpy as np
from scipy.stats import stats

concord_test = np.zeros(len(G))
discord_test = np.zeros(len(G))

for ii in H:
    for jj in H:
        if jj > ii:
            if H[ii] < H[jj]:
                concord_test[ii] = concord_test[ii] + 1
            elif H[ii] > H[jj]:
                discord_test[ii] = discord_test[ii] + 1
                
concord_total = np.sum(concord_test)
discord_total = np.sum(discord_test)

kendall_tau = (concord_total-discord_total)/(concord_total+discord_total)

print(kendall_tau)           
    
tau, p_value = stats.kendalltau(G,H)

print(tau)
    