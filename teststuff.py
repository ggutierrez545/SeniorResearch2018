# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 22:02:15 2018

@author: gjgut
"""

plt.plot(A, label='Max Transition Probability = 1')
plt.plot(B, label='Max Transition Probability = 0.1')
plt.show()
plt.legend(loc='upper right')
plt.xlabel('Node ID', fontsize=15)
plt.ylabel('Node Rank', fontsize=15)

A_top = sorted(range(len(A)), key=lambda i: A[i])[-10:]
B_top = sorted(range(len(B)), key=lambda i: B[i])[-10:]


from scipy.stats import stats

tau, p_value = stats.kendalltau(A,B)


Pagerank2 = []
for ii in range(len(Pagerank)):
    for jj in range(len(Pagerank[ii])):
        Pagerank2.append(Pagerank[ii][jj])
        
Algo_rank = sorted(range(len(Algorithm)), key=lambda i: Algorithm[i])[-5:]