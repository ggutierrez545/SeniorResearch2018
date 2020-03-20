# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 10:18:36 2017

@author: gjgut
"""

def InLists(AdjacencyMatrix):
    
    import numpy as np
   
    nrows, ncols = AdjacencyMatrix.shape
    if nrows != ncols:
        print(f"Error, not adjacencey matrix - has {nrows} rows vs {ncols} columns.") 
        return
    
    Inlist = [[] for j in range(ncols)] 
    InWeights = [[] for j in range(ncols)] 
    
    #index by columns first to get inward connections to node
    for col in range(ncols): 
        for row in range(nrows):
            if AdjacencyMatrix[row][col] > 0.0: #if there is a connection
                Inlist[col].append(row)
                InWeights[col].append(AdjacencyMatrix[row][col]) 
                
    InDeg = np.zeros(ncols)
    for ii in range(len(Inlist)):
        InDeg[ii] = np.count_nonzero(Inlist[ii])
                
    return Inlist, InWeights, InDeg, ncols


def OutLists(AdjacencyMatrix):
    
    import numpy as np
    
    nrows, ncols = AdjacencyMatrix.shape
    if nrows != ncols:
        print(f"Error, not adjacencey matrix - has {nrows} rows vs {ncols} columns.") 
        return
    
    Outlist = [[] for j in range(ncols)] 
    OutWeights = [[] for j in range(ncols)] 
    
    #index by rows first to get outward connections to node
    for row in range(nrows): 
        for col in range(ncols):
            if AdjacencyMatrix[row][col] > 0.0: #if there is a connection
                Outlist[row].append(col)
                OutWeights[row].append(AdjacencyMatrix[row][col])
                
    OutDeg = np.zeros(ncols)
    for ii in range(len(Outlist)):
        OutDeg[ii] = np.count_nonzero(Outlist[ii])
                
    return Outlist, OutWeights, OutDeg, ncols


def InList_AdjMat(InList,InWeights):
    
    import numpy as np
    nodes = len(InList)
    AdjMat = np.zeros([nodes,nodes])

    for col in range(nodes):
        for row in range(len(InList[col])):
            AdjMat[InList[col][row]][col] = InWeights[col][row]
            
    return AdjMat


def OutList_AdjMat(OutList,OutWeights):
    
    import numpy as np
    nodes = len(OutList)
    AdjMat = np.zeros([nodes,nodes])
    
    for row in range(nodes):
        for col in range(len(OutList[row])):
            AdjMat[row][OutList[row][col]] = OutWeights[row][col]
    
    return AdjMat


def GroundTrials(DirAdjMat,Trials,Beta):
    '''
    DirAdjMat is a directed, weighted adjacency matrix.
    Trials are the number of individual simulations to perform per node.
    Beta is a renormalization variable which sets the max weight in network to Beta.
    When a network has many cyclical connections and high probs on transmission,
    (weights close to one), it is necessary to renormalize the weights so
    as to mitigate runaway transmission probability propagation. Also, this Beta 
    should match the Beta used in the McAlgorithm function as that function's purpose
    is to approximate the results of this function. Set Beta to one if not using
    renormalization.
    '''
    import numpy as np
    import random as rd
    import matplotlib.pyplot as plt
    import time

    OutList, OutWeights, OutDeg, Nodes = OutLists(DirAdjMat)
    SimHist = np.zeros([Trials, Nodes])

    start = time.time()
    for infnode in range(Nodes):
        print(infnode)
        for trial in range(Trials):
    
            NowState = np.zeros(Nodes)
            NxtState = np.zeros(Nodes)
            NowState[infnode] = 1
            InfHist = np.zeros((Nodes + 1,Nodes))

            cont_sim = 1
            t = 0
            while cont_sim != 0:
                cont_sim = 0
                InfHist[t] = NowState #set row t InfectionHistory to the updated NowStates 
                for n in range(len(NowState)): #cycle through NowStates
                    if NowState[n] == 1: #if infected
                        cont_sim = 1 #as long as one node is infected, the trial continues
                        NxtState[n] = 2 #set infected node to recovered in NextStates 
                        for out_con in range(len(OutWeights[n])): # cycle through infnode's outgoing cons
                            if NowState[OutList[n][out_con]] == 0:
                                frand = rd.uniform(0,1) #roll the dice
                                if frand < Beta*OutWeights[n][out_con]: #if infection successful. Beta multiplies weights here
                                    NxtState[OutList[n][out_con]] = 1 #set suscept node to infected in next state list
                NowState = np.copy(NxtState) #set updated node state list as original node states for cycle
                t = t + 1
 
            SimHist[trial][infnode] = np.count_nonzero(InfHist[t-1]) - 1
            #keeps track of each node's propagation of infection for every trial
            #subtracts one at end to exclude the initial beginning of the infection at that node
    
    end = time.time()
    runtime = end - start
    
    NAvgInfs = SimHist.mean(axis=0)
    MeanInfs = NAvgInfs.mean()
    StdDev = NAvgInfs.std()
    Error = SimHist.std(axis=0)/np.sqrt(Trials)
    plt.errorbar(range(Nodes),NAvgInfs,yerr=Error,fmt='-',ecolor='tomato')  
    plt.title(f'Average Infections Caused by Seed Node (N={Nodes},Trials={Trials},Mean={MeanInfs:.3f},\u03c3={StdDev:.3f},Runtime={runtime:.3f} seconds)', fontsize=15)
    plt.xlabel('Node ID',fontsize=15)
    plt.ylabel('Infections',fontsize=15)

    return SimHist, NAvgInfs, runtime


def McAlgorithm(DirAdjMat,Beta,Termination,GroundAvgInfs,GroundSimHist,GroundRunTime):
    '''
    this algorithm operates using the inward bound connections, instead of the outward bound connections.
    Beta is a renormalization variable. Multiplying by Beta makes the max weight in network Beta
    this term dictates how long the sim will get propagate the smaller the term, the longer the 
    probabilities will be calculated
    '''
    import time
    import numpy as np
    import matplotlib.pyplot as plt
    
    InList, InWeights, InDeg, Nodes = InLists(DirAdjMat)
    AvgInfs = np.zeros(Nodes)

    start = time.time()
    for seed in range(Nodes):
        print(seed)
        prev_uninfected = np.ones(Nodes) #prob that node has remained uninfected until t-1
        uninfected = np.ones(Nodes) #prob that node has remained uninfected until t. starts at one for each node
        last_infected = np.zeros(Nodes) #probability that node was infected on last timestep
        cur_infected = np.zeros(Nodes) #probability that node was infected on current timestep
        last_infected[seed] = 1
        uninfected[seed] = 0
    
        t = 0
        while sum(prev_uninfected - uninfected) > Termination:
            prev_uninfected = np.copy(uninfected)
            
            for node in range(Nodes):
                prob_uninfected = 1 #set probability of specific node being uninfected to one
                
                for con in range(len(InList[node])): #look at each incoming connection
                    prob_uninfected = prob_uninfected*(1-(last_infected[InList[node][con]]*(Beta*InWeights[node][con])))
                    #prob that node remains uninfected decreases as we go through each connection.
                
                cur_infected[node] = (1-prob_uninfected)*uninfected[node]
                #prob of current infection is prob of being infected times the prob that node hasn't been infected yet
        
            for node in range(Nodes):
                last_infected[node] = cur_infected[node] 
                #update our last infected list for each timestep
                uninfected[node] = uninfected[node] - cur_infected[node]
                #new prob of being uninfected is old prob - the prob that node is currently infected
        
            t = t+1
            
        AvgInfs[seed] = sum(1 - uninfected) - 1 #dont want to include seed node infection in total
    
    end = time.time()
    runtime = end - start
    
    MeanInfs = AvgInfs.mean()
    GroundMeanInfs = GroundAvgInfs.mean()
    
    MeanPercDif = ((GroundMeanInfs - MeanInfs)/GroundMeanInfs)*100
    RunTimePercDif =((GroundRunTime - runtime)/GroundRunTime)*100 
    
    GroundError = GroundSimHist.std(axis=0)/np.sqrt(len(GroundSimHist)) 
    plt.plot(range(Nodes),AvgInfs,color='red',label='Probability Algorithm')
    plt.errorbar(range(Nodes),GroundAvgInfs,yerr=GroundError,fmt='-',ecolor='green',label='Ground Truth')
    plt.legend(loc='upper right')
    plt.title(f'Comparison b/w Alrogithm and Ground Truth Results (N={Nodes},'
              f' Mean and % Difference=({MeanInfs:.3f}, {MeanPercDif:.3f})'
              f' Runtime and % Difference=({runtime:.3f}, {RunTimePercDif:.3f})')
    plt.xlabel('Node ID',fontsize=15)
    plt.ylabel('Infections',fontsize=15)
               
    return AvgInfs


def KendallsTau(GroundAvgInfs,AlgoAvgInfs,PageRanks,OutDegree):
    
    import numpy as np
    from scipy.stats import stats
    
    Ground_rank = sorted(range(len(GroundAvgInfs)), key=lambda i: GroundAvgInfs[i])[-len(GroundAvgInfs):]
    Algo_rank = sorted(range(len(AlgoAvgInfs)), key=lambda i: AlgoAvgInfs[i])[-len(AlgoAvgInfs):]
    Pagerank_rank = sorted(range(len(PageRanks)), key=lambda i: PageRanks[i])[-len(PageRanks):]
    OutDeg_rank = sorted(range(len(OutDegree)), key=lambda i: OutDegree[i])[-len(OutDegree):]

    Ground_rank.reverse()
    Pagerank_rank.reverse()
    Algo_rank.reverse()
    OutDeg_rank.reverse()

    positions = list(range(453))

    dictionary = dict(zip(Ground_rank,positions))

    Ground_rankings = [dictionary[i] for i in Ground_rank]
    Page_rankings = [dictionary[i] for i in Pagerank_rank]
    Algo_rankings = [dictionary[i] for i in Algo_rank]
    OutDeg_rankings = [dictionary[i] for i in OutDeg_rank]

    '''
    Below is a full code of Kendall's tau just in case you want to know how it works
    concord = np.zeros(453)
    discord = np.zeros(453)

    for ii in Page_rankings:
        for jj in Page_rankings:
            if jj > ii:
                if Page_rankings[ii] < Page_rankings[jj]:
                    concord[ii] = concord[ii] + 1
                elif Page_rankings[ii] > Page_rankings[jj]:
                    discord[ii] = discord[ii] + 1
                
    concord_total = np.sum(concord)
    discord_total = np.sum(discord)
    
    kendall_tau = (concord_total-discord_total)/(concord_total+discord_total)
    '''
    AlgoTau, p_value = stats.kendalltau(Ground_rankings,Algo_rankings)
    PageTau, Pp_value = stats.kendalltau(Ground_rankings,Page_rankings)
    OutTau, Op_value = stats.kendalltau(Ground_rankings,OutDeg_rankings)
    
    return AlgoTau, PageTau, OutTau







