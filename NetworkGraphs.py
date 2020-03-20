# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 09:26:46 2017

@author: gjgut
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from AdlistCreate import NetCons
from InfectSim import Sim

AdjM = np.array([[0,0,0,0,0,0],[0.5,0,0,0,0,0],[0.4,0.3,0,0,0,0],[0,0.9,0,0,0,0],[0,0.2,0.1,0,0,0],[0,0,0,0.8,0.7,0]])

#%%
def SIR(AdjM,infnode,seed):

    Adlist, Weights, ncols = NetCons(AdjM)    
    InfHist = Sim(AdjM,infnode,seed)
    tsteps,nodes = np.shape(InfHist)

    G = nx.DiGraph() #make an empty graph
    for ni in range(len(Adlist)): #go through connections in Adlist
        for cons in Adlist[ni]: #find those connections
            G.add_edges_from([(ni,cons)]) #add the edges to those connected nodes
    pos = nx.circular_layout(G) #set a layout
    nx.draw_networkx_edges(G,pos)
    nx.draw_networkx_labels(G,pos)

    for step in range(tsteps): #go through timesteps in InfHist
        color_map = [] #make an empty list for the node colors
        for n in range(nodes): #check the nodes
            if InfHist[step][n] == 0: #if susceptible
                color_map.append('yellow') #make yeller
            elif InfHist[step][n] == 1: #if infected
                color_map.append('tomato') #make fruit
            elif InfHist[step][n] == 2: #if recovered
                color_map.append('turquoise') #make turtle
        nx.draw_networkx_nodes(G,pos,node_size=300,node_color=color_map)
        plt.axis('off')
        plt.title('t = %d'%step)
        plt.savefig('SavedGraphs/graph_t=%d.png'%step)
    
    return G
    
    
        
        
        