# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 11:35:52 2017

@author: gjgut
"""
import numpy as np
import random as rd
import networkx as nx
import matplotlib.pyplot as plt

Adlist = np.load('inList2.npy')
Weights = np.load('weights2.npy')
Nodes = len(Adlist)
Trials = 1
infnode = 27

    
Nstate = np.zeros(Nodes)
NxtState = np.zeros(Nodes)
Nstate[infnode] = 1
InfHist = np.zeros((Nodes + 1,Nodes))

cont_sim=1
t = 0
while cont_sim != 0:
    cont_sim=0
    InfHist[t,:] = Nstate #set row in history equal to the updated node state list
    for n in range(len(Nstate)): #cycle through node states
        if Nstate[n] == 1: #if infected
            cont_sim=1
            NxtState[n] = 2 #set infected node to recovered in next state list
            for ncon_out in range(len(Weights[n])): # cycle through infnode's connections
                if Nstate[Adlist[n][ncon_out]] == 0:
                    frand = rd.uniform(0,1) #roll the dice
                    if frand < Weights[n][ncon_out]: #if infection successful
                        NxtState[Adlist[n][ncon_out]] = 1 #set suscept node to infected in next state list
    Nstate = np.copy(NxtState) #set updated node state list as original node states for cycle
    t = t + 1

InfHist = InfHist[~(InfHist==0).all(1)] #delete rows after end of sim

tsteps,nodes = np.shape(InfHist)

G = nx.DiGraph() #make an empty graph
for ni in range(len(Adlist)): #go through connections in Adlist
    for cons in Adlist[ni]: #find those connections
        G.add_edges_from([(ni,cons)]) #add the edges to those connected nodes
pos = nx.random_layout(G) #set a layout
nx.draw_networkx_edges(G,pos,width=.1)

for step in range(tsteps): #go through timesteps in InfHist
    color_map = [] #make an empty list for the node colors
    for n in range(nodes): #check the nodes
        if InfHist[step][n] == 0: #if susceptible
            color_map.append('yellow') #make yeller
        elif InfHist[step][n] == 1: #if infected
            color_map.append('tomato') #make fruit
        elif InfHist[step][n] == 2: #if recovered
            color_map.append('turquoise') #make turtle
    nx.draw_networkx_nodes(G,pos,node_size=15,node_color=color_map)
    plt.axis('off')
    plt.title('t = %d'%step)
    plt.savefig('1000NodeGraph2/graph_t=%d.png'%step)