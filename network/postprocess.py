__author__ = 'SungJoonPark'

import pandas as pd
import numpy as np
import os
import networkx as nx

def make_network_from_geneset(network_file, geneset_file, thr=100):
    df = pd.read_csv(network_file,sep='\t',header=None, names =['Source','Target','Weight'])
    G = nx.DiGraph()
    G.add_weighted_edges_from([tuple(x) for x in df.values])

    #construct subgraph
    subgraph_dict = {}
    with open(geneset_file,'r') as r:
        for line in r:
            cluster_id = int(line.strip().split("\t")[0])
            node_list = map(lambda x : int(x),line.strip().split("\t")[2:])
            #reconstruct only geneset size is more than 100
            if len(node_list)>thr:
                SG = G.subgraph(node_list)
                SG_edge_list = SG.edges(data=True)
                subgraph_dict[cluster_id] = SG_edge_list

    #output the result
    for key in subgraph_dict.keys():
        output_file = os.path.split(geneset_file)[0]+"/"+"cluster_"+str(key)+".txt"
        w = open(output_file,'w')
        for edge_list in subgraph_dict[key]:
            source = int(edge_list[0])
            target = int(edge_list[1])
            weight = edge_list[2]['weight']
            w.write(str(source)+"\t"+str(target)+"\t"+str(weight)+"\n")
        w.close()