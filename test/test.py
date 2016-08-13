__author__ = 'SungJoonPark'

# import pandas as pd
#
# input_node_file = "Q:/COSSY+/source code/COSSY_source_code/data/networks/string/node-string.txt"
# input_edge_file = "Q:/COSSY+/source code/COSSY_source_code/data/networks/string/edge-string.txt"
# output_file = "Q:/COSSY+/data/network/string/string_stringid.sif"
#
# node_df =  pd.read_csv(input_node_file,sep='/t',header=None,names=['local_id','string_id'],index_col=0)['string_id']
#
# edge_df =  pd.read_csv(input_edge_file,header=None,sep='/t',names=['node1','node2'])
# edge_df['node1'] = node_df[edge_df['node1']].values
# edge_df['node2'] = node_df[edge_df['node2']].values
# edge_df.insert(1,'interaction',['pp']*len(edge_df))
#
# edge_df.to_csv(output_file,header=False,index=False,sep='/t')

import os
print os.path.splitext("Q:/COSSY+/tools/hotnet2/mutsigcv/data/IGCG_TCGA/BRCA/simple_somatic_mutation.open.BRCA-US.tsv")[0]