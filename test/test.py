__author__ = 'SungJoonPark'

import pandas as pd

input_node_file = "Q:/COSSY+/source code/COSSY_source_code/data/networks/string/node-string.txt"
input_edge_file = "Q:/COSSY+/source code/COSSY_source_code/data/networks/string/edge-string.txt"
output_file = "Q:/COSSY+/data/network/string/string_stringid.sif"

print pd.read_csv(input_node_file)