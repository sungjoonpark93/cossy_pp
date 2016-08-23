__author__ = 'SungJoonPark'

import pandas as pd
import numpy as np
import os
import networkx as nx
import sys


def make_gmt_file(inputfile=None,outputfile=None):
    w = open(outputfile,'w')
    with open(inputfile,'r') as r:
        for i,line_string in enumerate(r):
            w.write(str(i)+"\t"+"clusterONE"+"\t"+line_string)
    w.close()

if __name__ =='__main__':
    outputfile = "Q:/COSSY+/source code/COSSY_source_code/data/keggWhole_clusterONE.gmt"
    inputfile = "Q:/COSSY+/tools/clusterONE/output\/keggWhole_clusterONE_minsize5.txt"
    make_gmt_file(inputfile=inputfile,outputfile=outputfile)