__author__ = 'SungJoonPark'

import base.baseVariables as basevar
import pandas as pd
import os

def load_gmt_file(gmt_file):
    print "loading gmt file"
    #gmt file name should like "kegg.gmt". In this case kegg is used as a part of key of mis dit
    #read gmt file and make mis_dict
    with open(gmt_file,'r') as r:
        filename = os.path.split(gmt_file)[1]
        mis_source = filename.split('.')[0]
        mis_dict = {mis_source+str(i):line.strip().split('\t')[2:] for i,line in enumerate(r)}
    return mis_dict


def load_gct_file(gct_file):
    print "loading gct file"
    #read gct file and make as gene x patient dataframe
    exp_df = pd.read_csv(gct_file, sep='\t',skiprows=2).drop('Description',axis=1).set_index('Name')
    return exp_df


def load_cls_file(cls_file):
    print "loading cls file"
    with open(cls_file,'r') as r:
        label_list = [line.strip().split(' ') for i,line in enumerate(r) if i==2][0]
    return label_list


def load_mutation_file(mutation_file):
    print "loading mutation file"
    #mutation_file should be gene x patient matrix in which each element is binary enconded value of mutation information. e.g. 1 is mutated, 0 is not
    mut_df=pd.read_csv(mutation_file,index_col=0)
    return mut_df


def load_data(gct_file=None, cls_file=None, mutation_file=None, gmt_file=None):
    print "loading the gct,cls,mutation,gmt data"
    exp_df = load_gct_file(gct_file)
    label_list = load_cls_file(cls_file)
    mut_df = load_mutation_file(mutation_file)
    mis_dict = load_gmt_file(gmt_file)

    return exp_df, label_list, mut_df,mis_dict

if __name__ =='__main__':
    temp_gct_file = basevar.exp_gct_data['BRCA']
    temp_cls_file = basevar.exp_cls_data['BRCA']
    temp_gmt_file = basevar.gmt['kegg']
    temp_mutation_file = basevar.mut_data['BRCA']
    exp_df, label_list, mut_df , mis_dict = load_data(temp_gct_file,temp_cls_file,temp_mutation_file,temp_gmt_file)
