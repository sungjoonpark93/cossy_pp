__author__ = 'SungJoonPark'

import pandas as pd

root_dir = "Q:\COSSY+\data\mis\icossy/"
mapping_dir = "Q:\COSSY+\data\mapping_file/"

string_mapping_file = mapping_dir+"HGNC_Entrez_STRING_fix_excel_problem.csv"
kegg_mapping_file = mapping_dir+"HGNC_Entrez_fix_excel_problem.csv"


string_dict = pd.read_csv(string_mapping_file,index_col=1)['HGNC_Symbol'].to_dict()
kegg_dict = pd.read_csv(kegg_mapping_file,index_col=5)['Approved Symbol'].to_dict()


kegg_cossy = root_dir+"string_cossy.gmt"
w = open(root_dir+"string_cossy_symbol.gmt",'w')
with open(kegg_cossy,'r') as r:
    for line in r:
        line_list =  line.strip().split("\t")
        w.write(line_list[0]+"\t"+line_list[1]+"\t")
        for i,element in enumerate(line_list):
            if i>=2:
                if  element in string_dict:
                    w.write(string_dict[element])
                    w.write("\t")
        w.write("\n")

w.close()
