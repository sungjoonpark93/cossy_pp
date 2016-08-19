__author__ = 'SungJoonPark'

import os
import pandas as pd
import preprocess.icgc as icgc
import numpy as np



def ready_mutsigcv_for_generateHeat_from_mutsigcvoutput(mutsigcvoutput):
    #erase header, and filter Q value that is 0 to 1e-15
    df = pd.read_csv(mutsigcvoutput,sep="\t")
    df.loc[df['q']==0,'q'] = 1e-15
    df = df.fillna('NaN')
    outputfile = os.path.splitext(mutsigcvoutput)[0]+"_ready_for_generateHeat.txt"

    df.to_csv(outputfile,index=False,header=False,sep="\t")


def ready_mutation_for_generateHeat_from_icgc_mut(icgc_mutation_file):
    #make mutation file for input of generateHeat, we should prepare minimum cnafile and snv file
    def dict_to_file(dict_, outputfile):
        w = open(outputfile,'w')
        for sample in dict_.keys():
            w.write(sample)
            w.write("\t")
            for gene in dict_[sample]:
                w.write(gene)
                w.write("\t")
            w.write("\n")
        w.close()

    def icgc_to_cnafile(df,outputfile):
        df = df[np.logical_or( (df['mutation_type'] =='insertion of <=200bp') , (df['mutation_type'] =='deletion of <=200bp'))]

        cna_dict = {}
        for idx, row in df.iterrows():
            sample = row['submitted_sample_id']
            gene = row['Gene Symbol']
            
            if sample not in cna_dict:
                if row['mutation_type'] == 'insertion of <=200bp':
                    cna_dict[sample]=[gene+"(A)"]
                elif row['mutation_type'] == 'deletion of <=200bp':
                    cna_dict[sample]=[gene+"(D)"]
            else:
                if row['mutation_type'] == 'insertion of <=200bp':
                    cna_dict[sample].append(gene+"(A)")
                elif row['mutation_type'] == 'deletion of <=200bp':
                    cna_dict[sample].append(gene+"(D)")

        dict_to_file(cna_dict, outputfile)

    def icgc_to_snvfile(df,outputfile):
        df = df[df['mutation_type'] =='single base substitution']
        snv_dict = {}
        for idx, row in df.iterrows():
            sample = row['submitted_sample_id']
            gene = row['Gene Symbol']
            if sample not in snv_dict:
                snv_dict[sample]=[gene]
            else:
                snv_dict[sample].append(gene)

        dict_to_file(snv_dict, outputfile)

    df = icgc.get_gene_named_added_icgc_mut_df(icgc_mutation_file)

    snv_outputfile = os.path.splitext(icgc_mutation_file)[0]+"_for_generateHeat.snv"
    cna_outputfile = os.path.splitext(icgc_mutation_file)[0]+"_for_generateHeat.cna"
    icgc_to_snvfile(df,snv_outputfile)
    icgc_to_cnafile(df,cna_outputfile)


if __name__ == '__main__':
    icgc_mutation_file = "D:\hotnet2\hotnet2-master\data\ICGC_TCGA\BRCA\simple_somatic_mutation.open.BRCA-US.tsv"
    ready_mutation_for_generateHeat_from_icgc_mut(icgc_mutation_file)