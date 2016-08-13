__author__ = 'SungJoonPark'

import pandas as pd
import os

def make_maf_for_mutsigcv(inputfile):
    df = pd.read_csv(inputfile,sep='\t',low_memory=False)
    #there is no gene symbol in ICGC mutation file, only Ensemble ID, so wee need to make gene symbol from EnsembleID
    mapping_file = "Q:/COSSY+/data/mapping_file/HGNC_ApprovedSymbol_EnsembleID.txt"
    map_df = pd.read_csv(mapping_file,sep='\t')
    map_dict = map_df.set_index('Ensembl ID(supplied by Ensembl)')['Approved Symbol'].to_dict()
    #filter the ensemble id into the one in the mapping file.
    df = df[df['gene_affected'].isin(map_dict.keys())]
    df['Gene Symbol'] = [map_dict[ensemble_id] for ensemble_id in df['gene_affected']]
    df= df[['Gene Symbol','submitted_sample_id','consequence_type']]

    #select only below type of mutation type
    df = df[df['consequence_type'].isin(['disruptive_inframe_deletion','disruptive_inframe_insertion','frameshift_variant','inframe_deletion','inframe_insertion','missense_variant','stop_gained','exon_variant'])]
    df = df.replace(['disruptive_inframe_deletion','disruptive_inframe_insertion','frameshift_variant','inframe_deletion','inframe_insertion','missense_variant','stop_gained','exon_variant'],
               ['In_Frame_Del','In_Frame_Ins','Frame_Shift_Del','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','RNA'])


    df.columns = [['gene','patient','Variant_Classification']]

    outputfile = os.path.splitext(inputfile)[0]+"_for_mutsigcv.maf"
    df.to_csv(outputfile,sep='\t',index=False)
    return outputfile


if __name__== '__main__':
    #dataset_name = ['BRCA']
    dataset_name = ['BRCA','COAD','PRAD','LUSC','STAD']
    rootdir = "Q:/COSSY+/tools/hotnet2/mutsigcv/data/maf file/IGCG_TCGA"
    for one_dataset in dataset_name:
        inputfile = rootdir+"/"+one_dataset+"/simple_somatic_mutation.open."+one_dataset+"-US.tsv"
        make_maf_for_mutsigcv(inputfile)