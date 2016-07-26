__author__ = 'SungJoonPark'

import code.base.filter as filter
import pandas as pd
import os

def match(mut_preproceesed_filename=None,exp_preproceesed_filename=None,conserve_exp_normal_sample=True):
    '''
    match the sample and gene of mut_preprocessed_file, exp_preprocessed_file.
    Same gene and same patient with sample ordered prprocessed file is made
    output file in the input file directory with _patient_matched mark
    '''
    mut_df = pd.read_csv(mut_preproceesed_filename,index_col=0)
    exp_df = pd.read_csv(exp_preproceesed_filename,index_col=0)
    mut_df,exp_df = filter.match_gene(mut_preprocessed_df=mut_df,exp_preprocessed_df=exp_df)
    mut_df,exp_df = filter.match_sample(mut_df,exp_df,conserve_exp_normal_sample=conserve_exp_normal_sample)

    mut_output_filename = os.path.splitext(mut_preproceesed_filename)[0]+'_patient_matched.csv'
    exp_output_filename = os.path.splitext(exp_preproceesed_filename)[0]+'_patient_matched.csv'

    mut_df.to_csv(mut_output_filename)
    exp_df.to_csv(exp_output_filename)


if __name__ =='__main__':
    match(mut_preproceesed_filename="simple_somatic_mutation.open.BRCA-US.tsv_smoothing_string_preprocessed.csv",exp_preproceesed_filename="exp_seq.BRCA-US.tsv_preprocessed.csv")