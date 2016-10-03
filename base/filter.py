__author__ = 'SungJoonPark'
import pandas as pd
import numpy as np

#Before making gct file, need to order column as tumor sample and normal sample order, alo return dividing index.
def order_sample(df):
    normal_sample_list = [sample for sample in df.columns if int(sample[13:15])>=10]
    tumor_sample_list = [sample for sample in df.columns if int(sample[13:15])<10]
    #divide tumor sample and normal sampel and sort for each barcode order
    ordered_sample_list = sorted(tumor_sample_list)+sorted(normal_sample_list)
    ordered_df = df.loc[:,ordered_sample_list]
    tumor_normal_divide_index = len(tumor_sample_list)
    return ordered_df,tumor_normal_divide_index

def get_normal_sample_list(preprocessed_df):
    normal_sample_list = [sample for sample in preprocessed_df.columns if int(sample[13:15]) >= 10]
    return normal_sample_list

def get_tumor_sample_list(preprocessed_df):
    tumor_sample_list = [sample for sample in preprocessed_df.columns if int(sample[13:15]) < 10]
    return tumor_sample_list

def get_sample_list(preprocessed_df):
    sample_list =[sample for sample in preprocessed_df.columns]
    return sample_list

def get_sample_classes_and_lables(preprocessed_df):
    #conver sample list into classes and labels in which classes is [1,0,0] and lables is [tumor, normal]
    #tumor is 1 , normal is 0
    sample_list = get_sample_list(preprocessed_df)
    sample_label_list = ['normal' if int(sample[13:15])>=10 else "tumor" for sample in sample_list]
    labels = ('normal','tumor')

    classes = []
    for label in sample_label_list:
        if label == 'tumor':
            classes.append(1)
        elif label =='normal':
            classes.append(0)

    return {'classes':classes,'labels':labels}


def replace_normal_sample(preprocessed_df,normal_sample_list):
    '''
    replace the normal sample column of preprocessed_df with input's normal_sample_list having value 0
    '''
    preprocessed_df = preprocessed_df.drop(get_normal_sample_list(preprocessed_df),axis=1)
    data=np.zeros((len(preprocessed_df.index),len(normal_sample_list)))
    normal_sample_df = pd.DataFrame(data, index=preprocessed_df.index, columns =normal_sample_list)
    preprocessed_df = pd.concat([preprocessed_df,normal_sample_df],axis=1)
    return preprocessed_df


def add_exp_only_gene_to_mut(mut_preprocessed_df = None, exp_preprocessed_df=None):
    '''
    We have to add exp_only_gene mut having value 0 meaning no somatic mutation
    '''
    exp_only_gene = list(set(exp_preprocessed_df.index) -set(mut_preprocessed_df.index))
    exp_only_gene_df = pd.DataFrame(np.zeros((len(exp_only_gene),len(mut_preprocessed_df.columns))),index=exp_only_gene,columns=mut_preprocessed_df.columns)
    mut_preprocessed_df = pd.concat([mut_preprocessed_df,exp_only_gene_df],axis=0)
    return mut_preprocessed_df


def match_gene(mut_preprocessed_df=None , exp_preprocessed_df = None, conserve_exp_gene = True):
    '''
    conserve_exp_gene means add exp_only to mut having value 0
    '''
    if conserve_exp_gene == True:
        mut_preprocessed_df = add_exp_only_gene_to_mut(mut_preprocessed_df=mut_preprocessed_df, exp_preprocessed_df=exp_preprocessed_df)
    common_gene_list = list(mut_preprocessed_df.index.intersection(exp_preprocessed_df.index))
    return mut_preprocessed_df.loc[common_gene_list,:], exp_preprocessed_df.loc[common_gene_list,:]


def match_gene_with_network_data(mut_preprocessed_df=None,network_df=None, is_exp = False, exp_preprocessed_df = None, conserve_exp_gene=True):

    '''
    conserve_exp_gene flag is meaningful only when is_exp is true
    '''
    def get_overlapped_gene_list(mut_preprocessed_df=None, network_df=None, is_exp=False, exp_preprocessed_df=None):
        if is_exp==False:
            return list(mut_preprocessed_df.index.intersection(network_df.index))
        elif is_exp==True:
            patient_genes = mut_preprocessed_df.index.intersection(exp_preprocessed_df.index)
            return list(patient_genes.intersection(network_df.index))

    if is_exp==False:
        common_gene_list= get_overlapped_gene_list(mut_preprocessed_df=mut_preprocessed_df,network_df=network_df,is_exp=False)
        return mut_preprocessed_df.loc[common_gene_list,:],network_df.loc[common_gene_list,common_gene_list]

    elif is_exp==True:
        if conserve_exp_gene == True:
            mut_preprocessed_df = add_exp_only_gene_to_mut(mut_preprocessed_df=mut_preprocessed_df,exp_preprocessed_df=exp_preprocessed_df)
        common_gene_list = get_overlapped_gene_list(mut_preprocessed_df=mut_preprocessed_df,network_df=network_df,is_exp=True,exp_preprocessed_df=exp_preprocessed_df)


        return mut_preprocessed_df.loc[common_gene_list,:],network_df.loc[common_gene_list,common_gene_list],exp_preprocessed_df.loc[common_gene_list,:]



def match_sample(mut_preprocessed_df=None,exp_preprocessed_df=None, conserve_exp_normal_sample=False):
    '''
    :param mut_preprocessed_df:
    :param exp_preprocessed_df:
    :param conserve_exp_normal_sample: True if you want to conserve the normal sample of exp data.
    :return:
    '''

    #match the sample between mut_preprocessed_df, and exp_preprocessed_df, and make the order of sample same for each df.
    #make dict that key is sample and value is one sample barcode having the sample id.
    #there could be multiple sample barcode mapped into one sample id, we just asign one barcode to the sample id

    if conserve_exp_normal_sample ==True:
        mut_preprocessed_df = replace_normal_sample(mut_preprocessed_df,get_normal_sample_list(exp_preprocessed_df))

    mut_sample_barcode_dict = {}
    for barcode in mut_preprocessed_df.columns.tolist():
        sample = barcode[0:16]
        if sample not in mut_sample_barcode_dict:
            mut_sample_barcode_dict[sample]=barcode


    exp_sample_barcode_dict = {}
    for barcode in exp_preprocessed_df.columns.tolist():
        sample = barcode[0:16]
        if sample not in exp_sample_barcode_dict:
            exp_sample_barcode_dict[sample]=barcode

    overlap_sample_list = list(set(mut_sample_barcode_dict.keys()) & set(exp_sample_barcode_dict.keys()))
    mut_matched_barcode_list = pd.Series(mut_sample_barcode_dict)[overlap_sample_list].sort_values().values.tolist()
    exp_matched_barcode_list = pd.Series(exp_sample_barcode_dict)[overlap_sample_list].sort_values().values.tolist()

    sample_matched_mut_preprocessed_df=mut_preprocessed_df.loc[:,mut_matched_barcode_list]
    sample_matched_exp_preprocessed_df=exp_preprocessed_df.loc[:,exp_matched_barcode_list]

    #we have to make the sample order same for mut and exp
    sample_matched_mut_preprocessed_df,_ = order_sample(sample_matched_mut_preprocessed_df)
    sample_matched_exp_preprocessed_df,_ = order_sample(sample_matched_exp_preprocessed_df)

    return sample_matched_mut_preprocessed_df,sample_matched_exp_preprocessed_df


