__author__ = 'SungJoonPark'

import base.baseVariables as basevar
import pandas as pd
import os
import network.preprocess as network_preprocess
import preprocess.smoothing_for_icossy as smoothing
import base.filter as basefilter


def load_gmt_file(gmt_file):
    print "loading gmt file"
    #gmt file name should like "kegg.gmt". In this case kegg is used as a part of key of mis dit
    #read gmt file and make mis_dict
    with open(gmt_file,'r') as r:
        filename = os.path.split(gmt_file)[1]
        mis_source = filename.split('.')[0]
        mis_dict = {"mis_"+mis_source+"_"+str(i):line.strip().split('\t')[2:] for i,line in enumerate(r)}
    return mis_dict


def load_exp_data(gct_file,type='gct'):
    if type=='gct':
        #exp file is gct_format
        print "loading gct file"
        #read gct file and make as gene x patient dataframe
        exp_df = pd.read_csv(gct_file, sep='\t',skiprows=2).drop('Description',axis=1).set_index('Name')


    return exp_df


def load_cls_file(cls_file):
    print "loading cls file"
    with open(cls_file,'r') as r:
        label_list = [line.strip().split(' ') for i,line in enumerate(r) if i==2][0]
    return label_list


def load_class_data(class_realted_file=None,type='cls'):
    if type=='cls':
        label_list = load_cls_file(class_realted_file)
    elif type == 'preprocessed':
        label_list = basefilter.get_sample_list(class_realted_file)
    #there are only two labels
    labels = list(set(label_list))
    classes = [{labels[0]:1,labels[1]:0}[label] for label in label_list]
    return classes,labels



def load_mutation_data(mutation_file):
    print "loading mutation file"
    #mutation_file should be gene x patient matrix in which each element is binary enconded value of mutation information. e.g. 1 is mutated, 0 is not
    mut_df=pd.read_csv(mutation_file,index_col=0)
    return mut_df



def load_data(exp_file=None,  mutation_file=None,  gmt_file=None, network_file_for_smoothing=None ,type=None):
    #exp file is gct file, mutation_file is preprocessed file

    print "loading data"

    mis_dict = load_gmt_file(gmt_file)

    if type=='expression':
        exp_df = load_exp_data(exp_file, type='gct')
        profile = get_profile(exp_df=exp_df)
        profileData = {'profile':profile, 'classes': basefilter.get_sample_classes_and_lables(profile)['classes'], 'labels': basefilter.get_sample_classes_and_lables(profile)['labels']}
        return {'profileData':profileData ,'misList':mis_dict}

    elif type =='mutation':
        mut_df = load_mutation_data(mutation_file)
        network_df_for_smoothing = network_preprocess.get_network(network_file_for_smoothing)
        profile = get_profile(mut_df = mut_df ,network_df_for_smoothing=network_df_for_smoothing,type='mut')
        profileData = {'profile':profile, 'classes': basefilter.get_sample_classes_and_lables(profile)['classes'], 'labels': basefilter.get_sample_classes_and_lables(profile)['labels']}
        return {'profileData':profileData ,'misList':mis_dict}

    elif type =='mut_with_exp':
        mut_df = load_mutation_data(mutation_file)
        exp_df = load_exp_data(exp_file, type='gct')
        network_df_for_smoothing = network_preprocess.get_network(network_file_for_smoothing)
        profile = get_profile(mut_df=mut_df , exp_df=exp_df , network_df_for_smoothing=network_df_for_smoothing , type='mut_with_exp')
        profileData = {'profile':profile, 'classes': basefilter.get_sample_classes_and_lables(profile)['classes'], 'labels': basefilter.get_sample_classes_and_lables(profile)['labels']}
        return {'profileData':profileData ,'misList':mis_dict}

    else:
        print "type is not specified "


def get_profile(mut_df=None, exp_df=None, network_df_for_smoothing = None, type='exp'):
    if type=='exp':
        profile = exp_df
    elif type=='mut':
        #smoothing is needed
        profile = smoothing.smoothing(mut_df=mut_df, network_df=network_df_for_smoothing,type='mut')
    elif type== 'mut_with_exp':
        profile = smoothing.smoothing(mut_df=mut_df,exp_df=exp_df,network_df=network_df_for_smoothing,type='mut_with_exp',conserve_exp_normal_sample=True)
    return profile


if __name__ =='__main__':
    temp_gct_file = basevar.exp_gct_data['BRCA']
    temp_cls_file = basevar.exp_cls_data['BRCA']
    temp_gmt_file = basevar.gmt['kegg']
    temp_mutation_file = basevar.mut_data['BRCA']
    temp_network_for_smoothing_file = basevar.network['kegg']
    print load_data(mutation_file=temp_mutation_file , gmt_file=temp_gmt_file ,network_file_for_smoothing=temp_network_for_smoothing_file,type='mutation')