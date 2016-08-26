__author__ = 'SungJoonPark'
import pandas as pd
import numpy as np
import time
import base.filter as filter

def smoothing(exp_df=None, mut_df =None, network_df = None, type='mut',alpha=0.7, conserve_exp_normal_sample=True):
    #alpha is restart parameter for smoothing
    #conserve_exp_normal_sampe is only for type mut_with_exp


    #mutation
    if type=='mut':
        print("you are going to smoothing mutation data")
        s_time = time.time()
        n_iter = 0

       #reverse_mut_df 1 to 0 0 to 1
        mutation_value_reversed_mut_df = mut_df.replace([0, 1], [1, 0])
        mut_preprocessed_df, network_df = filter.match_gene_with_network_data(mut_preprocessed_df=mutation_value_reversed_mut_df, network_df=network_df)
        print str(len(mut_preprocessed_df.columns)) + " genes are matched with network data"
        patient_mutation_df = mut_preprocessed_df.T
        F_0 = patient_mutation_df.copy(deep=True)
        F = patient_mutation_df.copy(deep=True)

        while True:
            n_iter = n_iter + 1

            F_t = F.copy(deep=True)
            F = (alpha * np.matmul(F, network_df)) + ((1 - alpha) * F_0)

            if np.sum(np.asmatrix(np.square((F - F_t)))) < (1 / 1000000.0):
                print "network smoothing end. ", str(n_iter), "iterations, " , str((time.time() - s_time) / 60.0) , " minutes elapsed\n"
                break

        return F.T

    #mut_with_exp
    elif type=='mut_with_exp':
        s_time = time.time()
        n_iter = 0

        #reverse_mut_df 1 to 0 0 to 1
        mutation_value_reversed_mut_df = mut_df.replace([0, 1], [1, 0])
        mut_preprocessed_df, exp_preprocessed_df = filter.match_sample(mut_preprocessed_df=mutation_value_reversed_mut_df, exp_preprocessed_df=exp_df, conserve_exp_normal_sample=conserve_exp_normal_sample)
        mut_preprocessed_df, network_df, exp_preprocessed_df = filter.match_gene_with_network_data(mut_preprocessed_df=mut_preprocessed_df, network_df=network_df, is_exp=True, exp_preprocessed_df=exp_preprocessed_df)

        patient_mutation_df = mut_preprocessed_df.T
        patient_exp_df = exp_preprocessed_df.T


        if (patient_mutation_df.shape != patient_exp_df.shape):
            raise Exception("mutation and exp data shape are not matched")

        F_0 = patient_mutation_df.copy(deep=True)
        F = patient_mutation_df.copy(deep=True)
        while True:
            n_iter = n_iter + 1

            F_t = F.copy(deep=True)

            # get max by index but, when dividing operation meet, it divides by column so one method is to transpose.
            F = np.multiply(((alpha * np.matmul(F, network_df)) + ((1 - alpha) * F_0)), (patient_exp_df.T / np.max(np.abs(patient_exp_df), axis=1)).T)


            if np.sum(np.asmatrix(np.square((F - F_t)))) < (1 / 1000000.0):
                print "network smoothing end. ", str(n_iter), "iterations, " , str((time.time() - s_time) / 60.0) , " minutes elapsed\n"
                break

        return F.T