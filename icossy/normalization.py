__author__ = 'SungJoonPark'

import pandas as pd

def normalize_exp_data(exp_df, type='fuzzy'):
    if type=='fuzzy':
        total_gene_num = len(exp_df.index)
        for sample in exp_df.columns :
            sorted(exp_df[sample],reverse=True)

            print exp_df[sample]



if __name__ =='__main__':
    tempData = pd.DataFrame([range(1,11), range(100,110)],index=['sample1','sample2'], columns=['gene']*10).T
    normalize_exp_data(tempData)
    print tempData