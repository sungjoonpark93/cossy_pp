__author__ = 'SungJoonPark'

import pandas as pd
import numpy as np

def normalize_exp_data(exp_df, type='fuzzy'):
    if type=='fuzzy':
        for sample in exp_df.columns :
            column_data = exp_df[sample]
            top5percent_value = np.percentile(column_data, 95)
            low15percent_value = np.percentile(column_data, 15)
            exp_df[sample] = column_data.apply(lambda x:1 if x>top5percent_value else ((x-low15percent_value)/(top5percent_value-low15percent_value)) if (x<= top5percent_value and x>=low15percent_value) else 0)


    return exp_df




if __name__ =='__main__':
    tempData = pd.DataFrame([range(1,11), range(100,110) , [23,124,22,11,33,44,10,45,2,3]],index=['sample1','sample2','sample3'], columns=['gene']*10).T
    print tempData
    normalize_exp_data(tempData)

