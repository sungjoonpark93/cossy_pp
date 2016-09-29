__author__ = 'SungJoonPark'
#make dataset file (exp, mut) 10folds
import pandas as pd
from sklearn.cross_validation import KFold
import os
import base.baseVariables as data

def make_nfold_dataset(datasetfile, n_Folds = 10):
    #dataset file is preprocessed file
    basename = os.path.basename(os.path.splitext(datasetfile)[0])
    outputdir = os.path.dirname(datasetfile)+"/"+ basename+"_10folds"
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    df = pd.read_csv(datasetfile,index_col=0)
    kf = KFold(len(df.columns), n_folds=n_Folds,shuffle=True)
    fold_num = 1
    for train_index , test_index in kf:
        train_outputfile= outputdir + "/"+ basename+"_train"+str(fold_num)+".csv"
        test_outputfile =  outputdir +"/"+basename +"_test"+str(fold_num) +".csv"
        train_df , test_df = df.iloc[:,train_index], df.iloc[:,test_index]
        train_df.to_csv(train_outputfile)
        test_df.to_csv(test_outputfile)
        fold_num = fold_num +1



if __name__ =='__main__':
    for dataset in  data.dataset_names:
        print data.exp_data[dataset]
        make_nfold_dataset(data.exp_data[dataset])
        print data.mut_data[dataset]
        make_nfold_dataset(data.mut_data[dataset])