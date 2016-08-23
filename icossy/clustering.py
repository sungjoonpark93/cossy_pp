'''
Created on 2016. 8. 22.

@author: sun
'''

from sklearn import datasets
from sklearn import metrics
from sklearn.cluster import KMeans
from sklearn import cluster

import pandas
from pandas.core.frame import DataFrame

if __name__ == "__main__":
    dataset = datasets.load_iris()
    
    x = dataset.data
    y = dataset.target
    
    
    print x, type(x)
    print y, type(y)
    
    for a in dataset:
        print a
        
    frame = DataFrame(data=x, columns=dataset.feature_names)
    
    print frame
    
    km = KMeans(n_clusters=3, random_state=1).fit(frame)
    
    labels = km.labels_
    score = metrics.silhouette_score(x, labels, metric='euclidean')
    
    print score
    print labels

    