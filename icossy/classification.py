__author__ = 'SungJoonPark'
from sklearn import datasets
from scipy.stats import ttest_ind  # @UnresolvedImport

from pandas.core.frame import DataFrame
from operator import itemgetter
from sklearn.cluster.hierarchical import AgglomerativeClustering

import numpy as np
from scipy.spatial.distance import euclidean
import pandas

def classify(data, clustering, dist="euclidean"):

    if type(data) != "":
        print "data should be type of DataFrame!"
        print "please check the type"

    classification = {}

    if dist == "euclidean":
        dist = euclidean

    classificationList = {}


    print
    print
    for pid in data:
        classificationList[pid] = {}
        onerec = data[pid]
        for misid in clustering:
            clusters = clustering[misid]['c']
            repgenes = clustering[misid]['repGenes']

            datasubset = onerec.ix[repgenes].T

            distances = [(cid, clusters[cid]['label'], dist(datasubset, clusters[cid]['centroid'])) for cid in clusters]

            closestcid = min(distances, key=itemgetter(2))

            classificationList[pid][misid] = closestcid

    return classificationList

def fit(trainTopkClustering_result, testData):

    #data_to_be_predicted is dataframe which is gene x sample dataframe
    fittingResult = {}

    testData_profile = testData['profileData']['profile']
    classification = classify(testData_profile, trainTopkClustering_result)

    for pid in classification:
        patient = classification[pid]
        cls = {0:0, 1:0}
        for misid in patient:
            cls[ int(patient[misid][1] + 0.5) ] += 1

        if cls[0] > cls[1]:
            fittingResult[pid] = 0
        else:
            fittingResult[pid] = 1

    return fittingResult


