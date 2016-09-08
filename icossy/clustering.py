'''
Created on 2016. 8. 22.

@author: sun
'''

from sklearn import datasets
from scipy.stats import ttest_ind  # @UnresolvedImport

from pandas.core.frame import DataFrame
from operator import itemgetter
from sklearn.cluster.hierarchical import AgglomerativeClustering

import numpy as np
from scipy.spatial.distance import euclidean
import pandas


def ttestRepFunc(data):
    
    tf = data["profile"]
    cls = data["classes"]
    
    repScoreList = []
    
    for gene in tf.index:
        
        orilist = tf.ix[gene]
        
        list1 = [x for x, c in zip(orilist, cls) if c == 1]
        list2 = [x for x, c in zip(orilist, cls) if c == 0]
        
        tv, pv = ttest_ind(list1, list2)
        
        repScoreList.append((gene, pv))
    
    return sorted(repScoreList, key=itemgetter(1))

def representativeScore(data, scoreFunc="ttest"):
    if scoreFunc == "ttest":
        scoreFunc = ttestRepFunc
    
    repGeneList = scoreFunc(data)
    
    return [x[0] for x in repGeneList]


def clusteringInMIS(data, clusterNum, topKgene, misList, method="ttest"):
    #print "start clustering"

    result = {}
    
    clustering = AgglomerativeClustering(n_clusters=clusterNum)
    
    ranked = representativeScore(data, method)
    
    for misid in misList:
        misgene = misList[misid]
        
        '''
        representativeGenes = []
        
        for i in range(topKgene):
            representativeGenes.append(misgene[i])
        '''
        
        representativeGenes = [x for x in ranked if x in misgene][0:topKgene]
        
        if len(representativeGenes) < topKgene:
            #raise Exception("not enough genes")
            continue
        
        datasubset = data["profile"].ix[representativeGenes].T
        clusteringResult = clustering.fit(X=datasubset)
        clusterLabels = clusteringResult.labels_
        
        tripleList = zip(clusterLabels, data["classes"], data["profile"].columns.values.tolist())
        
        
        result[misid] = {'c':{cid:{"positive":[], "negative":[]} for cid in set(clusterLabels)}, 'repGenes':representativeGenes}
        
        for triple in tripleList:
            
            if triple[1] == 1:
                result[misid]['c'][triple[0]]["positive"].append(triple[2])
            else: 
                result[misid]['c'][triple[0]]["negative"].append(triple[2])
        
        for cid in result[misid]['c']:
            cluster = result[misid]['c'][cid]
            patients = cluster['positive'] + cluster['negative']
            
            centroid = np.array(np.mean(datasubset.T[patients].T))
            
            result[misid]['c'][cid]['centroid'] = centroid
            result[misid]['c'][cid]['label'] = float(len(result[misid]['c'][cid]['positive'])) / float(len(result[misid]['c'][cid]['positive']) + len(result[misid]['c'][cid]['negative']))
        '''
        result.append( {misid:[{"positive":[x[2] for x in tripleList if x[0] == clab and x[1] == 1], 
                                "negative":[x[2] for x in tripleList if x[0] == clab and x[1] == 0] } for clab in set(clusterLabels)] } )
        '''
                
    return result

# def classify(data, clustering, dist="euclidean"):
#
#     if type(data) != "":
#         print "data should be type of DataFrame!"
#         print "please check the type"
#
#     classification = {}
#
#     if dist == "euclidean":
#         dist = euclidean
#
#     classificationList = {}
#
#
#     print
#     print
#     for pid in data:
#         classificationList[pid] = {}
#         onerec = data[pid]
#         for misid in clustering:
#             clusters = clustering[misid]['c']
#             repgenes = clustering[misid]['repGenes']
#
#             datasubset = onerec.ix[repgenes].T
#
#             distances = [(cid, clusters[cid]['label'], dist(datasubset, clusters[cid]['centroid'])) for cid in clusters]
#
#             closestcid = min(distances, key=itemgetter(2))
#
#             classificationList[pid][misid] = closestcid
#
#     return classificationList

if __name__ == "__main__":
    
    import pprint
    
    tmpframe = [ [0.9, 0.9, 0.9, 0.8, 0.9, 0.8, 0.9, 0.0, 0.0, 0.1, 0.0, 0.1, 0.2 ],
                 [0.8, 0.7, 0.6, 0.6, 0.7, 0.4, 0.5, 0.3, 0.1, 0.2, 0.1, 0.1, 0.2 ],
                 [0.4, 0.3, 0.4, 0.4, 0.5, 0.2, 0.5, 0.5, 0.6, 0.7, 0.5, 0.6, 0.6 ],
                 [0.3, 0.5, 0.3, 0.4, 0.3, 0.3, 0.2, 0.3, 0.1, 0.1, 0.4, 0.1, 0.5 ],
                 [0.7, 0.5, 0.6, 0.6, 0.7, 0.8, 0.6, 0.2, 0.1, 0.2, 0.2, 0.6, 0.1 ],
                 [0.4, 0.1, 0.2, 0.1, 0.3, 0.2, 0.1, 0.7, 0.6, 0.7, 0.5, 0.6, 0.4 ],
                 [0.5, 0.3, 0.3, 0.4, 0.3, 0.1, 0.3, 0.2, 0.1, 0.0, 0.0, 0.1, 0.1 ],
                 [0.1, 0.1, 0.2, 0.1, 0.3, 0.1, 0.2, 0.4, 0.3, 0.6, 0.8, 0.6, 0.4 ],
                 [0.7, 0.7, 0.4, 0.7, 0.6, 0.5, 0.7, 0.4, 0.3, 0.5, 0.3, 0.5, 0.4 ],
                 [0.9, 0.8, 0.8, 0.9, 0.7, 0.8, 0.9, 0.7, 0.5, 0.6, 0.7, 0.6, 0.5 ]]
    tmppatients = ["p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10", "p11", "p12", "p13"]
    tmpgenes = ["g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10"]
    
    profileData = {"profile" : DataFrame(data=tmpframe, index=tmpgenes, columns=tmppatients), "classes": [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0], "labels":("tumor", "normal")}

    misList = {"mis1":["g1", "g3", "g4"], "mis2":["g2", "g5", "g8", "g9"], "mis3":["g6", "g7", "g10"]}
    
    pp = pprint.PrettyPrinter()
    
    repgenes = representativeScore(profileData)
    print repgenes
    
    print profileData["profile"].ix[misList["mis1"]].T
    
    me = np.mean(profileData["profile"].ix[misList["mis1"]].T)
    
    print me, type(me) 
    
    clusterLabs = [1,3,1,2,2,2,2,3,1,3]
    
    result = clusteringInMIS(data=profileData, clusterNum=3, topKgene=3, misList=misList)

    pp.pprint(result)
    
    test = profileData['profile'][['p1','p11']]
    print test
    
#    classes = classify(data=test, clustering=result)
    
    #pp.pprint(classes)
    
    import random
    
    f=5
    
    folds = [[] for x in range(f)]
    profile = profileData["profile"]
    classes = profileData["classes"]
        
    pairlist = [ (classes[x], profile.columns[x]) for x in range(len(classes))]
    random.shuffle(pairlist)

    pospairs = enumerate([ x for x in pairlist if x[0] == 1])
    negpairs = enumerate([ x for x in pairlist if x[0] == 0])
    
    
    for idx, v in pospairs:
        i = idx%f
        folds[i].append(v)
    
    for idx, v in negpairs:
        i = f - idx%f -1
        folds[i].append(v)
    
    foldedData = []
    
    for fold in folds:
        pids = [x[1] for x in fold]
        classes = [x[0] for x in fold]
        
        profileSubset = profile[pids]
        
        foldedData.append({"profile":profileSubset, "classes":classes, "labels":profileData["labels"]})
    
    print
    pp.pprint(pairlist)
    print
    pp.pprint(folds)
    print " -------------------- "
    pp.pprint(foldedData)
    
    print 
    print
    
    idxs = [0,1,2]
#    new = {"profile":pandas.concat([foldedData[x]["profile"] for x in idxs], axis=1), "classes": , "labels":foldedData[idx[0]]["labels"]}
    
#    pp.pprint(new)
