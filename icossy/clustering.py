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
    
    result = []
    
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
        
        datasubset = data["profile"].ix[representativeGenes].T
        clusteringResult = clustering.fit(X=datasubset)
        clusterLabels = clusteringResult.labels_
        
        
        tripleList = zip(clusterLabels, data["classes"], data["profile"].feature_names)
        
        
        result = {cid:{"positive":[], "negative":[]} for cid in set(clusterLabs)}
    
        for triple in tripleList:
            if triple[1] == 1:
                result[triple[0]]["positive"].append(triple[2])
            else: 
                result[triple[0]]["negative"].append(triple[2])
        
        for cid in result:
            cluster = result[cid]
            patients = cluster['positive'] + cluster['negative']
            centroid = np.array(np.mean(data["profile"][patients].T))
            
            result[cid]['centroid'] = centroid
        '''
        result.append( {misid:[{"positive":[x[2] for x in tripleList if x[0] == clab and x[1] == 1], 
                                "negative":[x[2] for x in tripleList if x[0] == clab and x[1] == 0] } for clab in set(clusterLabels)] } )
        '''
                
    return result



if __name__ == "__main__":
    
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
    
    
    repgenes = representativeScore(profileData)
    print repgenes
    
    print profileData["profile"].ix[misList["mis1"]].T
    
    me = np.mean(profileData["profile"].ix[misList["mis1"]].T)
    
    print me, type(me) 
    
    clusterLabs = [1,3,1,2,2,2,2,3,1,3]
    
    tripleList = zip(clusterLabs, profileData["classes"], tmppatients)
    
    result = {cid:{"positive":[], "negative":[]} for cid in set(clusterLabs)}
    
    for triple in tripleList:
        if triple[1] == 1:
            result[triple[0]]["positive"].append(triple[2])
        else: 
            result[triple[0]]["negative"].append(triple[2])
    
    print
    for cid in result:
        cluster = result[cid]
        patients = cluster['positive'] + cluster['negative']
        centroid = np.array(np.mean(profileData["profile"][patients].T))
        print centroid
           
    '''
    result = []
    result.append( [{"positive":[x[2] for x in tripleList if x[0] == clab and x[1] == 1], 
                     "negative":[x[2] for x in tripleList if x[0] == clab and x[1] == 0] } for clab in set(clusterLabs)] )
    '''
    print result
    
    '''
    dataset = datasets.load_iris()
    x = dataset.data
    y = dataset.target_names
    
    print y
    print dataset.feature_names
    frame = DataFrame(data=x, columns=dataset.feature_names)
    
#    km = KMeans(n_clusters=3, random_state=1).fit(frame[['sepal width (cm)', 'petal length (cm)']])
    km = AgglomerativeClustering(n_clusters=3).fit(frame[['sepal width (cm)', 'petal length (cm)']])
    
    print km
    print km.labels_
    '''