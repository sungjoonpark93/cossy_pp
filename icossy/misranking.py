'''
Created on 2016. 8. 22.

@author: sun
'''

from math import log
from operator import itemgetter

def log2(x):
    
    return log(x) / log(2)

def computeEntropy(clusterResults, misNum=31):
    print "calculating entropy"
    misRanking = []
    
    for misid in clusterResults:
        clusteringResult = clusterResults[misid]
        freqlist = {}
        
        clusters = clusteringResult['c']
        
        allpos = 0.0
        allneg = 0.0
            
        freqlist[misid] = []
        sumofEntropy = 0
        sumofNormalizedSampleNumber = 0
        
        for cid in clusters:
            cluster = clusters[cid]
            poslist = cluster['positive']
            neglist = cluster['negative']

            allpos += len(poslist)
            allneg += len(neglist)
            
            freqlist[misid].append((len(poslist), len(neglist)))
        
        
        for freq in freqlist[misid]:
            
            P = float(freq[0] / allpos)
            N = float(freq[1] / allneg)
            
            pi = P / (P+N)
            
            # for computational safety
            if pi == 0 or pi == 1:
                entropy = 0
            else :
                entropy = - pi * log2(pi) - (1-pi) * log2(1-pi)
            
            sumofEntropy += (P+N) * entropy
            sumofNormalizedSampleNumber = (P+N)
            
        misEntropy = sumofEntropy / sumofNormalizedSampleNumber
        
        misRanking.append( (misid, misEntropy) )
    
    return sorted(misRanking, key=itemgetter(1))[0:misNum]

if __name__ == "__main__":
    print "start"
    
    clusterResults = {'mis1': {'c': {0: {'centroid': [ 0.87142857,  0.38571429,  0.32857143],
                                         'label': 1.0,
                                         'negative': ['p5', 'p6', 'p7'],
                                         'positive': ['p1', 'p2', 'p3', 'p4']},
                                     1: {'centroid': [ 0.06666667,  0.53333333,  0.4       ],
                                         'label': 0.0,
                                         'negative': ['p8', 'p11'],
                                         'positive': ['p13']},
                                     2: {'centroid': [ 0.06666667,  0.63333333,  0.1       ],
                                         'label': 0.0,
                                         'negative': ['p9', 'p10'],
                                         'positive': ['p12']}},          
                                'repGenes': ['g1', 'g3', 'g4']},
                      'mis2': {'c': {0: {'centroid': [ 0.13333333,  0.33333333,  0.66666667],
                                         'label': 0.0,
                                         'negative': ['p10', 'p11'],
                                         'positive': ['p12']},
                                     1: {'centroid': [ 0.61428571,  0.64285714,  0.15714286],
                                         'label': 1.0,
                                         'negative': [],
                                         'positive': ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7']},
                                     2: {'centroid': [ 0.2       ,  0.13333333,  0.36666667],
                                         'label': 0.0,
                                         'negative': ['p8', 'p9', 'p13'],
                                         'positive': []}},
                               'repGenes': ['g2', 'g5', 'g8']},
                      'mis3': {'c': {0: {'centroid': [ 0.16666667,  0.81666667,  0.28333333],
                                         'label': 1.0,
                                         'negative': ['p6'],
                                         'positive': ['p2', 'p3', 'p4', 'p5', 'p7']},
                                     1: {'centroid': [ 0.58333333,  0.6       ,  0.08333333],
                                         'label': 0.0,
                                         'negative': ['p8', 'p9', 'p10', 'p11', 'p12', 'p13'],
                                         'positive': []},
                                     2: {'centroid': [ 0.4,  0.9,  0.5],
                                         'label': 1.0,
                                         'negative': [],
                                         'positive': ['p1']}},
                               'repGenes': ['g6', 'g10', 'g7']}}
    
    ent = computeEntropy(clusterResults=clusterResults, misNum=3)
    
    print ent
    
    
    print
    print
    
    ent2 = [('mis4', 0.6928760584965482), ('mis1', 3.5178502274200967), ('mis2', 4.236759932620531)]
    ent3 = [('mis7', 0.6928760584965482), ('mis6', 3.5178502274200967), ('mis5', 4.236759932620531)]
    ent4 = [('mis2', 0.6928760584965482), ('mis5', 3.5178502274200967), ('mis1', 4.236759932620531)]
    
    allent = [ent,ent2,ent3,ent4]
    allent2 = reduce(lambda x, y : x + y, allent)
    
    misidlist = [x[0] for x in allent2]
    miscounts = [(misid, misidlist.count(misid)) for misid in set(misidlist)]
    miscountsorted = sorted(miscounts, key=itemgetter(1), reverse=True)[0:3]
    
    print allent
    print
    print allent2
    print
    print misidlist
    print miscountsorted
    