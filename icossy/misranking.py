'''
Created on 2016. 8. 22.

@author: sun
'''

from math import log
from operator import itemgetter

def log2(x):
    
    return log(x) / log(2)

def computeEntropy(clusterResults):
    misRanking = []
    
    for clusteringResult in clusterResults:

        freqlist = {}
        for misid in clusteringResult:
            clusters = clusteringResult[misid]
            
            allpos = 0.0
            allneg = 0.0
                
            freqlist[misid] = []
            sumofEntropy = 0
            sumofNormalizedSampleNumber = 0
            
            for cluster in clusters:
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
                    entropy = -pi * log2(pi) - (1-pi) * log2(1-pi)
                
                sumofEntropy += (P+N) * entropy
                sumofNormalizedSampleNumber = (P+N)
                
            misEntropy = sumofEntropy / sumofNormalizedSampleNumber
            
            misRanking.append( (misid, misEntropy) )
    
    return sorted(misRanking, key=itemgetter(1))

if __name__ == "__main__":
    print "start"
    
    clusterResult = [{"MIS1":[{'positive': ['p1', 'p3'], 'negative': ['p9']}, {'positive': ['p4', 'p5', 'p6', 'p7'], 'negative': []}, {'positive': ['p2'], 'negative': ['p8', 'p10']}]}]
    
    ent = computeEntropy(clusterResults=clusterResult)
    
    print ent