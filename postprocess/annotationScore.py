'''
Created on 2016. 8. 16.

@author: sun
'''
from gprofiler import GProfiler
import json
from postprocess.goEnrichment import GOEnrichmentTester

def jaccard(list1, list2):
    intersect = [x for x in list1 if x in list2]
    union = [x for x in list1] + [x for x in list2 if x not in list1]
    
    return float(len(intersect)) / float(len(union))

def loadCOSMICAnnotation(fname):
    with open(fname,"r") as jsonfp:
        annotation = json.load(jsonfp)
    
    return annotation

def getscore(termlist, disease, annotation):
    diseaseAnnotation = annotation[disease]
    annotatedTerms_answer = [a["term"] for a in diseaseAnnotation]
    annotatedTerms_query = [a["term"] for a in termlist]
    
    
    return jaccard(annotatedTerms_answer, annotatedTerms_query)

