'''
Created on 2016. 8. 16.

@author: sun
'''
import json
import os
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

def getAllDescendantFiles(root):
    retfiles = []
    
    for (path, dir, files) in os.walk(root):
        for filename in files:
            retfiles.append ( (path,filename) )
    return retfiles

rootdir = "Q:/COSSY+/misResultCSVOnly/ICGC_TCGA/fortest/"
misrootdirs = {"BRCA":getAllDescendantFiles(rootdir + "BRCA/"), "COAD":getAllDescendantFiles(rootdir + "COAD/"), "LUSC":getAllDescendantFiles(rootdir + "LUSC"), "PRAD":getAllDescendantFiles(rootdir + "PRAD"), "STAD":getAllDescendantFiles(rootdir + "STAD")}

gotester = GOEnrichmentTester()

cosmicFilename = "G:/cossy++/cosmicGO.json"
cosmicAnswer = loadCOSMICAnnotation(cosmicFilename)


for dataset in misrootdirs:
    files = misrootdirs[dataset]
    
    for (path, filename) in files:
        print path
        print filename
        with open(path + "/" + filename,'r') as r:
            for line in r:
                tokens = [x.strip() for x in line.replace('\"','').split(",") if len(x.strip()) > 0]
                
                if tokens[0] == 'gene1':
                    continue;
                
                genelist = tokens[1:len(tokens)]
                
                goTerms = [x["term"] for x in gotester.getGoTerms(genelist)]
                cosmicTerms = [x["term"] for x in cosmicAnswer[dataset]]
                
                score = jaccard(goTerms, cosmicTerms)
                
                result = 
                print tokens
        break
    
    break
