'''
Created on 2016. 8. 16.

@author: sun
'''
import json
import os
import pprint
from postprocess.goEnrichment import GOEnrichmentTester

def jaccard(list1, list2):
    intersect = [x for x in list1 if x in list2]
    union = [x for x in list1] + [x for x in list2 if x not in list1]
    
    return float(len(intersect)) / float(len(union))

def loadCOSMICAnnotation(fname, mapfname):
    with open(fname, "r") as jsonfp:
        rawannotation = json.load(jsonfp)
    
    annotation = {}
    for d in rawannotation:
        annotation[d] = [x['term'] for x in rawannotation[d]]
        
    map = {}
    
    with open(mapfname, "r") as r:
        for line in r:
            toks = line.strip().split(",")
            if toks[1] != '-':
                map[toks[0]] = toks[1]
    
    for d in map:
        if map[d] not in annotation:
            annotation[map[d]] = []
        
        annotation[map[d]] = [x for x in annotation[d]] + [x for x in annotation[map[d]] if x not in annotation[d]]
        
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
            retfiles.append ((path, filename))
    return retfiles

if __name__ == '__main__':
    rootdir = "Q:/COSSY+/misResultCSVOnly/ICGC_TCGA/fortest/"
    resdir = "Q:/COSSY+/misResultCSVOnly/ICGC_TCGA/gprofile/"
    misrootdirs = {"BRCA":getAllDescendantFiles(rootdir + "BRCA/"),
                   "COAD":getAllDescendantFiles(rootdir + "COAD/"),
                   "LUSC":getAllDescendantFiles(rootdir + "LUSC/"),
                   "PRAD":getAllDescendantFiles(rootdir + "PRAD/"),
                   "STAD":getAllDescendantFiles(rootdir + "STAD/")}
    
    gotester = GOEnrichmentTester()
    pp = pprint.PrettyPrinter(indent=4)
    
    cosmicFilename = "G:/cossy++/cosmicGO.json"
    cosmicMapFilename = "G:/cossy++/cosmicDiseaseMapping.csv"
    cosmicAnswer = loadCOSMICAnnotation(cosmicFilename, cosmicMapFilename)
    
    tmp = cosmicAnswer["BRCA"]
    
    pp.pprint(tmp)
    print len(tmp)
    
    for dataset in misrootdirs:
        files = misrootdirs[dataset]
        
        for (path, filename) in files:
            print path
            print filename
            with open(path + "/" + filename, 'r') as r:
                print filename
                scorefilename = resdir + "/" + dataset + "/" + "score_" + filename
                w = open(scorefilename, 'w')
                w.write("mis#,queryGOterm,cosmicGOterm,jaccard,accuracy\n")
                for line in r:
                    tokens = [x.strip() for x in line.replace('\"', '').split(",") if len(x.strip()) > 0]
                    
                    if tokens[0] == 'gene1':
                        continue;
                    
                    genelist = tokens[1:len(tokens)]
                    goList = gotester.getGoTerms(genelist)
                    goTerms = [x["term"] for x in goList]
    #                cosmicTerms = [x["term"] for x in cosmicAnswer[dataset]]
                    cosmicTerms = cosmicAnswer[dataset]
                    
                    score = jaccard(goTerms, cosmicTerms)
                    if(len(goTerms) == 0) :
                        accScore = 0
                    else:
                        accScore = float(len([x for x in goTerms if x in cosmicTerms])) / float(len(goTerms))
                        
                    print tokens[0], ": ", len(goTerms), "/", len(cosmicTerms), score, accScore
                    w.write("%s,%d,%d,%.4f,%.4f\n" % (tokens[0], len(goTerms), len(cosmicTerms), score, accScore))
                
                w.close()
