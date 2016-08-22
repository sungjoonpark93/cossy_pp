'''
Created on 2016. 8. 18.

@author: sun
'''

import math
from postprocess.annotationScore import getAllDescendantFiles,\
    loadCOSMICAnnotation, jaccard
import pprint

class GOrilla():
    def __init__(self, filename):
        self.GOAnnotation = self.loadGOAnnotation(filename)
    
    def loadGOAnnotation(self, filename):
        GOAnnotation = {"genemap":{}, "gomap":{}, "go":{}}
        
        with open(filename, 'r') as r:
            for line in r:
                if line[0] == '!':
                    continue
                
                toks = line.split('\t')
                goid = toks[4]
                gostr = toks[9]
                gene = toks[2]
                
                if goid not in GOAnnotation["gomap"]:
                    GOAnnotation["gomap"][goid] = set()
                GOAnnotation["gomap"][goid].add(gene)
                
                if gene not in GOAnnotation["genemap"]:
                    GOAnnotation["genemap"][gene] = set()
                GOAnnotation["genemap"][gene].add(goid)
                
                GOAnnotation["go"][goid]=gostr
                GOAnnotation["go"][gostr]=goid
        
        return GOAnnotation
    
    def comb(self, n,r):
        if r == 0 or n==r:
            return 1
        
        lognf = reduce(lambda x, y: x+y, [math.log(x) for x in range(n+1)[1:]])
        logrf = reduce(lambda x, y: x+y, [math.log(x) for x in range(r+1)[1:]])
        lognrf = reduce(lambda x, y: x+y, [math.log(x) for x in range(n-r+1)[1:]])
        
        logncr = lognf - logrf - lognrf
        
        return math.floor(math.exp(logncr)+0.5)
        
    def HGT(self, b,N,B,n):
        sum = 0
        for i in xrange(b, min(n,B) +1):
            sum += self.comb(n=n,r=i) * self.comb(n=N-n, r=B-i) / self.comb(n=N, r=B)
            
        return sum
    
    def mHG(self, labellist, B):
        N = len(labellist)
        
        return min([self.HGT(b=sum(labellist[0:n+1]), N=N, B=B, n=n+1) for n in range(len(labellist))])
    
    def getLabels(self, genelist):
        labels = {}
        goterms = set()
        
        for gene in genelist:
            if gene in self.GOAnnotation:
                goterms |= self.GOAnnotation["genemap"][gene]
        
        for goid in goterms:
            annotatedGenes = self.GOAnnotation["gomap"][goid]
            
            label = list(map(lambda x : int(x in annotatedGenes), genelist))
            labels[goid] = label
            
        return labels
    
    def computeGOEnrichment(self, genelist):
        result = []
        
        labels = self.getLabels(genelist)
        for goid in labels:
            B = len(self.GOAnnotation["gomap"][goid] & set(genelist))
            
            enrichmentScore = self.mHG(labellist=labels[goid], B=B)
            result.append( {"goid":goid, "score":enrichmentScore, "goterm":self.GOAnnotation["go"][goid]} )
        
        return result

if __name__ == '__main__':
    
    gofile = "Q:/gene ontology/goa_human_2016-08-18.gaf"
    print 'loading go annotation...'
    gorilla = GOrilla(gofile)
    print 'end go annotation'
    
    '''
    genes = ["ABL1", "ABL2", "AKT1"]
    
    
    annoLabels = gorilla.getLabels(genelist=genes)
    print annoLabels
    
    en = gorilla.computeGOEnrichment(genelist=genes)
    
    termset = [x["goterm"] for x in en if x["score"] < 0.05]
    print termset
    print en
    print "end!"
    '''
    

    rootdir = "Q:/COSSY+/misResultCSVOnly/ICGC_TCGA/fortest/"
    resdir = "Q:/COSSY+/misResultCSVOnly/ICGC_TCGA/gprofile/"
    misrootdirs = {"BRCA":getAllDescendantFiles(rootdir + "BRCA/"),
                   "COAD":getAllDescendantFiles(rootdir + "COAD/"),
                   "LUSC":getAllDescendantFiles(rootdir + "LUSC/"),
                   "PRAD":getAllDescendantFiles(rootdir + "PRAD/"),
                   "STAD":getAllDescendantFiles(rootdir + "STAD/")}
    
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
#                    goList = gotester.getGoTerms(genelist)
                    en = gorilla.computeGOEnrichment(genelist=genelist)
                    goTerms = [x["goterm"] for x in en if x["score"] < 0.05]
                        
#                    goTerms = [x["term"] for x in goList]
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

    