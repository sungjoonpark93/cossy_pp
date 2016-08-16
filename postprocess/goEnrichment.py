'''
Created on 2016. 8. 12.

@author: sun
'''

from gprofiler import GProfiler
from operator import itemgetter

import numpy as np
import json

class GOEnrichmentTester():
    
    def __init__(self):
        self.gp = GProfiler("COSSY++/1.5")
    
    def getGoTerms(self, genelist):
        result = []
        res = self.gp.gprofile(query=genelist)
        
        for i in range(len(res)):
            pvalue = res[i][2]
            goid = res[i][8]
            gocat = res[i][9]
            goterm = res[i][11]
            
            '''
            if (gocat =="MF" or gocat == "CC" or gocat == "BP"):
                result.append({"pvalue":pvalue, "id":goid, "category":gocat, "term":goterm})
            '''
            result.append({"pvalue":pvalue, "id":goid, "category":gocat, "term":goterm})
        
        return result
    
    def readTSV(self, fname):
        records = []
        with open(fname) as reader:
            headers = []
            
            for line in reader:
                values = [x.replace("\"","") for x in line.split("\t")]
                if line.startswith("Gene Symbol"):
                    headers = values
                    continue
                
                rec = {headers[i] : values[i] for i in range(len(headers))}
                
                records.append(rec)
                
        return records
    
    def loadCOSMIC(self, fname):
        self.result = {"somatic":{}, "germline":{}}
        self.diseaseList = []
        
        records = self.readTSV(fname=fname)
        
        for rec in records:
            geneSymbol = rec["Gene Symbol"]
            somaticTumors = [x.strip() for x in rec["Tumour Types(Somatic)"].strip().split(",")]
            germlineTumors = [x.strip() for x in rec["Tumour Types(Germline)"].strip().split(",")]
            
            for tumorType in somaticTumors:
                if tumorType == "":
                    continue;
                
                if tumorType not in self.result["somatic"]:
                    self.result["somatic"][tumorType] = []
                self.result["somatic"][tumorType].append(geneSymbol)
                
                if tumorType not in self.diseaseList:
                    self.diseaseList.append(tumorType)
                
                
            for tumorType in germlineTumors:
                if tumorType == "":
                    continue;
                
                if tumorType not in self.result["germline"]:
                    self.result["germline"][tumorType] = []
                self.result["germline"][tumorType].append(geneSymbol)
                
                if tumorType not in self.diseaseList:
                    self.diseaseList.append(tumorType)
        
        self.makeGOList()

    def getGenes(self, disease):
        
        if disease in self.result["somatic"]:
            somaticGenes = self.result["somatic"][disease]
        else:
            somaticGenes = []
        
        if disease in self.result["germline"]:
            germlineGenes = self.result["germline"][disease]
        else:
            germlineGenes = []
        
        return somaticGenes + germlineGenes
    
    def makeGOList(self):
        self.GOList = {}
        
        for tumorType in self.diseaseList:
            print "."
            genes = self.getGenes(tumorType)
            goTerms = self.getGoTerms(genes)
            
            goTerms = sorted(goTerms, cmp=self.pvaluecomp)
            
            self.GOList[tumorType] = goTerms
    
    def writeCOSMICGO(self, fname):
        with open(fname, "w") as w:
            json.dump(self.GOList, w, indent=4)
    
    def corr(self, genes, disease):
        inputGO = sorted(self.getGoTerms(genes), cmp=self.pvaluecomp)
        inputGO_terms = [x["term"] for x in inputGO]
        
        answerGO = sorted([x for x in self.GOList[disease] if x["term"] in inputGO_terms], cmp=self.pvaluecomp)
        answerGO_terms = [x["term"] for x in answerGO]
        
        assert(len(inputGO_terms) != len(answerGO_terms))
        
        inputGO_ranks_pair = [(x,inputGO_terms.index(x)) for x in inputGO_terms]
        answerGO_ranks_pair = [(x,answerGO_terms.index(x)) for x in answerGO_terms]
        
        inputGO_ranks = [x[1] for x in sorted(inputGO_ranks_pair, key=itemgetter(0))]
        answerGO_ranks = [x[1] for x in sorted(answerGO_ranks_pair, key=itemgetter(0))]
        
        np.correlate(inputGO_ranks, answerGO_ranks, "same")
        
    def pvaluecomp(self, a,b):
        x = a['pvalue']
        y = b['pvalue']
        if x > y:
            return 1
        elif x < y:
            return -1
        else:
            return 0