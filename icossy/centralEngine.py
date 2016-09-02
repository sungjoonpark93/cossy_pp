'''
Created on 2016. 8. 25.

@author: sun
'''

import dataload as dl
import clustering as cl
import misranking as mr
import random
import pandas as pd
import parameter
from pandas.util.testing import makeFloatIndex
from operator import itemgetter

class cossyPlus():
    def __init__(self,param):
        
        self.analyze_type = param.analyze_type
        self.cluster_num = param.cluster_num
        self.representativeGene_num = param.representativeGene_num
        self.mis_num = param.mis_num
        self.exp_normalization_type = param.exp_normalization_type
        
        self.enhancedRobustness = param.enhancedRobustness

#        self.is_loading_class_file = param.is_loading_class_file
#        self.clustering_method = param.clustering_method



        self.exp_file = param.exp_file
        self.mutation_file = param.mutation_file
        self.smoothing_source_file = param.smoothing_source_file

        self.gmt_file = param.gmt_file

        self.misReulst_file = param.misResult_file

        self.run(self.enhancedRobustness)

    def run(self, enhancedRobustness=False):
        if enhancedRobustness == False:
            print "not doing enhancedRobustness"
            self.dataload_result = self.loadData()
            self.misList = self.dataload_result['misList']
            self.clustering_result = self.clustering(self.dataload_result)
            self.entropy_result = self.misranking(self.clustering_result)
            self.write_misResult_from_entropyResult(self.entropy_result,self.misReulst_file,self.misList)
        else:
            print "doing enhancedRobustness"
            self.dataload_result = self.loadData()
            allMISList = self.dataload_result['misList']
            
            numOfFolds = 10
            foldData = self.makeFolds(self.dataload_result['profileData'], numOfFolds)
            
            robustMisList = []
            
            for foldID in range(numOfFolds):
                
                curFoldIdx = range(numOfFolds)
                curFoldIdx.remove(foldID)
                
                curFold = self.merged(foldData, curFoldIdx)
                inData = {"profileData" : curFold, "misList":allMISList}
                
                curClustering_result = self.clustering(inData)
                curEntropy_result = self.misranking(curClustering_result)
#                self.write_misResult_from_entropyResult(curEntropy_result, self.misReulst_file, self.misList)
                
                robustMisList.append(curEntropy_result)
            
            candidateMISwithEntropy = reduce(lambda x, y : x + y, robustMisList)
    
            cnadidateMISids = [x[0] for x in candidateMISwithEntropy]
            miscounts = [(misid, cnadidateMISids.count(misid)) for misid in set(cnadidateMISids)]
            miscounts.sort(key=itemgetter(1), reverse=True)
            
            finalMISList = [x[0] for x in miscounts[0:self.mis_num]]
            
            self.misList = {x : allMISList[x] for x in allMISList if x in finalMISList}
            self.dataload_result['misList'] = self.misList
            self.clustering_result = self.clustering(self.dataload_result)
            self.entropy_result = self.misranking(self.clustering_result)
            self.write_misResult_from_entropyResult(self.entropy_result,self.misReulst_file,self.misList)
            
        return self.entropy_result

    def makeFolds(self, profileData, numOfFolds=10):
        
        folds = [[] for x in range(numOfFolds)]
        profile = profileData["profile"]
        classes = profileData["classes"]
            
        pairlist = [ (classes[x], profile.columns[x]) for x in range(len(classes))]
        random.shuffle(pairlist)
    
        pospairs = enumerate([ x for x in pairlist if x[0] == 1])
        negpairs = enumerate([ x for x in pairlist if x[0] == 0])
        
        
        for idx, v in pospairs:
            i = idx%numOfFolds
            folds[i].append(v)
        
        for idx, v in negpairs:
            i = numOfFolds - idx%numOfFolds -1
            folds[i].append(v)
        
        foldedData = []
    
        for fold in folds:
            pids = [x[1] for x in fold]
            classes = [x[0] for x in fold]
            
            profileSubset = profile[pids]
            
            foldedData.append({"profile":profileSubset, "classes":classes, "labels":profileData["labels"]})
        
        return foldedData
    
    def merged(self, foldedData, mergingIdx):
        #labels don't have to be merged just use one of foldData because it's same for all foldData's element. ('tumor', 'normal')
        profileMerged = pd.concat([foldedData[x]["profile"] for x in mergingIdx], axis=1)

        classesMerged = []
        for idx in mergingIdx:
            classesMerged.extend(foldedData[idx]['classes'])

        labelsMerged = foldedData[mergingIdx[0]]["labels"]
        return {"profile":profileMerged, "classes":classesMerged, "labels":labelsMerged}

    def loadData(self):
        print "start loading data.."

        if self.analyze_type =='expression':
            dataload_result = dl.load_data(exp_file=self.exp_file,  gmt_file=self.gmt_file, analyzing_type=self.analyze_type , exp_normalize_tpye= self.exp_normalization_type)
        elif self.analyze_type =='mutation':
            dataload_result = dl.load_data(mutation_file =self.mutation_file , gmt_file=self.gmt_file, analyzing_type=self.analyze_type , network_file_for_smoothing=self.smoothing_source_file)
        elif self.analyze_type =='mut_with_exp':
            dataload_result = dl.load_data(exp_file=self.exp_file, mutation_file =self.mutation_file , gmt_file=self.gmt_file, analyzing_type=self.analyze_type , network_file_for_smoothing=self.smoothing_source_file , exp_normalize_tpye= self.exp_normalization_type)
        else:
            raise Exception('unspecified analyzing type')
        return dataload_result
    
    def clustering(self, data):
        print "start clustering..."
        return cl.clusteringInMIS(data["profileData"], self.cluster_num,self.representativeGene_num, data['misList'])
    
    def misranking(self, clusternig_result):
        print "calcluating entropy and ranking mis"
        return mr.computeEntropy(clusternig_result, self.mis_num)
    
    # classification
    def fit(self, data):
        
        fittingResult = {}
        classification = cl.classify(data, self.clustering_result)
        
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
    
    # making classification model
    def makeModel(self):
        
        pass
    
    def write_misResult_from_entropyResult(self,entropy_result, outputfile, misList):
        if outputfile == None:
            raise Exception("you didn't specified the outputfile")

        w = open(outputfile,'w')
        for misid_result_tuple in entropy_result:

            misid = misid_result_tuple[0]
            entropy = misid_result_tuple[1]
            mis_genes = misList[misid]

            w.write(misid+"\t"+str(entropy)+"\t")
            for mis_gene in mis_genes:
                w.write(mis_gene)
                w.write("\t")
            w.write("\n")
        w.close()


if __name__ == "__main__":
    data_dir = "Q:\COSSY+\data\preprocessed\TCGA\ICGC/release21/BRCA/for_test/"
    misoutput_dir = "Q:\COSSY+\mis_result\iCOSSY\TCGA_ICGC/"

    p=parameter.Parameter()
    p.analyze_type = "expression"
    p.exp_file=data_dir+"exp.BRCA-US.tsv_preprocessed.gct"
    p.mutation_file = data_dir+"mut.BRCA-US.tsv_preprocessed.csv"
    p.misResult_file =misoutput_dir + "BRCA/temp_icossy_brca_result.txt"
    p.gmt_file = "Q:\COSSY+\data\mis\icossy/kegg_cossy_symbol.gmt"

    cossyPlus(p)
    pass


