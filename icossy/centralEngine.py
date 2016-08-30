'''
Created on 2016. 8. 25.

@author: sun
'''

import dataload as dl
import clustering as cl
import misranking as mr
import random

class cossyPlus():
    def __init__(self,param):
        self.analyze_type = param.analyze_type
        self.cluster_num = param.cluster_num
        self.representativeGene_num = param.representativeGene_num
        self.clustering_method = param.clustering_method
        self.exp_normalization_type = param.exp_normalization_type
        self.is_loading_class_file =param.is_loading_class_file
        self.do_tenFold = param.do_tenFold

        self.exp_file = param.exp_file
        self.mutation_file = param.mutation_file
        self.smoothing_source_file = param.smoothing_source_file

        self.gmt_file = param.gmt_file

        self.misReulst_file = param.misResult_file

        self.run()

    def run(self, enhancedRobustness=False):
        if enhancedRobustness == False:
            self.dataload_result = self.loadData()
            self.misList = self.dataload_result['misList']
            self.clustering_result = self.clustering(self.dataload_result)
            self.entropy_result = self.misranking(self.clustering_result)
            self.write_misResult_from_entropyResult(self.entropy_result,self.misReulst_file,self.misList)
        else:
            self.dataload_result = self.loadData()
            
            
            pass
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
        
        return {"profile":pandas.concat([foldedData[x]["profile"] for x in idxs], axis=1), "classes": , "labels":foldedData[idx[0]]["labels"]}
            
    
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
        return cl.clusteringInMIS(data["profileData"], self.cluster_num,self,self.representativeGene_num, data['misList'])
    
    def misranking(self, clusternig_result):
        print "calcluating entropy and ranking mis"
        return mr.computeEntropy(clusternig_result)
    
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
    pass