'''
Created on 2016. 8. 25.

@author: sun
'''

from dataload import load_data
import clustering as cl
import misranking as mr

class cossyPlus():
    def __init__(self,param):
        self.analyze_type = param.analyze_type
        self.cluster_num = param.cluster_num
        self.representativeGene_num = param.representativeGene_num
        self.clustering_method = param.clustering_method
        self.exp_normalization_type = param.exp_normalization_type

        self.exp_file = param.exp_file
        self.mutation_file = param.mutation_file
        self.smoothing_source_file = param.smoothing_source_file

        self.gmt_file = param.gmt_file


        self.run()

    def run(self):
        dataload_result = self.loadData()
        clustering_result = self.clustering(dataload_result)
        entropy_result = self.misranking(clustering_result)
        return entropy_result

    def loadData(self):
        print "start loading data.."

        if self.analyze_type =='expression':
            dataload_result = self.loadData(exp_file=self.exp_file,  gmt_file=self.gmt_file, analyzing_type=self.analyze_type , exp_normalize_tpye= self.exp_normalization_type)
        elif self.analyze_type =='mutation':
            dataload_result = self.loadData(mutation_file =self.mutation_file , gmt_file=self.gmt_file, analyzing_type=self.analyze_type , network_file_for_smoothing=self.smoothing_source_file)
        elif self.analyze_type =='mut_with_exp':
            dataload_result = self.loadData(exp_file=self.exp_file, mutation_file =self.mutation_file , gmt_file=self.gmt_file, analyzing_type=self.analyze_type , network_file_for_smoothing=self.smoothing_source_file , exp_normalize_tpye= self.exp_normalization_type)
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
    def fit(self):
        
        pass
    
    # making classification model
    def makeModel(self):
        
        pass
    
    
if __name__ == "__main__":
    print "start"