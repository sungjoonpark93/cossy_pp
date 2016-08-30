__author__ = 'SungJoonPark'


class Parameter:
    def __init__(self):
        self.analyze_type = "exp"
        self.cluster_num = 5
        self.representativeGene_num = 5
        self.clustering_method = ""
        self.exp_normalization_type =""

        self.exp_file = None
        self.mutation_file = None
        self.smoothing_source_file =None

        self.gmt_file = None



if __name__ =='__main__':
    p = Parameter()
    p.analyze_type = "exp_mut"
    print p.analyze_type