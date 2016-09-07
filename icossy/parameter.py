__author__ = 'SungJoonPark'


class Parameter:
    def __init__(self):
        self.analyze_type = "expression"
        self.cluster_num = 8
        self.representativeGene_num = 5
        self.exp_normalization_type ="fuzzy"

        self.enhancedRobustness = True
        self.mis_num =31

        self.exp_file = None
        self.mutation_file = None
        self.smoothing_source_file =None

        self.gmt_file = None
        self.misResult_file = None


        self.doTenfolds = False

if __name__ =='__main__':
    p = Parameter()
    p.analyze_type = "exp_mut"
    print p.analyze_type