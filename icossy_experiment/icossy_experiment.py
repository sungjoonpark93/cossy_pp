'''
Created on 2016. 9. 2.

@author: sun
'''
from icossy import parameter
from icossy.centralEngine import cossyPlus

if __name__ == '__main__':
    
    datasets = ["BRCA", "COAD", "LUSC", "PRAD", "STAD"]
    
    #analysis_types = ["expression", "mut_with_exp", "mutation"]
    analysis_types = ["expression"]
    misFiles = {"kegg":"kegg_cossy_symbol.gmt", "keggClusterOne":"keggWhole_clusterONE_symbol.gmt", "stringClusterOne":"string_clusterONE_symbol.gmt", "string":"string_cossy_symbol.gmt"}
    
    for dataset in datasets:
        for aType in analysis_types:
            for mistype in misFiles:
                misFile = misFiles[mistype]
                exp_data_dir = "Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/" + dataset + "/for_test/"    
                misoutput_dir = "Q:/COSSY+/mis_result/iCOSSY/TCGA_ICGC/"
            
                exp_p = parameter.Parameter()
    #            exp_p.analyze_type = "expression"
                exp_p.analyze_type = aType
                exp_p.exp_file = exp_data_dir + "exp."+dataset+"-US.tsv_preprocessed.gct"
                exp_p.mutation_file = exp_data_dir + "mut."+dataset+"-US.tsv_preprocessed.csv"
                exp_p.misResult_file = misoutput_dir + dataset+"/"+aType+"_icossy_"+dataset+"_"+mistype+"_result.txt"
                # exp_p.gmt_file = "Q:/COSSY+/data/mis/icossy/kegg_cossy_symbol.gmt"
#                exp_p.smoothing_source_file = ""
                exp_p.gmt_file = "Q:/COSSY+/data/mis/icossy/"+misFile
                
                cossyPlus(exp_p)
    
                print "done! : " + dataset + " / " + aType + " / " + mistype