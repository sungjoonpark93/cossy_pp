__author__ = 'SungJoonPark'
__author__ = 'SungJoonPark'


from icossy import parameter
from icossy.centralEngine import cossyPlus

if __name__ == '__main__':

    datasets = ['BRCA','COAD','PRAD','LUSC' ,'STAD']

    #smoothingSourceFiles = {"string":"string_fix_excel_problem_tab_seperator.sif","reactome":"PathwayCommons.8.reactome.BINARY_SIF_DRUG_deleted.hgnc.txt.sif" , "kegg":"KEGG_matrix_GeneSymbol_0_1.sif"}
    output_dir = "Q:/COSSY+/icossy_result/TCGA/ICGC/"

    for dataset in datasets:
        misFile = "Q:\COSSY+\data\mis\icossy/hotnet2_string_"+dataset+".gmt"
        exp_data_dir = "Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/" + dataset + "/for_test/"
        output_dir = "Q:/COSSY+/icossy_result/TCGA/ICGC/"
        p = parameter.Parameter()
        p.analyze_type ='expression'
        p.exp_file = exp_data_dir + "exp."+dataset+"-US.tsv_preprocessed.gct"
        outputfile = output_dir + dataset +"/for_test/expression/" + "hotnet2_string_"+dataset+"_10foldsAccuracy.txt"
        p.gmt_file = misFile
        p.doTenfolds=True
        p.enhancedRobustness=True
        cossy_plus = cossyPlus(p)
        accuracy = cossy_plus.run_CV()
        w = open(outputfile,'w')
        w.write(str(accuracy))
        print "done! : " + dataset



