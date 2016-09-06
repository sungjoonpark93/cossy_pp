__author__ = 'SungJoonPark'


from icossy import parameter
from icossy.centralEngine import cossyPlus

if __name__ == '__main__':

    datasets = ["BRCA", "COAD", "LUSC", "PRAD", "STAD"]

    #analysis_types = ["expression", "mut_with_exp", "mutation"]
    analysis_types = ["expression","mutation","mut_with_exp"]
    misFiles = {"kegg":"kegg_cossy_symbol.gmt", "keggClusterOne":"keggWhole_clusterONE_symbol.gmt", "stringClusterOne":"string_clusterONE_symbol.gmt", "string":"string_cossy_symbol.gmt"}
    smoothingSourceFiles = {"string":"string_fix_excel_problem_tab_seperator.sif","reactome":"PathwayCommons.8.reactome.BINARY_SIF_DRUG_deleted.hgnc.txt.sif" , "keggWhole":"KEGG_matrix_GeneSymbol_0_1.sif"}


    for dataset in datasets:
        for aType in analysis_types:
            for mistype in misFiles:
                for smoothingSource in smoothingSourceFiles.keys():
                    misFile = misFiles[mistype]
                    exp_data_dir = "Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/" + dataset + "/for_test/"
                    misoutput_dir = "Q:/COSSY+/mis_result/iCOSSY/TCGA_ICGC/"
                    smoothingSourceDir = "Q:/COSSY+/data/network/"+smoothingSource+"/"

                    p = parameter.Parameter()

                    p.analyze_type = aType
                    if aType =='exp':
                        p.exp_file = exp_data_dir + "exp."+dataset+"-US.tsv_preprocessed.gct"
                    elif aType =='mutation':
                        p.mutation_file = exp_data_dir + "mut."+dataset+"-US.tsv_preprocessed.csv"
                        p.smoothing_source_file = smoothingSourceDir + smoothingSourceFiles[smoothingSource]
                    elif aType =='mut_with_exp':
                        p.exp_file = exp_data_dir + "exp."+dataset+"-US.tsv_preprocessed.gct"
                        p.mutation_file = exp_data_dir + "mut."+dataset+"-US.tsv_preprocessed.csv"
                        p.smoothing_source_file = smoothingSourceDir + smoothingSourceFiles[smoothingSource]

                    p.misResult_file = misoutput_dir + dataset+"/"+aType+"_icossy_"+dataset+"_"+mistype+"_result.txt"
                    p.gmt_file = "Q:/COSSY+/data/mis/icossy/"+misFile


                    cossyPlus(p)

                    print "done! : " + dataset + " / " + aType + " / " + mistype +"/" + smoothingSource +"/"