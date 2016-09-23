__author__ = 'SungJoonPark'


from icossy import parameter, centralEngine, classification
from icossy.centralEngine import cossyPlus

if __name__ == '__main__':
    experimentGuy = 'sunwon'
    if experimentGuy =='sungjoon':
        datasets = ["BRCA", "COAD"]
    elif experimentGuy =='sunwon':
        datasets = ["LUSC", "PRAD", "STAD"]
    else:
        raise Exception("wrong experiment guy")

    #analysis_types = ["expression", "mut_with_exp", "mutation"]
    analysis_types = ["mut_with_exp","mutation"]
    #analysis_types = ['mut_with_exp']
    misFiles = {"kegg":"kegg_cossy_symbol.gmt", "keggClusterOne":"keggWhole_clusterONE_symbol.gmt", "stringClusterOne":"string_clusterONE_symbol.gmt", "string":"string_cossy_symbol.gmt"}
    #smoothingSourceFiles = {"string":"string_fix_excel_problem_tab_seperator.sif","reactome":"PathwayCommons.8.reactome.BINARY_SIF_DRUG_deleted.hgnc.txt.sif" , "kegg":"KEGG_matrix_GeneSymbol_0_1.sif"}
    smoothingSourceFiles = {"reactome":"PathwayCommons.8.reactome.BINARY_SIF_DRUG_deleted.hgnc.txt.sif" , "kegg":"KEGG_matrix_GeneSymbol_0_1.sif"}
    output_dir = "Q:/COSSY+/icossy_result/TCGA/ICGC/"

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
                    if aType =='expression':
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

                    p.doTenfolds=True
                    p.enhancedRobustness=True

                    cossy_plus = cossyPlus(p)
                    
                    numOfFolds = 10
                    num_of_correct = 0
                    
                    foldData = cossy_plus.makeFolds(cossy_plus.dataload_result['profileData'], numOfFolds)
                    for foldID in range(numOfFolds):

                        trainFoldIdx = range(numOfFolds)
                        trainFoldIdx.remove(foldID)
            
                        trainFolds = cossy_plus.merged(foldData, trainFoldIdx)
                        trainData = {"profileData" : trainFolds, "misList":cossy_plus.dataload_result['misList']}
            
                        testFold = foldData[foldID]
                        testData = {"profileData" : testFold, "misList":cossy_plus.dataload_result['misList']}
                        
                        cossy_plus.run()
                        
                        trainEntropy_result = cossy_plus.entropy_result
                        trainTopkClustering_result = {mis_entropy[0] : cossy_plus.clustering_result[mis_entropy[0]] for mis_entropy in trainEntropy_result}
                        predict_dict = classification.fit(trainTopkClustering_result,testData)
                        predict_dict = classification.fit(cossy_plus.clustering_result,testData)
                        obs_dict = {patient:testData['profileData']['classes'][i] for i,patient in enumerate(testData['profileData']['profile'].columns)}
            
                        for patient in predict_dict.keys():
                            if predict_dict[patient] == obs_dict[patient]:
                                num_of_correct = num_of_correct+1
                    
                    #accuracy = cossy_plus.run_CV()
                    num_of_total = len(cossy_plus.dataload_result['profileData']['profile'].columns)
                    accuracy = float(num_of_correct) / float(num_of_total)
                    
                    outputfile = output_dir + dataset + "/" + "for_test/" + aType +"/"+"mis_"+mistype+"_smoothing_"+smoothingSource+"_10foldsAccuracy.txt"
                    w = open(outputfile,'w')
                    w.write(str(accuracy))
                    w.close()

                    print "done! : " + dataset + " / " + aType + " / " + mistype +"/" + smoothingSource +"/" +str(accuracy)
                    
else:
    experimentGuy = 'sunwon'
    if experimentGuy =='sungjoon':
        datasets = ["BRCA", "COAD"]
    elif experimentGuy =='sunwon':
        datasets = ["LUSC", "PRAD", "STAD"]
    else:
        raise Exception("wrong experiment guy")

    #analysis_types = ["expression", "mut_with_exp", "mutation"]
    analysis_types = ["mut_with_exp","mutation"]
    #analysis_types = ['mut_with_exp']
    misFiles = {"kegg":"kegg_cossy_symbol.gmt", "keggClusterOne":"keggWhole_clusterONE_symbol.gmt", "stringClusterOne":"string_clusterONE_symbol.gmt", "string":"string_cossy_symbol.gmt"}
    #smoothingSourceFiles = {"string":"string_fix_excel_problem_tab_seperator.sif","reactome":"PathwayCommons.8.reactome.BINARY_SIF_DRUG_deleted.hgnc.txt.sif" , "kegg":"KEGG_matrix_GeneSymbol_0_1.sif"}
    smoothingSourceFiles = {"reactome":"PathwayCommons.8.reactome.BINARY_SIF_DRUG_deleted.hgnc.txt.sif" , "kegg":"KEGG_matrix_GeneSymbol_0_1.sif"}
    output_dir = "Q:/COSSY+/icossy_result/TCGA/ICGC/"

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
                    if aType =='expression':
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

                    p.doTenfolds=True
                    p.enhancedRobustness=True

                    cossy_plus = cossyPlus(p)
                    accuracy = cossy_plus.run_CV()

                    outputfile = output_dir + dataset + "/" + "for_test/" + aType +"/"+"mis_"+mistype+"_smoothing_"+smoothingSource+"_10foldsAccuracy.txt"
                    w = open(outputfile,'w')
                    w.write(str(accuracy))
                    w.close()

                    print "done! : " + dataset + " / " + aType + " / " + mistype +"/" + smoothingSource +"/" +str(accuracy)