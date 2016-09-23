'''
Created on 2016. 9. 12.

@author: sun
'''

import os

if __name__ == '__main__':
    output_dir = "Q:/COSSY+/icossy_result/TCGA/ICGC/"
    datasets = ["BRCA", "COAD", "LUSC", "PRAD", "STAD"]
    analysis_types = ["expression", "mut_with_exp", "mutation"]
    misFiles = ["kegg", "keggClusterOne", "stringClusterOne", "string"]
    smoothingSourceFiles = ["string", "reactome", "kegg"]
    
    outfname = output_dir + "summary.txt"
    outf = open(outfname, 'w')
    
    for dataset in datasets:
        for aType in analysis_types:
            for mistype in misFiles:
                
                if aType == 'expression':
                    fname = output_dir + dataset + "/" + "for_test/" + aType + "/" + "mis_" + mistype + "_10foldsAccuracy.txt"
                        
                    r = open(fname, 'r')
                    
                    content = dataset + '\t' + aType + '\t' + mistype + '\t' + 'none' + '\t' + r.read()
                    print content
                    outf.write(content + '\n')
                    
                    r.close()
                else:
                    for smoothingSource in smoothingSourceFiles:
                        fname = output_dir + dataset + "/" + "for_test/" + aType + "/" + "mis_" + mistype + "_smoothing_" + smoothingSource + "_10foldsAccuracy.txt"
                        
                        r = open(fname, 'r')
                        
                        content = dataset + '\t' + aType + '\t' + mistype + '\t' + smoothingSource + '\t' + r.read()
                        print content
                        outf.write(content + '\n')
                        
                        r.close()
    outf.close()
