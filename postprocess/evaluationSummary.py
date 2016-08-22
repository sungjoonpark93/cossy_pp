'''
Created on 2016. 8. 18.

@author: sun
'''

from postprocess.annotationScore import getAllDescendantFiles
from gorillaImpl import GOrilla

import numpy as np

if __name__ == '__main__':
    summarydir = "Q:/COSSY+/misResultCSVOnly/ICGC_TCGA/gprofile/"
    
   
    datasets = ["BRCA", "COAD", "LUSC", "PRAD", "STAD"]
    
    for dataset in datasets:
        resdir = summarydir + dataset + "/"
        
        resfiles = getAllDescendantFiles(resdir)
        
        with open(summarydir + dataset + "_result(jaccard+accuracy).csv",'w') as w:
            w.write("condition,avg.jaccard,avg.accuracy\n")
            for (path, resfile) in resfiles:
                print "opended!" + dataset
                with open(path + "/" + resfile, 'r') as r:
                    print "-- " + resfile
                    jacscores = []
                    accscores = []
                    for line in r:
                        (misid, querygo, cosmicgo, jacscore, accscore) = line.split(",")
                        if misid == 'mis#':
                            continue
                        
                        jacscores.append(float(jacscore))
                        accscores.append(float(accscore))
                
                    meanjac = np.mean(jacscores)
                    meanacc = np.mean(accscores)
                    
                    #w.write(resfile+"\t"+ str(meanjac)+"\t"+ str(meanacc)+"\n")
                    w.write("%s,%.4f,%.4f\n" % (resfile, meanjac, meanacc))
                
                print "."