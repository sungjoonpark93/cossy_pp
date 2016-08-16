'''
Created on 2016. 8. 12.

@author: sun
'''


from gprofiler import GProfiler
from postprocess import goEnrichment as en

import pprint

import postprocess.goEnrichment
from postprocess.goEnrichment import GOEnrichmentTester

gp = GProfiler("COSSY++/1.5")

pp = pprint.PrettyPrinter(indent=4)
genelist = ["AKT1", "ABL1", "ABL2", "BCR"]

'''
res = gp.gprofile(genelist)

for i in range(len(res[0])):
    print str(i) + ": " + str(res[0][i])

pvalue = res[0][2]
goid = res[0][8]
gocat = res[0][9]
geterm = res[0][11]

print ""



res = en.getGoTerms(genelist)

pp.pprint(res)

'''
tester = GOEnrichmentTester()

cosmicfile = "G:/cossy++/Census_allWed Jul 20 07-40-20 2016.tsv"

tester.loadCOSMIC(cosmicfile)

cmlgene1 = tester.result["germline"]["breast"]
cmlgene2 = tester.result["somatic"]["breast"]

pp.pprint(cmlgene1)
pp.pprint(cmlgene2)

genes = tester.getGenes("breast")

pp.pprint(genes)

res = tester.getGoTerms(genes)

outfname = "G:/cossy++/cosmicGO.json"
tester.writeCOSMICGO(outfname)

pp.pprint(res)