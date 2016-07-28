'''
Created on 2016. 7. 27.

@author: sun
'''

network = {"reactome":"Q:/COSSY+/data/network/reactome/PathwayCommons.8.reactome.BINARY_SIF_DRUG_deleted.hgnc.txt.sif", "string":"Q:/COSSY+/data/network/string/string_fix_excel_problem_tab_seperator.sif"}

rootdir = "Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/"

dataset_names = ["brca", "coad", "lusc", "prad", "stad"]

mut_src_data = {"brca":rootdir + "BRCA/patient_not_matched/mut/simple_somatic_mutation.open.BRCA-US.tsv_preprocessed.csv",
                "coad":rootdir + "COAD/patient_not_matched/mut/simple_somatic_mutation.open.COAD-US.tsv_preprocessed.csv",
                "lusc":rootdir + "LUSC/patient_not_matched/mut/simple_somatic_mutation.open.LUSC-US.tsv_preprocessed.csv",
                "prad":rootdir + "PRAD/patient_not_matched/mut/simple_somatic_mutation.open.PRAD-US.tsv_preprocessed.csv",
                "stad":rootdir + "STAD/patient_not_matched/mut/simple_somatic_mutation.open.STAD-US.tsv_preprocessed.csv"}

exp_src_data = {"brca":rootdir + "BRCA/patient_not_matched/exp/exp_seq.BRCA-US.tsv_preprocessed.csv",
                "coad":rootdir + "COAD/patient_not_matched/exp/exp_seq.COAD-US.tsv_preprocessed.csv",
                "lusc":rootdir + "LUSC/patient_not_matched/exp/exp_seq.LUSC-US.tsv_preprocessed.csv",
                "prad":rootdir + "PRAD/patient_not_matched/exp/exp_seq.PRAD-US.tsv_preprocessed.csv",
                "stad":rootdir + "STAD/patient_not_matched/exp/exp_seq.STAD-US.tsv_preprocessed.csv"}

mut_data = {"brca":rootdir + "BRCA/for_test/mut.BRCA-US.tsv_preprocessed.csv",
            "coad":rootdir + "COAD/for_test/mut.COAD-US.tsv_preprocessed.csv",
            "lusc":rootdir + "LUSC/for_test/mut.LUSC-US.tsv_preprocessed.csv",
            "prad":rootdir + "PRAD/for_test/mut.PRAD-US.tsv_preprocessed.csv",
            "stad":rootdir + "STAD/for_test/mut.STAD-US.tsv_preprocessed.csv"}

exp_data = {"brca":rootdir + "BRCA/for_test/exp.BRCA-US.tsv_preprocessed.csv",
            "coad":rootdir + "COAD/for_test/exp.COAD-US.tsv_preprocessed.csv",
            "lusc":rootdir + "LUSC/for_test/exp.LUSC-US.tsv_preprocessed.csv",
            "prad":rootdir + "PRAD/for_test/exp.PRAD-US.tsv_preprocessed.csv",
            "stad":rootdir + "STAD/for_test/exp.STAD-US.tsv_preprocessed.csv"}

smoothed_dir = {"brca":rootdir + "BRCA/for_test/smoothed/",
                "coad":rootdir + "COAD/for_test/smoothed/",
                "lusc":rootdir + "LUSC/for_test/smoothed/",
                "prad":rootdir + "PRAD/for_test/smoothed/",
                "stad":rootdir + "STAD/for_test/smoothed/"}