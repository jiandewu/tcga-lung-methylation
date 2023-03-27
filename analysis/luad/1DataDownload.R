#______________loading packages__________________#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
# ERROR: dependencies ‘EDASeq’, ‘genefilter’, ‘sva’ are not available for package ‘TCGAbiolinks’
BiocManager::install("Bioconductor")
BiocManager::install("sva")
BiocManager::install("genefilter")
BiocManager::install("EDASeq")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("missMethyl")
BiocManager::install("artyfarty")




#increase memory size
#memory.limit(size = 28000) # do this on your own risk!

setwd(Sys.getenv("LUAD_DATA_DIR", unset = "."))

#_________ DNA methylation Download____________#
# DNA methylation aligned to hg19

library(TCGAbiolinks)

query_met <- GDCquery(project= "TCGA-LUAD", 
                      data.category = "DNA methylation", 
                      platform = "Illumina Human Methylation 450", 
                      legacy = TRUE)
GDCdownload(query_met)
#putting files togathers
data.met <- GDCprepare(query_met)
#saving the met object
saveRDS(object = data.met,
        file = "data.met.RDS",
        compress = FALSE)
