

library(TCGAbiolinks)
library(IlluminaHumanMethylation450kmanifest)

library(missMethyl)
library(edgeR)

setwd(Sys.getenv("LUSC_DATA_DIR", unset = "."))
# loading saved session: Once you saved your data, you can load it into your environment: 
data.met = readRDS(file = "data.met.RDS")
# met matrix
met <- as.data.frame(SummarizedExperiment::assay(data.met))
# clinical data
clinical <- data.frame(data.met@colData)

#___________inspectiong methylation data_______________#

# get the 450k annotation data
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## remove probes with NA
# probe.na <- rowSums(is.na(met)
probe.na <- is.na(rowSums(met))
table(probe.na == 0)
# FALSE   TRUE 
# 110772 374805 

# chose those has not NA values in rows
probe <- probe.na[probe.na == 0]
met <- met[row.names(met) %in% names(probe), ]

## remove probes that match to chromosome  X and Y 
keep <- !(row.names(met) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep)
met <- met[keep, ]
rm(keep) # remove no further needed probes.

## remove SNPs overlapped probe
table(is.na(ann450k$Probe_rs))
# probes without snp
no.snp.probe <- ann450k$Name[is.na(ann450k$Probe_rs)]

snp.probe <- ann450k[!is.na(ann450k$Probe_rs), ]
#snps with maf <= 0.05
snp5.probe <- snp.probe$Name[snp.probe$Probe_maf <= 0.05]

# filter met
met <- met[row.names(met) %in% c(no.snp.probe, snp5.probe), ]

#remove no-further needed dataset
rm(no.snp.probe, probe, probe.na, snp.probe, snp5.probe)

## Removing probes that have been demonstrated to map to multiple places in the genome.
# list adapted from https://www.tandfonline.com/doi/full/10.4161/epi.23470

crs.reac <- read.csv(Sys.getenv("CROSS_REACTIVE_PROBES_FILE", unset = "cross_reactive_probe.chen2013.csv"))
crs.reac <- crs.reac$TargetID[-1]

# filtre met
met <- met[ -which(row.names(met) %in% crs.reac), ]
bval <- met

## converting beta values to m_values
## m = log2(beta/1-beta)
mval <- t(apply(met, 1, function(x) log2(x/(1-x))))
#______________saving/loading_____________________#
# save data sets
saveRDS(mval, file = "mval.RDS", compress = FALSE)
saveRDS (bval, file = "bval.RDS", compress = FALSE)
