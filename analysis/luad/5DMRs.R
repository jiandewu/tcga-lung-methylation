library(RColorBrewer)
library(Gviz)


setwd(Sys.getenv("LUSC_DATA_DIR", unset = "."))
mval <- readRDS("mval_filtered.RDS")
bval <- readRDS("bval_filtered.RDS")
clinical = readRDS(file = "clinical.RDS")

design <- model.matrix(~0+sample_type, data = clinical)
colnames(design) <- gsub("sample_type", "", colnames(design))
library(limma)
contr.matrix <- makeContrasts(
  # This will test for whether the log-expression of group Lean is equal to 
  # the average of groups ObeseI, ObeseII and ObeseIII.
  Normal.vs.Tumor = Normal-Tumor, 
  levels = colnames(design))
contr.matrix

# setting some annotation


library(DMRcate)
myAnnotation <- cpg.annotate(object = mval, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contr.matrix, 
                             coef = "Normal.vs.Tumor", arraytype = "450K")
str(myAnnotation)

#endif /* NEWSTUFF */
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges
write.table(results.ranges, file = "DMR_genes2.txt", sep="\t", quote = F, row.names=T)

# set up the grouping variables and colours
#install.packages('riskyr')
library(riskyr)
groups <- pal[1:length(unique(clinical$sample_type))]
names(groups) <- levels(factor(clinical$sample_type))
cols <- groups[as.character(factor(clinical$sample_type))]

# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 2, CpGs = bval, phen.col = cols, what = "Beta", arraytype = "450K", genome = "hg19")
