library(minfi)

setwd(Sys.getenv("LUAD_DATA_DIR", unset = "."))
# loading saved session: Once you saved your data, you can load it into your environment: 
data.met = readRDS(file = "data.met.RDS")
mval <- readRDS("mval.RDS")
bval <- readRDS("bval.RDS")

# clinical data
clinical <- data.frame(data.met@colData)

table(clinical$sample_type, clinical$race)

# filtering and grouping
clinical <- clinical[, c("barcode", "sample_type", "race")]
clinical <- na.omit(clinical)
#clinical <- clinical[-which(clinical$sample_type == "Metastatic"), ]
clinical <- clinical[-which(clinical$sample_type == "Recurrent Tumor"), ]
barcode <- clinical$barcode

# removing samples from meth matrixes
bval <- bval[, colnames(bval) %in% barcode]
mval <- mval[, colnames(mval) %in% barcode]
saveRDS(mval, file = "mval_filtered.RDS", compress = FALSE)
saveRDS(bval, file = "bval_filtered.RDS", compress = FALSE)

# Making sure about samples in clinical and matrixes and their order
table(colnames(mval) %in% row.names(clinical))
table(colnames(bval) %in% row.names(clinical))
#
all(row.names(clinical) == colnames(bval))
all(row.names(clinical) == colnames(mval))

clinical$sample_type <- gsub("Solid Tissue Normal","Normal",clinical$sample_type)
clinical$sample_type <- gsub("Primary Tumor","Tumor",clinical$sample_type)

#Making grouping variable
clinical$sample_type <- as.factor(clinical$sample_type)
#levels(clinical$sample_type)
clinical$sample_type <- relevel(clinical$sample_type, ref = "Normal")

saveRDS(clinical, file = "clinical.RDS", compress = FALSE)

#_____________ DMC analysis________________#
design <- model.matrix(~0+sample_type, data = clinical)
colnames(design) <- gsub("sample_type", "", colnames(design))

library(limma)

contr.matrix <- makeContrasts(
  # This will test for whether the log-expression of group Lean is equal to 
  # the average of groups ObeseI, ObeseII and ObeseIII.
  Normal.vs.Tumor = Normal-Tumor, 
  levels = colnames(design))
contr.matrix

# fit the linear model 
fit <- lmFit(mval, design)
fit <- contrasts.fit(fit, contrasts=contr.matrix)
fit2 <- eBayes(fit)
summary(decideTests(fit2))


# extracting significantly methylated probes
#deff.meth = topTable(fit2, coef=ncol(design), sort.by="p",number = nrow(mval), adjust.method = "BY")
#write.table(deff.meth,file="DMC_sig_probes.txt",sep="\t",row.names=T,quote=F)

# get the 450k annotation data
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(mval),ann450k$Name),c(1:4,12:19,24:ncol(ann450k))]

deff.meth <- topTreat(fit2, coef=1, n=Inf, genelist=ann450kSub)
write.table(deff.meth[deff.meth$adj.P.Val<0.05,], file = "DMC_sig_probes.txt", sep="\t", quote = F, row.names=T)

# Visualization
# plot the top 10 most significantly differentially methylated CpGs 
par(mfrow=c(2,5))
sapply(rownames(deff.meth)[1:10], function(cpg){
  plotCpg(bval, cpg=cpg, pheno=clinical$sample_type, ylab = "Beta values")
})

# making a volcano plot
# #making dataset
# dat <- data.frame(foldchange = fit[["coefficients"]][,2], logPvalue =  -log10(fit2[["p.value"]][,2]))
# dat$threshold <- as.factor(abs(dat$foldchange) < 2.0)
# 
# #Visualization
# cols <- c("TRUE" = "grey", "FALSE" = "blue")
# ggplot(data=dat, aes(x=foldchange, y = logPvalue, color=threshold)) +
#   geom_point(alpha=.6, size=1.2) +
#   scale_colour_manual(values = cols) +
#   geom_vline(xintercept = 0.4, colour="#990000", linetype="dashed") + 
#   geom_vline(xintercept = - 0.4, colour="#990000", linetype="dashed") +
#   theme(legend.position="none") +
#   xlab("Fold Change") +
#   ylab("-log10 p value") +
#   theme_bw() +
#   theme(legend.position = "none")

logFC.cutoff <- 2
adj.P.Val.cutoff <- 0.05
Volcano.color <- as.factor(ifelse(deff.meth$adj.P.Val<adj.P.Val.cutoff & 
                                    abs(deff.meth$logFC)>=logFC.cutoff,
                                  ifelse(deff.meth$logFC > logFC.cutoff,'Up','Down'),'Not'))

library(ggplot2)

p <- ggplot(data=deff.meth, aes(x=logFC, y=-log10(adj.P.Val), col=Volcano.color)) +
  geom_point(alpha=0.8, size = 1.8) +
  theme_bw(base_size = 15) +
  labs(x="log2 (fold change)",y="-log10 (adj.p-value)",title = "Normal vs. Tumor") +
  geom_hline(yintercept = -log10(adj.P.Val.cutoff), lty=4,col="black",lwd=0.6) +
  geom_vline(xintercept = c(-logFC.cutoff, logFC.cutoff), lty=4,col="black",lwd=0.6) +
  scale_color_manual(name = "", values = c("red", "green", "gray"), limits = c("Up", "Down", "Not")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
plot(p)

# Heatmap draw
library(gplots)
up.gene <- row.names(deff.meth[deff.meth$adj.P.Val<0.001&deff.meth$logFC>2,])
down.gene <- row.names(deff.meth[deff.meth$adj.P.Val<0.001&deff.meth$logFC<(-2),])
# Lean.vs.ObeseIII.topgenes <- row.names(Lean.vs.ObeseIII)[1:length(de.ObeseIII)]
Ctrl.vs.Treat.topgenes <- c(up.gene,down.gene)
i <- which(row.names(mval) %in% Ctrl.vs.Treat.topgenes)
j <- c(which(clinical$sample_type=="Normal"), which(clinical$sample_type=="Tumor"))
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(data.matrix(bval[i,j]), scale="row", Rowv = TRUE, Colv=FALSE, dendrogram="none",
          key = TRUE, keysize = 1.5, density.info="none",
          labRow=row.names(mval)[i], labCol=rep(c("N","T"),c(sum(clinical$sample_type=="Normal"),sum(clinical$sample_type=="Tumor"))), 
          col=mycol, trace="none", 
          margin=c(5,10), lhei=c(2,10), main = "Normal vs Tumor")
