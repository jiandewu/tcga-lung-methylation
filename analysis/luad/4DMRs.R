

setwd(Sys.getenv("LUAD_DATA_DIR", unset = "."))
mval <- readRDS("mval_filtered.RDS")
bval <- readRDS("bval_filtered.RDS")
clinical = readRDS(file = "clinical.RDS")

design <- model.matrix(~sample_type, data = clinical)
colnames(design) <- gsub("sample_type", "", colnames(design))

#BiocManager::install("DMRcate")
library(DMRcate)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# setting some annotation
myAnnotation <- cpg.annotate(object = mval, datatype = "array", 
                             what = "M", 
                             analysis.type = "differential", 
                             design = design, 
                             contrasts = FALSE, 
                             coef = "Tumor", 
                             arraytype = "450K",
                             fdr = 0.001)
str(myAnnotation)

# DMR analysis
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges
write.table(results.ranges, file = "DMR_genes.txt", sep="\t", quote = F, row.names=T)

# visualization
dmr.table <- data.frame(results.ranges)

# setting up a variable for grouping and color

library(RColorBrewer)

pal <- brewer.pal(8,"Dark2")
groups <- pal[1:length(unique(clinical$sample_type))]
names(groups) <- levels(factor(clinical$sample_type))

#setting up the genomic region 
gen <- "hg19"
# the index of the DMR that we will plot 
dmrIndex <- 2
# coordinates are stored under results.ranges[dmrIndex]

chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <- as.numeric(start(results.ranges[dmrIndex]))
end <- as.numeric(end(results.ranges[dmrIndex]))

# add 25% extra space to plot
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))


# defining CpG islands track
# download cpgislands for chromosome number 6 from ucsc
chr.cpg <- read.table(Sys.getenv("CPG_ISLAND_BED", unset = "cpgIslandExt.hg19.bed"),header = TRUE)
chr.cpg <- chr.cpg[chr.cpg$chrom==chrom,]

islandData <- GRanges(seqnames=Rle(chr.cpg[,1]), 
                      ranges=IRanges(start=chr.cpg[,2],
                                     end=chr.cpg[,3]),
                      strand=Rle(strand(rep("*",nrow(chr.cpg)))))

# DNAseI hypersensitive sites track
#downloaded from ucsc
chr.dnase <- read.table(Sys.getenv("DNASE_SITES_BED", unset = "DnaseSites.bed"),header = TRUE)
chr.dnase <- chr.dnase[chr.dnase$chrom==chrom,]

dnaseData <- GRanges(seqnames=chr.dnase[,1],
                     ranges=IRanges(start=chr.dnase[,2], end=chr.dnase[,3]),
                     strand=Rle(rep("*",nrow(chr.dnase))),
                     data=chr.dnase[,5])

#Setting up the ideogram, genome, and RefSeq tracks 
library(Gviz)

iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name=paste0(chrom))
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
                    from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                    rstarts="exonStarts", rends="exonEnds", gene="name", 
                    symbol="name2", transcript="name", strand="strand", 
                    fill="darkblue",stacking="squish", name="RefSeq", 
                    showId=TRUE, geneSymbol=TRUE)

#Ensure that the methylation data is ordered by chromosome and base position.

# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kOrd <- ann450k[order(ann450k$chr,ann450k$pos),]
bvalOrd <- bval[match(ann450kOrd$Name,rownames(bval)),]

#Create the data tracks:
# create genomic ranges object from methylation data
cpgData <- GRanges(seqnames=Rle(ann450kOrd$chr),
                   ranges=IRanges(start=ann450kOrd$pos, end=ann450kOrd$pos),
                   strand=Rle(rep("*",nrow(ann450kOrd))),
                   betas=bvalOrd)

# methylation data track
methTrack <- DataTrack(range=cpgData, 
                       groups=clinical$sample_type, # change this if your groups are diffrent
                       genome = gen,
                       chromosome=chrom,
                       ylim=c(-0.05,1.05),
                       col=pal,
                       type=c("a","p"), 
                       name="DNA Meth.\n(beta value)",
                       background.panel="white", 
                       legend=TRUE, 
                       cex.title=0.8,
                       cex.axis=0.8, 
                       cex.legend=0.8)

# CpG island track
islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG Is.", 
                               chromosome=chrom,fill="darkgreen")

# DNaseI hypersensitive site data track
dnaseTrack <- DataTrack(range=dnaseData, genome=gen, name="DNAseI", 
                        type="gradient", chromosome=chrom)

# DMR position data track
dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                            chromosome=chrom,fill="darkred")


# Set up the tracklist and indicate the relative sizes of the different tracks. 
# Finally, draw the plot using the plotTracks function
tracks <- list(iTrack, gTrack, methTrack, dmrTrack, islandTrack, dnaseTrack,
               rTrack)
sizes <- c(2,2,5,2,2,2,3) # set up the relative sizes of the tracks
# to save figure and scaling graphic device
tiff(filename = "dmr.tiff", width = 15, height = 10, units = "in", res = 400)
plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))
dev.off()


myAnnotation <- cpg.annotate(object = mval, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "naive - rTreg", arraytype = "450K")