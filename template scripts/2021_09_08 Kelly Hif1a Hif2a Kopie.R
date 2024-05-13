update.packages()

BiocManager::install("tximport")
n
BiocManager::install("EnsDb.Hsapiens.v101")
n
BiocManager::install("tximeta")
n
BiocManager::install("Rcpp")
n
BiocManager::install("dplyr")
n
BiocManager::install("htmltools")
BiocManager::install("httpuv")
n
BiocManager::install("AnnotationHub")
n
BiocManager::install("biomaRt")
n
BiocManager::install("apeglm")
n
BiocManager::install("ashr")
n
BiocManager::install("ReportingTools")
n
BiocManager::install("org.Hs.eg.db")

BiocManager::install("pheatmap")
n
BiocManager::install("ReportingTools")
n
BiocManager::install("Glimma")
n
BiocManager::install("edgeR")
n
BiocManager::install("pathview")

BiocManager::install("AnnotationHub")
n
BiocManager::install("EnhancedVolcano")
n
BiocManager::install("goseq")
n
BiocManager::install("airway")
n
BiocManager::install("EnrichmentBrowser")
n
BiocManager::install("regioneR")
n
BiocManager::install("ALL")
n
BiocManager::install("regioneR")
n
BiocManager::install("GO.db")
n
BiocManager::install("clusterProfiler")
n
BiocManager::install("pathview")
n
BiocManager::install("enrichplot")
n
BiocManager::install("ComplexHeatmap")
n
BiocManager::install("circlize")
n
BiocManager::install("tidyr")
n
BiocManager::install("data.table")
n
BiocManager::install()


BiocManager::version()

install.packages("stringr")
install.packages("R.utils")
install.packages("ggplot2")
install.packages("ALL")


library(AnnotationHub)
library(limma)
library(Glimma)
library(edgeR)
# library(Mus.musculus)\n
# library(Homo.sapiens)\n
# library(RColorBrewer)\n
# library(RNAseq123)\n
# library(xlsx)\n
# library(writexl)\n
library(RColorBrewer)
library(EnsDb.Hsapiens.v86)
library(tximport)
library(tximeta)
library(stringr)
library(gplots)
library(R.utils)
library(biomaRt)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(DESeq2)
library(apeglm)
library(ashr)
library("ggplot2")
library(ReportingTools)
library(pheatmap)
library(tidyverse)
library(xlsx)
library(writexl)
library("org.Hs.eg.db")

library(EnhancedVolcano)

library(pathview)
library(airway)
library(goseq)
library(airway)
library(regioneR)
library(EnrichmentBrowser)
library(ALL)
data(ALL)
require(DOSE) 
require(ComplexHeatmap)
library(DOSE)

library(regioneR)
library(EnrichmentBrowser)
library(airway)
library("GO.db")

library()

# get error info: "options(error=recover)"; options(error=NULL)
options(error=NULL)

###############################################################################
# 19.01.2021
# tximeta
###############################################################################

browseVignettes("tximeta")
# http://127.0.0.1:23876/library/tximeta/doc/tximeta.html\n
# use list of dirs of quants = files


###############################################################################
# Prepare files:
###############################################################################

setwd("C:/Users/kelterbs/OneDrive - Charit? - Universit?tsmedizin Berlin/RNA-seq/Kelly CRISPR HIF KO")
# mac setwd("~/OneDrive - Charité - Universitätsmedizin Berlin/RNA-seq/Kelly CRISPR HIF KO")
getwd()
# first file:
example.quant <- "./quants/CH_HS_KK_077_S33_quant/quant.sf"
example.table <- read.table(example.quant, header = T)
head(example.table)
basename(example.quant)
dirname(example.quant)
basename(dirname(example.quant))

# generate file list
files = {
}
for (i in list.dirs(path = dirname(dirname(example.quant)), full.names = TRUE, recursive = FALSE)) {
	files <- c(files, paste(i, "/quant.sf", sep = ""))
	print(basename(i))
	print(head(read.table(files[length(files)], header = T)))
	print(files[length(files)])
}

# remove unwanted folders/files
files
files <- files[1:16]
files

# Sample names:
# define Samples
group <- factor(rep(c("LV_Nx", "LV_Hx", "Hif1a_KO", "Hif2a_KO"), 4),level=c("LV_Nx", "LV_Hx", "Hif1a_KO", "Hif2a_KO"))
group
summary(group)

Passage <- factor(c(rep(c(rep("P7", 2), rep("P13", 2)), 2), rep(c(rep("P8", 2), rep("P14",2)), 2)),levels = c("P7", "P8", "P13", "P14"))
summary(Passage)
Passage

batch <- factor(c(rep(1,8), rep(2,8)))
replicate <- factor(c(rep(1:4,2), rep(5:8,2)))

rep <- {}
for (i in 1:4){
  print(i)
  rep <- c(rep,rep(i,4))
}
repetition <- factor(rep)

# define Column names
basename(dirname(files))

col.names <- {
}
for (i in 33:48) {
	print(i)
	Si <- paste("Kelly_", i + 44, "_S", sprintf("%02d", i), sep = "")
	print(Si)
	col.names <- c(col.names, Si)
	print(col.names)
}
length(col.names)
samples <- col.names

file.exists(files)
# names as working sample names\n
coldata <- data.frame(files, fname = basename(dirname(files)), names = samples, condition = group, Passage = Passage, batch = batch, replicate = replicate, repetition = repetition, stringsAsFactors = FALSE)
coldata

###############################################################################
# Load into tximeta
###############################################################################

se <- tximeta(coldata)
se

colData(se)
assayNames(se)
rowRanges(se)
seqinfo(se)
rownames(se)

gse <- se
head(assays(gse)[["counts"]])
head(assays(gse)[["abundance"]])
head(assays(gse)[["length"]])

gse <- summarizeToGene(se)
rowRanges(gse)
length(rownames(mcols(gse)["gene_id"]))


col.names <- colnames(mcols(gse))
col.names_o <- col.names
# col.names <- c(col.names[1:3],"SYMBOL_list",col.names[5])
colnames(mcols(gse)) <- col.names

# add other metadata
columns(retrieveDb(se))
txidb <- retrieveDb(se)
columns(txidb)

columns(org.Hs.eg.db)
gse <- addIds(gse, "REFSEQ", gene=TRUE, multiVals = "list")
rowRanges(gse)
gse <- addIds(gse, "ENTREZID")
rowRanges(gse)
gse <- addIds(gse, "SYMBOL")
colnames(mcols(gse)["SYMBOL"])<-"Symbol_list"
gse <- addIds(gse, "ENSEMBL")

### change rownames

# add from ENSEMBLE
# library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
ah = AnnotationHub()
query(hub, c("EnsDb", "sapiens"))
edb <- ah[["AH95744"]]

# edb <- EnsDb.Hsapiens.v86
columns(edb)


mapIds(edb, keys = "ENSG00000000005", column = "GENENAME", keytype = "GENEID")
mapIds(edb, keys = "ENSG00000000005", column = "SYMBOL", keytype = "GENEID")
mapIds(edb, keys = "ENSG00000000005", column = "SEQNAME", keytype = "GENEID")
mapIds(edb, keys = "ENSG00000000005", column = "TXNAME", keytype = "GENEID")
mapIds(edb, keys = "ENSG00000000005", column = "UNIPROTID", keytype = "GENEID")

mcols(gse) <- cbind(mcols(gse),
                    Symbol_ens=mapIds(edb, keys = mcols(gse)[,"ENSEMBL"], column = "SYMBOL", keytype = "GENEID"))
mcols(gse) <- cbind(mcols(gse),
                    Symbol_ens=mapIds(edb, keys = mcols(gse)[,"ENSEMBL"], column = "UNIPROTID", keytype = "GENEID"))

rowRanges(gse)[,"SYMBOL"]
length(mcols(gse)[,"SYMBOL"])
sum(is.na(mcols(gse)[,"SYMBOL"]))
sum(is.na(mcols(gse)[,"GENENAME"]))
sum(is.na(mcols(gse)[,"SYMBOL_entrez"]))

###########################################
# save/load gse file

saveRDS(gse, file = "C:/Users/kelterbs/OneDrive - Charit? - Universit?tsmedizin Berlin/RNA-seq/Kelly CRISPR HIF KO/2021_09_09_gse_txi.RDS") 
gse <- readRDS("2021_09_09_gse_txi.RDS")
list.files()

# gse <- gse_load

# factor(goi$Genes, levels = goi$Genes)

group <- factor(rep(c("LV_Nx", "LV_Hx", "Hif1a_KO", "Hif2a_KO"), 4),level=c("LV_Nx", "LV_Hx", "Hif1a_KO", "Hif2a_KO"))
group
summary(group)

Passage <- factor(c(rep(c(rep("P7", 2), rep("P13", 2)), 2), rep(c(rep("P8", 2), rep("P14",2)), 2)),levels = c("P7", "P8", "P13", "P14"))
summary(Passage)
Passage

batch <- factor(c(rep(1,8), rep(2,8)))
replicate <- factor(c(rep(1:4,2), rep(5:8,2)))

rep <- {}
for (i in 1:4){
  print(i)
  rep <- c(rep,rep(i,4))
}
repetition <- factor(rep)

##########################################\n
# DESeq2 Tutorial\n
##########################################\n
# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# condition already loaded above
colData(gse)$condition <- group
colData(gse)$Passage <- Passage
colData(gse)$batch <- batch
colData(gse)$replicate <- replicate
colData(gse)$repetition <- repetition
samples <- rownames(colData(gse))

colData(gse)

# setwd(xy)
getwd()
# saveRDS(gse, file = "2021_09_24_gse_txi.RDS") 
gse <- readRDS("2021_09_24_gse_txi.RDS")


# Load gse file for DESeq2
dds <- DESeqDataSet(gse, ~condition)

##########################################\n
# filter all rows with rowsum = 0

keep0 <- rowSums(counts(dds) == 0) == length(colData(dds)[, 1])
head(keep0)
summary(keep0)
dds <- dds[!keep0,]
dds
head(assays(dds)[["counts"]])
length(rownames(counts(dds)))

# or filter all lower as 5 counts
dds <- estimateSizeFactors(dds)

head(counts(dds))
counts(dds)["ENSG00000000003.15",]
keep5 <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 4  # List of genes with more or equal 4 samples, that show norm. counts greater than 5
head(keep5)
summary(keep5)

dds
head(assays(dds)[["counts"]])
length(rownames(counts(dds)))

dds <- dds[keep5,]
# dds.k0 <- dds

samplenames <- samples
samplenames
sample <- {}
for (i in c(1:4, 9:12)) {
	s <- paste(group[i], "_", Passage[i], sep = "")
	sample <- c(sample, s)
}
sample

sample.16 <- {}
for (i in c(1:16)) {
	s <- paste(group[i], "_", Passage[i], sep = "")
	sample.16 <- c(sample.16, s)
}
sample.16

dds$sample <- factor(sample.16, levels = c("LV_Nx_P7","LV_Nx_P8","LV_Hx_P7","LV_Hx_P8","Hif1a_KO_P13","Hif1a_KO_P14","Hif2a_KO_P13","Hif2a_KO_P14"))


## all counts, without replicates...
# dds <- DESeq(dds)
# saveRDS(dds, file = "C:/Users/kelterbs/OneDrive - Charit? - Universit?tsmedizin Berlin/RNA-seq/Kelly CRISPR HIF KO/2021_09_14_dds_allcounts.RDS") 


dds$run <- samplenames
dds <- collapseReplicates(dds, dds$sample,dds$run)
colData(dds)
dds$runsCollapsed
colnames(dds)
matchFirstLevel <- dds$sample == levels(dds$sample)[1]
stopifnot(all(rowSums(counts(dds[,matchFirstLevel])) == counts(dds[,1])))

# run DESeq ############
dds <- DESeq(dds)

# creating results
res <- results(dds)
res
dim(res)
plotMA(res)

head(mcols(dds))
head(mcols(dds)["SYMBOL"])
head(assays(dds)[["counts"]])

# Plot overview (first "result": LV Nx vs Hif1a KO)
par(mfrow=c(1,1))
plotMA(res)
title(main="Hif2a KO vs LV Nx")

# gene-wise comparision:
res.LV_Nx.vs.LV_Hx <- results(dds, contrast= c("condition","LV_Hx","LV_Nx"))
res.LV_Nx.vs.LV_Hx$gene_id <- mcols(dds)["gene_id"]
res.LV_Nx.vs.LV_Hx$SYMBOL <- mcols(dds)["SYMBOL"]
res.LV_Nx.vs.LV_Hx$ENSEMBL <- mcols(dds)["ENSEMBL"]

par(mfrow=c(1,1))
plotMA(res.LV_Nx.vs.LV_Hx)
title(main="Ctrl: Hx vs. Nx")

# LV_Hx vs. Hif1a_KO
res.LV_Hx_vs_Hif1a_KO <- results(dds, contrast= c("condition","Hif1a_KO","LV_Hx"))
res.LV_Hx_vs_Hif1a_KO$gene_id <- mcols(dds)["gene_id"]
res.LV_Hx_vs_Hif1a_KO$SYMBOL <- mcols(dds)["SYMBOL"]

plotMA(res.LV_Hx_vs_Hif1a_KO)
title(main="Hif1a vs. Ctrl (Hx)")

# LV_Hx vs. Hif2a_KO
res.LV_Hx_vs_Hif2a_KO <- results(dds, contrast= c("condition","Hif2a_KO","LV_Hx"))
res.LV_Hx_vs_Hif2a_KO$gene_id <- mcols(dds)["gene_id"]
res.LV_Hx_vs_Hif2a_KO$SYMBOL <- mcols(dds)["SYMBOL"]

plotMA(res.LV_Hx_vs_Hif2a_KO)
title(main="Hif2a vs. Ctrl (Hx)")

### all
# pdf("2021_09_08 MA plots2.pdf") 
par(mfrow=c(1,3))
plotMA(res.LV_Hx_vs_Hif1a_KO, colSig = "#BF6D39", ylim = c(-6,6))
title(main="HIF-1a vs. Ctrl (Hx)")
plotMA(res.LV_Nx.vs.LV_Hx, colSig = "#0066CC", ylim = c(-6,6))
title(main="Hypoxia vs. Normoxia")
plotMA(res.LV_Hx_vs_Hif2a_KO, colSig = "#00CC99", ylim = c(-6,6))
title(main="HIF-2a vs. Ctrl (Hx)")
dev.off()
resultsNames(dds)

# Shrinkage

# ashr
# install ashr
library(AnnotationHub)
ah = AnnotationHub()
query(ah, c("EnsDb", "sapiens"))
edb <- ah[["AH95744"]]


columns(edb)
resLFC.LV_Nx.vs.LV.Hx <- lfcShrink(dds, contrast= c("condition","LV_Hx","LV_Nx"), type="ashr")
resLFC.LV_Nx.vs.LV.Hx$gene_id <- mcols(dds)[["gene_id"]]
resLFC.LV_Nx.vs.LV.Hx$SYMBOL <- mcols(dds)[["SYMBOL"]]
resLFC.LV_Nx.vs.LV.Hx$ENTREZ <- mcols(dds)[["ENTREZID"]]
resLFC.LV_Nx.vs.LV.Hx$ENTREZ2 <- as.character(mapIds(edb,keys = rownames(resLFC.LV_Nx.vs.LV.Hx), column = "ENTREZID", keytype = "GENEID", multiVals="first"))
resLFC.LV_Nx.vs.LV.Hx$SYMBOL2 <- as.character(mapIds(edb,keys = rownames(resLFC.LV_Nx.vs.LV.Hx), column = "SYMBOL", keytype = "GENEID", multiVals="first"))
resLFC.LV_Nx.vs.LV.Hx$DESCRIPTION <- as.character(mapIds(edb,keys = rownames(resLFC.LV_Nx.vs.LV.Hx), column = "DESCRIPTION", keytype = "GENEID", multiVals="first"))

as.character(mapIds(edb,keys = rownames(resLFC.LV_Nx.vs.LV.Hx)[1], column="UNIPROTMAPPINGTYPE", keytype = "GENEID", multiVals="first"))

summary(mcols(dds)[["ENTREZID"]])
summary(resLFC.LV_Nx.vs.LV.Hx$ENTREZ)

summary(resLFC.LV_Nx.vs.LV.Hx$ENTREZ==resLFC.LV_Nx.vs.LV.Hx$ENTREZ2)

resLFC.LV_Hx.vs.Hif1a_KO <- lfcShrink(dds, contrast= c("condition","Hif1a_KO","LV_Hx"), type="ashr")
resLFC.LV_Hx.vs.Hif1a_KO$gene_id <- mcols(dds)[["gene_id"]]
resLFC.LV_Hx.vs.Hif1a_KO$SYMBOL <- mcols(dds)[["SYMBOL"]]
resLFC.LV_Hx.vs.Hif1a_KO$ENTREZ <- mcols(dds)[["ENTREZID"]]
resLFC.LV_Hx.vs.Hif1a_KO$ENTREZ2 <- as.character(mapIds(edb,keys = rownames(resLFC.LV_Hx.vs.Hif1a_KO), column = "ENTREZID", keytype = "GENEID", multiVals="first"))
resLFC.LV_Hx.vs.Hif1a_KO$SYMBOL2 <- as.character(mapIds(edb,keys = rownames(resLFC.LV_Hx.vs.Hif1a_KO), column = "SYMBOL", keytype = "GENEID", multiVals="first"))
resLFC.LV_Hx.vs.Hif1a_KO$DESCRIPTION <- as.character(mapIds(edb,keys = rownames(resLFC.LV_Hx.vs.Hif1a_KO), column = "DESCRIPTION", keytype = "GENEID", multiVals="first"))


resLFC.LV_Hx.vs.Hif2a_KO <- lfcShrink(dds, contrast= c("condition","Hif2a_KO","LV_Hx"), type="ashr")
resLFC.LV_Hx.vs.Hif2a_KO$gene_id <- mcols(dds)[["gene_id"]]
resLFC.LV_Hx.vs.Hif2a_KO$SYMBOL <- mcols(dds)[["SYMBOL"]]
resLFC.LV_Hx.vs.Hif2a_KO$ENTREZ <- mcols(dds)[["ENTREZID"]]
resLFC.LV_Hx.vs.Hif2a_KO$ENTREZ2 <- as.character(mapIds(edb,keys = rownames(resLFC.LV_Hx.vs.Hif2a_KO), column = "ENTREZID", keytype = "GENEID", multiVals="first"))
resLFC.LV_Hx.vs.Hif2a_KO$SYMBOL2 <- as.character(mapIds(edb,keys = rownames(resLFC.LV_Hx.vs.Hif2a_KO), column = "SYMBOL", keytype = "GENEID", multiVals="first"))
resLFC.LV_Hx.vs.Hif2a_KO$DESCRIPTION <- as.character(mapIds(edb,keys = rownames(resLFC.LV_Hx.vs.Hif2a_KO), column = "DESCRIPTION", keytype = "GENEID", multiVals="first"))



resLFC.LV_Nx.vs.Hif1a_KO <- lfcShrink(dds, contrast= c("condition","Hif1a_KO","LV_Nx"), type="ashr")
resLFC.LV_Nx.vs.Hif1a_KO$gene_id <- mcols(dds)[["gene_id"]]
resLFC.LV_Nx.vs.Hif1a_KO$SYMBOL <- mcols(dds)[["SYMBOL"]]
resLFC.LV_Nx.vs.Hif1a_KO$ENTREZ <- mcols(dds)[["ENTREZID"]]
resLFC.LV_Nx.vs.Hif1a_KO$ENTREZ2 <- as.character(mapIds(edb,keys = rownames(resLFC.LV_Nx.vs.Hif1a_KO), column = "ENTREZID", keytype = "GENEID", multiVals="first"))

resLFC.LV_Nx.vs.Hif2a_KO <- lfcShrink(dds, contrast= c("condition","Hif2a_KO","LV_Nx"), type="ashr")
resLFC.LV_Nx.vs.Hif2a_KO$gene_id <- mcols(dds)[["gene_id"]]
resLFC.LV_Nx.vs.Hif2a_KO$SYMBOL <- mcols(dds)[["SYMBOL"]]
resLFC.LV_Nx.vs.Hif2a_KO$ENTREZ <- mcols(dds)[["ENTREZID"]]
resLFC.LV_Nx.vs.Hif2a_KO$ENTREZ2 <- as.character(mapIds(edb,keys = rownames(resLFC.LV_Nx.vs.Hif2a_KO), column = "ENTREZID", keytype = "GENEID", multiVals="first"))

resLV.Nx_Hx <- resLFC.LV_Nx.vs.LV.Hx
resHx.LV_Hif1a <- resLFC.LV_Hx.vs.Hif1a_KO
resHx.LV_Hif2a <- resLFC.LV_Hx.vs.Hif2a_KO
resNx_Hif1a <- resLFC.LV_Nx.vs.Hif1a_KO
resNx_Hif2a <- resLFC.LV_Nx.vs.Hif2a_KO


# pdf("2021_09_09 Normalization2.pdf") 
par(mfrow=c(2,3))

plotMA(res.LV_Nx.vs.LV_Hx, ylim=c(-6,6))
title(main="Ctrl: Hx vs. Nx")
plotMA(res.LV_Hx_vs_Hif1a_KO, ylim=c(-6,6))
title(main="Hif1a vs. Ctrl (Hx)")
plotMA(res.LV_Hx_vs_Hif2a_KO, ylim=c(-6,6))
title(main="Hif2a vs. Ctrl (Hx)")

plotMA(resLV.Nx_Hx, ylim=c(-6,6))
title(main="with ASHR shrinkage")
plotMA(resHx.LV_Hif1a, ylim=c(-6,6))
title(main="with ASHR shrinkages")
plotMA(resHx.LV_Hif2a, ylim=c(-6,6))
title(main="with ASHR shrinkages")
dev.off() 

par(mfrow=c(1,1))
plot(res.LV_Hx_vs_Hif1a_KO$log2FoldChange,resLV.Nx_Hx$log2FoldChange)


length()

resLV.Nx_Hx
resHx.LV_Hif1a
resHx.LV_Hif2a
resNx_Hif1a
resNx_Hif2a


summary(resLV.Nx_Hx)
summary(resHx.LV_Hif1a)
summary(resHx.LV_Hif2a)
summary(resNx_Hif1a)
summary(resNx_Hif2a)

resLV.Nx_Hx %>% subset(padj<0.05) %>% rownames() %>% length()
resHx.LV_Hif1a %>% subset(padj<0.05) %>% rownames() %>% length()
resHx.LV_Hif2a %>% subset(padj<0.05) %>% rownames() %>% length()
resNx_Hif1a %>% subset(padj<0.05) %>% rownames() %>% length()
resNx_Hif2a %>% subset(padj<0.05) %>% rownames() %>% length()




## Dispersion plot 
par(mfrow=c(1,1))
plotDispEsts(dds)

##############################################
# goi

library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)

columns(org.Hs.eg.db) # Entrez Gene identifiers
columns(EnsDb.Hsapiens.v86) # ENSEMBL Gene identifiers
columns(retrieveDb(se)) #  Taximeta database

# choose database:
edb <- EnsDb.Hsapiens.v86

# not in dataset: "NBPF23", "DDX4", "IL31RA", "CPZ", "MMP20" 

goi_N <- data.frame(Genes=c("PHOX2B", "ALK", "KIF1B", "CASC15", "NBAT1", "BARD1", "SDHB", "LMO1", "DUSP12", "HSD17B12", "HACE1", "LIN28B", "SPAG16", "EZH2", "CHEK2", "BARD1", "PALB2", "PINK1", "NEFL", "TP53", "HRAS", "PTPN11", "NF1", "APC", "BRCA2", "MLF1", "CDKN1B", "SMARCA4", "KIF15", "MYCN"))

goi <- data.frame(Genes=c("PHOX2B", "ALK", "KIF1B", "CASC15", "NBAT1", "BARD1", "SDHB", "LMO1", "DUSP12", "HSD17B12", "HACE1", "LIN28B", "SPAG16", "EZH2", "CHEK2", "BARD1", "PALB2", "PINK1", "NEFL", "TP53", "HRAS", "PTPN11", "NF1", "APC", "BRCA2", "MLF1", "CDKN1B", "SMARCA4", "KIF15", "MYCN"))
goi

goi$ENSEMBL <- mapIds(edb,
                      keys=goi$Genes,
                      column="GENEID",
                      keytype="SYMBOL",
                      multiVals="first")

# goi$ENSEMBL[7] <- "ENSG00000179571"

goi$GENENAME <- mapIds(edb,
                       keys=goi$Genes,
                       column="GENENAME",
                       keytype="SYMBOL",
                       multiVals="first")
goi
summary(goi$ENSEMBL %in% rownames(res))
setdiff(goi$ENSEMBL,rownames(res))

subset(goi, ENSEMBL %in% c("ENSG00000152670", "ENSG00000164509", "ENSG00000109625", "ENSG00000137674"))$Genes

goi_H <- data.frame(Genes = c("HIF1A", "EPAS1", "HIF3A", "ARNT"))
goi_H

goi_H$ENSEMBL <- mapIds(edb,
                      keys=goi_H$Genes,
                      column="GENEID",
                      keytype="SYMBOL",
                      multiVals="first")

goi_H$GENENAME <- mapIds(edb,
                       keys=goi_H$Genes,
                       column="GENENAME",
                       keytype="SYMBOL",
                       multiVals="first")

rownames(goi_H) <- goi_H$ENSEMBL
goi_H <- mutate(goi_H,GENENAME = str_replace(goi_H$GENENAME,pattern="EPAS1","HIF2A"))
goi_H <- mutate(goi_H,GENENAME = str_replace(goi_H$GENENAME,pattern="ARNT","HIF1B"))
goi_H

summary(goi_H$ENSEMBL %in% rownames(res))
setdiff(goi_H$ENSEMBL,rownames(res))
subset(goi_H, ENSEMBL %in% setdiff(goi_H$ENSEMBL,rownames(res)))$Genes
intersect(goi_H$ENSEMBL,rownames(res))


################# goi Lina

goi_L <- data.frame(Genes = c("VEGFA","SLC2A1", "CA9", "NTRK2","EGLN1","EGLN3","EDN1","FGFR4", "EPO", "WT1", "CXXC5"))

goi_L <- data.frame(Genes = c("WT1","CA9","SULF2"))

goi_L$ENSEMBL <- mapIds(edb,
                        keys=goi_L$Genes,
                        column="GENEID",
                        keytype="SYMBOL",
                        multiVals="first")

goi_L$GENENAME <- mapIds(edb,
                         keys=goi_L$Genes,
                         column="GENENAME",
                         keytype="SYMBOL",
                         multiVals="first")

rownames(goi_L) <- goi_L$ENSEMBL

summary(goi_L$ENSEMBL %in% rownames(res))
setdiff(goi_L$ENSEMBL,rownames(res))
subset(goi_L, ENSEMBL %in% setdiff(goi_L$ENSEMBL,rownames(res)))$Genes
intersect(goi_L$ENSEMBL,rownames(res))

res


# Plot counts
# Color
# dds file with all counts:
# saveRDS(dds, file = "C:/Users/kelterbs/OneDrive - Charit? - Universit?tsmedizin Berlin/RNA-seq/Kelly CRISPR HIF KO/2021_09_14_dds_allcounts.RDS") 
# dds <- readRDS("C:/Users/kelterbs/OneDrive - Charit? - Universit?tsmedizin Berlin/RNA-seq/Kelly CRISPR HIF KO/2021_09_14_dds_allcounts.RDS")



library(ggpubr)
dds$condition

plots <- {}
l <- length(goi$Genes)
n <- 1
for (n in 1:l){
  print(n)
  e <- goi$ENSEMBL[n]
  print(e)
  data <- plotCounts(dds, gene=e, returnData = TRUE)
  data$condition <- colData(dds)$condition
  i <- goi$GENENAME[n]
  print(i)
  y <- sym(i)
  g <- ggplot(data = data, aes(x = condition, y =count, fill = condition, colour=condition, xlab="", ylab="", las=2)) +
    # stat_summary(fun=mean, geom="point", shape=18, size=5) +
    # geom_violin(adjust = .3,alpha = .4,position=position_dodge(0.8)) +
    geom_boxplot(alpha = .1, position=position_dodge(0.8), colour="black") +
    # geom_line(group=c(
    #  subset(crat_data,time1 == levels(crat_data$time1)[1])[,"patient1"],
    #  subset(crat_data, time1 == levels(crat_data$time1)[2])[,"patient1"]),
    #  position=position_dodge(0.8),alpha = .3,size=1) +
    geom_dotplot(binaxis="y",stackdir='center', dotsize=1.3,position=position_dodge(0.8)) +
    xlab("")+
    ylab("expr.")+
    ggtitle(paste(i,sep=" "))
  g
  l <- length(goi$Genes)
  assign(paste("p",n,sep=""),g)
}

ggarrange(p1+rremove("x.text"),
          p2+rremove("x.text")+rremove("ylab"),
          p3+rremove("ylab")+rremove("x.text"),
          p4+rremove("ylab")+rremove("x.text"),
          p5+rremove("ylab")+rremove("x.text"),
          p6+rremove("ylab")+rremove("x.text"),
          p7+rremove("x.text"),
          p8+rremove("ylab")+rremove("x.text"),
          p9+rremove("ylab")+rremove("x.text"),
          p10+rremove("ylab")+rremove("x.text"),
          p11+rremove("ylab")+rremove("x.text"),
          p12+rremove("ylab")+rremove("x.text"),
          p13+rremove("x.text"),
          p14+rremove("ylab")+rremove("x.text"),
          p15+rremove("ylab")+rremove("x.text"),
          p16+rremove("ylab")+rremove("x.text"),
          p17+rremove("ylab")+rremove("x.text"),
          p18+rremove("ylab")+rremove("x.text"),
          p19+rremove("x.text"),
          p20+rremove("ylab")+rremove("x.text"),
          p21+rremove("ylab")+rremove("x.text"),
          p22+rremove("ylab")+rremove("x.text"),
          p23+rremove("ylab")+rremove("x.text"),
          p24+rremove("ylab")+rremove("x.text"),
          p25+ggpubr::rotate_x_text(),
          p26+rremove("ylab")+ggpubr::rotate_x_text(),
          p27+rremove("ylab")+ggpubr::rotate_x_text(),
          p28+rremove("ylab")+ggpubr::rotate_x_text(),
          p29+rremove("ylab")+ggpubr::rotate_x_text(),
          p30+rremove("ylab")+ggpubr::rotate_x_text(),
          common.legend = TRUE, legend = "right", nrow = 5, ncol= 6)



####################### Hifs
goi_H
intersect(goi_H$ENSEMBL,rownames(res))
plots <- {}
l <- length(intersect(goi_H$ENSEMBL,rownames(res)))
n <- 1

for (n in 1:l){
  print(n)
  e <- intersect(goi_H$ENSEMBL,rownames(res))[n]
  print(e)
  data <- plotCounts(dds, gene=e, returnData = TRUE)
  data$condition <- colData(dds)$condition
  i <- goi_H[e,"GENENAME"]
  print(i)
  y <- sym(i)
  g <- ggplot(data = data, aes(x = condition, y =count, fill = condition, colour=condition, xlab="", ylab="", las=2)) +
    # stat_summary(fun=mean, geom="point", shape=18, size=5) +
    # geom_violin(adjust = .3,alpha = .4,position=position_dodge(0.8)) +
    geom_boxplot(alpha = .1, position=position_dodge(0.8), colour="black") +
    # geom_line(group=c(
    #  subset(crat_data,time1 == levels(crat_data$time1)[1])[,"patient1"],
    #  subset(crat_data, time1 == levels(crat_data$time1)[2])[,"patient1"]),
    #  position=position_dodge(0.8),alpha = .3,size=1) +
    geom_dotplot(binaxis="y",stackdir='center', dotsize=1.3,position=position_dodge(0.8)) +
    xlab("")+
    ylab("transcript counts")+
    ggtitle(paste(i,sep=" "))
  g
  l <- length(intersect(goi_H$ENSEMBL,rownames(res)))
  assign(paste("p",n,sep=""),g)
}

ggarrange(p1+ggpubr::rotate_x_text(),
          p2+rremove("ylab")+ggpubr::rotate_x_text(),
          p3+ggpubr::rotate_x_text(),
          p4+rremove("ylab")+ggpubr::rotate_x_text(),
          common.legend = TRUE, legend = "right", nrow = 2, ncol= 2)

####################### Lina
goi_L
intersect(goi_L$ENSEMBL,rownames(res))
plots <- {}
l <- length(intersect(goi_L$ENSEMBL,rownames(res)))
n <- 1

goi_L

colData(dds)
colData(dds)$repetition <- repetition
goi_L
batch3 <- factor(rep(c("brown", "brown", "orange", "orange"),4))
colbox <- factor(levels(data$condition))

display.brewer.all()
brewer.pal(4, "Set1")

col_cond <- data$condition
levels(col_cond) <- brewer.pal(4, "Set1")
col_cond <- as.vector(col_cond)

col_cond <- (brewer.pal(4, "Set1") -> levels(data$condition)) %>%
  as.vector()

  
  data$condition
brewer.pal(4, "Set1")
col_cond <- as.vector(col_cond)

colbox <- brewer.pal(4, "Set1")
colpoint <- 

  l <- length(intersect(goi_L$ENSEMBL,rownames(res)))
n <- 1
for (n in 1:l){
  print(n)
  e <- intersect(goi_L$ENSEMBL,rownames(res))[n]
  print(e)
  data <- plotCounts(dds, gene=e, returnData = TRUE)
  data$condition <- colData(dds)$condition
  i <- goi_L[e,"Genes"]
  print(i)
  y <- sym(i)
  g <- ggplot(data = data, aes(x = condition, y =count, colour=condition, fill= condition, xlab="", ylab="", las=2)) +
    # stat_summary(fun=mean, geom="point", shape=18, size=5) +
    # geom_violin(adjust = .3,alpha = .4,position=position_dodge(0.8)) +
    geom_boxplot(alpha = .2, position=position_dodge(0.8), lwd=0.1, colour=colbox, fill=colbox) +
    # geom_line(group=c(
    #  subset(crat_data,time1 == levels(crat_data$time1)[1])[,"patient1"],
    #  subset(crat_data, time1 == levels(crat_data$time1)[2])[,"patient1"]),
    #  position=position_dodge(0.8),alpha = .3,size=1) +
    geom_dotplot(binaxis="y",stackdir='center', dotsize=1.3,position=position_dodge(0.8), fill=batch3, colour=batch3) +
    xlab("")+
    ylab("norm. txn counts")+
    ggtitle(paste(i,sep=" "))
  g
l <- length(intersect(goi_L$ENSEMBL,rownames(res)))
  assign(paste("p",n,sep=""),g)
}

ggarrange(p1+ggpubr::rotate_x_text(),
          p2+rremove("ylab")+ggpubr::rotate_x_text(),
          p3+rremove("ylab")+ggpubr::rotate_x_text(),
          p4+rremove("ylab")+ggpubr::rotate_x_text(),
          p5+ggpubr::rotate_x_text(),
          p6+rremove("ylab")+ggpubr::rotate_x_text(),
          p7+rremove("ylab")+ggpubr::rotate_x_text(),
          p8+rremove("ylab")+ggpubr::rotate_x_text(),
          p9+ggpubr::rotate_x_text(),
          p10+rremove("ylab")+ggpubr::rotate_x_text(),
          p11+rremove("ylab")+ggpubr::rotate_x_text(),
          common.legend = TRUE, legend = "right", nrow = 3, ncol= 4)

###################################################
### log2F

goi_L$Genes

pdf(paste("2021_09_17 Kelly Lina_log2FC_",".pdf",sep="")) 
ggplot(goi_L) + 
  geom_col(aes(x = Genes, y = log2FoldChange, fill=Genes)) +
  geom_label(x = goi_L$Genes, y = goi_L$log2FoldChange, label = paste(round(goi_L$log2FoldChange,digits = 1),goi_L$signif,sep=" ,"), col=factor(goi_L$signif)) +
  # coord_cartesian(ylim = c(-1, 1)) +
  labs(title=paste(sample," (",paste("Hypoxia vs. Normoxia (Ctrl)", collapse=", "),")",sep=""))
dev.off()


pdf(paste("2021_07_26 Polyamines_log2FC_",sample,"_all",".pdf",sep="")) 
ggplot(goi_L, aes(fill=results, y=log2FoldChange, x=Genes)) + 
  geom_label(x = goi_r$Genes, y = goi_r$log2FoldChange, label = goi_r$signif) +
  geom_bar(position="dodge", stat="identity")
dev.off()


###################
## data export
length(rownames(counts(dds)))
length(rownames(resLV.Nx_Hx))
summary(rownames(resLV.Nx_Hx)==resLV.Nx_Hx$gene_id[,1])

resLV.Nx_Hx$GENENAME <- mapIds(edb,
                           keys=rownames(resLV.Nx_Hx),
                           column="GENENAME",
                           keytype="GENEID",
                           multiVals="first")
resLV.Nx_Hx$SYMBOL2 <- mapIds(edb,
                         keys=rownames(resLV.Nx_Hx),
                         column="SYMBOL",
                         keytype="GENEID",
                         multiVals="first")

res_Hx_all <- cbind(data.frame(ENSEMBL=rownames(resLV.Nx_Hx),resLV.Nx_Hx[,c("GENENAME","baseMean","log2FoldChange","padj")],counts(dds)[rownames(resLV.Nx_Hx),]))

l <- length(colnames(res_Hx_all))
names <- paste(str_split(colData(dds)$run,pattern = "_", simplify=TRUE)[,3],colData(dds)$sample,sep="_")
summary(str_split(colnames(res_Hx_all)[(l-15):l],pattern = "_", simplify=TRUE)[,3]==str_split(names,pattern = "_", simplify=TRUE)[,1])
colnames(res_Hx_all)[(l-15):l] <- names

res_Hx_all
tail(res_Hx_all)

write_xlsx(res_Hx_all,"2021_09_17 Kelly Hx all.xlsx")

# Lina ##################
## 



resLV.Nx_Hx[goi_L$ENSEMBL,]
resLV.Nx_Hx[goi_L$ENSEMBL,c(1,2,5)]
counts(dds)[goi_L$ENSEMBL,]

goi_L2 <- cbind(goi_L[,1:2],resLV.Nx_Hx[goi_L$ENSEMBL,c(1,2,5)])
goi_L2$signif <- ifelse(is.na(goi_L2$padj),"NA",ifelse(goi_L2$padj < 0.001, "***", ifelse(goi_L2$padj < 0.01, "**", ifelse(goi_L2$padj < 0.05, "*", "ns"))))

goi_L2 <- cbind(goi_L2,counts(dds)[goi_L$ENSEMBL,])
goi_L2


l <- length(colnames(goi_L2))
names <- paste(str_split(colData(dds)$run,pattern = "_", simplify=TRUE)[,3],colData(dds)$sample,sep="_")
summary(str_split(colnames(goi_L2)[(l-15):l],pattern = "_", simplify=TRUE)[,3]==str_split(names,pattern = "_", simplify=TRUE)[,1])
colnames(goi_L2)[(l-15):l] <- names

goi_L2
goi_L <- goi_L2
class(goi_L$Genes)
write_xlsx(goi_L,"2021_09_17 Kelly Lina.xlsx")



library(RColorBrewer)
col.P1 <- brewer.pal(4, "Set1")
col.P2 <- brewer.pal(4, "Pastel1")
col.eigth <- c(rep(col.P1,2), rep(col.P2,2))
col.eigth
col.eigth <- as.character(col.eigth)
col <- col.eigth

par(mfrow=c(1,1))
plotCounts(dds, gene=which.min(resLV.Nx_Hx$padj), intgroup="condition", col=col)
title(sub=paste(mcols(dds)[which.min(resLV.Nx_Hx$padj),"SYMBOL"]))

mytitle = "Counts:"
mysubtitle = paste(mcols(dds)[which.min(resLV.Nx_Hx$padj),"SYMBOL"],": ",mcols(dds)[which.min(resLV.Nx_Hx$padj),"GENENAME"], sep="")
# mtext(side=3, line=2, at=0.5, adj=0, cex=1.5, mytitle)
mtext(side=3, line=0.25, at=2, adj=0.1, cex=1.25, mysubtitle)

# Plot Top 16 LV_Hx vs LV_Nx
pdf("2020_02_01 TOP16 LV_Hx_vs_LV_Nx.pdf") 
par(mfrow=c(4,4))
for (i in 1:16){
plotCounts(dds.k0, gene=rownames(res.dds.k0.LV_Nx.vs.LV_Hx.Ordered[i,]), intgroup="condition", col=col, main=paste(mcols(dds.k0)[rownames(res.dds.k0.LV_Nx.vs.LV_Hx.Ordered[i,]),"SYMBOL"]))
mysubtitle = paste(mcols(dds.k0)[rownames(res.dds.k0.LV_Nx.vs.LV_Hx.Ordered[i,]),"GENENAME"], sep="")
# mtext(side=3, line=2, at=0.5, adj=0, cex=1.5, mytitle)
mtext(side=3, line=0.25, at=0, adj=0.1, cex=0.7, mysubtitle)
}
dev.off() 

# Plot Top 16 LV_Hx vs Hif1a_KO
pdf("2020_02_01 TOP16 LV_Hx_vs_Hif1a_KO.pdf") 
par(mfrow=c(4,4))
for (i in 1:16){
plotCounts(dds.k0, gene=rownames(res.dds.k0.LV_Hx_vs_Hif1a_KO.Ordered[i,]), intgroup="condition", col=col, main=paste(mcols(dds.k0)[rownames(res.dds.k0.LV_Hx_vs_Hif1a_KO.Ordered[i,]),"SYMBOL"]))
mysubtitle = paste(mcols(dds.k0)[rownames(res.dds.k0.LV_Hx_vs_Hif1a_KO.Ordered[i,]),"GENENAME"], sep="")
# mtext(side=3, line=2, at=0.5, adj=0, cex=1.5, mytitle)
mtext(side=3, line=0.25, at=0, adj=0.1, cex=0.7, mysubtitle)
}
dev.off() 

# Plot Top 16 LV_Hx vs Hif2a_KO
pdf("2020_02_01 TOP16 LV_Hx_vs_Hif2a_KO.pdf") 
par(mfrow=c(4,4))
for (i in 1:16){
plotCounts(dds.k0, gene=rownames(res.dds.k0.LV_Hx_vs_Hif2a_KO.Ordered[i,]), intgroup="condition", col=col, main=paste(mcols(dds.k0)[rownames(res.dds.k0.LV_Hx_vs_Hif2a_KO.Ordered[i,]),"SYMBOL"]))
mysubtitle = paste(mcols(dds.k0)[rownames(res.dds.k0.LV_Hx_vs_Hif2a_KO.Ordered[i,]),"GENENAME"], sep="")
# mtext(side=3, line=2, at=0.5, adj=0, cex=1.5, mytitle)
mtext(side=3, line=0.25, at=0, adj=0.1, cex=0.7, mysubtitle)
}
dev.off()



resLV.Nx_Hx
resHx.LV_Hif1a
resHx.LV_Hif2a
resNx_Hif1a
resNx_Hif2a

subset(resLV.Nx_Hx, SYMBOL == "DPM1")

# Plot Top 16 LV_Nx vs Hif1a_KO
pdf("2020_02_01 TOP16 LV_Nx_vs_Hif1a_KO.pdf") 
par(mfrow=c(4,4))
for (i in 1:16){
  plotCounts(dds.k0, gene=rownames(resNx_Hif1a_o[i,]), intgroup="condition", col=col, main=paste(mcols(dds.k0)[rownames(resNx_Hif1a_o[i,]),"SYMBOL"]))
  mysubtitle = paste(mcols(dds.k0)[rownames(resNx_Hif1a_o[i,]),"GENENAME"], sep="")
  # mtext(side=3, line=2, at=0.5, adj=0, cex=1.5, mytitle)
  mtext(side=3, line=0.25, at=0, adj=0.1, cex=0.7, mysubtitle)
}
dev.off()

# Plot Top 16 LV_Nx vs Hif2a_KO
pdf("2020_02_01 TOP16 LV_Nx_vs_Hif2a_KO.pdf") 
par(mfrow=c(4,4))
for (i in 1:16){
  plotCounts(dds.k0, gene=rownames(resNx_Hif2a_o[i,]), intgroup="condition", col=col, main=paste(mcols(dds.k0)[rownames(resNx_Hif2a_o[i,]),"SYMBOL"]))
  mysubtitle = paste(mcols(dds.k0)[rownames(resNx_Hif2a_o[i,]),"GENENAME"], sep="")
  # mtext(side=3, line=2, at=0.5, adj=0, cex=1.5, mytitle)
  mtext(side=3, line=0.25, at=0, adj=0.1, cex=0.7, mysubtitle)
}
dev.off()


resLFC.ashr.Coll.k0.LV_Nx.vs.Hif2a_KO



# Advanced plotting with ggblot
# did not work
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


# Advances results visualizaiton in HTML
browseVignettes("ReportingTools")

desReport <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq',title = 'RNA-seq analysis of differential expression using DESeq',reportDirectory = "./reports")

des2Report <- HTMLReport(shortName = 'Kelly_LV_Hx_vs._Hif1a', title = 'RNA-seq analysis of differential expression using DESeq2', reportDirectory = "./reports")
publish(dds.k0,des2Report, pvalueCutoff=0.2, factor = colData(dds.k0)$condition, reportDir="./reports")
finish(des2Report)
# no Gene symbol...

df <- rownames(dds.k0)
add.anns <- function(df, mart, ...)
{
    nm <- rownames(df)
    anns <- getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol",
"description"),
        filters = "ensembl_gene_id", values = nm, mart = mart)
    anns <- anns[match(nm, anns[, 1]), ]
    colnames(anns) <- c("ID", "Gene Symbol", "Gene Description")
    df <- cbind(anns, df[, 2:ncol(df)])
    rownames(df) <- nm
    df
}

desReport <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq',title = 'RNA-seq analysis of differential expression using DESeq',reportDirectory = "./reports")

des2Report <- HTMLReport(shortName = 'Kelly_LV_Hx_vs._Hif1a', title = 'RNA-seq analysis of differential expression using DESeq2', reportDirectory = "./reports")
publish(dds.k0,des2Report, pvalueCutoff=0.1, factor = colData(dds.k0)$condition, .modifyDF = list(add.anns,modifyReportDF), mart = mart, reportDir="./reports")
finish(des2Report)

## heat-map
vsd <- vst(dds)
ntd <- normTransform(dds)

mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, vsd$batch)
assay(vsd) <- mat
plotPCA(vsd)


library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("treatment","genotype")])

par(mfrow=c(2,1))

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE,show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)



sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

### PCA
pdf("2021_07_20 Kelly PCA.pdf")
plotPCA(vsd.k0, intgroup=c("condition", "batch"))
dev.off() 
# removed batch effect

vsd.nobatch <- vsd
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, vsd$batch)
assay(vsd.nobatch) <- mat
plotPCA(vsd.nobatch)


##### Check for overlappin genes:

resLV.Nx_Hx_s
resHx.LV_Hif1a
resHx.LV_Hif2a
resNx_Hif1a
resNx_Hif2a

df <- subset(resLV.Nx_Hx, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
df <- rbind(df,subset(resHx.LV_Hif1a, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)))
df <- rbind(df,subset(resHx.LV_Hif2a, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)))



resLV.Nx_Hx_s <- subset(resLV.Nx_Hx, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
resHx.LV_Hif1a_s <- subset(resHx.LV_Hif1a, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
resHx.LV_Hif2a_s <- subset(resHx.LV_Hif2a, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
# resNx_Hif1a_s <- subset(resNx_Hif1a, padj<0.05 & baseMean > 10 & (log2FoldChange > 0.5 | log2FoldChange < -0.5))
# resNx_Hif2a_s <- subset(resNx_Hif2a, padj<0.05 & baseMean > 10 & (log2FoldChange > 0.5 | log2FoldChange < -0.5))

length(rownames(resLV.Nx_Hx_s))
length(rownames(resHx.LV_Hif1a_s))
length(rownames(resHx.LV_Hif2a_s))
length(rownames(resNx_Hif1a_s))
length(rownames(resNx_Hif2a_s))

mcols(dds.k0)["ENSG00000161958.11","SYMBOL"]


ol_LV_Hif1a_Hx <- merge(data.frame(resHx.LV_Hif1a_s),data.frame(resLV.Nx_Hx_s)[2:5], by = 'row.names')
length(rownames(ol_LV_Hif1a_Hx))

ol_LV_Hif2a_Hx <- merge(data.frame(resHx.LV_Hif2a_s),data.frame(resLV.Nx_Hx_s)[2:5], by = 'row.names')
length(rownames(ol_LV_Hif2a_Hx))


ol_LV_Hif1a_Nx <- merge(data.frame(resNx_Hif1a_s),data.frame(resLV.Nx_Hx_s)[2:5], by = 'row.names')
length(rownames(ol_LV_Hif1a_Nx))

ol_LV_Hif2a_Nx <- merge(data.frame(resNx_Hif2a_s),data.frame(resLV.Nx_Hx_s)[2:5], by = 'row.names')
length(rownames(ol_LV_Hif2a_Nx))



ol_Hif1a_Hif2a_Nx <- merge(data.frame(resNx_Hif1a_s),data.frame(resNx_Hif2a_s)[2:5], by = 'row.names')
length(rownames(ol_Hif1a_Hif2a_Nx))

ol_Hif1a_Hif2a_Hx <- merge(data.frame(resHx.LV_Hif1a_s),data.frame(resHx.LV_Hif2a_s)[2:5], by = 'row.names')
length(rownames(ol_Hif1a_Hif2a_Hx))

# examples
# lostgenes <- subset(merged_base10,padj.y<0.05 & baseMean.y > 10 & log2FoldChange.y > 0.5)
# merged.res <- merge(data.frame(res), omiqa[,c(1,5:7,16,32)], by = 'GENEID', all.x = TRUE)
# merged.omiqa <- merge(data.frame(res), omiqa[,c(1,5:7,16,32)], by = 'GENEID', all.y = TRUE)
# merged.ol <- merged.omiqa <- merge(data.frame(res), omiqa[,c(1,5:7,16,32)], by = 'GENEID')



# Venn Diagramm ############################################
# Chart Venn Diagramm
library(VennDiagram)
library(RColorBrewer)
install.packages("svglite")
library(tidyverse)
library(cowplot)

length(rownames(resLV.Nx_Hx_s))
length(rownames(resHx.LV_Hif1a_s))
length(rownames(resHx.LV_Hif2a_s))

set1 <- rownames(resLV.Nx_Hx_s)
set2 <- rownames(resHx.LV_Hif1a_s)
set3 <- rownames(resHx.LV_Hif2a_s)



myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("LV" , "Hif1a" , "Hif2a"),
  filename = 'venn_diagramm_HIFs_Nx.svg',
  output=TRUE,
  
  # Output features
  imagetype="svg" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .3,
  # fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  # cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


## venn.diagram ##########################



set1 <- rownames(resLV.Nx_Hx_s)
set2 <- rownames(resHx.LV_Hif1a_s)
set3 <- rownames(resHx.LV_Hif2a_s)

myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("LV" , "Hif1a" , "Hif2a"),
  filename = 'venn_diagramm_HIFs_Hx.png',
  output=TRUE,

  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .3,
  # fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  # cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

# Calcualte numbers #########
# 
vn <- calculate.overlap(x = list(set1, set2, set3))
a1 <- data.frame(SYMBOL=mcols(dds.k0)[vn$a1,"SYMBOL"],
                 baseMean=mcols(dds.k0)[vn$a1,"baseMean"],
                 counts(dds.k0)[vn$a1,]
                 )
write_xlsx(as.data.frame(a1),"a1.xlsx")

a2 <- data.frame(SYMBOL=mcols(dds.k0)[vn$a2,"SYMBOL"],
                 baseMean=mcols(dds.k0)[vn$a2,"baseMean"],
                 counts(dds.k0)[vn$a2,]
)
write_xlsx(as.data.frame(a2),"a2.xlsx")

a3 <- data.frame(SYMBOL=mcols(dds.k0)[vn$a3,"SYMBOL"],
                 baseMean=mcols(dds.k0)[vn$a3,"baseMean"],
                 counts(dds.k0)[vn$a3,]
)
write_xlsx(as.data.frame(a3),"a3.xlsx")

a4 <- data.frame(SYMBOL=mcols(dds.k0)[vn$a4,"SYMBOL"],
                 baseMean=mcols(dds.k0)[vn$a4,"baseMean"],
                 counts(dds.k0)[vn$a4,]
)
write_xlsx(as.data.frame(a4),"a4.xlsx")

a5 <- data.frame(SYMBOL=mcols(dds.k0)[vn$a5,"SYMBOL"],
                 baseMean=mcols(dds.k0)[vn$a5,"baseMean"],
                 counts(dds.k0)[vn$a5,]
)
write_xlsx(as.data.frame(a5),"a5.xlsx")

a6 <- data.frame(SYMBOL=mcols(dds.k0)[vn$a6,"SYMBOL"],
                 baseMean=mcols(dds.k0)[vn$a6,"baseMean"],
                 
                 counts(dds.k0)[vn$a6,]
)
write_xlsx(as.data.frame(a6),"a6.xlsx")

a7 <- data.frame(SYMBOL=mcols(dds.k0)[vn$a7,"SYMBOL"],
                 baseMean=mcols(dds.k0)[vn$a7,"baseMean"],
                 counts(dds.k0)[vn$a7,]
)
write_xlsx(as.data.frame(a7),"a7.xlsx")

# Plot genes ######################################################################
# 

pdf("2021_07_21 TOP16 a1.pdf") 
par(mfrow=c(4,4))
for (i in 1:16){
  plotCounts(dds.k0, gene=rownames(resNx_Hif2a_o[i,]), intgroup="condition", col=col, main=paste(mcols(dds.k0)[rownames(resNx_Hif2a_o[i,]),"SYMBOL"]))
  mysubtitle = paste(mcols(dds.k0)[rownames(resNx_Hif2a_o[i,]),"GENENAME"], sep="")
  # mtext(side=3, line=2, at=0.5, adj=0, cex=1.5, mytitle)
  mtext(side=3, line=0.25, at=0, adj=0.1, cex=0.7, mysubtitle)
}
dev.off()




# Export R Data ########################################################################
# 
saveRDS(dds.k0, file = "C:/Users/kelterbs/OneDrive - Charit? - Universit?tsmedizin Berlin/RNA-seq/Kelly CRISPR HIF KO/2021_07_21_dds.k0.RDS") 
saveRDS(dds, file = "C:/Users/kelterbs/OneDrive - Charit? - Universit?tsmedizin Berlin/RNA-seq/Kelly CRISPR HIF KO/2021_09_14_dds_allcounts.RDS") 

dds.k0_load <- readRDS("C:/Users/kelterbs/OneDrive - Charit? - Universit?tsmedizin Berlin/RNA-seq/Kelly CRISPR HIF KO/2021_07_21_dds.k0.RDS")


# Export to CSV ########################################################################
# 
write.csv2(as.data.frame(res.dds.k0.LV_Hx_vs_Hif1a_KO),file="LV_Hx_vs_Hif1a_KO.csv")

# Export to XLS
write_xlsx(as.data.frame(res.dds.k0.LV_Hx_vs_Hif1a_KO),"Kelly_results_combined.xlsx")

# Combined output
# LV_Hx vs. Hif1a_KO
res.dds.k0.LV_Nx.vs.LV_Hx
res.ddsColl.k0.LV_Nx.vs.LV_Hx

# LV_Hx vs. Hif1a_KO
res.dds.k0.LV_Hx_vs_Hif1a_KO
res.ddsColl.k0.LV_Hx_vs_Hif1a_KO

# LV_Hx vs. Hif2a_KO
res.dds.k0.LV_Hx_vs_Hif2a_KO
res.ddsColl.k0.LV_Hx_vs_Hif2a_KO





mcols(dds.k0)

results_combined <- data.frame(
  mcols(dds.k0)["SYMBOL"],
  mcols(dds.k0)["gene_id"],
  mcols(dds.k0)["GENENAME"],
  mcols(dds.k0)["baseMean"],
  res.ddsColl.k0.LV_Nx.vs.LV_Hx[,c(2,5:6)],
  res.ddsColl.k0.LV_Hx_vs_Hif1a_KO[,c(2,5:6)],
  res.ddsColl.k0.LV_Hx_vs_Hif2a_KO[,c(2,5:6)],
  unlist(mcols(dds.k0)["tx_ids"]),
  mcols(dds.k0)["REFSEQ"],
  counts(dds.k0)
)

write_xlsx(results_combined,"Kelly_results_combined20210720.xlsx")






head(results_combined)
head(res.ddsColl.k0.LV_Nx.vs.LV_Hx[,c(2,5:6)])
head(counts(dds.k0))
length(rownames(counts(dds.k0)))
head(mcols(dds.k0)["gene_id"])
head(mcols(dds.k0)["GENENAME"])
head(unlist(mcols(dds.k0)["GENENAME"]))
head(as.character(mcols(dds.k0)["GENENAME"]))
head(unlist(mcols(dds.k0)["GENENAME"], use.names = FALSE)[])

as.integer(mcols(dds.k0)["GENENAME"])
purrr::flatten_dbl(mcols(dds.k0)["GENENAME"])
rlang::last_error()

head(paste(mcols(dds.k0)["gene_id"]))
head(as.data.frame(mcols(dds.k0)["gene_id"])[,1])
mcols(res.dds.k0.LV_Hx_vs_Hif1a_KO)

write_xlsx(as.data.frame(results_combined),"Kelly_results_combined.xlsx")


# poster #####################################################################################################################################################################
### 

### MA plots

getwd()
setwd("/Users/foxios/OneDrive - Charité - Universitätsmedizin Berlin/RNA-seq/Kelly CRISPR HIF KO/poster")
getwd()

pdf("2021_09_23 MA plots_all.pdf", width = 9, height = 4) 
par(mfrow=c(1,3))
plotMA(resHx.LV_Hif1a, alpha = 0.05, colSig = "#BF6D39", ylim = c(-6,6))
title(main="HIF-1a vs. Ctrl (Hx)")
plotMA(resLV.Nx_Hx, alpha = 0.05, colSig = "#0066CC", ylim = c(-6,6))
title(main="Hypoxia vs. Normoxia")
plotMA(resHx.LV_Hif2a, alpha = 0.05, colSig = "#00CC99", ylim = c(-6,6))
title(main="HIF-2a vs. Ctrl (Hx)")
dev.off()
resultsNames(dds)


summary(resLV.Nx_Hx)
summary(resHx.LV_Hif1a)
summary(resHx.LV_Hif2a)
summary(resNx_Hif1a)
summary(resNx_Hif2a)
## uses 0.1

length(rownames(subset(resLV.Nx_Hx,padj < 0.05 & log2FoldChange > 1)))
length(rownames(subset(resLV.Nx_Hx,padj < 0.05 & log2FoldChange < -1)))
length(rownames(subset(resHx.LV_Hif1a,padj < 0.05 & log2FoldChange > 1)))
length(rownames(subset(resHx.LV_Hif1a,padj < 0.05 & log2FoldChange < -1)))
length(rownames(subset(resHx.LV_Hif2a,padj < 0.05 & log2FoldChange > 1)))
length(rownames(subset(resHx.LV_Hif2a,padj < 0.05 & log2FoldChange < -1)))

length(rownames(subset(resHx.LV_Hif1a,padj < 0.05 & log2FoldChange < 0)))
length(rownames(subset(resHx.LV_Hif1a,padj < 0.05 & log2FoldChange < 0)))

## GO terms ###############################################
## 
data(gse16873.d)
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",species = "hsa", out.suffix = "gse16873")
i <- 1
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],species = "hsa", out.suffix = "gse16873", kegg.native = T)
str(pv.out)
resLV.Nx_Hx$log2FoldChange
resLV.Nx_Hx[,2]
pvdata <- data.frame(log2F = resLV.Nx_Hx[,2])
rownames(pvdata) <- resLV.Nx_Hx[,8]
pv.out <- pathview(gene.data = pvdata, pathway.id = "04110",species = "hsa", out.suffix = "resLV.Nx_Hx")

resHx <- subset(resLV.Nx_Hx, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
genelistUp <- factor(as.integer(resHx$padj < .05 & resHx$log2FoldChange > 0 ) )
names(genelistUp) <- rownames(resHx)

library("topGO")
myGOdata <- new( "topGOdata",
                 ontology = "BP",
                 allGenes = genelistUp,
                 nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="ensembl" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
GenTable( myGOdata, goTestResults )

data(genes)
pwf <- nullp(genes,'hg19','ensGene')
pvals <- goseq(pwf,'hg19','ensGene')
head(pvals)

library(goseq)
supportedGenomes()
vignette(goseq)
supportedGeneIDs()

class(resLV.Nx_Hx)


#####################################
# https://bioconductor.github.io/BiocWorkshops/functional-enrichment-analysis-of-high-throughput-omics-data.html#workshop-description-4

library("goseq")
library("geneLenDataBase")
library("org.Hs.eg.db")

library(EnrichmentBrowser)
library("GO.db")
db = c("go", "kegg", "msigdb", "enrichr")
showAvailableSpecies(db = c("go", "kegg", "msigdb", "enrichr"), cache = TRUE)
showAvailableCollections(
  org,
  db = c("go", "kegg", "msigdb", "enrichr"),
  cache = TRUE
)
kegg.gs <- getGenesets(org="hsa", db="kegg")
go.gs <- getGenesets(org="hsa", db="go")
# geht nicht: go.gs <- getGenesets(org="hsa", db="go", go.onto="BP", go.mode="GO.db")

ora.all <- sbea(method="ora", se=resLV.Nx_Hx, gs=hsa.gs, perm=0, alpha=0.2)

gsRanking(ora.all)
# example
library(ALL)
data(ALL)
####################################
# https://wikis.utexas.edu/display/bioiteam/GO+Enrichment+using+goseq
library("goseq")
library("geneLenDataBase")
library("org.Hs.eg.db")

ALL <- resLV.Nx_Hx$gene_id
DEG <- subset(resLV.Nx_Hx, padj < 0.05)$gene_id
DEG <- subset(resHx.LV_Hif1a, padj < 0.05)$gene_id
DEG <- subset(resHx.LV_Hif2a, padj < 0.05)$gene_id


ALL.vector <- c(t(ALL))
DEG.vector <- c(t(DEG))

gene.vector=as.integer(ALL.vector%in%DEG.vector)
names(gene.vector)=ALL.vector 
#lets explore this new vector a bit
head(gene.vector)
tail(gene.vector)

supportedGenomes()
pwf <- nullp(gene.vector,"hg19","ensGene")
GO.wall <- goseq(pwf,"hg19","ensGene")

#How many enriched GO terms do we have
class(GO.wall)
head(GO.wall)
nrow(GO.wall)

go <- GO.wall[1:10,]

goseq_Hx <- GO.wall
goseq_HIF1a <- GO.wall
goseq_HIF2a <- GO.wall

goseq_Hx$condition <- rep("Hx",length(rownames(goseq_Hx)))
goseq_Hx_bp <- subset(goseq_Hx,ontology=="BP")
goseq_Hx_mf <- subset(goseq_Hx,ontology=="MF")
goseq_Hx_cc <- subset(goseq_Hx,ontology=="CC")
goseq_Hx_mix <- rbind(goseq_Hx_bp[1:5,],goseq_Hx_mf[1:5,],goseq_Hx_cc[1:5,])
goseq_Hx_mix$order <- c(15:1)

goseq_HIF1a$condition <- rep("HIF-1a",length(rownames(goseq_HIF1a)))
goseq_HIF1a_bp <- subset(goseq_HIF1a,ontology=="BP")
goseq_HIF1a_mf <- subset(goseq_HIF1a,ontology=="MF")
goseq_HIF1a_cc <- subset(goseq_HIF1a,ontology=="CC")
goseq_HIF1a_mix <- rbind(goseq_HIF1a_bp[1:5,],goseq_HIF1a_mf[1:5,],goseq_HIF1a_cc[1:5,])
goseq_HIF1a_mix$order <- c(15:1)

goseq_HIF2a$condition <- rep("HIF-2a",length(rownames(goseq_HIF2a)))
goseq_HIF2a_bp <- subset(goseq_HIF2a,ontology=="BP")
goseq_HIF2a_mf <- subset(goseq_HIF2a,ontology=="MF")
goseq_HIF2a_cc <- subset(goseq_HIF2a,ontology=="CC")
goseq_HIF2a_mix <- rbind(goseq_HIF2a_bp[1:5,],goseq_HIF2a_mf[1:5,],goseq_HIF2a_cc[1:5,])
goseq_HIF2a_mix$order <- c(15:1)

gos_all <- rbind(goseq_Hx,goseq_HIF1a,goseq_HIF2a)
gos_all <- gos_all[order(gos_all$over_represented_pvalue),]

ggplot(data=goseq_HIF1a_mix, aes(x=reorder(term,-over_represented_pvalue), y=-log10(over_represented_pvalue), fill=ontology)) +
         geom_bar(stat="identity") +
          coord_flip()+
  facet_wrap(~.+ontology)

ggplot(data=goseq_Hx_mix, aes(x=term, y=-log10(over_represented_pvalue), fill=ontology)) +
  geom_bar(stat="identity") +
  coord_flip()
ggplot(data=goseq_Hx_mix, aes(x=reorder(term,order), y=-log10(over_represented_pvalue), fill=ontology)) +
  geom_bar(stat="identity") +
  coord_flip()+
  geom_hline(yintercept = 2)
ggplot(data=goseq_HIF1a_mix, aes(x=reorder(term,order), y=-log10(over_represented_pvalue), fill=ontology)) +
  geom_bar(stat="identity") +
  coord_flip()+
  geom_hline(yintercept = 2)
ggplot(data=goseq_HIF2a_mix, aes(x=reorder(term,order), y=-log10(over_represented_pvalue), fill=ontology)) +
  geom_bar(stat="identity") +
  coord_flip() +
geom_hline(yintercept = 2)

-log10(0.01)

-log10(7.538164e-28)

enriched.GO=GO.wall$category[GO.wall$over_represented_pvalue<.05]
#NOTE: They recommend using a more stringent multiple testing corrected p value here

#How many GO terms do we have now?
class(enriched.GO)
head(enriched.GO)
length(enriched.GO)

library(GO.db)
capture.output(for(go in enriched.GO[1:258]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file="SigGo.txt")

less "SigGo.txt"

# KEGG ##########################################
#
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)


BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(DOSE)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
require(clusterProfiler)
data(geneList)


# reading in data from deseq2
df <- subset(resLV.Nx_Hx, padj < 0.05)
df <- subset(resHx.LV_Hif1a, padj < 0.05)
df <- subset(resHx.LV_Hif2a, padj < 0.05)

# df <- subset(resLV.Nx_Hx, padj < 0.05 & (log2FoldChange > 0))

# list of genes with 6524
# df <- read.csv("drosphila_example_de.csv", header=TRUE)
length(rownames(df))

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$gene_id

# omit any NA values 
# gene_list<-na.omit(original_gene_list)
gene_list<-(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
length(gene_list)

ego <- enrichGO(gene          = names(gene_list),
                keyType = "ENSEMBL",
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

head(as.data.frame(ego))
dotplot(ego)

df["ENSG00000168477",]

gse <- gseGO(geneList=gene_list, 
             nPermSimple = 100000,
             ont ="ALL", 
             keyType = "ENSEMBL",
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             minGSSize    = 10,
             maxGSSize    = 500,
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")
sgse <- as.data.frame(gse)
# sgse <- sgse[order(sgse$rank),]
head(sgse)
head(sgse)[,1:10]
subset(sgse, Description == "cell cycle checkpoint signaling")[,1:10]

# gse_Hx <- gse
sgse_Hx <- as.data.frame(gse_Hx)
write_xlsx(as.data.frame(sgse_Hx),"2021_09_27 GO_Hx.xlsx")

# gse_HIF1a <- gse
sgse_HIF1a <- as.data.frame(gse_Hx)
write_xlsx(as.data.frame(sgse_HIF1a),"2021_09_27 GO_HIF1a.xlsx")

gse_HIF2a <- gse
sgse_HIF2a <- as.data.frame(gse_Hx)
write_xlsx(as.data.frame(sgse_HIF2a),"2021_09_27 GO_HIF2a.xlsx")

dotplot(gse, showCategory=30, split=".sign") + facet_grid(.~.sign)


plotGOgraph(gse)

write_xlsx(as.data.frame(df),"2021_09_27 gene list.xlsx")
gseaplot(gse, gse$ID[50])


# Get Go term children ##############################

as.list(GOMFCHILDREN['GO:0019900'])
as.list(GOMFOFFSPRING['GO:0019900']) 
length(as.list(GOMFOFFSPRING['GO:0019900'])[[1]])


# compare different clusters... #####################

dfhx <- subset(resLV.Nx_Hx,padj < 0.05)
dfhx <- dfhx[order(dfhx$log2FoldChange, decreasing = TRUE),]
dfhx <- dfhx$ENTREZ

dfhif1 <- subset(resHx.LV_Hif1a,padj < 0.05)
dfhif1 <- dfhif1[order(dfhif1$log2FoldChange, decreasing = TRUE),]
dfhif1 <- dfhif1$ENTREZ

dfhif2 <- subset(resHx.LV_Hif2a,padj < 0.05)
dfhif2 <- dfhif2[order(dfhif2$log2FoldChange, decreasing = TRUE),]
dfhif2 <- dfhif2$ENTREZ

dfhif1
dfhif2
df <- list(HIF1a = dfhif1,
           Hx = dfhx,
           HIF2a = dfhif2)

ck <- compareCluster(geneCluster = df, fun = "enrichGO", OrgDb = org.Hs.eg.db)
ck <- compareCluster(geneCluster = df, fun = "enrichKEGG", organism="hsa", pvalueCutoff=0.05)

ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck)[1:7]
dotplot(ck)
dotplot(ck, showCategory=20)

# ck_kegg <- ck




# heat map of top genes ###################################

require(RColorBrewer)
require(ComplexHeatmap)
require(circlize)
require(digest)
require(cluster)



# collapse replicates
dds$run <- samplenames
dds_heat <- collapseReplicates(dds, dds$condition,dds$run)
colData(dds_heat)
dds_heat$runsCollapsed
colnames(dds_heat)
matchFirstLevel <- dds_heat$sample == levels(dds_heat$sample)[1]
stopifnot(all(rowSums(counts(dds_heat[,matchFirstLevel])) == counts(dds_heat[,1])))

df <- subset(resLV.Nx_Hx, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
df <- rbind(df,subset(resHx.LV_Hif1a, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)))
df <- rbind(df,subset(resHx.LV_Hif2a, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)))

colnames(resLV.Nx_Hx)
colnames(resHx.LV_Hif1a)
colnames(resHx.LV_Hif2a)

summary(duplicated(rownames(df)))
summary(unique(rownames(df)))

df <- df[!duplicated(rownames(df)),]
df <- subset(df,!SYMBOL2 =="")
summary(df$SYMBOL2)
summary(is.na(df$SYMBOL2))
summary(df$SYMBOL2=="")

mat <- counts(dds_heat)
length(rownames(counts(dds_heat)))
length(rownames(df))
rownames(mat)<- res$SYMBOL2
length(rownames(mat))
summary(subset(df, padj < 0.05)$SYMBOL2 %in% rownames(mat))

# mat <- mat[subset(df, padj < 0.05)$SYMBOL2,]
# mat <- mat[rownames(mat) %in% subset(df, padj < 0.05)$SYMBOL2,]
mat <- mat[df$SYMBOL2,]

"WT1" %in% rownames(mat)
goi_L2
summary(goi_L[[1]] %in% rownames(mat))
summary(goi_N[[1]] %in% rownames(mat))


#### TOP genes

hx <- resLV.Nx_Hx[order(resLV.Nx_Hx$log2FoldChange),] %>% subset(padj < 0.05 & baseMean > 200 & !SYMBOL2 == "")
# hxup <- resLV.Nx_Hx[order(resLV.Nx_Hx$log2FoldChange, decreasing = TRUE),10]
h1 <- resLV.Nx_Hx[order(resHx.LV_Hif1a$log2FoldChange),] %>% subset(padj < 0.05 & baseMean > 200 & !SYMBOL2 == "")
# h1up <- resLV.Nx_Hx[order(resHx.LV_Hif1a$log2FoldChange, decreasing = TRUE),10]
h2 <- resLV.Nx_Hx[order(resHx.LV_Hif2a$log2FoldChange),] %>% subset(padj < 0.05 & baseMean > 200 & !SYMBOL2 == "")
#h2up <- resLV.Nx_Hx[order(resHx.LV_Hif2a$log2FoldChange, decreasing = TRUE),10]

colnames(resLV.Nx_Hx)
colnames(resHx.LV_Hif1a)
colnames(resHx.LV_Hif2a)

goi_N[[1]]

# genes <- subset(goi_L[[1]],goi_L[[1]] %in% rownames(mat))
# genes <- c(genes,subset(goi_N[[1]],goi_N[[1]] %in% rownames(mat)))
genes <- c(goi_L[[1]],
           goi_N[[1]],
           hx[1:10,10],tail(hx[,10],n=10),
           h1[1:10,10],tail(h1[,10],n=10),
           h2[1:10,10],tail(h2[,10],n=10))
length(genes)
genes <- subset(rownames(mat),rownames(mat) %in% genes)
length(genes)
genes <- unique(genes)
length(genes)

position <- mat
l <- length(rownames(mat))
pos <- data.frame(genes = rownames(mat),pos = 1:l)
rownames(pos) <- rownames(mat)
pos <- pos[genes,2]

# mat <- cbind(res,counts(dds_heat))
# mat <- subset(mat,!SYMBOL2=="" & SYMBOL2 %in% df$SYMBOL2)
# rownames(mat) <- mat$SYMBOL2
# mat <- mat[,c(12:15)]
heat <- t(scale(t(mat)))
summary(heat)
myCol <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
myBreaks <- seq(-1.5, 1.5, length.out = 100)

ann <- data.frame(condition = colData(dds_heat)[,c(3)],
  stringsAsFactors = FALSE)

colAnn <- HeatmapAnnotation(
  df = ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  # col = colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'))

boxplotCol <- HeatmapAnnotation(
  boxplot = anno_boxplot(
    heat,
    border = FALSE,
    gp = gpar(fill = '#CCCCCC'),
    pch = '.',
    size = unit(2, 'mm'),
    axis = TRUE,
    axis_param = list(
      gp = gpar(fontsize = 12),
      side = 'left')),
  annotation_width = unit(c(2.0), 'cm'),
  which = 'col')

boxplotRow <- HeatmapAnnotation(
  boxplot = row_anno_boxplot(
    heat,
    border = FALSE,
    gp = gpar(fill = '#CCCCCC'),
    pch = '.',
    size = unit(2, 'mm'),
    axis = TRUE,
    axis_param = list(
      gp = gpar(fontsize = 12),
      side = 'top')),
  annotation_width = unit(c(2.0), 'cm'),
  which = 'row')

genelabels <- rowAnnotation(
  Genes = anno_mark(
    at = pos,
    labels = genes,
    labels_gp = gpar(fontsize = 15, fontface = 'bold'),
    padding = 0.75),
  width = unit(2.0, 'cm') +
    
    max_text_width(
      rownames(heat)[seq(1, nrow(heat), 40)],
      gp = gpar(fontsize = 10,  fontface = 'bold')))

# Clusters
# pamClusters <- cluster::pam(heat, k = 4) # pre-select k = 4 centers
pamClusters <- cluster::pam(heat, k= 4) # pre-select k = 4 centers

pamClusters$clustering <- paste0('Cluster ', pamClusters$clustering)
pamClusters$clustering <- factor(pamClusters$clustering,
                                 levels = c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4','Cluster 5', 'Cluster 6'))
# Heatmap(heat)
hmap <- Heatmap(heat,
                
                # split the genes / rows according to the PAM clusters
                split = pamClusters$clustering,
                cluster_row_slices = FALSE,
                
                name = 'Gene\nZ-\nscore',
                
                col = colorRamp2(myBreaks, myCol),
                
                # parameters for the colour-bar that represents gradient of expression
                heatmap_legend_param = list(
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(8, 'cm'),
                  legend_height = unit(5.0, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 12, fontface = 'bold'),
                  labels_gp=gpar(fontsize = 12, fontface = 'bold')),
                
                # row (gene) parameters
                cluster_rows = TRUE,
                show_row_dend = TRUE,
                #row_title = 'Statistically significant genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
                row_title_rot = 90,
                show_row_names = FALSE,
                row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                row_names_side = 'left',
                row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = TRUE,
                show_column_dend = TRUE,
                column_title = '',
                column_title_side = 'bottom',
                column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                column_title_rot = 0,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                column_names_max_height = unit(10, 'cm'),
                column_dend_height = unit(25,'mm'),
                
                # cluster methods for rows and columns
                clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                clustering_method_columns = 'ward.D2',
                clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                clustering_method_rows = 'ward.D2',
                
                # specify top and bottom annotations
                top_annotation = colAnn,
                bottom_annotation = boxplotCol)

######################
#####

hmap <- Heatmap(heat,
                
                # split the genes / rows according to the PAM clusters
                split = pamClusters$clustering,
                cluster_row_slices = TRUE,
                
                name = 'Gene\nZ-\nscore',

                col = colorRamp2(myBreaks, myCol),
                
                # parameters for the colour-bar that represents gradient of expression
                heatmap_legend_param = list(
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(8, 'cm'),
                  legend_height = unit(5.0, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 12, fontface = 'bold'),
                  labels_gp=gpar(fontsize = 12, fontface = 'bold')),
                
                # row (gene) parameters
                cluster_rows = TRUE,
                show_row_dend = TRUE,
                row_title = 'Statistically significant genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
                row_title_rot = 90,
                show_row_names = FALSE,
                row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                row_names_side = 'left',
                row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = TRUE,
                show_column_dend = TRUE,
                column_title = '',
                column_title_side = 'bottom',
                column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                column_title_rot = 0,
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                column_names_max_height = unit(10, 'cm'),
                column_dend_height = unit(25,'mm'),
                
                # cluster methods for rows and columns
                clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                clustering_method_columns = 'ward.D2',
                clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                clustering_method_rows = 'ward.D2',
                
                # specify top and bottom annotations
                top_annotation = colAnn)
                # bottom_annotation = boxplotCol)


pdf("Heatmap2.pdf", height = 18, width = 6)
draw(hmap + genelabels,
     heatmap_legend_side = 'left',
     annotation_legend_side = 'right',
     row_sub_title_side = 'left')
dev.off()



# Venn Diagramm Hypoxia, Hif1a, Hif2a
####################################################################

set1 <- rownames(resLV.Nx_Hx_s)
set2 <- rownames(resHx.LV_Hif1a_s)
set3 <- rownames(resHx.LV_Hif2a_s)

venn <- venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("LV" , "Hif1a" , "Hif2a"),
  filename = NULL,
  # output=TRUE,
  
  # Output features
  # imagetype="svg" ,
  # height = 480 , 
  # width = 480 , 
  # resolution = 300,
  # compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .3,
  # fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  # cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

ggsave(venn, file="Venn.svg", device = "svg")


####################################################################
## Neuroblastoma genes

nb <- read.table("./neuroblastoma dataset/Childhood_Neuroblastoma_from_DisGeNET.gmt", header=TRUE)
head(nb)
class(nb)

files <- list.files("./neuroblastoma dataset/", patter=".gmt")
nb_all <- {}
for (i in files){
  print(i)
  nb <- read.table(paste("./neuroblastoma dataset/",i,sep=""),header=TRUE)
  nb <- colnames(nb)
  nb_all <- c(nb_all,nb)
  }
nb_all

nb_all <- unique(nb_all)
nb_all

set1 <- resLV.Nx_Hx_s$SYMBOL2
set2 <- resHx.LV_Hif1a_s$SYMBOL2
set3 <- resHx.LV_Hif2a_s$SYMBOL2
set23 <- unique(c(set2,set3))

nb_all <- nb_all[nb_all %in% res$SYMBOL2]
set4 <- nb_all

venn <- venn.diagram(
  x = list(set1, set23, set4),
  category.names = c("LV" , "Hif1a/Hif2a" ,"Neuroblastoma"),
  filename = NULL,
  output=TRUE,
  
  # Output features
  # imagetype="svg" ,
  # height = 480 , 
  # width = 480 , 
  # resolution = 300,
  # compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .3,
  # fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  # cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

ggsave(venn, file="Venn_Neuro.svg", device = "svg")

## Plot counts
#####################################
####################### Lina

nb_i<-intersect(nb_all,set1)
nb_i<-intersect(nb_i,set23)
length(nb_i)

nb_i <- subset(res,SYMBOL2 %in% nb_i)
nb_i <- nb_i[order(nb_i$log2FoldChange),]


goi_L <- data.frame(Genes = c("WT1","CA9","SULF2"))

goi_L$ENSEMBL <- mapIds(edb,
                        keys=goi_L$Genes,
                        column="GENEID",
                        keytype="SYMBOL",
                        multiVals="first")

goi_L$GENENAME <- mapIds(edb,
                         keys=goi_L$Genes,
                         column="GENENAME",
                         keytype="SYMBOL",
                         multiVals="first")

rownames(goi_L) <- goi_L$ENSEMBL

goi_L
intersect(goi_L$ENSEMBL,rownames(res))
plots <- {}
l <- length(intersect(goi_L$ENSEMBL,rownames(res)))
n <- 1

goi_L

colData(dds)
colData(dds)$repetition <- repetition
goi_L
batch3 <- factor(rep(c("black", "black", "darkgrey", "darkgrey"),4))
colbox <- factor(levels(data$condition))

display.brewer.all()
brewer.pal(4, "Set1")

col_cond <- data$condition
levels(col_cond) <- brewer.pal(4, "Set1")
col_cond <- as.vector(col_cond)

col_cond <- (brewer.pal(4, "Set1") -> levels(data$condition)) %>%
  as.vector()


data$condition
brewer.pal(4, "Set1")
col_cond <- as.vector(col_cond)

# colbox <- brewer.pal(4, "Set1")
colbox <- c("#ED1C24","#0066CC", "#BF6D39", "#00CC99")

  l <- length(intersect(goi_L$ENSEMBL,rownames(res)))
n <- 1
for (n in 1:l){
  print(n)
  e <- intersect(goi_L$ENSEMBL,rownames(res))[n]
  print(e)
  data <- plotCounts(dds, gene=e, returnData = TRUE)
  data$condition <- colData(dds)$condition
  i <- goi_L[e,"Genes"]
  print(i)
  y <- sym(i)
  g <- ggplot(data = data, aes(x = condition, y =count, colour=condition, fill= condition, xlab="", ylab="", las=2)) +
    # stat_summary(fun=mean, geom="point", shape=18, size=5) +
    # geom_violin(adjust = .3,alpha = .4,position=position_dodge(0.8)) +
    geom_boxplot(alpha = .2, position=position_dodge(0.8), lwd=0.1, colour=colbox, fill=colbox) +
    # geom_line(group=c(
    #  subset(crat_data,time1 == levels(crat_data$time1)[1])[,"patient1"],
    #  subset(crat_data, time1 == levels(crat_data$time1)[2])[,"patient1"]),
    #  position=position_dodge(0.8),alpha = .3,size=1) +
    geom_dotplot(binaxis="y",stackdir='center', dotsize=0.8,position=position_dodge(0.8), fill=batch3, colour=batch3) +
    xlab("")+
    ylab("norm. txn counts")+
    ggtitle(paste(i,sep=" "))
  g
  l <- length(intersect(goi_L$ENSEMBL,rownames(res)))
  assign(paste("p",n,sep=""),g)
}

ggarrange(p1+ggpubr::rotate_x_text(),
          p2+rremove("ylab")+ggpubr::rotate_x_text(),
          p3+rremove("ylab")+ggpubr::rotate_x_text(),
          common.legend = TRUE, legend = "right", nrow = 1, ncol= 3)


### Volcanos:







EnhancedVolcano(table_AKI,
                title = 'AKI: zero-hour vs. post_TX',
                subtitle="",
                lab = table_AKI$SYMBOL,
                selectLab = goi$SYMBOL,
                xlim = c(-2, 2),
                ylim = c(0, 2),
                FCcutoff = 0.5,
                # pCutoff = 0.05,
                x = 'logFC',
                boxedLabels = TRUE,
                y = 'AveExpr')

volcano_names <- ifelse(abs(palmieri_fit_PBx$coefficients)>=1, 
                        palmieri_fit_PBx$genes$SYMBOL, NA)

colnames(table_AKI)
is.numeric(table_AKI$logFC)
is.numeric(table_AKI$AveExpr_AKI)
head(1/(table_AKI$AveExpr))
-log(0.1)
table_AKI$SYMBOL
table_AKI$neg_AveExpr <- 10^(-table_AKI$AveExpr)

pdf(paste("2021_08_31 volcano_AKI2",sample,".pdf",sep="")) 
EnhancedVolcano(table_AKI,
                title = 'AKI: zero-hour vs. post_TX',
                subtitle="",
                lab = table_AKI$SYMBOL,
                selectLab = goi$SYMBOL,
                xlim = c(-2, 2),
                ylim = c(0, 2),
                FCcutoff = 0.5,
                # pCutoff = 0.05,
                x = 'logFC',
                boxedLabels = TRUE,
                y = 'AveExpr')

plot(table_AKI$logFC~table_AKI$AveExpr, main="AKI")

dev.off()



######## annotation fix
ah <- AnnotationHub()
query(ah,"EnsDb.Hsapiens.v99")
rappdirs::user_cache_dir(appname="AnnotationHub")
tools::R_user_dir("AnnotationHub", which="cache")
moveFiles<-function(package){
  olddir <- path.expand(rappdirs::user_cache_dir(appname=package))
  newdir <- tools::R_user_dir(package, which="cache")
  dir.create(path=newdir, recursive=TRUE)
  files <- list.files(olddir, full.names =TRUE)
  moveres <- vapply(files,
                    FUN=function(fl){
                      filename = basename(fl)
                      newname = file.path(newdir, filename)
                      file.rename(fl, newname)
                    },
                    FUN.VALUE = logical(1))
  if(all(moveres)) unlink(olddir, recursive=TRUE)
}
package="AnnotationHub"
moveFiles(package)

query(ah,"EnsDb.Hsapiens.v99")
query(hub, c("EnsDb", "sapiens", "99"))
query(hub, c("EnsDb", "sapiens"))
query(hub, c("EnsDb", "sapiens","104"))

AnnotationHub("AH95744")
?AnnotationHub
