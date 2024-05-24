.libPaths()

# Example dataset
library(DESeq2)
sessionInfo()
set.seed(123)
dds <- makeExampleDESeqDataSet(n = 10000, m = 100)
dds$group <- factor(rep(rep(c("X","Y"),each=10),5))
save(dds,file="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/data/BLAS/dds_sample_L.dds")

# Run for each BLAS
# restart R
sessionInfo()$BLAS
sessionInfo()$LAPACK
blas <- basename(dirname(sessionInfo()$BLAS))
library(DESeq2)
start.time <- Sys.time()
load(file="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/data/BLAS/dds_sample_L.dds")
design(dds) <-  ~condition+group+condition:group
dds <- DESeq(dds)
summary(results(dds))
end.time <- Sys.time()
time.taken <- end.time - start.time
cat(blas,time.taken)
metadata(dds)$BLAS <- sessionInfo()$BLAS
metadata(dds)$LAPACK <- sessionInfo()$LAPACK

save(dds,file=paste("/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/data/BLAS/dds_sample_L_",blas,".dds",sep=""))


# L Dataset
## Times:
# mkl: 11.90276 sec
# openBLAS: 16.71817 sec
# ATLAS: 9.826937 sec
# BLAS: 12.04894 sec





# Standard Dataset
## Times:
# mkl: 0.9552917 sec
# openBLAS: 1.005392 sec
# ATLAS: 0.8985503 sec
# BLAS: 0.9572477 sec

load(file="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/data/BLAS/dds_sample_blas.dds")




save(dds_l_list,file="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/data/BLAS/dds_sample_l_list")



# Variant 1 (Standard)

data <- "/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/data"
design = ~genotype+treatment+genotype:treatment
dds <- DESeqDataSet(gse, design = design)
keep.sn <- rowSums(counts(dds)) >= 10
dds <- dds[keep.sn,]
dds <- DESeq(dds)
summary(results(dds))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# dds_atlas <- dds
dds_BLAS <- dds
save(dds_BLAS,file="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/data/BLAS/dds_BLAS.dds")

## Times:
# mkl: 1.06 min
# openBLAS: 5.69 min
# ATLAS: 1.23 min
# BLAS: 1.37 min

load(file="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/data/BLAS/dds_openBLAS.dds")





...Relevent codes...



# Variant 2 (+Replicates)
data <- "/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/data"
load(file=paste(data,"tximeta.txm", sep="/"))
design = ~genotype+treatment+genotype:treatment
dds <- DESeqDataSet(gse, design = design)
keep.sn <- rowSums(counts(dds)) >= 10
dds <- dds[keep.sn,]
dds <- collapseReplicates(dds, dds$samplename, dds$names)
dds <- DESeq(dds)
summary(results(dds))

# Variant 3 (+Experiment)
data <- "/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/data"
load(file=paste(data,"tximeta.txm", sep="/"))
design = ~experiment+genotype+treatment+genotype:treatment
dds <- DESeqDataSet(gse, design = design)
keep.sn <- rowSums(counts(dds)) >= 10
dds <- dds[keep.sn,]
dds <- collapseReplicates(dds, dds$samplename, dds$names)
dds <- DESeq(dds)
summary(results(dds))

# Variant 1
# -- replacing outliers and refitting for 84 genes
# outliers [1]       : 11, 0.032%
# low counts [2]     : 7891, 23%

# Variant 2 (+Replicates)
# -- replacing outliers and refitting for 66 genes
# outliers [1]       : 9, 0.027%
# low counts [2]     : 7235, 21%

# Variant 3 (+Experiment)
# -- replacing outliers and refitting for 18326 genes
# outliers [1]       : 28, 0.083%
# low counts [2]     : 6572, 19%

dds
gse

resultsNames(dds)
counts(dds, normalized=TRUE)[1,1:5]
mcols(dds)[15000,]
colData(dds)[20]
packageVersion("DESeq2")




