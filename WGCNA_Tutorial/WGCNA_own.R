---
  title: "RNA-Seq Kelly analysis"
author: "Kelterborn"
date: "2024-03-07"
output:
  md_document:
  variant: gfm
toc: true
always_allow_html: true
# editor_options: 
  chunk_output_type: console
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = T,
                      error=FALSE, warning=FALSE, message=FALSE)
```

# 0. Load
##  - R 
BiocManager::install()

BiocManager::install("CorLevelPlot") 
```{r librarys, include=FALSE}

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

ifelse(Sys.info()["sysname"]== "Linux",
       s <- "/mnt/s",
       s <- "S:")
dir <- paste(s,"AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all",sep="/")
gitdir <- paste(dir,"git_RNAseq_Kelly_Hx",sep="/")
data <- paste(dir,"data",sep="/")

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 40)

```

# 1. WGCNA
## -load dds

load(file=paste(data,"deseq2.dds", sep="/"))
s75 <- (nrow(colData(dds))*0.75) %>% round()
dds75 <- dds[rowSums(counts(dds) >= 15) >= 66,]
nrow(dds75) # 13284 genes

vsd <- vst(dds75, blind = FALSE) #transform while accounting for design 
# counts <- counts(dds, normalized=TRUE)
# mcols(dds)$SYMBOL %>% head()
# input_mat = t(counts)
# input_mat[1:10,1:10]

colData <- colData(dds75)
norm.counts <- assay(vsd) %>% 
  t()

## outliners?
gsg <- goodSamplesGenes(norm.counts)
summary(gsg)
summary(gsg$goodGenes)

gsg$allOK

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(norm.counts), method = "average")
plot(htree) # S46, S50, S58?

# PCA 
pca <- prcomp(norm.counts)
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


## -network construction

# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
n.cores <- 20
my.cluster <- parallel::makeCluster(n.cores,type = "FORK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

sft <- pickSoftThreshold(norm.counts[,1:300],
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

summary(norm.counts)

parallel::stopCluster(my.cluster)


sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)
# power 10 or 26

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 26
temp_cor <- cor
cor <- WGCNA::cor

# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 30000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

cor <- temp_cor

save(bwnet,file="data_analysis/bwnet.wgcna")

```

## -Module Eigengenes 
```{r}
load(file="data_analysis/bwnet.wgcna")

module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)

# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module



# 6A. Relate modules to traits --------------------------------------------------
# module trait associations



# create traits file - binarize categorical variables
traits <- colData %>% 
  mutate(disease_state_bin = ifelse(grepl('COVID', disease_state), 1, 0)) %>% 
  select(8)

traits <- colData %>% 
  mutate(treatment_bin = ifelse(grepl('Hx', treatment), 1, 0)) %>% 
  select("treatment_bin")

# binarize categorical variables

colData$genotype %>% levels()

genotype_bin <- binarizeCategoricalColumns(colData$genotype,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)

condition_bin <- binarizeCategoricalColumns(colData$condition,
                                            includePairwise = FALSE,
                                            includeLevelVsAll = TRUE,
                                            minCount = 1)

traits <- cbind(traits, genotype_bin,condition_bin)


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[18:22],
             y = names(heatmap.data)[1:17],
             col = c("blue1", "skyblue", "white", "pink", "red"))

```





