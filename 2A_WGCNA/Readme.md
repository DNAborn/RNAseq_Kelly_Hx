WGCNA
================
Kelterborn
2024-03-07

- [0. Load](#0-load)
  - [- R](#--r)
- [1. WGCNA](#1-wgcna)
  - [-load dds](#-load-dds)
  - [(-pickSoftThreshold extern)](#-picksoftthreshold-extern)
  - [-pickSoftThreshold](#-picksoftthreshold)
  - [-network construction](#-network-construction)
  - [-Module Eigengenes](#-module-eigengenes)
  - [-Intramodular analysis: Identifying driver
    genes](#-intramodular-analysis-identifying-driver-genes)
  - [-TS Analysis](#-ts-analysis)
  - [-GO terms enrichment](#-go-terms-enrichment)
  - [-module to sample](#-module-to-sample)
- [Export into dds](#export-into-dds)

# 0. Load

## - R

BiocManager::install()

BiocManager::install(“CorLevelPlot”)

# 1. WGCNA

## -load dds

``` r
load(file=paste(data,"deseq2.dds", sep="/"))
s75 <- (nrow(colData(dds))*0.75) %>% round()
dds75 <- dds[rowSums(counts(dds) >= 15) >= 66,]
nrow(dds75) # 15687 genes left
```

    ## [1] 15687

``` r
vsd <- vst(dds75, blind = FALSE) #transform while accounting for design 
# counts <- counts(dds, normalized=TRUE)
# mcols(dds)$SYMBOL %>% head()
# input_mat = t(counts)
# input_mat[1:10,1:10]

colData <- colData(dds75)
norm.counts <- assay(vsd) %>% 
  t()
dim(norm.counts)
```

    ## [1]    88 15687

``` r
## outliners?
gsg <- goodSamplesGenes(norm.counts)
```

    ##  Flagging genes and samples with too many missing values...
    ##   ..step 1

``` r
summary(gsg)
```

    ##             Length Class  Mode   
    ## goodGenes   15687  -none- logical
    ## goodSamples    88  -none- logical
    ## allOK           1  -none- logical

``` r
summary(gsg$goodGenes)
```

    ##    Mode    TRUE 
    ## logical   15687

``` r
gsg$allOK
```

    ## [1] TRUE

``` r
# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(norm.counts), method = "average")
plot(htree) # S46, S50, S58?
```

<img src="README_files/figure-gfm/load_dds-1.png" width="100%" />

``` r
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
```

<img src="README_files/figure-gfm/load_dds-2.png" width="100%" />

## (-pickSoftThreshold extern)

## -pickSoftThreshold

    ## pickSoftThreshold: will use block size 2851.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 2851 of 15687
    ##    ..working on genes 2852 through 5702 of 15687
    ##    ..working on genes 5703 through 8553 of 15687
    ##    ..working on genes 8554 through 11404 of 15687
    ##    ..working on genes 11405 through 14255 of 15687
    ##    ..working on genes 14256 through 15687 of 15687
    ##    Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1  0.09950 77.5000          0.829  7850.0   7850.00   7950
    ## 2      2  0.18000  3.3600          0.907  4720.0   4750.00   5370
    ## 3      3  0.17500  1.3100          0.876  3160.0   3200.00   4130
    ## 4      4  0.10100  0.5910          0.841  2260.0   2310.00   3370
    ## 5      5  0.01900  0.1840          0.802  1700.0   1730.00   2850
    ## 6      6  0.00614 -0.0852          0.788  1320.0   1340.00   2460
    ## 7      7  0.07700 -0.2790          0.794  1050.0   1060.00   2170
    ## 8      8  0.20700 -0.4430          0.822   858.0    853.00   1930
    ## 9      9  0.32900 -0.5700          0.853   710.0    696.00   1730
    ## 10    10  0.42500 -0.6760          0.878   597.0    574.00   1570
    ## 11    12  0.56500 -0.8560          0.914   435.0    403.00   1310
    ## 12    14  0.64600 -0.9870          0.937   329.0    291.00   1120
    ## 13    16  0.69400 -1.0900          0.950   255.0    215.00    969
    ## 14    18  0.73100 -1.1700          0.962   202.0    162.00    849
    ## 15    20  0.76000 -1.2400          0.972   163.0    123.00    751
    ## 16    22  0.78500 -1.2900          0.981   134.0     95.00    670
    ## 17    24  0.79800 -1.3300          0.983   111.0     74.10    602
    ## 18    25  0.80000 -1.3600          0.981   101.0     65.50    571
    ## 19    26  0.80500 -1.3800          0.983    93.0     58.40    543
    ## 20    28  0.81200 -1.4200          0.981    78.7     46.40    493
    ## 21    30  0.82200 -1.4600          0.981    67.2     37.30    450
    ## 22    32  0.83100 -1.4900          0.984    57.8     30.00    412
    ## 23    34  0.83700 -1.5000          0.986    50.0     24.30    379
    ## 24    36  0.84700 -1.5200          0.988    43.5     19.90    349
    ## 25    38  0.85500 -1.5400          0.990    38.1     16.40    323
    ## 26    40  0.86300 -1.5600          0.991    33.5     13.50    299
    ## 27    42  0.86600 -1.5700          0.990    29.6     11.30    278
    ## 28    44  0.87000 -1.5900          0.991    26.3      9.39    259
    ## 29    46  0.86800 -1.6000          0.989    23.4      7.88    242
    ## 30    48  0.86600 -1.6200          0.986    21.0      6.59    226
    ## 31    50  0.87200 -1.6200          0.988    18.8      5.56    212

![](README_files/figure-gfm/pickSoftThreshold-1.png)<!-- -->

    ##    Power  SFT.R.sq     slope truncated.R.sq   mean.k. median.k.   max.k.
    ## 18    25 0.7998721 -1.358033      0.9810872 101.47754  65.53304 571.3261
    ## 19    26 0.8052342 -1.378630      0.9825375  93.02945  58.38016 543.3971

## -network construction

``` r
# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 26
temp_cor <- cor
cor <- WGCNA::cor

# memory estimate w.r.t blocksize
# bwnet <- blockwiseModules(norm.counts,
#                  maxBlockSize = 15000,
#                  TOMType = "signed",
#                  power = soft_power,
#                  mergeCutHeight = 0.25,
#                  numericLabels = FALSE,
#                  randomSeed = 1234,
#                  verbose = 3)
# 
# cor <- temp_cor
# 
# save(bwnet,file=paste(data,"bwnet.RDS", sep="/"))

# TS
cor <- WGCNA::cor
bwnet <- blockwiseModules(norm.counts,               

                          # == Adjacency Function ==
                          power = soft_power,                
                          networkType = "signed",

                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 40000,

                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,

                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = paste(data,"blockwiseTOM", sep="/"),

                          # == Output Options
                          numericLabels = F,
                          verbose = 3)
```

    ##  Calculating module eigengenes block-wise from all genes
    ##    Flagging genes and samples with too many missing values...
    ##     ..step 1
    ##  ..Working on block 1 .
    ##     TOM calculation: adjacency..
    ##     ..will use 40 parallel threads.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 1 into file /mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/data/blockwiseTOM-block.1.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 1 genes from module 2 because their KME is too low.
    ##      ..removing 1 genes from module 4 because their KME is too low.
    ##      ..removing 1 genes from module 12 because their KME is too low.
    ##      ..removing 1 genes from module 24 because their KME is too low.
    ##      ..removing 1 genes from module 49 because their KME is too low.
    ##  ..merging modules that are too close..
    ##      mergeCloseModules: Merging modules whose distance is less than 0.25
    ##        Calculating new MEs...

``` r
cor <- temp_cor

save(bwnet,file=paste(data,"bwnet_TS.RDS", sep="/"))
```

## -Module Eigengenes

``` r
load(file=paste(data,"bwnet_TS.RDS", sep="/"))

module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:400px; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEgreenyellow
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEblack
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEblue
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEpink
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEgreen
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEpurple
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEbrown
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEturquoise
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEmagenta
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEred
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEyellow
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEgrey
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RNA_P2041_S37
</td>
<td style="text-align:right;">
-0.0225829
</td>
<td style="text-align:right;">
-0.0724136
</td>
<td style="text-align:right;">
-0.1042940
</td>
<td style="text-align:right;">
-0.0228416
</td>
<td style="text-align:right;">
0.0930260
</td>
<td style="text-align:right;">
-0.0461527
</td>
<td style="text-align:right;">
0.0633819
</td>
<td style="text-align:right;">
0.1254729
</td>
<td style="text-align:right;">
-0.0611301
</td>
<td style="text-align:right;">
-0.0418531
</td>
<td style="text-align:right;">
-0.1245557
</td>
<td style="text-align:right;">
-0.0326363
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S38
</td>
<td style="text-align:right;">
0.0330771
</td>
<td style="text-align:right;">
-0.0860900
</td>
<td style="text-align:right;">
-0.0729152
</td>
<td style="text-align:right;">
-0.0712031
</td>
<td style="text-align:right;">
-0.0487715
</td>
<td style="text-align:right;">
-0.0840368
</td>
<td style="text-align:right;">
0.0949757
</td>
<td style="text-align:right;">
0.0481258
</td>
<td style="text-align:right;">
0.0039414
</td>
<td style="text-align:right;">
0.0733570
</td>
<td style="text-align:right;">
-0.0041574
</td>
<td style="text-align:right;">
0.1533559
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S39
</td>
<td style="text-align:right;">
0.0696383
</td>
<td style="text-align:right;">
-0.0952882
</td>
<td style="text-align:right;">
-0.1201861
</td>
<td style="text-align:right;">
-0.0170000
</td>
<td style="text-align:right;">
0.1070351
</td>
<td style="text-align:right;">
-0.0252364
</td>
<td style="text-align:right;">
0.0986152
</td>
<td style="text-align:right;">
0.1619557
</td>
<td style="text-align:right;">
-0.0488475
</td>
<td style="text-align:right;">
-0.0425104
</td>
<td style="text-align:right;">
-0.1558429
</td>
<td style="text-align:right;">
0.0381960
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S40
</td>
<td style="text-align:right;">
0.0423960
</td>
<td style="text-align:right;">
-0.0881383
</td>
<td style="text-align:right;">
-0.1043251
</td>
<td style="text-align:right;">
-0.0398706
</td>
<td style="text-align:right;">
0.0547433
</td>
<td style="text-align:right;">
-0.0354322
</td>
<td style="text-align:right;">
0.0901087
</td>
<td style="text-align:right;">
0.1156110
</td>
<td style="text-align:right;">
-0.0217501
</td>
<td style="text-align:right;">
-0.0036398
</td>
<td style="text-align:right;">
-0.0994400
</td>
<td style="text-align:right;">
0.0191356
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S41
</td>
<td style="text-align:right;">
0.1509036
</td>
<td style="text-align:right;">
0.0496434
</td>
<td style="text-align:right;">
-0.0938437
</td>
<td style="text-align:right;">
0.0644966
</td>
<td style="text-align:right;">
0.1691126
</td>
<td style="text-align:right;">
-0.0449969
</td>
<td style="text-align:right;">
-0.0258649
</td>
<td style="text-align:right;">
0.1298936
</td>
<td style="text-align:right;">
-0.0194218
</td>
<td style="text-align:right;">
-0.2368260
</td>
<td style="text-align:right;">
-0.2280376
</td>
<td style="text-align:right;">
0.0329747
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S42
</td>
<td style="text-align:right;">
0.0115133
</td>
<td style="text-align:right;">
-0.1101813
</td>
<td style="text-align:right;">
-0.0938943
</td>
<td style="text-align:right;">
-0.0781219
</td>
<td style="text-align:right;">
-0.0014160
</td>
<td style="text-align:right;">
-0.0632929
</td>
<td style="text-align:right;">
0.1084979
</td>
<td style="text-align:right;">
0.0938583
</td>
<td style="text-align:right;">
-0.0259494
</td>
<td style="text-align:right;">
0.0666340
</td>
<td style="text-align:right;">
-0.0461265
</td>
<td style="text-align:right;">
0.0081877
</td>
</tr>
</tbody>
</table>

</div>

``` r
# get number of genes for each module
table(bwnet$colors)
```

    ## 
    ##       black        blue       brown       green greenyellow        grey 
    ##        1195        2695        2424        1697          98         425 
    ##     magenta        pink      purple         red   turquoise      yellow 
    ##         167         870         135        1417        2804        1760

``` r
# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)
```

<img src="README_files/figure-gfm/eigengenes-1.png" width="100%" />

``` r
# grey module = all genes that doesn't fall into other modules were assigned to the grey module



# 6A. Relate modules to traits --------------------------------------------------
# module trait associations



# create traits file - binarize categorical variables
traits <- colData$treatment_bin <- ifelse(grepl('Hx', colData$treatment), 1, 0)

# binarize categorical variables

colData$genotype %>% levels()
```

    ## [1] "Kelly" "HIF1A" "HIF2A" "HIF1B"

``` r
genotype_bin <- binarizeCategoricalColumns(colData$genotype,
                           includePairwise = FALSE,
                           includeLevelVsAll = TRUE,
                           dropFirstLevelVsAll = FALSE,
                           minCount = 1)
colnames(genotype_bin) <- levels(colData$genotype)

condition_bin <- binarizeCategoricalColumns(colData$condition,
                           includePairwise = FALSE,
                           includeLevelVsAll = TRUE,
                           dropFirstLevelVsAll = FALSE,
                           minCount = 1)
colnames(condition_bin) <- levels(colData$condition)

traits <- cbind(traits, genotype_bin,condition_bin)
rownames(traits) <- rownames(colData)
dim(traits)
```

    ## [1] 88 13

``` r
orig.colnames <- colnames(traits)
colnames(traits)[1] <- c("Hypoxia")

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:400px; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
Row.names
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEgreenyellow
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEblack
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEblue
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEpink
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEgreen
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEpurple
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEbrown
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEturquoise
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEmagenta
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEred
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEyellow
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEgrey
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Hypoxia
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Kelly
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
HIF1A
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
HIF2A
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
HIF1B
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Kelly_Nx
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Kelly_Hx
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
HIF1A_Nx
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
HIF1A_Hx
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
HIF2A_Nx
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
HIF2A_Hx
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
HIF1B_Nx
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
HIF1B_Hx
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RNA_P2041_S37
</td>
<td style="text-align:right;">
-0.0225829
</td>
<td style="text-align:right;">
-0.0724136
</td>
<td style="text-align:right;">
-0.1042940
</td>
<td style="text-align:right;">
-0.0228416
</td>
<td style="text-align:right;">
0.0930260
</td>
<td style="text-align:right;">
-0.0461527
</td>
<td style="text-align:right;">
0.0633819
</td>
<td style="text-align:right;">
0.1254729
</td>
<td style="text-align:right;">
-0.0611301
</td>
<td style="text-align:right;">
-0.0418531
</td>
<td style="text-align:right;">
-0.1245557
</td>
<td style="text-align:right;">
-0.0326363
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S38
</td>
<td style="text-align:right;">
0.0330771
</td>
<td style="text-align:right;">
-0.0860900
</td>
<td style="text-align:right;">
-0.0729152
</td>
<td style="text-align:right;">
-0.0712031
</td>
<td style="text-align:right;">
-0.0487715
</td>
<td style="text-align:right;">
-0.0840368
</td>
<td style="text-align:right;">
0.0949757
</td>
<td style="text-align:right;">
0.0481258
</td>
<td style="text-align:right;">
0.0039414
</td>
<td style="text-align:right;">
0.0733570
</td>
<td style="text-align:right;">
-0.0041574
</td>
<td style="text-align:right;">
0.1533559
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S39
</td>
<td style="text-align:right;">
0.0696383
</td>
<td style="text-align:right;">
-0.0952882
</td>
<td style="text-align:right;">
-0.1201861
</td>
<td style="text-align:right;">
-0.0170000
</td>
<td style="text-align:right;">
0.1070351
</td>
<td style="text-align:right;">
-0.0252364
</td>
<td style="text-align:right;">
0.0986152
</td>
<td style="text-align:right;">
0.1619557
</td>
<td style="text-align:right;">
-0.0488475
</td>
<td style="text-align:right;">
-0.0425104
</td>
<td style="text-align:right;">
-0.1558429
</td>
<td style="text-align:right;">
0.0381960
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S40
</td>
<td style="text-align:right;">
0.0423960
</td>
<td style="text-align:right;">
-0.0881383
</td>
<td style="text-align:right;">
-0.1043251
</td>
<td style="text-align:right;">
-0.0398706
</td>
<td style="text-align:right;">
0.0547433
</td>
<td style="text-align:right;">
-0.0354322
</td>
<td style="text-align:right;">
0.0901087
</td>
<td style="text-align:right;">
0.1156110
</td>
<td style="text-align:right;">
-0.0217501
</td>
<td style="text-align:right;">
-0.0036398
</td>
<td style="text-align:right;">
-0.0994400
</td>
<td style="text-align:right;">
0.0191356
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S41
</td>
<td style="text-align:right;">
0.1509036
</td>
<td style="text-align:right;">
0.0496434
</td>
<td style="text-align:right;">
-0.0938437
</td>
<td style="text-align:right;">
0.0644966
</td>
<td style="text-align:right;">
0.1691126
</td>
<td style="text-align:right;">
-0.0449969
</td>
<td style="text-align:right;">
-0.0258649
</td>
<td style="text-align:right;">
0.1298936
</td>
<td style="text-align:right;">
-0.0194218
</td>
<td style="text-align:right;">
-0.2368260
</td>
<td style="text-align:right;">
-0.2280376
</td>
<td style="text-align:right;">
0.0329747
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S42
</td>
<td style="text-align:right;">
0.0115133
</td>
<td style="text-align:right;">
-0.1101813
</td>
<td style="text-align:right;">
-0.0938943
</td>
<td style="text-align:right;">
-0.0781219
</td>
<td style="text-align:right;">
-0.0014160
</td>
<td style="text-align:right;">
-0.0632929
</td>
<td style="text-align:right;">
0.1084979
</td>
<td style="text-align:right;">
0.0938583
</td>
<td style="text-align:right;">
-0.0259494
</td>
<td style="text-align:right;">
0.0666340
</td>
<td style="text-align:right;">
-0.0461265
</td>
<td style="text-align:right;">
0.0081877
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>

</div>

``` r
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')
dim(heatmap.data)
```

    ## [1] 88 25

``` r
MEs <- heatmap.data %>% colnames() %>% str_detect(pattern="ME") %>% sum()
max <- heatmap.data %>% ncol()
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[(MEs+1):max],
             y = names(heatmap.data)[1:MEs],
             col = viridis(100, option='plasma'))
```

<img src="README_files/figure-gfm/eigengenes-2.png" width="100%" />

``` r
             # col = c("blue1", "skyblue", "white", "pink", "red"))

module.gene.mapping <- as.data.frame(bwnet$colors)

# Genes related to Hypoxia
module.gene.mapping %>% 
  dplyr::filter(`bwnet$colors` == 'turquoise') %>% 
  rownames() %>% head() %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:400px; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
ENSG00000001497
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000001630
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000002549
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000002919
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000003056
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000003147
</td>
</tr>
</tbody>
</table>

</div>

## -Intramodular analysis: Identifying driver genes

``` r
module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure[1:10,1:10] %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:400px; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000000003
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000000419
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000000457
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000000460
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000001084
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000001167
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000001460
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000001461
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000001497
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000001617
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
MEgreenyellow
</td>
<td style="text-align:right;">
-0.1159951
</td>
<td style="text-align:right;">
-0.5703501
</td>
<td style="text-align:right;">
-0.0813194
</td>
<td style="text-align:right;">
-0.1325352
</td>
<td style="text-align:right;">
-0.2492833
</td>
<td style="text-align:right;">
-0.0509276
</td>
<td style="text-align:right;">
0.2931047
</td>
<td style="text-align:right;">
0.2975899
</td>
<td style="text-align:right;">
0.2282183
</td>
<td style="text-align:right;">
0.3583264
</td>
</tr>
<tr>
<td style="text-align:left;">
MEblack
</td>
<td style="text-align:right;">
0.1250298
</td>
<td style="text-align:right;">
-0.3929736
</td>
<td style="text-align:right;">
0.4127026
</td>
<td style="text-align:right;">
-0.7173010
</td>
<td style="text-align:right;">
-0.0947735
</td>
<td style="text-align:right;">
-0.5742677
</td>
<td style="text-align:right;">
-0.0906712
</td>
<td style="text-align:right;">
0.9262625
</td>
<td style="text-align:right;">
-0.6178117
</td>
<td style="text-align:right;">
0.6491368
</td>
</tr>
<tr>
<td style="text-align:left;">
MEblue
</td>
<td style="text-align:right;">
0.0916580
</td>
<td style="text-align:right;">
-0.1653780
</td>
<td style="text-align:right;">
0.8184299
</td>
<td style="text-align:right;">
-0.9142004
</td>
<td style="text-align:right;">
-0.4423330
</td>
<td style="text-align:right;">
-0.4335178
</td>
<td style="text-align:right;">
0.2458620
</td>
<td style="text-align:right;">
0.7183112
</td>
<td style="text-align:right;">
-0.9166825
</td>
<td style="text-align:right;">
-0.0919951
</td>
</tr>
<tr>
<td style="text-align:left;">
MEpink
</td>
<td style="text-align:right;">
0.4644931
</td>
<td style="text-align:right;">
-0.7127120
</td>
<td style="text-align:right;">
0.7001243
</td>
<td style="text-align:right;">
-0.7540650
</td>
<td style="text-align:right;">
-0.7814765
</td>
<td style="text-align:right;">
-0.7914143
</td>
<td style="text-align:right;">
0.5761266
</td>
<td style="text-align:right;">
0.7938665
</td>
<td style="text-align:right;">
-0.6481853
</td>
<td style="text-align:right;">
0.1757105
</td>
</tr>
<tr>
<td style="text-align:left;">
MEgreen
</td>
<td style="text-align:right;">
0.5999831
</td>
<td style="text-align:right;">
-0.7466202
</td>
<td style="text-align:right;">
0.0322159
</td>
<td style="text-align:right;">
-0.0078138
</td>
<td style="text-align:right;">
-0.4649663
</td>
<td style="text-align:right;">
-0.6623493
</td>
<td style="text-align:right;">
0.3818358
</td>
<td style="text-align:right;">
0.3498492
</td>
<td style="text-align:right;">
0.1059467
</td>
<td style="text-align:right;">
0.4565627
</td>
</tr>
<tr>
<td style="text-align:left;">
MEpurple
</td>
<td style="text-align:right;">
0.6751690
</td>
<td style="text-align:right;">
-0.1445500
</td>
<td style="text-align:right;">
0.3694312
</td>
<td style="text-align:right;">
-0.1006530
</td>
<td style="text-align:right;">
-0.2127415
</td>
<td style="text-align:right;">
-0.5287967
</td>
<td style="text-align:right;">
0.0372571
</td>
<td style="text-align:right;">
0.2223280
</td>
<td style="text-align:right;">
-0.3682379
</td>
<td style="text-align:right;">
0.0551029
</td>
</tr>
<tr>
<td style="text-align:left;">
MEbrown
</td>
<td style="text-align:right;">
-0.2745423
</td>
<td style="text-align:right;">
0.4925550
</td>
<td style="text-align:right;">
-0.7471740
</td>
<td style="text-align:right;">
0.9157300
</td>
<td style="text-align:right;">
0.5410580
</td>
<td style="text-align:right;">
0.6977581
</td>
<td style="text-align:right;">
-0.3108125
</td>
<td style="text-align:right;">
-0.9145102
</td>
<td style="text-align:right;">
0.8355361
</td>
<td style="text-align:right;">
-0.2477169
</td>
</tr>
<tr>
<td style="text-align:left;">
MEturquoise
</td>
<td style="text-align:right;">
0.1533269
</td>
<td style="text-align:right;">
-0.1579351
</td>
<td style="text-align:right;">
-0.6522148
</td>
<td style="text-align:right;">
0.7960361
</td>
<td style="text-align:right;">
0.0949327
</td>
<td style="text-align:right;">
0.1627130
</td>
<td style="text-align:right;">
0.0635655
</td>
<td style="text-align:right;">
-0.5612057
</td>
<td style="text-align:right;">
0.8505339
</td>
<td style="text-align:right;">
0.1059646
</td>
</tr>
<tr>
<td style="text-align:left;">
MEmagenta
</td>
<td style="text-align:right;">
-0.3870149
</td>
<td style="text-align:right;">
0.3452500
</td>
<td style="text-align:right;">
-0.4214512
</td>
<td style="text-align:right;">
0.1424331
</td>
<td style="text-align:right;">
0.8185590
</td>
<td style="text-align:right;">
0.2844966
</td>
<td style="text-align:right;">
-0.7668347
</td>
<td style="text-align:right;">
0.0794656
</td>
<td style="text-align:right;">
0.1568495
</td>
<td style="text-align:right;">
0.5614501
</td>
</tr>
<tr>
<td style="text-align:left;">
MEred
</td>
<td style="text-align:right;">
-0.4964423
</td>
<td style="text-align:right;">
0.8627144
</td>
<td style="text-align:right;">
-0.2984007
</td>
<td style="text-align:right;">
0.3971874
</td>
<td style="text-align:right;">
0.6210349
</td>
<td style="text-align:right;">
0.7628490
</td>
<td style="text-align:right;">
-0.5013851
</td>
<td style="text-align:right;">
-0.6580135
</td>
<td style="text-align:right;">
0.2170335
</td>
<td style="text-align:right;">
-0.4850426
</td>
</tr>
</tbody>
</table>

</div>

``` r
module.membership.measure.pvals[1:10,1:10] %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:400px; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000000003
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000000419
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000000457
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000000460
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000001084
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000001167
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000001460
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000001461
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000001497
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
ENSG00000001617
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
MEgreenyellow
</td>
<td style="text-align:right;">
0.2818335
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.4513386
</td>
<td style="text-align:right;">
0.2183409
</td>
<td style="text-align:right;">
0.0191729
</td>
<td style="text-align:right;">
0.6374845
</td>
<td style="text-align:right;">
0.0055816
</td>
<td style="text-align:right;">
0.0048653
</td>
<td style="text-align:right;">
0.0324701
</td>
<td style="text-align:right;">
0.0006082
</td>
</tr>
<tr>
<td style="text-align:left;">
MEblack
</td>
<td style="text-align:right;">
0.2457745
</td>
<td style="text-align:right;">
0.0001523
</td>
<td style="text-align:right;">
0.0000645
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.3797691
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.4008289
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
</tr>
<tr>
<td style="text-align:left;">
MEblue
</td>
<td style="text-align:right;">
0.3957009
</td>
<td style="text-align:right;">
0.1236031
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000160
</td>
<td style="text-align:right;">
0.0000245
</td>
<td style="text-align:right;">
0.0209449
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.3939582
</td>
</tr>
<tr>
<td style="text-align:left;">
MEpink
</td>
<td style="text-align:right;">
0.0000051
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.1015240
</td>
</tr>
<tr>
<td style="text-align:left;">
MEgreen
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.7657273
</td>
<td style="text-align:right;">
0.9424009
</td>
<td style="text-align:right;">
0.0000050
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0002418
</td>
<td style="text-align:right;">
0.0008342
</td>
<td style="text-align:right;">
0.3258902
</td>
<td style="text-align:right;">
0.0000078
</td>
</tr>
<tr>
<td style="text-align:left;">
MEpurple
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.1790521
</td>
<td style="text-align:right;">
0.0003968
</td>
<td style="text-align:right;">
0.3507781
</td>
<td style="text-align:right;">
0.0465909
</td>
<td style="text-align:right;">
0.0000001
</td>
<td style="text-align:right;">
0.7303766
</td>
<td style="text-align:right;">
0.0373485
</td>
<td style="text-align:right;">
0.0004157
</td>
<td style="text-align:right;">
0.6101158
</td>
</tr>
<tr>
<td style="text-align:left;">
MEbrown
</td>
<td style="text-align:right;">
0.0096380
</td>
<td style="text-align:right;">
0.0000011
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000001
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0032050
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0199675
</td>
</tr>
<tr>
<td style="text-align:left;">
MEturquoise
</td>
<td style="text-align:right;">
0.1538058
</td>
<td style="text-align:right;">
0.1416648
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.3789657
</td>
<td style="text-align:right;">
0.1298555
</td>
<td style="text-align:right;">
0.5562859
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.3258077
</td>
</tr>
<tr>
<td style="text-align:left;">
MEmagenta
</td>
<td style="text-align:right;">
0.0001954
</td>
<td style="text-align:right;">
0.0009865
</td>
<td style="text-align:right;">
0.0000433
</td>
<td style="text-align:right;">
0.1855705
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0072225
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.4617550
</td>
<td style="text-align:right;">
0.1444578
</td>
<td style="text-align:right;">
0.0000000
</td>
</tr>
<tr>
<td style="text-align:left;">
MEred
</td>
<td style="text-align:right;">
0.0000009
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0047449
</td>
<td style="text-align:right;">
0.0001273
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0422426
</td>
<td style="text-align:right;">
0.0000017
</td>
</tr>
</tbody>
</table>

</div>

``` r
# Genes correlating with Hypoxia
gene.hypoxia.corr <- cor(norm.counts, traits$Hypoxia, use = 'p')
gene.hypoxia.corr.pvals <- corPvalueStudent(gene.hypoxia.corr, nSamples)
gene.hypoxia.corr %>% head() %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:400px; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
ENSG00000000003
</td>
<td style="text-align:right;">
0.1803307
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000000419
</td>
<td style="text-align:right;">
-0.3782877
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000000457
</td>
<td style="text-align:right;">
0.5623292
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000000460
</td>
<td style="text-align:right;">
-0.7737649
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000001084
</td>
<td style="text-align:right;">
-0.3085498
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000001167
</td>
<td style="text-align:right;">
-0.5922289
</td>
</tr>
</tbody>
</table>

</div>

``` r
# TOP 10 (pval) genes correlating with Hypoxia
gene.hypoxia.corr.pvals %>%
  log(base = 10) %>% abs() %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  arrange(desc(V1)) %>%
  head(10) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:400px; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
gene_id
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
V1
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
ENSG00000114268
</td>
<td style="text-align:right;">
63.25669
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000185633
</td>
<td style="text-align:right;">
62.61945
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000186918
</td>
<td style="text-align:right;">
58.82411
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000148926
</td>
<td style="text-align:right;">
56.08801
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000134107
</td>
<td style="text-align:right;">
55.26945
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000196968
</td>
<td style="text-align:right;">
53.40425
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000182379
</td>
<td style="text-align:right;">
53.31721
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000122884
</td>
<td style="text-align:right;">
53.26596
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000107819
</td>
<td style="text-align:right;">
51.72925
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000112715
</td>
<td style="text-align:right;">
51.27487
</td>
</tr>
</tbody>
</table>

</div>

## -TS Analysis

``` r
moduleLabelsAutomatic20 <- bwnet$colors   
moduleColorsAutomatic20 <- labels2colors(moduleLabelsAutomatic20)

# Dies ist die korrekte Funktion!!!
plotDendroAndColors(bwnet$dendrograms[[1]], moduleColorsAutomatic20,
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Cluster Dendrogram")

module_df <- data.frame(
  gene_id = names(bwnet$colors),
  colors = labels2colors(bwnet$colors)
)

length(moduleLabelsAutomatic20)
```

    ## [1] 15687

``` r
ncol(norm.counts)
```

    ## [1] 15687

``` r
MEs0 <- moduleEigengenes(norm.counts, moduleLabelsAutomatic20)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

plotEigengeneNetworks(MEs0, "", marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 
    4, 1, 2), cex.lab = 0.8, xLabelsAngle = 90)

# Add treatment names
MEs0$treatment <- colData$condition
# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
```

![](README_files/figure-gfm/TS_analysis-1.png)![](README_files/figure-gfm/TS_analysis-2.png)![](README_files/figure-gfm/TS_analysis-3.png)

## -GO terms enrichment

``` r
expr_universe <- rownames(dds)

# Get GO terms of top colors
modcols <- c("yellow","red")
module_go <- module_df[module_df$colors %in% modcols,] # Hypoxia
module_go_ens <- module_go$gene_id
go_enrich_test <- enrichGO(gene = module_go_ens,
                      universe = expr_universe,
                      OrgDb = "org.Hs.eg.db", 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.1, 
                      qvalueCutoff = 0.5)
go1 <- dotplot(clusterProfiler::simplify(go_enrich_test))+labs(title = paste(modcols, collapse=" & "))

modcols <- c("magenta","turquoise","brown")
module_go <- module_df[module_df$colors %in% modcols,] # Hypoxia
module_go_ens <- module_go$gene_id
go_enrich_test <- enrichGO(gene = module_go_ens,
                      universe = expr_universe,
                      OrgDb = "org.Hs.eg.db", 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.1, 
                      qvalueCutoff = 0.5)
go2 <- dotplot(clusterProfiler::simplify(go_enrich_test))+labs(title = paste(modcols, collapse=" & "))

modcols <- c("purple","green")
module_go <- module_df[module_df$colors %in% modcols,] # Hypoxia
module_go_ens <- module_go$gene_id
go_enrich_test <- enrichGO(gene = module_go_ens,
                      universe = expr_universe,
                      OrgDb = "org.Hs.eg.db", 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.1, 
                      qvalueCutoff = 0.5)
go3 <- dotplot(clusterProfiler::simplify(go_enrich_test))+labs(title = paste(modcols, collapse=" & "))

modcols <- c("pink","blue","black")
module_go <- module_df[module_df$colors %in% modcols,] # Hypoxia
module_go_ens <- module_go$gene_id
go_enrich_test <- enrichGO(gene = module_go_ens,
                      universe = expr_universe,
                      OrgDb = "org.Hs.eg.db", 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.1, 
                      qvalueCutoff = 0.5)
go4 <- dotplot(clusterProfiler::simplify(go_enrich_test))+labs(title = paste(modcols, collapse=" & "))

(go1 + go2) / (go3 + go4) +plot_layout(guides = "collect", axis_titles="collect")
```

<img src="README_files/figure-gfm/goa-1.png" width="100%" />

## -module to sample

``` r
# module Sample correlation

MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )
```

    ## # A tibble: 1,056 × 3
    ##    treatment name          value
    ##    <fct>     <chr>         <dbl>
    ##  1 Kelly_Nx  greenyellow -0.0226
    ##  2 Kelly_Nx  black       -0.0724
    ##  3 Kelly_Nx  blue        -0.104 
    ##  4 Kelly_Nx  pink        -0.0228
    ##  5 Kelly_Nx  green        0.0930
    ##  6 Kelly_Nx  purple      -0.0462
    ##  7 Kelly_Nx  brown        0.0634
    ##  8 Kelly_Nx  turquoise    0.125 
    ##  9 Kelly_Nx  magenta     -0.0611
    ## 10 Kelly_Nx  red         -0.0419
    ## # ℹ 1,046 more rows

``` r
mydata <- mtcars[, c(1,3,4,5,6,7)]
head(mydata)
```

    ##                    mpg disp  hp drat    wt  qsec
    ## Mazda RX4         21.0  160 110 3.90 2.620 16.46
    ## Mazda RX4 Wag     21.0  160 110 3.90 2.875 17.02
    ## Datsun 710        22.8  108  93 3.85 2.320 18.61
    ## Hornet 4 Drive    21.4  258 110 3.08 3.215 19.44
    ## Hornet Sportabout 18.7  360 175 3.15 3.440 17.02
    ## Valiant           18.1  225 105 2.76 3.460 20.22

``` r
cormat <- round(cor(mydata),2)
head(cormat)
```

    ##        mpg  disp    hp  drat    wt  qsec
    ## mpg   1.00 -0.85 -0.78  0.68 -0.87  0.42
    ## disp -0.85  1.00  0.79 -0.71  0.89 -0.43
    ## hp   -0.78  0.79  1.00 -0.45  0.66 -0.71
    ## drat  0.68 -0.71 -0.45  1.00 -0.71  0.09
    ## wt   -0.87  0.89  0.66 -0.71  1.00 -0.17
    ## qsec  0.42 -0.43 -0.71  0.09 -0.17  1.00

``` r
class(cormat)
```

    ## [1] "matrix" "array"

``` r
melted_cormat <- melt(cormat)
head(melted_cormat)
```

    ##   Var1 Var2 value
    ## 1  mpg  mpg  1.00
    ## 2 disp  mpg -0.85
    ## 3   hp  mpg -0.78
    ## 4 drat  mpg  0.68
    ## 5   wt  mpg -0.87
    ## 6 qsec  mpg  0.42

``` r
bwnet$MEs %>% class()
```

    ## [1] "data.frame"

``` r
bwnet$MEs %>% data.matrix() %>% head() %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:400px; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEgreenyellow
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEblack
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEblue
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEpink
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEgreen
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEpurple
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEbrown
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEturquoise
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEmagenta
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEred
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEyellow
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEgrey
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RNA_P2041_S37
</td>
<td style="text-align:right;">
-0.0225829
</td>
<td style="text-align:right;">
-0.0724136
</td>
<td style="text-align:right;">
-0.1042940
</td>
<td style="text-align:right;">
-0.0228416
</td>
<td style="text-align:right;">
0.0930260
</td>
<td style="text-align:right;">
-0.0461527
</td>
<td style="text-align:right;">
0.0633819
</td>
<td style="text-align:right;">
0.1254729
</td>
<td style="text-align:right;">
-0.0611301
</td>
<td style="text-align:right;">
-0.0418531
</td>
<td style="text-align:right;">
-0.1245557
</td>
<td style="text-align:right;">
-0.0326363
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S38
</td>
<td style="text-align:right;">
0.0330771
</td>
<td style="text-align:right;">
-0.0860900
</td>
<td style="text-align:right;">
-0.0729152
</td>
<td style="text-align:right;">
-0.0712031
</td>
<td style="text-align:right;">
-0.0487715
</td>
<td style="text-align:right;">
-0.0840368
</td>
<td style="text-align:right;">
0.0949757
</td>
<td style="text-align:right;">
0.0481258
</td>
<td style="text-align:right;">
0.0039414
</td>
<td style="text-align:right;">
0.0733570
</td>
<td style="text-align:right;">
-0.0041574
</td>
<td style="text-align:right;">
0.1533559
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S39
</td>
<td style="text-align:right;">
0.0696383
</td>
<td style="text-align:right;">
-0.0952882
</td>
<td style="text-align:right;">
-0.1201861
</td>
<td style="text-align:right;">
-0.0170000
</td>
<td style="text-align:right;">
0.1070351
</td>
<td style="text-align:right;">
-0.0252364
</td>
<td style="text-align:right;">
0.0986152
</td>
<td style="text-align:right;">
0.1619557
</td>
<td style="text-align:right;">
-0.0488475
</td>
<td style="text-align:right;">
-0.0425104
</td>
<td style="text-align:right;">
-0.1558429
</td>
<td style="text-align:right;">
0.0381960
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S40
</td>
<td style="text-align:right;">
0.0423960
</td>
<td style="text-align:right;">
-0.0881383
</td>
<td style="text-align:right;">
-0.1043251
</td>
<td style="text-align:right;">
-0.0398706
</td>
<td style="text-align:right;">
0.0547433
</td>
<td style="text-align:right;">
-0.0354322
</td>
<td style="text-align:right;">
0.0901087
</td>
<td style="text-align:right;">
0.1156110
</td>
<td style="text-align:right;">
-0.0217501
</td>
<td style="text-align:right;">
-0.0036398
</td>
<td style="text-align:right;">
-0.0994400
</td>
<td style="text-align:right;">
0.0191356
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S41
</td>
<td style="text-align:right;">
0.1509036
</td>
<td style="text-align:right;">
0.0496434
</td>
<td style="text-align:right;">
-0.0938437
</td>
<td style="text-align:right;">
0.0644966
</td>
<td style="text-align:right;">
0.1691126
</td>
<td style="text-align:right;">
-0.0449969
</td>
<td style="text-align:right;">
-0.0258649
</td>
<td style="text-align:right;">
0.1298936
</td>
<td style="text-align:right;">
-0.0194218
</td>
<td style="text-align:right;">
-0.2368260
</td>
<td style="text-align:right;">
-0.2280376
</td>
<td style="text-align:right;">
0.0329747
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S42
</td>
<td style="text-align:right;">
0.0115133
</td>
<td style="text-align:right;">
-0.1101813
</td>
<td style="text-align:right;">
-0.0938943
</td>
<td style="text-align:right;">
-0.0781219
</td>
<td style="text-align:right;">
-0.0014160
</td>
<td style="text-align:right;">
-0.0632929
</td>
<td style="text-align:right;">
0.1084979
</td>
<td style="text-align:right;">
0.0938583
</td>
<td style="text-align:right;">
-0.0259494
</td>
<td style="text-align:right;">
0.0666340
</td>
<td style="text-align:right;">
-0.0461265
</td>
<td style="text-align:right;">
0.0081877
</td>
</tr>
</tbody>
</table>

</div>

``` r
ME.heatmap <- bwnet$MEs %>% data.matrix() %>% melt()

ME.heatmap %>% ggplot(., aes(x=Var2, y=Var1, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
```

![](README_files/figure-gfm/module_sample-1.png)<!-- -->

``` r
pheatmap(bwnet$MEs)
```

![](README_files/figure-gfm/module_sample-2.png)<!-- -->

# Export into dds

``` r
load(file=paste(data,"deseq2.dds", sep="/"))
load(file=paste(data,"bwnet_TS.RDS", sep="/"))

# color per gene
mcols(dds)$colors <- bwnet$colors[match(mcols(dds)$gene_id, names(bwnet$colors))] %>% 
  factor() %>%
  relevel(ref="grey")

# color per sample
module.cols <- bind_cols(as.data.frame(colData(dds)),
                     as.data.frame(bwnet$MEs))
module.cols <- module.cols[,str_detect(colnames(module.cols),pattern="ME")]
colnames(module.cols) <- module.cols %>% colnames() %>% str_remove(pattern="ME")
i <- "greenyellow"
for (i in colnames(module.cols)){
colData(dds)$tmp <- module.cols[,i]
names(colData(dds))[names(colData(dds)) == 'tmp'] <- i
}

save(dds, file=paste(data,"deseq2_wgcna.dds", sep="/"))
```
