WGCNA
================
Kelterborn
2024-03-07

- [0. Load](#0-load)
  - [- Load R librarys](#--load-r-librarys)
- [1. WGCNA](#1-wgcna)
  - [-load dds](#-load-dds)
    - [-plot sample dist.](#-plot-sample-dist)
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

## - Load R librarys

# 1. WGCNA

## -load dds

### -plot sample dist.

    ##             Length Class  Mode   
    ## goodGenes   16268  -none- logical
    ## goodSamples    88  -none- logical
    ## allOK           1  -none- logical

    ## [1] TRUE

![](README_files/figure-gfm/unnamed-chunk-1-1.png)![](README_files/figure-gfm/unnamed-chunk-1-2.png)

## (-pickSoftThreshold extern)

## -pickSoftThreshold

    ## pickSoftThreshold: will use block size 2750.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 2750 of 16268
    ##    ..working on genes 2751 through 5500 of 16268
    ##    ..working on genes 5501 through 8250 of 16268
    ##    ..working on genes 8251 through 11000 of 16268
    ##    ..working on genes 11001 through 13750 of 16268
    ##    ..working on genes 13751 through 16268 of 16268
    ##    Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1  0.00640 16.8000          0.931  8140.0   8140.00   8240
    ## 2      2  0.15500  3.3900          0.824  4880.0   4910.00   5540
    ## 3      3  0.33300  1.8300          0.886  3250.0   3300.00   4190
    ## 4      4  0.18800  0.7890          0.838  2320.0   2360.00   3400
    ## 5      5  0.03790  0.2530          0.764  1740.0   1770.00   2870
    ## 6      6  0.00269 -0.0554          0.728  1340.0   1360.00   2480
    ## 7      7  0.07170 -0.2660          0.740  1070.0   1080.00   2180
    ## 8      8  0.20400 -0.4360          0.781   869.0    864.00   1940
    ## 9      9  0.32900 -0.5660          0.821   718.0    703.00   1740
    ## 10    10  0.42700 -0.6760          0.853   602.0    578.00   1580
    ## 11    12  0.56400 -0.8610          0.896   437.0    403.00   1320
    ## 12    14  0.64500 -0.9890          0.925   329.0    289.00   1120
    ## 13    16  0.69300 -1.0900          0.943   255.0    212.00    973
    ## 14    18  0.72900 -1.1700          0.956   202.0    159.00    851
    ## 15    20  0.75700 -1.2400          0.968   163.0    121.00    752
    ## 16    22  0.78300 -1.2800          0.978   133.0     92.60    671
    ## 17    24  0.79400 -1.3300          0.980   110.0     71.90    602
    ## 18    25  0.79700 -1.3600          0.979   101.0     63.60    572
    ## 19    26  0.80300 -1.3800          0.981    92.1     56.40    544
    ## 20    28  0.81000 -1.4200          0.980    77.9     44.70    494
    ## 21    30  0.81900 -1.4600          0.981    66.4     35.60    451
    ## 22    32  0.82800 -1.4800          0.984    57.0     28.70    413
    ## 23    34  0.83500 -1.5000          0.985    49.3     23.30    379
    ## 24    36  0.84500 -1.5200          0.988    42.9     19.00    350
    ## 25    38  0.85300 -1.5400          0.990    37.5     15.60    323
    ## 26    40  0.86000 -1.5500          0.991    33.0     12.90    300
    ## 27    42  0.86400 -1.5700          0.990    29.1     10.60    279
    ## 28    44  0.86800 -1.5800          0.991    25.8      8.84    260
    ## 29    46  0.86600 -1.6000          0.989    23.0      7.40    242
    ## 30    48  0.86400 -1.6200          0.986    20.6      6.17    227
    ## 31    50  0.87100 -1.6200          0.989    18.4      5.16    212

![](README_files/figure-gfm/pickSoftThreshold-1.png)<!-- -->

    ##    Power  SFT.R.sq     slope truncated.R.sq   mean.k. median.k.   max.k.
    ## 18    25 0.7972861 -1.356494      0.9793366 100.54602  63.64718 572.1999
    ## 19    26 0.8029531 -1.377792      0.9806661  92.11529  56.42164 544.1937

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
    ##      ..removing 1 genes from module 29 because their KME is too low.
    ##      ..removing 1 genes from module 42 because their KME is too low.
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
# head(module_eigengenes) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")

# get number of genes for each module
table(bwnet$colors)
```

    ## 
    ##       black        blue       brown        cyan       green greenyellow 
    ##        1019        2499        2426          76        1510         154 
    ##        grey     magenta        pink      purple         red      salmon 
    ##         450         496         962         177        1022         131 
    ##         tan   turquoise      yellow 
    ##         141        3124        2081

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

# colData$genotype %>% levels()

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
# dim(traits)
orig.colnames <- colnames(traits)
colnames(traits)[1] <- c("Hypoxia")

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

# head(heatmap.data) %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')
# dim(heatmap.data)
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
# module.gene.mapping %>% 
#  dplyr::filter(`bwnet$colors` == 'turquoise') %>% 
#  rownames() %>% head() %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
```

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
MEbrown
</td>
<td style="text-align:right;">
-0.2548163
</td>
<td style="text-align:right;">
0.3835913
</td>
<td style="text-align:right;">
-0.7964666
</td>
<td style="text-align:right;">
0.9328351
</td>
<td style="text-align:right;">
0.5130433
</td>
<td style="text-align:right;">
0.6542792
</td>
<td style="text-align:right;">
-0.2781453
</td>
<td style="text-align:right;">
-0.8718433
</td>
<td style="text-align:right;">
0.8947011
</td>
<td style="text-align:right;">
-0.1471881
</td>
</tr>
<tr>
<td style="text-align:left;">
MEgreen
</td>
<td style="text-align:right;">
-0.4786721
</td>
<td style="text-align:right;">
0.8533198
</td>
<td style="text-align:right;">
-0.4799353
</td>
<td style="text-align:right;">
0.5795083
</td>
<td style="text-align:right;">
0.7074275
</td>
<td style="text-align:right;">
0.7822774
</td>
<td style="text-align:right;">
-0.5503786
</td>
<td style="text-align:right;">
-0.7567056
</td>
<td style="text-align:right;">
0.4096604
</td>
<td style="text-align:right;">
-0.3875060
</td>
</tr>
<tr>
<td style="text-align:left;">
MEsalmon
</td>
<td style="text-align:right;">
0.6786161
</td>
<td style="text-align:right;">
-0.1598559
</td>
<td style="text-align:right;">
0.3695897
</td>
<td style="text-align:right;">
-0.1072858
</td>
<td style="text-align:right;">
-0.2279364
</td>
<td style="text-align:right;">
-0.5420069
</td>
<td style="text-align:right;">
0.0499154
</td>
<td style="text-align:right;">
0.2344174
</td>
<td style="text-align:right;">
-0.3685222
</td>
<td style="text-align:right;">
0.0636799
</td>
</tr>
<tr>
<td style="text-align:left;">
MEpink
</td>
<td style="text-align:right;">
0.6536159
</td>
<td style="text-align:right;">
-0.7876569
</td>
<td style="text-align:right;">
0.2664018
</td>
<td style="text-align:right;">
-0.2162103
</td>
<td style="text-align:right;">
-0.6413139
</td>
<td style="text-align:right;">
-0.7623657
</td>
<td style="text-align:right;">
0.5163500
</td>
<td style="text-align:right;">
0.4716704
</td>
<td style="text-align:right;">
-0.1160821
</td>
<td style="text-align:right;">
0.3402889
</td>
</tr>
<tr>
<td style="text-align:left;">
MEturquoise
</td>
<td style="text-align:right;">
0.2674448
</td>
<td style="text-align:right;">
-0.2728222
</td>
<td style="text-align:right;">
-0.5563340
</td>
<td style="text-align:right;">
0.7093303
</td>
<td style="text-align:right;">
-0.0208291
</td>
<td style="text-align:right;">
0.0162399
</td>
<td style="text-align:right;">
0.1436616
</td>
<td style="text-align:right;">
-0.4520571
</td>
<td style="text-align:right;">
0.7673475
</td>
<td style="text-align:right;">
0.1439595
</td>
</tr>
<tr>
<td style="text-align:left;">
MEblack
</td>
<td style="text-align:right;">
-0.3343860
</td>
<td style="text-align:right;">
0.7012574
</td>
<td style="text-align:right;">
0.2821298
</td>
<td style="text-align:right;">
-0.2218551
</td>
<td style="text-align:right;">
0.2983879
</td>
<td style="text-align:right;">
0.3783966
</td>
<td style="text-align:right;">
-0.3347173
</td>
<td style="text-align:right;">
-0.1310596
</td>
<td style="text-align:right;">
-0.4155377
</td>
<td style="text-align:right;">
-0.4931315
</td>
</tr>
<tr>
<td style="text-align:left;">
MEred
</td>
<td style="text-align:right;">
-0.5373699
</td>
<td style="text-align:right;">
0.3541931
</td>
<td style="text-align:right;">
0.1908512
</td>
<td style="text-align:right;">
-0.4405965
</td>
<td style="text-align:right;">
0.2413935
</td>
<td style="text-align:right;">
0.3167674
</td>
<td style="text-align:right;">
-0.2990997
</td>
<td style="text-align:right;">
0.2349852
</td>
<td style="text-align:right;">
-0.4046701
</td>
<td style="text-align:right;">
-0.0590650
</td>
</tr>
<tr>
<td style="text-align:left;">
MEpurple
</td>
<td style="text-align:right;">
0.1988043
</td>
<td style="text-align:right;">
-0.0976320
</td>
<td style="text-align:right;">
0.6615223
</td>
<td style="text-align:right;">
-0.4716665
</td>
<td style="text-align:right;">
-0.6735426
</td>
<td style="text-align:right;">
-0.1922573
</td>
<td style="text-align:right;">
0.5835958
</td>
<td style="text-align:right;">
0.0843186
</td>
<td style="text-align:right;">
-0.5384024
</td>
<td style="text-align:right;">
-0.6913231
</td>
</tr>
<tr>
<td style="text-align:left;">
MEblue
</td>
<td style="text-align:right;">
0.0208187
</td>
<td style="text-align:right;">
-0.0514062
</td>
<td style="text-align:right;">
0.7928315
</td>
<td style="text-align:right;">
-0.8853670
</td>
<td style="text-align:right;">
-0.3571309
</td>
<td style="text-align:right;">
-0.3452664
</td>
<td style="text-align:right;">
0.1699360
</td>
<td style="text-align:right;">
0.6555042
</td>
<td style="text-align:right;">
-0.9109243
</td>
<td style="text-align:right;">
-0.1383717
</td>
</tr>
<tr>
<td style="text-align:left;">
MEyellow
</td>
<td style="text-align:right;">
0.3415250
</td>
<td style="text-align:right;">
-0.6022128
</td>
<td style="text-align:right;">
0.6801317
</td>
<td style="text-align:right;">
-0.8508983
</td>
<td style="text-align:right;">
-0.5593854
</td>
<td style="text-align:right;">
-0.7636993
</td>
<td style="text-align:right;">
0.3353450
</td>
<td style="text-align:right;">
0.9423342
</td>
<td style="text-align:right;">
-0.7505036
</td>
<td style="text-align:right;">
0.3591655
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
MEbrown
</td>
<td style="text-align:right;">
0.0165794
</td>
<td style="text-align:right;">
0.0002250
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000003
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0086924
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.1711647
</td>
</tr>
<tr>
<td style="text-align:left;">
MEgreen
</td>
<td style="text-align:right;">
0.0000024
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000022
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
0.0000739
</td>
<td style="text-align:right;">
0.0001915
</td>
</tr>
<tr>
<td style="text-align:left;">
MEsalmon
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.1368235
</td>
<td style="text-align:right;">
0.0003944
</td>
<td style="text-align:right;">
0.3197777
</td>
<td style="text-align:right;">
0.0326907
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.6441952
</td>
<td style="text-align:right;">
0.0279266
</td>
<td style="text-align:right;">
0.0004112
</td>
<td style="text-align:right;">
0.5555741
</td>
</tr>
<tr>
<td style="text-align:left;">
MEpink
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0121124
</td>
<td style="text-align:right;">
0.0430494
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000003
</td>
<td style="text-align:right;">
0.0000035
</td>
<td style="text-align:right;">
0.2814702
</td>
<td style="text-align:right;">
0.0011789
</td>
</tr>
<tr>
<td style="text-align:left;">
MEturquoise
</td>
<td style="text-align:right;">
0.0117672
</td>
<td style="text-align:right;">
0.0101204
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.8472558
</td>
<td style="text-align:right;">
0.8806263
</td>
<td style="text-align:right;">
0.1817669
</td>
<td style="text-align:right;">
0.0000098
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.1808534
</td>
</tr>
<tr>
<td style="text-align:left;">
MEblack
</td>
<td style="text-align:right;">
0.0014518
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0077424
</td>
<td style="text-align:right;">
0.0377654
</td>
<td style="text-align:right;">
0.0047468
</td>
<td style="text-align:right;">
0.0002779
</td>
<td style="text-align:right;">
0.0014351
</td>
<td style="text-align:right;">
0.2235558
</td>
<td style="text-align:right;">
0.0000567
</td>
<td style="text-align:right;">
0.0000011
</td>
</tr>
<tr>
<td style="text-align:left;">
MEred
</td>
<td style="text-align:right;">
0.0000001
</td>
<td style="text-align:right;">
0.0007103
</td>
<td style="text-align:right;">
0.0748872
</td>
<td style="text-align:right;">
0.0000174
</td>
<td style="text-align:right;">
0.0234691
</td>
<td style="text-align:right;">
0.0026393
</td>
<td style="text-align:right;">
0.0046432
</td>
<td style="text-align:right;">
0.0275387
</td>
<td style="text-align:right;">
0.0000921
</td>
<td style="text-align:right;">
0.5846303
</td>
</tr>
<tr>
<td style="text-align:left;">
MEpurple
</td>
<td style="text-align:right;">
0.0633304
</td>
<td style="text-align:right;">
0.3654981
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000035
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0727291
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.4347667
</td>
<td style="text-align:right;">
0.0000001
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
0.8473312
</td>
<td style="text-align:right;">
0.6343220
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0006363
</td>
<td style="text-align:right;">
0.0009859
</td>
<td style="text-align:right;">
0.1134480
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.1985560
</td>
</tr>
<tr>
<td style="text-align:left;">
MEyellow
</td>
<td style="text-align:right;">
0.0011280
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
0.0014039
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
0.0005892
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
0.1802528
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000000419
</td>
<td style="text-align:right;">
-0.3779426
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000000457
</td>
<td style="text-align:right;">
0.5622608
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000000460
</td>
<td style="text-align:right;">
-0.7735153
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000001084
</td>
<td style="text-align:right;">
-0.3086314
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000001167
</td>
<td style="text-align:right;">
-0.5921564
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
63.19177
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000185633
</td>
<td style="text-align:right;">
63.01567
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000186918
</td>
<td style="text-align:right;">
58.83606
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000148926
</td>
<td style="text-align:right;">
56.50996
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000134107
</td>
<td style="text-align:right;">
55.12680
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000196968
</td>
<td style="text-align:right;">
53.48278
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000182379
</td>
<td style="text-align:right;">
53.36257
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000122884
</td>
<td style="text-align:right;">
53.09573
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000107819
</td>
<td style="text-align:right;">
51.79683
</td>
</tr>
<tr>
<td style="text-align:left;">
ENSG00000112715
</td>
<td style="text-align:right;">
51.27317
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

    ## [1] 16268

``` r
ncol(norm.counts)
```

    ## [1] 16268

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

    ## # A tibble: 1,320 × 3
    ##    treatment name         value
    ##    <fct>     <chr>        <dbl>
    ##  1 Kelly_Nx  brown      0.0704 
    ##  2 Kelly_Nx  green     -0.00321
    ##  3 Kelly_Nx  salmon    -0.0397 
    ##  4 Kelly_Nx  pink       0.0747 
    ##  5 Kelly_Nx  turquoise  0.131  
    ##  6 Kelly_Nx  black     -0.108  
    ##  7 Kelly_Nx  red       -0.124  
    ##  8 Kelly_Nx  purple    -0.0253 
    ##  9 Kelly_Nx  blue      -0.113  
    ## 10 Kelly_Nx  yellow    -0.0507 
    ## # ℹ 1,310 more rows

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
MEbrown
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEgreen
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEsalmon
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEpink
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEturquoise
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEblack
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEred
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEpurple
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEblue
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEyellow
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEtan
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEcyan
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEgreenyellow
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
MEmagenta
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
0.0704479
</td>
<td style="text-align:right;">
-0.0032104
</td>
<td style="text-align:right;">
-0.0396951
</td>
<td style="text-align:right;">
0.0746683
</td>
<td style="text-align:right;">
0.1311820
</td>
<td style="text-align:right;">
-0.1078184
</td>
<td style="text-align:right;">
-0.1242764
</td>
<td style="text-align:right;">
-0.0252541
</td>
<td style="text-align:right;">
-0.1131964
</td>
<td style="text-align:right;">
-0.0506865
</td>
<td style="text-align:right;">
-0.0298488
</td>
<td style="text-align:right;">
0.0178199
</td>
<td style="text-align:right;">
-0.0559553
</td>
<td style="text-align:right;">
-0.0841484
</td>
<td style="text-align:right;">
0.0154085
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S38
</td>
<td style="text-align:right;">
0.0935854
</td>
<td style="text-align:right;">
0.0803486
</td>
<td style="text-align:right;">
-0.0858448
</td>
<td style="text-align:right;">
-0.0583535
</td>
<td style="text-align:right;">
0.0320930
</td>
<td style="text-align:right;">
0.0226513
</td>
<td style="text-align:right;">
-0.0158542
</td>
<td style="text-align:right;">
-0.0420966
</td>
<td style="text-align:right;">
-0.0627371
</td>
<td style="text-align:right;">
-0.0923639
</td>
<td style="text-align:right;">
0.0305658
</td>
<td style="text-align:right;">
-0.0187309
</td>
<td style="text-align:right;">
-0.0041504
</td>
<td style="text-align:right;">
-0.0641040
</td>
<td style="text-align:right;">
0.2016300
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S39
</td>
<td style="text-align:right;">
0.1098184
</td>
<td style="text-align:right;">
-0.0146235
</td>
<td style="text-align:right;">
-0.0309220
</td>
<td style="text-align:right;">
0.0839746
</td>
<td style="text-align:right;">
0.1662897
</td>
<td style="text-align:right;">
-0.1046556
</td>
<td style="text-align:right;">
-0.1602565
</td>
<td style="text-align:right;">
-0.0279465
</td>
<td style="text-align:right;">
-0.1325435
</td>
<td style="text-align:right;">
-0.0688365
</td>
<td style="text-align:right;">
0.0492993
</td>
<td style="text-align:right;">
0.0078569
</td>
<td style="text-align:right;">
-0.0461766
</td>
<td style="text-align:right;">
-0.1080058
</td>
<td style="text-align:right;">
0.0394809
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S40
</td>
<td style="text-align:right;">
0.0986984
</td>
<td style="text-align:right;">
0.0218010
</td>
<td style="text-align:right;">
-0.0406139
</td>
<td style="text-align:right;">
0.0343487
</td>
<td style="text-align:right;">
0.1123585
</td>
<td style="text-align:right;">
-0.0669471
</td>
<td style="text-align:right;">
-0.1051205
</td>
<td style="text-align:right;">
-0.0468795
</td>
<td style="text-align:right;">
-0.1065608
</td>
<td style="text-align:right;">
-0.0739194
</td>
<td style="text-align:right;">
0.0296701
</td>
<td style="text-align:right;">
-0.0039843
</td>
<td style="text-align:right;">
-0.0187718
</td>
<td style="text-align:right;">
-0.0900087
</td>
<td style="text-align:right;">
0.0477808
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S41
</td>
<td style="text-align:right;">
0.0207382
</td>
<td style="text-align:right;">
-0.1990075
</td>
<td style="text-align:right;">
-0.0392581
</td>
<td style="text-align:right;">
0.1469398
</td>
<td style="text-align:right;">
0.1368386
</td>
<td style="text-align:right;">
-0.3083045
</td>
<td style="text-align:right;">
-0.1075287
</td>
<td style="text-align:right;">
-0.1255792
</td>
<td style="text-align:right;">
-0.1268923
</td>
<td style="text-align:right;">
0.0558872
</td>
<td style="text-align:right;">
0.1469921
</td>
<td style="text-align:right;">
0.1629028
</td>
<td style="text-align:right;">
-0.0152721
</td>
<td style="text-align:right;">
0.0089755
</td>
<td style="text-align:right;">
0.0057531
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA_P2041_S42
</td>
<td style="text-align:right;">
0.1054921
</td>
<td style="text-align:right;">
0.0779958
</td>
<td style="text-align:right;">
-0.0708915
</td>
<td style="text-align:right;">
-0.0189553
</td>
<td style="text-align:right;">
0.0859868
</td>
<td style="text-align:right;">
0.0054402
</td>
<td style="text-align:right;">
-0.0708119
</td>
<td style="text-align:right;">
-0.0069345
</td>
<td style="text-align:right;">
-0.0897694
</td>
<td style="text-align:right;">
-0.1027825
</td>
<td style="text-align:right;">
0.0015171
</td>
<td style="text-align:right;">
-0.0393284
</td>
<td style="text-align:right;">
-0.0298222
</td>
<td style="text-align:right;">
-0.1022959
</td>
<td style="text-align:right;">
-0.0083317
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
