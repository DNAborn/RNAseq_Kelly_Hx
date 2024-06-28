RNA-Seq Kelly Hx all data processing
================
Kelterborn
2024-02-21

- [0. Load](#0-load)
  - [- \#R_libraries](#--r_libraries)
  - [- \#Linux](#--linux)
  - [- \#Download Data](#--download-data)
- [1. Prepare Data](#1-prepare-data)
  - [- \#Make Index](#--make-index)
  - [- \#Mapping (Salmon)](#--mapping-salmon)
  - [- Sample names](#--sample-names)
- [2. Process](#2-process)
  - [- Mapping Rates](#--mapping-rates)
  - [- Tximeta](#--tximeta)
  - [- DESeq2](#--deseq2)
- [Run Deseq2](#run-deseq2)
- [3. Pre-Analysis](#3-pre-analysis)

# 0. Load

## - \#R_libraries

BiocManager::install(“gifski”, force=TRUE)

BiocManager::install(“mgcv”)

## - \#Linux

## - \#Download Data

# 1. Prepare Data

## - \#Make Index

## - \#Mapping (Salmon)

## - Sample names

### Extract filenames from quants

### P3302

### P2041\*

### P557

### combine lists

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:400px; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
experiment
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
RNAs
</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
conditions
</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
date
</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
seq_id
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Seq_runs
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
Control
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:left;">
Kelly_Nx Kelly_Hx HIF1A_Hx HIF2A_Hx
</td>
<td style="text-align:left;">
2018-09-13 2018-09-14
</td>
<td style="text-align:left;">
P557
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
Simon
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:left;">
Kelly_Nx Kelly_Hx HIF1A_Nx HIF1A_Hx HIF1B_Nx HIF1B_Hx
</td>
<td style="text-align:left;">
2017-05-04 2021-06-16 2021-08-25 2021-08-27
</td>
<td style="text-align:left;">
P2041
</td>
<td style="text-align:right;">
22
</td>
</tr>
<tr>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
Ulrike
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:left;">
Kelly_Nx Kelly_Hx HIF1A_Nx HIF1A_Hx HIF2A_Nx HIF2A_Hx
</td>
<td style="text-align:left;">
2023-06-02 2023-06-08 2023-06-15 2023-06-28
</td>
<td style="text-align:left;">
P3302
</td>
<td style="text-align:right;">
150
</td>
</tr>
</tbody>
</table>

</div>

# 2. Process

## - Mapping Rates

### Plot mapping rates

![](Readme_files/figure-gfm/plot_mappingr-1.png)<!-- -->

### clear BFC

## - Tximeta

### add gene symbols

## - DESeq2

### - remove & collapse samples

``` r
colData(dds)$condition %>% table()
```

    ## .
    ## Kelly_Nx Kelly_Hx HIF1A_Nx HIF1A_Hx HIF2A_Nx HIF2A_Hx HIF1B_Nx HIF1B_Hx 
    ##       34       33       25       30       24       28        8        6

``` r
# collapse technical replicates
dds <- collapseReplicates(dds, dds$samplename, dds$names)

colData(dds)$condition %>% table()
```

    ## .
    ## Kelly_Nx Kelly_Hx HIF1A_Nx HIF1A_Hx HIF2A_Nx HIF2A_Hx HIF1B_Nx HIF1B_Hx 
    ##       16       15        9       14        8       12        8        6

``` r
# # collapse repetitions Ulrike
# colData(dds)$Ulrike_rep <- colData(dds)$names
# colData(dds)$Ulrike_rep[str_detect(colData(dds)$Ulrike_rep,pattern="P3302")] <- colData(dds)$exp_rep[str_detect(colData(dds)$exp_rep, pattern="Ulrike")]
# colData(dds)$collapse_rep <- paste(colData(dds)$Ulrike_rep,colData(dds)$condition, sep="_")
# 
# colData(dds)$collapse_rep %>% factor()
# colData(dds)[colData(dds)$experiment=="Simon","condition"] %>% table()
# colData(dds)[colData(dds)$experiment=="Katharina","condition"] %>% table()
# colData(dds)[colData(dds)$experiment=="Ulrike","condition"] %>% table()
# 
# dds <- collapseReplicates(dds, dds$collapse_rep, dds$names)
# colnames(dds) <- dds$names

# Remove uncomplete samples

# Simon HIF1A
colData(dds)[colData(dds)$experiment=="Simon","condition"] %>% table()
```

    ## .
    ## Kelly_Nx Kelly_Hx HIF1A_Nx HIF1A_Hx HIF2A_Nx HIF2A_Hx HIF1B_Nx HIF1B_Hx 
    ##        2        2        1        2        0        0        8        6

``` r
simon_hif1a <- subset(colData(dds), experiment=="Simon" & genotype=="HIF1A")[,c("condition","samplename")]
simon_hif1a
```

    ## DataFrame with 3 rows and 2 columns
    ##               condition    samplename
    ##                <factor>   <character>
    ## RNA_P2041_S44  HIF1A_Nx RNA_P2041_S44
    ## RNA_P2041_S45  HIF1A_Hx RNA_P2041_S45
    ## RNA_P2041_S46  HIF1A_Hx RNA_P2041_S46

``` r
samplestoexclude <- simon_hif1a %>% rownames()

# Katharina HIF1A & HIF2A
colData(dds)[colData(dds)$experiment=="Katharina","condition"] %>% table()
```

    ## .
    ## Kelly_Nx Kelly_Hx HIF1A_Nx HIF1A_Hx HIF2A_Nx HIF2A_Hx HIF1B_Nx HIF1B_Hx 
    ##        4        4        0        4        0        4        0        0

``` r
katharina_hif1a_2a <- subset(colData(dds), experiment=="Katharina" & (genotype=="HIF1A" | genotype=="HIF2A"))[,c("condition","samplename")]
katharina_hif1a_2a
```

    ## DataFrame with 8 rows and 2 columns
    ##              condition   samplename
    ##               <factor>  <character>
    ## RNA_P557_S35  HIF1A_Hx RNA_P557_S35
    ## RNA_P557_S36  HIF2A_Hx RNA_P557_S36
    ## RNA_P557_S39  HIF1A_Hx RNA_P557_S39
    ## RNA_P557_S40  HIF2A_Hx RNA_P557_S40
    ## RNA_P557_S43  HIF1A_Hx RNA_P557_S43
    ## RNA_P557_S44  HIF2A_Hx RNA_P557_S44
    ## RNA_P557_S47  HIF1A_Hx RNA_P557_S47
    ## RNA_P557_S48  HIF2A_Hx RNA_P557_S48

``` r
samplestoexclude <- c(samplestoexclude,katharina_hif1a_2a %>% rownames())
samplestoexclude
```

    ##  [1] "RNA_P2041_S44" "RNA_P2041_S45" "RNA_P2041_S46" "RNA_P557_S35" 
    ##  [5] "RNA_P557_S36"  "RNA_P557_S39"  "RNA_P557_S40"  "RNA_P557_S43" 
    ##  [9] "RNA_P557_S44"  "RNA_P557_S47"  "RNA_P557_S48"

``` r
samplestoexclude <- colData(dds)$samplename %in% samplestoexclude
dds <- dds[,!samplestoexclude]

colData(dds)$condition %>% table()
```

    ## .
    ## Kelly_Nx Kelly_Hx HIF1A_Nx HIF1A_Hx HIF2A_Nx HIF2A_Hx HIF1B_Nx HIF1B_Hx 
    ##       16       15        8        8        8        8        8        6

### (- remove outlier)

DEseq2 with all samples, design = ~experiment + genotype + treatment +
genotype:treatment -\> outlier: 2480, dispOutlier: 681 DEseq2 with all
samples, design = ~genotype + treatment + genotype:treatment -\>
outlier: 27, dispOutlier: 288

# Run Deseq2

``` r
sample.number <- colData(dds)$condition %>% table() %>% min()
keep.sn <- rowSums(counts(dds) > 10) >= sample.number # keep genes with at least x samples with a count of 10 or higher
keep.sn %>% summary()
```

    ##    Mode   FALSE    TRUE 
    ## logical   43597   21583

``` r
dds <- dds[keep.sn,]

dds <- DESeq(dds)

mcols(dds)[!mcols(dds)$betaConv,]
```

    ## DataFrame with 0 rows and 64 columns

``` r
betaconv <- mcols(dds)$betaConv
betaconv %>% is.na() %>% table()
```

    ## .
    ## FALSE 
    ## 21583

``` r
# pheatmap(assays(dds)[["cooks"]])
max(assays(dds)[["cooks"]])
```

    ## [1] 26.35415

``` r
stats <- data.frame("mean" = colMeans(assays(dds)[["cooks"]]),
                    "g.mean" = colMedians(assays(dds)[["cooks"]]),
                    "min" = colMins(assays(dds)[["cooks"]]),
                    "max" = colMaxs(assays(dds)[["cooks"]]))
pheatmap(log(stats))
```

![](Readme_files/figure-gfm/run_deseq2-1.png)<!-- -->

``` r
par(mar=c(10,3,2,2)+.1)
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, color=colData(dds)$experiment)
```

![](Readme_files/figure-gfm/run_deseq2-2.png)<!-- -->

``` r
resultsNames(dds)
```

    ## [1] "Intercept"                 "genotype_HIF1A_vs_Kelly"  
    ## [3] "genotype_HIF2A_vs_Kelly"   "genotype_HIF1B_vs_Kelly"  
    ## [5] "treatment_Hx_vs_Nx"        "genotypeHIF1A.treatmentHx"
    ## [7] "genotypeHIF2A.treatmentHx" "genotypeHIF1B.treatmentHx"

``` r
summary(results(dds, name="treatment_Hx_vs_Nx"))
```

    ## 
    ## out of 21583 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 9740, 45%
    ## LFC < 0 (down)     : 7433, 34%
    ## outliers [1]       : 2, 0.0093%
    ## low counts [2]     : 0, 0%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
design(dds)
```

    ## ~genotype + treatment + genotype:treatment
    ## <environment: 0x55cefc5fd3c0>

``` r
summary(results(dds, name="treatment_Hx_vs_Nx"))
```

    ## 
    ## out of 21583 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 9740, 45%
    ## LFC < 0 (down)     : 7433, 34%
    ## outliers [1]       : 2, 0.0093%
    ## low counts [2]     : 0, 0%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
plotMA(dds)
plotDispEsts(dds)
resultsNames(dds)
```

    ## [1] "Intercept"                 "genotype_HIF1A_vs_Kelly"  
    ## [3] "genotype_HIF2A_vs_Kelly"   "genotype_HIF1B_vs_Kelly"  
    ## [5] "treatment_Hx_vs_Nx"        "genotypeHIF1A.treatmentHx"
    ## [7] "genotypeHIF2A.treatmentHx" "genotypeHIF1B.treatmentHx"

``` r
mcols(dds)$dispOutlier %>% table()
```

    ## .
    ## FALSE  TRUE 
    ## 21025   558

``` r
outlier <- mcols(dds)$dispOutlier %>% mcols(dds)[.,"gene_id"]


plot(x=log(mcols(dds)$baseMean),y=log(mcols(dds)$dispGeneEst))

mcols(dds)$dispGeneEst[log(mcols(dds)$dispGeneEst) < -5] %>% length()
```

    ## [1] 681

<img src="Readme_files/figure-gfm/dds_design-1.png" width="50%" /><img src="Readme_files/figure-gfm/dds_design-2.png" width="50%" /><img src="Readme_files/figure-gfm/dds_design-3.png" width="50%" />

``` r
# Remove Outlier samples: "Control", "P2041_HIF1A_Hx_42"
# dds <- dds[!colData(dds)$experiment=="Control",]
dds <- dds[,!colnames(dds)=="P2041_Kelly_Nx_51"]


dds <- DESeq(dds)

stats <- data.frame("mean" = colMeans(assays(dds)[["cooks"]]),
                    "g.mean" = colMedians(assays(dds)[["cooks"]]),
                    "min" = colMins(assays(dds)[["cooks"]]),
                    "max" = colMaxs(assays(dds)[["cooks"]]))
pheatmap(log(stats))

par(mar=c(10,3,2,2)+.1)
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, color=colData(dds)$experiment)

resultsNames(dds)
summary(results(dds, name="treatment_Hx_vs_Nx"))

getwd()
design(dds)
# save(dds,file=paste(data,"deseq2.dds", sep="/"))
# dds <- 1
# load(file=paste(data,"deseq2.dds", sep="/"))
dds

subset(colData(dds), condition=="HIF1A_Hx" & sequencing =="P2041") %>% kable() %>% kable_styling("striped", full_width = T) %>% scroll_box(height = "400px")
```

# 3. Pre-Analysis

### - Data transformations

#### -#rlog

``` r
load(file=paste(data,"rlog.rld", sep="/"))
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))
```

<img src="Readme_files/figure-gfm/pre_trans_fig, figures-side-1.png" width="33%" /><img src="Readme_files/figure-gfm/pre_trans_fig, figures-side-2.png" width="33%" /><img src="Readme_files/figure-gfm/pre_trans_fig, figures-side-3.png" width="33%" />

### - Check sample distance

<img src="Readme_files/figure-gfm/pre_sample_dist-1.png" width="100%" />

### - Perform principal component analysis

<img src="Readme_files/figure-gfm/pca-1.png" width="80%" />

###### – Advanced PCA

    ## PC5 
    ##   5

<img src="Readme_files/figure-gfm/pca_advanced-1.png" width="80%" /><img src="Readme_files/figure-gfm/pca_advanced-2.png" width="80%" /><img src="Readme_files/figure-gfm/pca_advanced-3.png" width="80%" /><img src="Readme_files/figure-gfm/pca_advanced-4.png" width="80%" />

###### – \#PCA gif

###### – \#PCA gif2

<a href="pca.gif" height="100%," width="100%">PCA Gif</a>

<a href="pca2.gif" height="100%," width="100%">PCA Gif</a>

### - Plot example counts

    ## [1] 21583

    ## [1] 21583

    ## [1] 17277

    ## [1] 17276

<img src="Readme_files/figure-gfm/example_counts-1.png" width="50%" /><img src="Readme_files/figure-gfm/example_counts-2.png" width="50%" />

``` r
sessionInfo()
```

    ## R version 4.4.1 (2024-06-14)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 22.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /opt/intel/oneapi/mkl/2024.0/lib/libmkl_rt.so.2;  LAPACK version 3.10.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=de_DE.UTF-8    
    ##  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8   
    ##  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Berlin
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] ensembldb_2.28.0            AnnotationFilter_1.28.0    
    ##  [3] GenomicFeatures_1.56.0      animation_2.7              
    ##  [5] viridis_0.6.5               viridisLite_0.4.2          
    ##  [7] writexl_1.5.0               knitr_1.47                 
    ##  [9] kableExtra_1.4.0            R.utils_2.12.3             
    ## [11] R.oo_1.26.0                 R.methodsS3_1.8.2          
    ## [13] curl_5.2.1                  data.table_1.15.4          
    ## [15] sessioninfo_1.2.2           VennDiagram_1.7.3          
    ## [17] futile.logger_1.4.3         readxl_1.4.3               
    ## [19] patchwork_1.2.0             gridExtra_2.3              
    ## [21] EnhancedVolcano_1.22.0      cowplot_1.1.3              
    ## [23] ggalt_0.4.0                 PCAtools_2.16.0            
    ## [25] ggrepel_0.9.5               pheatmap_1.0.12            
    ## [27] GOSemSim_2.30.0             biomaRt_2.60.0             
    ## [29] clusterProfiler_4.12.0      vsn_3.72.0                 
    ## [31] AnnotationHub_3.12.0        org.Mm.eg.db_3.19.1        
    ## [33] AnnotationDbi_1.66.0        RColorBrewer_1.1-3         
    ## [35] DESeq2_1.44.0               SummarizedExperiment_1.34.0
    ## [37] Biobase_2.64.0              MatrixGenerics_1.16.0      
    ## [39] matrixStats_1.3.0           GenomicRanges_1.56.1       
    ## [41] GenomeInfoDb_1.40.1         IRanges_2.38.0             
    ## [43] S4Vectors_0.42.0            BiocGenerics_0.50.0        
    ## [45] tximport_1.32.0             tximeta_1.22.1             
    ## [47] stringi_1.8.4               plyr_1.8.9                 
    ## [49] lubridate_1.9.3             forcats_1.0.0              
    ## [51] stringr_1.5.1               dplyr_1.1.4                
    ## [53] purrr_1.0.2                 readr_2.1.5                
    ## [55] tidyr_1.3.1                 tibble_3.2.1               
    ## [57] ggplot2_3.5.1               tidyverse_2.0.0            
    ## [59] BiocFileCache_2.12.0        dbplyr_2.5.0               
    ## [61] devtools_2.4.5              usethis_2.2.3              
    ## [63] BiocManager_1.30.23        
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] fs_1.6.4                  ProtGenerics_1.36.0      
    ##   [3] bitops_1.0-7              enrichplot_1.24.0        
    ##   [5] HDO.db_0.99.1             httr_1.4.7               
    ##   [7] ash_1.0-15                profvis_0.3.8            
    ##   [9] tools_4.4.1               utf8_1.2.4               
    ##  [11] R6_2.5.1                  lazyeval_0.2.2           
    ##  [13] urlchecker_1.0.1          withr_3.0.0              
    ##  [15] prettyunits_1.2.0         preprocessCore_1.66.0    
    ##  [17] cli_3.6.3                 formatR_1.14             
    ##  [19] scatterpie_0.2.3          labeling_0.4.3           
    ##  [21] systemfonts_1.1.0         Rsamtools_2.20.0         
    ##  [23] yulab.utils_0.1.4         gson_0.1.0               
    ##  [25] txdbmaker_1.0.1           svglite_2.1.3            
    ##  [27] DOSE_3.30.1               maps_3.4.2               
    ##  [29] limma_3.60.3              rstudioapi_0.16.0        
    ##  [31] RSQLite_2.3.7             generics_0.1.3           
    ##  [33] gridGraphics_0.5-1        BiocIO_1.14.0            
    ##  [35] vroom_1.6.5               GO.db_3.19.1             
    ##  [37] Matrix_1.7-0              fansi_1.0.6              
    ##  [39] abind_1.4-5               lifecycle_1.0.4          
    ##  [41] yaml_2.3.8                qvalue_2.36.0            
    ##  [43] SparseArray_1.4.8         blob_1.2.4               
    ##  [45] promises_1.3.0            dqrng_0.4.1              
    ##  [47] crayon_1.5.3              miniUI_0.1.1.1           
    ##  [49] lattice_0.22-6            beachmat_2.20.0          
    ##  [51] KEGGREST_1.44.1           pillar_1.9.0             
    ##  [53] fgsea_1.30.0              rjson_0.2.21             
    ##  [55] codetools_0.2-20          fastmatch_1.1-4          
    ##  [57] glue_1.7.0                ggfun_0.1.5              
    ##  [59] remotes_2.5.0             vctrs_0.6.5              
    ##  [61] png_0.1-8                 treeio_1.28.0            
    ##  [63] cellranger_1.1.0          gtable_0.3.5             
    ##  [65] cachem_1.1.0              xfun_0.45                
    ##  [67] S4Arrays_1.4.1            mime_0.12                
    ##  [69] tidygraph_1.3.1           statmod_1.5.0            
    ##  [71] ellipsis_0.3.2            nlme_3.1-165             
    ##  [73] ggtree_3.12.0             bit64_4.0.5              
    ##  [75] progress_1.2.3            filelock_1.0.3           
    ##  [77] affyio_1.74.0             irlba_2.3.5.1            
    ##  [79] KernSmooth_2.23-24        colorspace_2.1-0         
    ##  [81] DBI_1.2.3                 tidyselect_1.2.1         
    ##  [83] bit_4.0.5                 compiler_4.4.1           
    ##  [85] extrafontdb_1.0           httr2_1.0.1              
    ##  [87] xml2_1.3.6                DelayedArray_0.30.1      
    ##  [89] shadowtext_0.1.3          rtracklayer_1.64.0       
    ##  [91] scales_1.3.0              hexbin_1.28.3            
    ##  [93] proj4_1.0-14              affy_1.82.0              
    ##  [95] rappdirs_0.3.3            digest_0.6.36            
    ##  [97] rmarkdown_2.27            XVector_0.44.0           
    ##  [99] htmltools_0.5.8.1         pkgconfig_2.0.3          
    ## [101] extrafont_0.19            sparseMatrixStats_1.16.0 
    ## [103] highr_0.11                fastmap_1.2.0            
    ## [105] rlang_1.1.4               htmlwidgets_1.6.4        
    ## [107] UCSC.utils_1.0.0          shiny_1.8.1.1            
    ## [109] DelayedMatrixStats_1.26.0 farver_2.1.2             
    ## [111] jsonlite_1.8.8            BiocParallel_1.38.0      
    ## [113] BiocSingular_1.20.0       RCurl_1.98-1.14          
    ## [115] magrittr_2.0.3            GenomeInfoDbData_1.2.12  
    ## [117] ggplotify_0.1.2           munsell_0.5.1            
    ## [119] Rcpp_1.0.12               ape_5.8                  
    ## [121] ggraph_2.2.1              zlibbioc_1.50.0          
    ## [123] MASS_7.3-61               pkgbuild_1.4.4           
    ## [125] parallel_4.4.1            Biostrings_2.72.1        
    ## [127] graphlayouts_1.1.1        splines_4.4.1            
    ## [129] hms_1.1.3                 locfit_1.5-9.10          
    ## [131] igraph_2.0.3              reshape2_1.4.4           
    ## [133] ScaledMatrix_1.12.0       pkgload_1.3.4            
    ## [135] futile.options_1.0.1      BiocVersion_3.19.1       
    ## [137] XML_3.99-0.17             evaluate_0.24.0          
    ## [139] lambda.r_1.2.4            tzdb_0.4.0               
    ## [141] tweenr_2.0.3              httpuv_1.6.15            
    ## [143] Rttf2pt1_1.3.12           polyclip_1.10-6          
    ## [145] ggforce_0.4.2             rsvd_1.0.5               
    ## [147] xtable_1.8-4              restfulr_0.0.15          
    ## [149] tidytree_0.4.6            later_1.3.2              
    ## [151] aplot_0.2.3               memoise_2.0.1            
    ## [153] GenomicAlignments_1.40.0  timechange_0.3.0
