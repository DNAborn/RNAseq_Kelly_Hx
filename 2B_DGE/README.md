DGE
================
Kelterborn
2024-03-20

- [0. Load](#0-load)
  - [- Load R librarys](#--load-r-librarys)
  - [- Load dds](#--load-dds)
  - [- functions](#--functions)
- [1. Make results](#1-make-results)
  - [-Plot example counts](#-plot-example-counts)
- [2. Data Dive](#2-data-dive)
  - [Colour sheme](#colour-sheme)
  - [Volcanos](#volcanos)
  - [Overlaps (Venn)](#overlaps-venn)
  - [Heatmaps](#heatmaps)

# 0. Load

## - Load R librarys

## - Load dds

## - functions

# 1. Make results

#### (Advanced results troubleshooting)

<figure>
<img src="Contrasts.png" alt="Contrasts_overview" />
<figcaption aria-hidden="true">Contrasts_overview</figcaption>
</figure>

### -Generate toplist

    ## design

    ## ~experiment + genotype + treatment + genotype:treatment

    ## cutoffs
    ## differential expressed: p=0.05,bM=10,l2FC=1
    ## top genes:              p=0.01,bM=100,l2FC=2

|                            | all.DEGs | top.DEGs |
|:---------------------------|---------:|---------:|
| deg_Hif1a.Hx.vs.Nx         |     6166 |      983 |
| deg_Hif2a.Hx.vs.Nx         |     3394 |      542 |
| deg_Hif1b.Hx.vs.Nx         |     2051 |      310 |
| deg_Kelly.Hx.vs.Nx         |     5301 |      866 |
| deg_Nx.Hif1a.vs.Kelly      |      223 |       21 |
| deg_Nx.Hif2a.vs.Kelly      |      596 |       75 |
| deg_Nx.Hif1b.vs.Kelly      |      721 |       58 |
| deg_Hx.Hif1a.vs.Kelly      |     1013 |       77 |
| deg_Hx.Hif2a.vs.Kelly      |     2532 |      348 |
| deg_Hx.Hif1b.vs.Kelly      |     4208 |      531 |
| deg_Hx.Hif2a.vs.Hif1a      |     4174 |      607 |
| deg_Hx.Hif1b.vs.Hif1a      |     4972 |      643 |
| deg_Hx.Hif1b.vs.Hif2a      |     2043 |      249 |
| deg_Hif1aHxNx.vs.KellyHxNx |      686 |       54 |
| deg_Hif2aHxNx.vs.KellyHxNx |     2238 |      281 |
| deg_Hif1bHxNx.vs.KellyHxNx |     3887 |      431 |
| deg_Hx.Hif1b.vs.Hif12a     |     2741 |      224 |
| deg_Hx.Kelly.vs.allHIFs    |     1502 |      148 |
| deg_Hx.vs.Nx               |     3599 |      495 |

## -Plot example counts

![](Readme_files/figure-gfm/results_counts-1.png)<!-- -->![](Readme_files/figure-gfm/results_counts-2.png)<!-- -->![](Readme_files/figure-gfm/results_counts-3.png)<!-- -->![](Readme_files/figure-gfm/results_counts-4.png)<!-- -->![](Readme_files/figure-gfm/results_counts-5.png)<!-- -->![](Readme_files/figure-gfm/results_counts-6.png)<!-- -->![](Readme_files/figure-gfm/results_counts-7.png)<!-- -->![](Readme_files/figure-gfm/results_counts-8.png)<!-- -->![](Readme_files/figure-gfm/results_counts-9.png)<!-- -->![](Readme_files/figure-gfm/results_counts-10.png)<!-- -->![](Readme_files/figure-gfm/results_counts-11.png)<!-- -->![](Readme_files/figure-gfm/results_counts-12.png)<!-- -->![](Readme_files/figure-gfm/results_counts-13.png)<!-- -->![](Readme_files/figure-gfm/results_counts-14.png)<!-- -->![](Readme_files/figure-gfm/results_counts-15.png)<!-- -->![](Readme_files/figure-gfm/results_counts-16.png)<!-- -->![](Readme_files/figure-gfm/results_counts-17.png)<!-- -->![](Readme_files/figure-gfm/results_counts-18.png)<!-- -->![](Readme_files/figure-gfm/results_counts-19.png)<!-- -->

|                 | symbol  | baseMean | log2FoldChange |    lfcSE |       stat |    pvalue |      padj |
|:----------------|:--------|---------:|---------------:|---------:|-----------:|----------:|----------:|
| ENSG00000234964 | FABP5P7 | 163.7803 |       -1.07385 | 1.094697 | -0.9809562 | 0.3266143 | 0.6176694 |

![](Readme_files/figure-gfm/results_counts-20.png)<!-- -->

# 2. Data Dive

## Colour sheme

## Volcanos

### Draw Vulcanos

![](Readme_files/figure-gfm/draw%20vulcano-1.png)![](Readme_files/figure-gfm/draw%20vulcano-2.png)

### (continuous Vulcanos)

``` r
# gradient is fixed to padj = y-axis

# Vulcano
lcol="grey20"
xlim=10
ylim=300
n <- "Kelly.Hx.vs.Nx"
res <- results_list[[n]]
l <- length(res)

res_shrink <- lfcShrink(dds, res=res, type="ashr")
res_shrink$symbol <- res$symbol

# remove nas
res <- res[!is.na(res$padj),]
res <- res[!is.na(res$log2FoldChange),]

# rename genes
rownames(res) <- res$symbol

# change shape of outliers
shape <- ifelse(abs(res$log2FoldChange) > xlim, 18,
                ifelse(res$padj < 10^-ylim,18,16))
summary(is.na(shape))

# shape[is.na(shape)] <- 2
names(shape)[shape == 18] <- 'out of scale'
names(shape)[shape == 16] <- 'in range'

# move outliers to coord. max.
res$log2FoldChange[res$log2FoldChange > xlim] <- xlim
res$log2FoldChange[res$log2FoldChange < -xlim] <- -xlim
res$padj[res$padj < 10^-ylim] <- 10^-ylim
summary(res$padj < 10^-ylim)

 p1 <- EnhancedVolcano(res,
    lab = res$symbol,
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 10^(-50),
    FCcutoff = 2,
    xlim = c(-xlim, xlim),
    pointSize = c(ifelse(res$log2FoldChange>2, 8, 1)),
    labSize = 6.0,
    shape = c(6, 6, 19, 16),
    title = "DESeq2 results",
    subtitle = "Differential expression",
    caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 10e-4"),
    legendPosition = "right",
    legendLabSize = 14,
    colAlpha = 0.9,
    colGradient = c('red3', 'royalblue'),
    drawConnectors = TRUE,
    hline = c(10e-8),
    widthConnectors = 0.5)

  p1

ev_f <- EnhancedVolcano(res,
    x = 'log2FoldChange',
    y = 'padj',
    lab = res$symbol,
    labSize = 1.5,
    drawConnectors = TRUE,
    boxedLabels = TRUE,
    widthConnectors = 0.5,
    colConnectors = lcol,
    max.overlaps = 17,
    colGradient = c('red3', 'royalblue'),
    xlim = c(-xlim, xlim),
    ylim = c(0, ylim),
    ylab = "Padj (-log10)",
    title = n,
    subtitle = paste("DE genes:",l),
    # sub = "SVF",

    FCcutoff = 2,
    # pointSize = c(ifelse(rownames(res_WT_D_vs.WT_BL) %in% rownames(top_WT_BL_vs.pcry_BL), 8, 1)),
    legendLabels=c('Not sig.','|L2F| > 2.5','p-adj < 0.05',
                   'p-adj & L2F'),
    legendPosition = 'bottom',
    legendLabSize = 8,
    legendIconSize = 2.0,
    axisLabSize = 8,
    titleLabSize = 8,
    subtitleLabSize = 8,
    captionLabSize = 8,
    caption = {}
   )

ev_f
```

### (prepare data)

### (simple volcano (full))

#### (check cutoff)

## Overlaps (Venn)

### - Hif1a

    ## [1] "Element=7 (702) --> a3(702)"

<img src="Readme_files/figure-gfm/venn_hif1a-1.png" width="100%" /><img src="Readme_files/figure-gfm/venn_hif1a-2.png" width="100%" /><img src="Readme_files/figure-gfm/venn_hif1a-3.png" width="100%" />

    ## Kelly.Hx.vs.Nx

|                 |  baseMean | log2FoldChange |     lfcSE |      stat | pvalue | padj | symbol  |
|:----------------|----------:|---------------:|----------:|----------:|-------:|-----:|:--------|
| ENSG00000073060 | 11958.596 |       1.912261 | 0.0564790 |  33.85793 |      0 |    0 | SCARB1  |
| ENSG00000132382 |  6908.702 |      -1.840556 | 0.0836373 | -22.00641 |      0 |    0 | MYBBP1A |
| ENSG00000186469 |  8085.827 |       1.871876 | 0.0904267 |  20.70048 |      0 |    0 | GNG2    |
| ENSG00000117016 |  6871.537 |      -1.758324 | 0.0956286 | -18.38701 |      0 |    0 | RIMS3   |
| ENSG00000189241 | 11177.408 |       1.731705 | 0.0524954 |  32.98773 |      0 |    0 | TSPYL1  |
| ENSG00000184185 |  6118.126 |       1.972228 | 0.0767743 |  25.68864 |      0 |    0 | KCNJ12  |
| ENSG00000100285 | 39849.439 |      -1.267082 | 0.0781810 | -16.20705 |      0 |    0 | NEFH    |
| ENSG00000179222 | 22237.096 |       1.642876 | 0.0794722 |  20.67233 |      0 |    0 | MAGED1  |
| ENSG00000168348 | 18593.856 |      -1.528296 | 0.0792187 | -19.29212 |      0 |    0 | INSM2   |

    ## Hif1a.Hx.vs.Nx

|                 |  baseMean | log2FoldChange |     lfcSE |      stat | pvalue | padj | symbol  |
|:----------------|----------:|---------------:|----------:|----------:|-------:|-----:|:--------|
| ENSG00000073060 | 11958.596 |       2.513715 | 0.0674164 |  37.28643 |      0 |    0 | SCARB1  |
| ENSG00000132382 |  6908.702 |      -2.577358 | 0.1000207 | -25.76825 |      0 |    0 | MYBBP1A |
| ENSG00000186469 |  8085.827 |       2.579589 | 0.1081358 |  23.85510 |      0 |    0 | GNG2    |
| ENSG00000117016 |  6871.537 |      -2.641920 | 0.1143741 | -23.09895 |      0 |    0 | RIMS3   |
| ENSG00000189241 | 11177.408 |       2.279385 | 0.0626101 |  36.40604 |      0 |    0 | TSPYL1  |
| ENSG00000184185 |  6118.126 |       2.406483 | 0.0916057 |  26.27001 |      0 |    0 | KCNJ12  |
| ENSG00000100285 | 39849.439 |      -2.283594 | 0.0934795 | -24.42883 |      0 |    0 | NEFH    |
| ENSG00000179222 | 22237.096 |       2.267305 | 0.0949827 |  23.87071 |      0 |    0 | MAGED1  |
| ENSG00000168348 | 18593.856 |      -2.263527 | 0.0946962 | -23.90304 |      0 |    0 | INSM2   |

<img src="Readme_files/figure-gfm/venn_hif1a-4.png" width="100%" />

    ## [1] "Element=12 (31) --> a4(31)"
    ## [1] "Element=14 (702) --> a6(702)"
    ## [1] "Element=16 (15) --> a8(15)"
    ## [1] "Element=17 (703) --> a9(703)"
    ## [1] "Element=19 (266) --> a11(266)"
    ## [1] "Element=20 (133) --> a12(133)"
    ## [1] "Element=22 (1584) --> a14(1584)"
    ## [1] "Element=23 (3466) --> a15(3466)"

<img src="Readme_files/figure-gfm/venn_hif1a-5.png" width="100%" /><img src="Readme_files/figure-gfm/venn_hif1a-6.png" width="100%" />

    ## compare results with contrast vsvs (Hif1a Hx vs. Nx  VS.  Kelly Hx vs. Nx

<img src="Readme_files/figure-gfm/venn_hif1a-7.png" width="100%" />

    ## [1] "Element=7 (4) --> a1(4)"
    ## [1] "Element=8 (25) --> a2(25)"
    ## [1] "Element=9 (607) --> a3(607)"
    ## [1] "Element=10 (17) --> a5(17)"
    ## [1] "Element=11 (37) --> a6(37)"

<img src="Readme_files/figure-gfm/venn_hif1a-8.png" width="100%" /><img src="Readme_files/figure-gfm/venn_hif1a-9.png" width="100%" /><img src="Readme_files/figure-gfm/venn_hif1a-10.png" width="100%" /><img src="Readme_files/figure-gfm/venn_hif1a-11.png" width="100%" /><img src="Readme_files/figure-gfm/venn_hif1a-12.png" width="100%" />

|                 | baseMean | log2FoldChange |     lfcSE |     stat |    pvalue |      padj | symbol |
|:----------------|---------:|---------------:|----------:|---------:|----------:|----------:|:-------|
| ENSG00000105880 | 618.9299 |       0.190761 | 0.1545126 | 1.234598 | 0.2169801 | 0.3009223 | DLX5   |

|                 | baseMean | log2FoldChange |     lfcSE |     stat | pvalue | padj | symbol |
|:----------------|---------:|---------------:|----------:|---------:|-------:|-----:|:-------|
| ENSG00000105880 | 618.9299 |       1.981755 | 0.1830812 | 10.82446 |      0 |    0 | DLX5   |

<img src="Readme_files/figure-gfm/venn_hif1a-13.png" width="100%" />

    ## [1] "Element=7 (1897) --> a1(1897)"
    ## [1] "Element=8 (131) --> a2(131)"
    ## [1] "Element=9 (455) --> a3(455)"
    ## [1] "Element=10 (62) --> a4(62)"
    ## [1] "Element=11 (243) --> a5(243)"
    ## [1] "Element=12 (282) --> a6(282)"
    ## [1] "Element=13 (99) --> a7(99)"

### (- Hif2a)

    ## Hif2a

![](Readme_files/figure-gfm/venn_hif2a-1.png)<!-- -->

    ## [1] "Element=12 (188) --> a4(188)"
    ## [1] "Element=14 (417) --> a6(417)"
    ## [1] "Element=16 (35) --> a8(35)"
    ## [1] "Element=17 (2511) --> a9(2511)"
    ## [1] "Element=19 (90) --> a11(90)"
    ## [1] "Element=20 (261) --> a12(261)"
    ## [1] "Element=22 (757) --> a14(757)"
    ## [1] "Element=23 (1834) --> a15(1834)"

![](Readme_files/figure-gfm/venn_hif2a-2.png)<!-- -->![](Readme_files/figure-gfm/venn_hif2a-3.png)<!-- -->![](Readme_files/figure-gfm/venn_hif2a-4.png)<!-- -->

### (- Hif1b)

    ## Hif1b

![](Readme_files/figure-gfm/venn_hif1b-1.png)<!-- -->

    ## [1] "Element=12 (374) --> a4(374)"
    ## [1] "Element=14 (239) --> a6(239)"
    ## [1] "Element=16 (24) --> a8(24)"
    ## [1] "Element=17 (3446) --> a9(3446)"
    ## [1] "Element=19 (47) --> a11(47)"
    ## [1] "Element=20 (253) --> a12(253)"
    ## [1] "Element=22 (546) --> a14(546)"
    ## [1] "Element=23 (942) --> a15(942)"

![](Readme_files/figure-gfm/venn_hif1b-2.png)<!-- -->![](Readme_files/figure-gfm/venn_hif1b-3.png)<!-- -->![](Readme_files/figure-gfm/venn_hif1b-4.png)<!-- -->

### (- overlap)

    ## overlap of overlaps

![](Readme_files/figure-gfm/venn_overlap_res1-1.png)<!-- -->

    ## [1] "Element=7 (33) --> a1(33)"
    ## [1] "Element=8 (2) --> a2(2)"
    ## [1] "Element=9 (71) --> a3(71)"
    ## [1] "Element=10 (7) --> a4(7)"
    ## [1] "Element=11 (4) --> a5(4)"
    ## [1] "Element=12 (146) --> a6(146)"
    ## [1] "Element=13 (241) --> a7(241)"

    ##   overlap    gene
    ## 1      a1   TIGD5
    ## 2      a2   SPRY4
    ## 3      a3     KDR
    ## 4      a4 MT-RNR1
    ## 5      a5  H2AC19
    ## 6      a6 CYP26B1
    ## 7      a7   PCGF2

![](Readme_files/figure-gfm/venn_overlap_res1-2.png)<!-- -->

### (- other)

### Compare Results 1 2 3

#### Vulcano lists

![](Readme_files/figure-gfm/vulcano_lists-1.png)<!-- -->![](Readme_files/figure-gfm/vulcano_lists-2.png)<!-- -->![](Readme_files/figure-gfm/vulcano_lists-3.png)<!-- -->

|                 |    baseMean | log2FoldChange |     lfcSE |      stat | pvalue | padj | symbol |
|:----------------|------------:|---------------:|----------:|----------:|-------:|-----:|:-------|
| ENSG00000071564 |  9467.75008 |     -0.8376478 | 0.0379823 | -22.05366 |      0 |    0 | TCF3   |
| ENSG00000106003 |   928.09823 |     -3.0186727 | 0.1254119 | -24.07007 |      0 |    0 | LFNG   |
| ENSG00000112769 |   367.38256 |     -2.9390195 | 0.1661520 | -17.68874 |      0 |    0 | LAMA4  |
| ENSG00000135100 |   504.06456 |     10.5449933 | 0.3781407 |  27.88643 |      0 |    0 | HNF1A  |
| ENSG00000152256 | 12032.53494 |      2.3886419 | 0.0897902 |  26.60247 |      0 |    0 | PDK1   |
| ENSG00000237686 |    43.53667 |      4.2463036 | 0.2722903 |  15.59477 |      0 |    0 | SCIRT  |

|                 |    baseMean | log2FoldChange |     lfcSE |       stat | pvalue | padj | symbol |
|:----------------|------------:|---------------:|----------:|-----------:|-------:|-----:|:-------|
| ENSG00000071564 |  9467.75008 |     -1.2049716 | 0.0452621 | -26.622069 |      0 |    0 | TCF3   |
| ENSG00000106003 |   928.09823 |     -4.0890506 | 0.1510384 | -27.072928 |      0 |    0 | LFNG   |
| ENSG00000112769 |   367.38256 |     -1.7922770 | 0.1939710 |  -9.239923 |      0 |    0 | LAMA4  |
| ENSG00000135100 |   504.06456 |     11.0577374 | 0.6902547 |  16.019793 |      0 |    0 | HNF1A  |
| ENSG00000152256 | 12032.53494 |      0.7461557 | 0.1076933 |   6.928522 |      0 |    0 | PDK1   |
| ENSG00000237686 |    43.53667 |      5.5992032 | 0.3501121 |  15.992602 |      0 |    0 | SCIRT  |

|                 |    baseMean | log2FoldChange |     lfcSE |      stat |    pvalue |      padj | symbol |
|:----------------|------------:|---------------:|----------:|----------:|----------:|----------:|:-------|
| ENSG00000071564 |  9467.75008 |      0.0650920 | 0.0443193 |  1.468704 | 0.1419131 | 0.6665731 | TCF3   |
| ENSG00000106003 |   928.09823 |      0.2291789 | 0.1441730 |  1.589611 | 0.1119226 | 0.6088678 | LFNG   |
| ENSG00000112769 |   367.38256 |     -1.0438988 | 0.1851363 | -5.638544 | 0.0000000 | 0.0000028 | LAMA4  |
| ENSG00000135100 |   504.06456 |     -1.8502210 | 0.7711634 | -2.399259 | 0.0164283 | 0.2265094 | HNF1A  |
| ENSG00000152256 | 12032.53494 |     -0.9768009 | 0.1059285 | -9.221322 | 0.0000000 | 0.0000000 | PDK1   |
| ENSG00000237686 |    43.53667 |     -0.6947345 | 0.3908038 | -1.777707 | 0.0754521 | 0.5144705 | SCIRT  |

|                 |    baseMean | log2FoldChange |     lfcSE |       stat |    pvalue |      padj | symbol |
|:----------------|------------:|---------------:|----------:|-----------:|----------:|----------:|:-------|
| ENSG00000071564 |  9467.75008 |     -0.3022319 | 0.0389838 |  -7.752763 | 0.0000000 | 0.0000000 | TCF3   |
| ENSG00000106003 |   928.09823 |     -0.8411990 | 0.1332788 |  -6.311576 | 0.0000000 | 0.0000000 | LFNG   |
| ENSG00000112769 |   367.38256 |      0.1028437 | 0.1769762 |   0.581116 | 0.5611623 | 0.7263930 | LAMA4  |
| ENSG00000135100 |   504.06456 |     -1.3374769 | 0.1528782 |  -8.748644 | 0.0000000 | 0.0000000 | HNF1A  |
| ENSG00000152256 | 12032.53494 |     -2.6192870 | 0.0915034 | -28.625009 | 0.0000000 | 0.0000000 | PDK1   |
| ENSG00000237686 |    43.53667 |      0.6581651 | 0.2038346 |   3.228918 | 0.0012426 | 0.0056239 | SCIRT  |

|                 |    baseMean | log2FoldChange |     lfcSE |        stat |    pvalue |      padj | symbol |
|:----------------|------------:|---------------:|----------:|------------:|----------:|----------:|:-------|
| ENSG00000071564 |  9467.75008 |     -0.3673238 | 0.0591161 |  -6.2135963 | 0.0000000 | 0.0000000 | TCF3   |
| ENSG00000106003 |   928.09823 |     -1.0703779 | 0.1965294 |  -5.4464020 | 0.0000001 | 0.0000022 | LFNG   |
| ENSG00000112769 |   367.38256 |      1.1467425 | 0.2566489 |   4.4681375 | 0.0000079 | 0.0001808 | LAMA4  |
| ENSG00000135100 |   504.06456 |      0.5127441 | 0.7862990 |   0.6520982 | 0.5143379 | 0.7404637 | HNF1A  |
| ENSG00000152256 | 12032.53494 |     -1.6424862 | 0.1401892 | -11.7162070 | 0.0000000 | 0.0000000 | PDK1   |
| ENSG00000237686 |    43.53667 |      1.3528996 | 0.4411382 |   3.0668383 | 0.0021634 | 0.0174763 | SCIRT  |

#### WGCNA RES1,2,3

<img src="Readme_files/figure-gfm/wgcna_res123-1.png" width="33%" /><img src="Readme_files/figure-gfm/wgcna_res123-2.png" width="33%" /><img src="Readme_files/figure-gfm/wgcna_res123-3.png" width="33%" />

### Compare KO

``` r
cat("Results 3 of Hif1a, Hif2a, Hif1b")
```

    ## Results 3 of Hif1a, Hif2a, Hif1b

``` r
names(deg_genes_list)
```

    ##  [1] "deg_Hif1a.Hx.vs.Nx"         "deg_Hif2a.Hx.vs.Nx"        
    ##  [3] "deg_Hif1b.Hx.vs.Nx"         "deg_Kelly.Hx.vs.Nx"        
    ##  [5] "deg_Nx.Hif1a.vs.Kelly"      "deg_Nx.Hif2a.vs.Kelly"     
    ##  [7] "deg_Nx.Hif1b.vs.Kelly"      "deg_Hx.Hif1a.vs.Kelly"     
    ##  [9] "deg_Hx.Hif2a.vs.Kelly"      "deg_Hx.Hif1b.vs.Kelly"     
    ## [11] "deg_Hx.Hif2a.vs.Hif1a"      "deg_Hx.Hif1b.vs.Hif1a"     
    ## [13] "deg_Hx.Hif1b.vs.Hif2a"      "deg_Hif1aHxNx.vs.KellyHxNx"
    ## [15] "deg_Hif2aHxNx.vs.KellyHxNx" "deg_Hif1bHxNx.vs.KellyHxNx"
    ## [17] "deg_Hx.Hif1b.vs.Hif12a"     "deg_Hx.Kelly.vs.allHIFs"   
    ## [19] "deg_Hx.vs.Nx"

``` r
input_list <- deg_genes_list[c("deg_Hif1aHxNx.vs.KellyHxNx","deg_Hif2aHxNx.vs.KellyHxNx","deg_Hif1bHxNx.vs.KellyHxNx")]
plt <- venn.diagram(
    x = input_list,
    category.names = paste(names(input_list),"\n(",input_list %>% summary() %>% .[c(1:length(input_list))],")",sep=""),
    force.unique = TRUE, na = "remove",
    filename = NULL,
    main = "Hif1a Hx-Nx vs. Kelly Hx-Nx", main.fontface = "bold",
    lwd = 2,
    lty = 'blank',
    fill = colors[c(3,5,7)],
    cat.col=c(colors[c(4,6,8)]),
    cat.fontface = "bold")

grid.newpage()
grid.draw(plt)
```

![](Readme_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# plot example counts
overlaps <- calculate.overlap(input_list)
overlaps <- overlaps[order(names(overlaps))]
names(overlaps)
```

    ## [1] "a1" "a2" "a3" "a4" "a5" "a6" "a7"

``` r
getVennElements(plt)
```

    ## [1] "Element=7 (227) --> a1(227)"
    ## [1] "Element=8 (110) --> a2(110)"
    ## [1] "Element=9 (568) --> a3(568)"
    ## [1] "Element=10 (166) --> a4(166)"
    ## [1] "Element=11 (183) --> a5(183)"
    ## [1] "Element=12 (1377) --> a6(1377)"
    ## [1] "Element=13 (2161) --> a7(2161)"

``` r
# plot each 1 gene
goi <- sapply(overlaps,"[[",1) %>% .[order(names(.))]
 
data.frame(overlap = names(goi),
           gene = mcols(dds)[goi,"symbol"])
```

    ##   overlap    gene
    ## 1      a1   BNIP3
    ## 2      a2 FAM162A
    ## 3      a3  CLSTN2
    ## 4      a4  PFKFB3
    ## 5      a5 BHLHE41
    ## 6      a6 ARHGEF6
    ## 7      a7    TNXB

``` r
plotCounts_SK(goi)
```

![](Readme_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

#### WGCNA KOs

<img src="Readme_files/figure-gfm/wgcna_ko-1.png" width="33%" /><img src="Readme_files/figure-gfm/wgcna_ko-2.png" width="33%" /><img src="Readme_files/figure-gfm/wgcna_ko-3.png" width="33%" />

#### WGCNA overlapped KOs

    ## [1] "a1" "a2" "a3" "a4" "a5" "a6" "a7"

<img src="Readme_files/figure-gfm/wgcna_ko_overlaps-1.png" width="33%" /><img src="Readme_files/figure-gfm/wgcna_ko_overlaps-2.png" width="33%" /><img src="Readme_files/figure-gfm/wgcna_ko_overlaps-3.png" width="33%" /><img src="Readme_files/figure-gfm/wgcna_ko_overlaps-4.png" width="33%" /><img src="Readme_files/figure-gfm/wgcna_ko_overlaps-5.png" width="33%" /><img src="Readme_files/figure-gfm/wgcna_ko_overlaps-6.png" width="33%" /><img src="Readme_files/figure-gfm/wgcna_ko_overlaps-7.png" width="33%" />

#### -Remove log files

## Heatmaps

### -Compare results

``` r
# https://slowkow.com/notes/pheatmap-tutorial/
# Complex heatmap https://github.com/jokergoo/ComplexHeatmap/

# combined results
pick_genes <- str_detect(names(results_list),pattern="Hif1aHxNx")
pick_results <- c(4,1,5,8,14)
names(results_list)[pick_results]
pick_genes <- topgenes_list[pick_genes] %>% unlist() %>% unname() %>% unique()
pick_genes <- lapply(topgenes_list,'[',1:10) %>% unlist() %>% unname() %>% unique()
pick_genes <- c(lapply(res_hif1a,'[',1:5),
                lapply(res_1_ab,'[',1:5),
                lapply(res_2_ab,'[',1:5)) %>% unlist() %>% unname() %>% unique()
pick_genes <- res_hif1a[[1]][1:30]
pick_genes <- lapply(res3_list,'[',1:10) %>% unlist() %>% unname() %>% unique()
pick_genes <- res3_list[[1]][1:30]

res_comb <- res.Kelly.Hx.vs.Nx[pick_genes,c(7,1)] %>% data.frame(.)
res_comb <- cbind(res_comb,lapply(results_list[pick_results],function(i) i[pick_genes,2]) %>% do.call(cbind,.) %>% data.frame(.))
res_comb_matrix <- as.matrix(res_comb[,c(-1,-2)])
res_comb_matrix[res_comb_matrix<1 & res_comb_matrix>-1] <- 0
rownames(res_comb_matrix) <- res_comb$symbol

# adapt colors to uniform breaks
mat_breaks <- quantile_breaks(res_comb_matrix, n = 20)
vir_cols <- viridis(length(mat_breaks))
vir_cols[9] <- "white"
hm_cols <- colorRamp2(mat_breaks,vir_cols)

hm <- Heatmap(res_comb_matrix,
        col = hm_cols,
        column_title = "Compare results",
        na_col = "black",
        row_names_gp = gpar(fontsize = 10)
        ) 
hm

plotCounts_SK(goi=pick_genes[1:3])

patchwork::wrap_elements((c_graphic))
```

### -Top genes

Complex Heatmap:
<https://jokergoo.github.io/ComplexHeatmap-reference/book/>

``` r
ht_opt$fast_hclust = TRUE

# Choose genes
pick_genes <- lapply(res3_list,'[',1:10) %>% unlist() %>% unname() %>% unique()

# Get counts, with summarized replicates
dds_heat <- collapseReplicates(dds, dds$condition,dds$names)
vsd <- vst(dds_heat, blind=TRUE) #Variance stabilized transformation
ntd <- normTransform(dds_heat)
rld <- rlog(dds_heat)
mat <- assay(vsd)
# mat <- assay(ntd)
# mat <- assay(rld)

# reduce to picked genes and convert to matrix
mat <- mat[pick_genes,c(1,3,5,7,2,4,6,8)] %>% as.matrix()

# Get WGCNA colors
WGCNA <- mcols(dds)[rownames(mat),"colors"] %>% as.character()

# WGCNA[is.na(WGCNA)] <- 'grey'

rownames(mat) <- mcols(dds)[pick_genes,"symbol"]
names(WGCNA) <- rownames(mat)

# adapt colors to uniform breaks
mat_breaks <- quantile_breaks(mat, n = 20)
vir_cols <- viridis(length(mat_breaks))
# vir_cols[9] <- "white"
hm_cols <- colorRamp2(mat_breaks,vir_cols)


hm <- Heatmap(mat,
  ## heatmap colors
      #   col = hm_cols,
      na_col = "black",
      
  ## columns
      column_title = "TOP genes",
      cluster_columns = FALSE,
      # column_km = 4,
      column_split = c(rep("Hx", 4),rep("Nx", 4)),
      # column_names_gp = gpar(col = c("lightcoral","skyblue1"), fontsize = c(10)),
      top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("lightcoral","skyblue1")),
        labels = c("Nx", "Hx"), 
        labels_gp = gpar(col = "white", fontsize = 10))),
  ## rows
      row_names_gp = gpar(fontsize = 10),
      cluster_rows = TRUE,
      clustering_method_rows = "centroid",
      # clustering_distance_rows = "kendall",
      row_km = 6,
  right_annotation = rowAnnotation(WGCNA = WGCNA, col=list(WGCNA=WGCNA))
       )
hm


list(bar = c("a" = "red", "b" = "green", "c" = "blue"))

col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))
ha = HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))

anno_col <- as.data.frame(colData(dds_heat)[,c("treatment","genotype")])
anno_colors <- list(treatment = c("lightcoral","skyblue1"),
                    genotype = c("grey","seagreen3","turquoise3","tan2"))

names(anno_colors$treatment) <- levels(anno_col$treatment)
names(anno_colors$genotype) <- levels(anno_col$genotype)



hm <- Heatmap(mat,
        col = hm_cols,
        column_title = "Compare results",
        na_col = "black",
        row_names_gp = gpar(fontsize = 10)
        ) 
hm

pheatmap(mat,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         annotation_col=anno_col,
         annotation_colors = anno_colors,
         show_colnames     = FALSE,
         col=viridis(20),
         cutree_rows = 8,
         cutree_cols = 8,
         fontsize_row = 5)
```
