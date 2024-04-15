DGE
================
Kelterborn
2024-03-20

- [0. Load](#0-load)
  - [- Load R librarys](#--load-r-librarys)
  - [- dds](#--dds)
  - [- functions](#--functions)
- [1. Make results](#1-make-results)
  - [Advanced results
    troubleshooting](#advanced-results-troubleshooting)
- [2. Data Dive](#2-data-dive)
  - [Results](#results)
  - [Volcanos](#volcanos)
    - [prepare data](#prepare-data)
    - [simple volcano (full)](#simple-volcano-full)
    - [Draw Vulcanos](#draw-vulcanos)

# 0. Load

## - Load R librarys

## - dds

``` r
# load(file=paste(data,"deseq2.dds", sep="/"))

load(file=paste(data,"deseq2_treatment.dds", sep="/"))
dds_t <- dds
load(file=paste(data,"deseq2_condition.dds", sep="/"))
dds_c <- dds

load(file=paste(data,"deseq2_wgcna.dds", sep="/"))
```

## - functions

# 1. Make results

``` r
design(dds)
```

    ## ~genotype + treatment + genotype:treatment
    ## <environment: 0x562b8ffeb498>

``` r
names(results_list)
```

    ##  [1] "Hif1a.Hx.vs.Nx"      "Hif2a.Hx.vs.Nx"      "Hif1b.Hx.vs.Nx"     
    ##  [4] "Kelly.Hx.vs.Nx"      "Nx.Hif1a.vs.Kelly"   "Nx.Hif2a.vs.Kelly"  
    ##  [7] "Nx.Hif1b.vs.Kelly"   "Hx.Hif1a.vs.Kelly"   "Hx.Hif2a.vs.Kelly"  
    ## [10] "Hx.Hif1b.vs.Kelly"   "Hx.Hif2a.vs.Hif1a"   "Hx.Hif1b.vs.Hif1a"  
    ## [13] "Hx.Hif1b.vs.Hif2a"   "Hx.Hif1b.vs.Hif12a"  "Hx.Kelly.vs.Hif12ab"

### Advanced results troubleshooting

# 2. Data Dive

## Results

## Volcanos

### prepare data

### simple volcano (full)

#### check cutoff

### Draw Vulcanos

``` r
# Input

# results name
n <- "Hif1a.Hx.vs.Nx"

# colours
topcol <- "royalblue4"
hscol <- "royalblue1"
lcol <- "grey20"

# limits
xlim <- 10
ylim <- 300

# number of top genes
ntop <- 100

###################

# names(results_list)

ev_kelly <- eVukcano_SK(n="Kelly.Hx.vs.Nx", # results name
            ntop=200, # number of top genes
            topcol="royalblue4", # color top genes
            hscol="royalblue1", # color highly significant genes
            lcol="grey20",
            xlim=10,
            ylim=300)

ev_hif1a <- eVukcano_SK(n <- "Hif1a.Hx.vs.Nx",
            ntop=200,
            topcol="orangered3",
            hscol="salmon1",
            lcol="grey20")

ev_hif2a <- eVukcano_SK(n <- "Hif2a.Hx.vs.Nx",
            ntop=200,
            topcol="hotpink4",
            hscol="hotpink1",
            lcol="grey20")

ev_hif1b <- eVukcano_SK(n <- "Hif1b.Hx.vs.Nx",
            ntop=200,
            topcol="darkseagreen4",
            hscol="darkseagreen1",
            lcol="grey20")

( ev_kelly + ev_hif1b ) + plot_layout(guides = "collect", axis_titles="collect")
( ev_hif1a + ev_hif2a) + plot_layout(guides = "collect", axis_titles="collect")
```

![](Readme_files/figure-gfm/draw%20vulcano-1.png)![](Readme_files/figure-gfm/draw%20vulcano-2.png)
