DGE
================
Kelterborn
2024-03-20

- [0. Load](#0-load)
  - [- Load R librarys](#--load-r-librarys)
  - [- dds](#--dds)
  - [- functions](#--functions)
- [1. Make results](#1-make-results)
  - [Plot example counts](#plot-example-counts)
- [2. Data Dive](#2-data-dive)
  - [Colour sheme](#colour-sheme)
  - [Overlaps (Venn)](#overlaps-venn)
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
    ## <environment: 0x55c9bb476870>

``` r
names(results_list)
```

    ##  [1] "Hif1a.Hx.vs.Nx"      "Hif2a.Hx.vs.Nx"      "Hif1b.Hx.vs.Nx"     
    ##  [4] "Kelly.Hx.vs.Nx"      "Nx.Hif1a.vs.Kelly"   "Nx.Hif2a.vs.Kelly"  
    ##  [7] "Nx.Hif1b.vs.Kelly"   "Hx.Hif1a.vs.Kelly"   "Hx.Hif2a.vs.Kelly"  
    ## [10] "Hx.Hif1b.vs.Kelly"   "Hx.Hif2a.vs.Hif1a"   "Hx.Hif1b.vs.Hif1a"  
    ## [13] "Hx.Hif1b.vs.Hif2a"   "Hx.Hif1b.vs.Hif12a"  "Hx.Kelly.vs.allHIFs"
    ## [16] "Hx.vs.Nx"

#### (Advanced results troubleshooting)

## Plot example counts

``` r
cols = brewer.pal(n=8,name = 'Paired')

# get 3 sample genes
topgenes_list <- lapply(results_list,topgenes) %>%  lapply(.,rownames) 
sample_plots <- list()
sample_plots_list <- {}
plot_list <- list()
li <- results_list %>% length()
for (i in 1:li){
  n <- results_list[i] %>% names()
  goi <- c(topgenes_list[[i]][1],sample(topgenes_list[[i]],size=2))
  l <- length(goi)
       for (ig in 1:l){
  s <- mcols(dds)[goi[ig],"symbol"]
  if (s ==""){s <- goi[ig]}
    d <- plotCounts(dds, gene=goi[ig], intgroup=c("condition","experiment","genotype","treatment"), main=s,returnData=TRUE)

  gcounts <- ggplot(d, aes(x = genotype, y = count, fill=treatment, color=treatment)) +
    geom_boxplot(color="black") +
    geom_point(shape=21,color="black",aes(fill=treatment),position=position_dodge(width=0.75), alpha=1) +
    scale_fill_manual(values=cols[c(1,5)]) +
    scale_color_manual(values=cols[c(1,5)]) +
    scale_y_continuous(trans = "log2") +
    labs(title = paste(s,"(",goi[ig],")",sep=" "))
  assign(paste("gcounts_",ig,sep=""),gcounts)
  
  plot_list[[paste(n,s,sep="_")]] <- gcounts
       }
  
  plot <- gcounts_1 + gcounts_2 + gcounts_3 +  
    plot_annotation(title = n) +
    plot_layout(guides = "collect", axis_titles="collect") & theme(legend.position = 'bottom',plot.title = element_text(size=10))
  print(plot)
}
```

![](Readme_files/figure-gfm/results_counts-1.png)<!-- -->![](Readme_files/figure-gfm/results_counts-2.png)<!-- -->![](Readme_files/figure-gfm/results_counts-3.png)<!-- -->![](Readme_files/figure-gfm/results_counts-4.png)<!-- -->![](Readme_files/figure-gfm/results_counts-5.png)<!-- -->![](Readme_files/figure-gfm/results_counts-6.png)<!-- -->![](Readme_files/figure-gfm/results_counts-7.png)<!-- -->![](Readme_files/figure-gfm/results_counts-8.png)<!-- -->![](Readme_files/figure-gfm/results_counts-9.png)<!-- -->![](Readme_files/figure-gfm/results_counts-10.png)<!-- -->![](Readme_files/figure-gfm/results_counts-11.png)<!-- -->![](Readme_files/figure-gfm/results_counts-12.png)<!-- -->![](Readme_files/figure-gfm/results_counts-13.png)<!-- -->![](Readme_files/figure-gfm/results_counts-14.png)<!-- -->![](Readme_files/figure-gfm/results_counts-15.png)<!-- -->![](Readme_files/figure-gfm/results_counts-16.png)<!-- -->

``` r
for (n in names(results_list)){
pn <-   str_detect(names(plot_list), pattern=n) %>% names(plot_list)[.]
plot.3 <- patchwork::wrap_plots(plot_list[pn],ncol = 3) + plot_layout(guides = "collect", axis_titles="collect") + plot_annotation(title = str_split(pn[[1]],pattern="_",simplify =TRUE)[1]) & theme(legend.position = 'bottom',plot.title = element_text(size=10))
        print(plot.3)
        }            
```

![](Readme_files/figure-gfm/results_counts-17.png)<!-- -->![](Readme_files/figure-gfm/results_counts-18.png)<!-- -->![](Readme_files/figure-gfm/results_counts-19.png)<!-- -->![](Readme_files/figure-gfm/results_counts-20.png)<!-- -->![](Readme_files/figure-gfm/results_counts-21.png)<!-- -->![](Readme_files/figure-gfm/results_counts-22.png)<!-- -->![](Readme_files/figure-gfm/results_counts-23.png)<!-- -->![](Readme_files/figure-gfm/results_counts-24.png)<!-- -->![](Readme_files/figure-gfm/results_counts-25.png)<!-- -->![](Readme_files/figure-gfm/results_counts-26.png)<!-- -->![](Readme_files/figure-gfm/results_counts-27.png)<!-- -->![](Readme_files/figure-gfm/results_counts-28.png)<!-- -->![](Readme_files/figure-gfm/results_counts-29.png)<!-- -->![](Readme_files/figure-gfm/results_counts-30.png)<!-- -->![](Readme_files/figure-gfm/results_counts-31.png)<!-- -->![](Readme_files/figure-gfm/results_counts-32.png)<!-- -->

# 2. Data Dive

## Colour sheme

``` r
colors_paired <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00")
colors_v <- c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF")
colors <- c("lavenderblush3","lavenderblush4","#90caf9","#1976d2", "#82e0aa", "#239b56", "#f8c471", "#b9770e") 
colors4 <- colors[c(1,3,5,7)]
```

## Overlaps (Venn)

``` r
# Hx vs. Nx
dev.off()
```

    ## null device 
    ##           1

``` r
plt <- venn.diagram(
    x = topgenes_list[c(4,3,1,2)],
    category.names = paste(names(topgenes_list[c(4,3,1,2)]),"\n(",topgenes_list %>% summary() %>% .[c(4,3,1,2)],")",sep=""),
    force.unique = TRUE, na = "remove",
    filename = NULL,
    main = "Hx vs. Nx", main.fontface = "bold",
    lwd = 2,
    lty = 'blank',
    fill = colors4[c(1,4,2,3)])
grid.draw(plt)

# Hx vs. Hx
dev.off()
```

    ## null device 
    ##           1

``` r
plt <- venn.diagram(
    x = topgenes_list[c(8,9,10)],
    category.names = paste(names(topgenes_list[c(8,9,10)]),"\n(",topgenes_list %>% summary() %>% .[c(8,9,10)],")",sep=""),
    force.unique = TRUE, na = "remove",
    filename = NULL,
    main = "Hx vs. Hx", main.fontface = "bold",
    lwd = 2,
    lty = 'blank',
    fill = colors4[c(2,3,4)])
grid.draw(plt)

# Hif1a
dev.off()
```

    ## null device 
    ##           1

``` r
plt <- venn.diagram(
    x = topgenes_list[c(4,1,5,8)],
    category.names = paste(names(topgenes_list[c(4,1,5,8)]),"\n(",topgenes_list %>% summary() %>% .[c(4,1,5,8)],")",sep=""),
    force.unique = TRUE, na = "remove",
    filename = NULL,
    main = "Hif1a", main.fontface = "bold",
    lwd = 2,
    lty = 'blank',
    fill = colors[c(2,4,1,3)])
grid.draw(plt)

# Hif2a
dev.off()
```

    ## null device 
    ##           1

``` r
plt <- venn.diagram(
    x = topgenes_list[c(4,2,6,9)],
    category.names = paste(names(topgenes_list[c(4,2,6,9)]),"\n(",topgenes_list %>% summary() %>% .[c(4,2,6,9)],")",sep=""),
    force.unique = TRUE, na = "remove",
    filename = NULL,
    main = "Hif2a", main.fontface = "bold",
    lwd = 2,
    lty = 'blank',
    fill = colors[c(2,6,1,5)])
grid.draw(plt)

# Hif1b
dev.off()
```

    ## null device 
    ##           1

``` r
plt <- venn.diagram(
    x = topgenes_list[c(4,3,7,10)],
    category.names = paste(names(topgenes_list[c(4,3,7,10)]),"\n(",topgenes_list %>% summary() %>% .[c(4,3,7,10)],")",sep=""),
    force.unique = TRUE, na = "remove",
    filename = NULL,
    main = "Hif2a", main.fontface = "bold",
    lwd = 2,
    lty = 'blank',
    fill = colors[c(2,8,1,7)])
grid.draw(plt)
```

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

( ev_kelly + ev_hif1b ) + plot_layout(guides = "collect", axes="collect", axis_titles="collect") & theme(legend.position = 'bottom', axis.title=element_text(size=8))
( ev_hif1a + ev_hif2a) + plot_layout(guides = "collect", axis_titles="collect") & theme(legend.position = 'bottom')
```

![](Readme_files/figure-gfm/draw%20vulcano-1.png)![](Readme_files/figure-gfm/draw%20vulcano-2.png)
