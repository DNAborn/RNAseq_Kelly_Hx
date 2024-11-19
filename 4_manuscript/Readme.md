Figures
================
Kelterborn
2024-07-18

- [0. Load](#0-load)
  - [- libraries, folders, R_utils](#--libraries-folders-r_utils)
  - [- Load dds](#--load-dds)
  - [- Colour sheme](#--colour-sheme)
  - [-Prepare Results](#-prepare-results)
- [Figure 1: Samples QC](#figure-1-samples-qc)
  - [PCA](#pca)
- [Figure 2: Differential expressed
  genes](#figure-2-differential-expressed-genes)
  - [-Volcano_function](#-volcano_function)
  - [-Volcano_function2](#-volcano_function2)
  - [-Plot Vulcanos](#-plot-vulcanos)
  - [-Venn](#-venn)
- [Figure 3: Gene Cluster](#figure-3-gene-cluster)
  - [Gene Cluster](#gene-cluster)
- [Table 1: Gene List](#table-1-gene-list)
- [Figure 4: Gene Set enrichment](#figure-4-gene-set-enrichment)
  - [GO Analysis](#go-analysis)
  - [Cluster GO terms](#cluster-go-terms)
  - [KEGG](#kegg)
- [Figure 5: Compare with ChIP-Seq](#figure-5-compare-with-chip-seq)
  - [Load datasets](#load-datasets)
  - [ChIP Venns](#chip-venns)
- [](#section)
- [\###Old code](#old-code)
  - [Enhanced volcano](#enhanced-volcano)
  - [Volcanos](#volcanos)
  - [Cluster genes](#cluster-genes)

# 0. Load

## - libraries, folders, R_utils

Load R libraries. If package is missing, install with
‘BiocManager::install(“PackageName”)’

## - Load dds

## - Colour sheme

## -Prepare Results

``` r
deg_genes_list <- lapply(results_list,topgenes_f) %>%  lapply(.,rownames) 
names(deg_genes_list) <- paste("deg",names(deg_genes_list),sep="_")

main_degs <- c(list("Kelly: Hx.vs.Nx" = deg_genes_list[["deg_Kelly.Hx.vs.Nx"]],
                     "Hif1b" = deg_genes_list[["deg_Hif1bHxNx.vs.KellyHxNx"]],
                     "Hif1a" = deg_genes_list[["deg_Hif1aHxNx.vs.KellyHxNx"]],
                     "Hif2a" = deg_genes_list[["deg_Hif2aHxNx.vs.KellyHxNx"]] ))


# Select genes
hif1a_2a_genes <- c(deg_genes_list[["deg_Hif1aHxNx.vs.KellyHxNx"]],
                     deg_genes_list[["deg_Hif2aHxNx.vs.KellyHxNx"]],
                    deg_genes_list[["deg_Hif2aHxNx.vs.Hif1aHxNx"]]) %>%
                  unique()

hif1a_2a_genes %>% length()
```

    ## [1] 5374

``` r
# Filter results
res_names <- names(results_list)
res_final <- results_list[c("Kelly.Hx.vs.Nx","Hif1a.Hx.vs.Nx","Hif2a.Hx.vs.Nx", "Hif1aHxNx.vs.KellyHxNx","Hif2aHxNx.vs.KellyHxNx","Hif1bHxNx.vs.KellyHxNx","Hif2aHxNx.vs.Hif1aHxNx")] 

# create table with all results
res_table <- lapply(res_final,data.frame) %>% lapply(.,"[", , c("log2FoldChange","padj"))
res_table <- do.call('cbind',res_table)
res_table_final <- res_final[[1]][,c("ENSEMBL","ENTREZ","symbol","baseMean")] %>% data.frame()
res_table_final <- cbind(res_table_final,res_table)
res_hif1a_2a <- res_table_final[hif1a_2a_genes,]
colnames(res_hif1a_2a)
```

    ##  [1] "ENSEMBL"                              
    ##  [2] "ENTREZ"                               
    ##  [3] "symbol"                               
    ##  [4] "baseMean"                             
    ##  [5] "Kelly.Hx.vs.Nx.log2FoldChange"        
    ##  [6] "Kelly.Hx.vs.Nx.padj"                  
    ##  [7] "Hif1a.Hx.vs.Nx.log2FoldChange"        
    ##  [8] "Hif1a.Hx.vs.Nx.padj"                  
    ##  [9] "Hif2a.Hx.vs.Nx.log2FoldChange"        
    ## [10] "Hif2a.Hx.vs.Nx.padj"                  
    ## [11] "Hif1aHxNx.vs.KellyHxNx.log2FoldChange"
    ## [12] "Hif1aHxNx.vs.KellyHxNx.padj"          
    ## [13] "Hif2aHxNx.vs.KellyHxNx.log2FoldChange"
    ## [14] "Hif2aHxNx.vs.KellyHxNx.padj"          
    ## [15] "Hif1bHxNx.vs.KellyHxNx.log2FoldChange"
    ## [16] "Hif1bHxNx.vs.KellyHxNx.padj"          
    ## [17] "Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange"
    ## [18] "Hif2aHxNx.vs.Hif1aHxNx.padj"

``` r
# create table with all shrinked results
# Filter results
res_shrink <- res_shrink_list[c("Kelly.Hx.vs.Nx","Hif1a.Hx.vs.Nx","Hif2a.Hx.vs.Nx", "Hif1aHxNx.vs.KellyHxNx","Hif2aHxNx.vs.KellyHxNx","Hif1bHxNx.vs.KellyHxNx","Hif2aHxNx.vs.Hif1aHxNx")] 

# create table with all results
res_shrink <- lapply(res_shrink,data.frame) %>% lapply(.,"[", , c("log2FoldChange","padj"))
res_shrink <- do.call('cbind',res_shrink)
res_shrink_final <- res_final[[1]][,c("symbol","baseMean")] %>% data.frame()
res_shrink_final <- cbind(res_shrink_final,res_shrink)
res_shrink_hif1a_2a <- res_shrink_final[hif1a_2a_genes,]


# Gene universe
# Expression min
results_list[[1]]$baseMean %>% hist(breaks=100000, xlim = c(0,100)) 

# in results
res_hif1a_2a$baseMean %>% min()
```

    ## [1] 1.333859

``` r
res_hif1a_2a$baseMean %>% hist(breaks=100000, xlim = c(0,100))

# use all genes in dataset
genes_univ_ens <- results_list[[1]]$ENSEMBL
genes_univ_ens %>% length()
```

    ## [1] 21583

``` r
genes_univ <- results_list[[1]]$ENTREZ
genes_univ %>% length()
```

    ## [1] 21583

``` r
genes_univ %>% is.na() %>% summary()
```

    ##    Mode   FALSE    TRUE 
    ## logical   16750    4833

<img src="Readme_files/figure-gfm/pre_results-1.png" width="50%" /><img src="Readme_files/figure-gfm/pre_results-2.png" width="50%" />

# Figure 1: Samples QC

## PCA

``` r
vst_dat <- assay(vst(dds))

p <- pca(vst_dat, metadata = colData(dds), removeVar = 0.99)
pca_table <- cbind(p$rotated,p$metadata)
pca1 <- ggplot(pca_table, aes(PC2, PC1, color=genotype, shape=treatment)) +
  geom_hline(yintercept = 0, linewidth = 0.1) + 
  geom_vline(xintercept = 0, linewidth = 0.1) +
  geom_point(size=4, alpha=0.5, stroke=1) +
  labs(title = "top 1% variable genes") +
  ylab(paste0("PC1: ",p$variance["PC1"] %>% round(digits = 1),"% variance")) +
  xlab(paste0("PC2: ",p$variance["PC2"] %>% round(digits = 1),"% variance")) +
  scale_color_manual(values=colors[c(2,4,6,8)]) +
  scale_shape_manual(values = c(21,16)) + 
  scale_fill_manual(values=c(colors[1],"white",colors[2],"white")) +
  theme_bw() +
  scale_y_reverse() +
  removeGrid(x=T, y=T)

p <- pca(vst_dat, metadata = colData(dds), removeVar = 0.95)
pca_table <- cbind(p$rotated,p$metadata)
pca5 <- ggplot(pca_table, aes(PC2, PC1, color=genotype, shape=treatment)) +
  geom_hline(yintercept = 0, linewidth = 0.1) + 
  geom_vline(xintercept = 0, linewidth = 0.1) +
  geom_point(size=4, alpha=0.5, stroke=1) +
  labs(title = "top 5% variable genes") +
  ylab(paste0("PC1: ",p$variance["PC1"] %>% round(digits = 1),"% variance")) +
  xlab(paste0("PC2: ",p$variance["PC2"] %>% round(digits = 1),"% variance")) +
  scale_color_manual(values=colors[c(2,4,6,8)]) +
  scale_shape_manual(values = c(21,16)) + 
  scale_fill_manual(values=c(colors[1],"white",colors[2],"white")) +
  theme_bw() +
  scale_y_reverse() +
  removeGrid(x=T, y=T)

p <- pca(vst_dat, metadata = colData(dds), removeVar = 0.90)
pca_table <- cbind(p$rotated,p$metadata)
pca10 <- ggplot(pca_table, aes(PC2, PC1, color=genotype, shape=treatment)) +
  geom_hline(yintercept = 0, linewidth = 0.1) + 
  geom_vline(xintercept = 0, linewidth = 0.1) +
  geom_point(size=4, alpha=0.5, stroke=1) +
  labs(title = "top 10% variable genes") +
  ylab(paste0("PC1: ",p$variance["PC1"] %>% round(digits = 1),"% variance")) +
  xlab(paste0("PC2: ",p$variance["PC2"] %>% round(digits = 1),"% variance")) +
  scale_color_manual(values=colors[c(2,4,6,8)]) +
  scale_shape_manual(values = c(21,16)) + 
  scale_fill_manual(values=c(colors[1],"white",colors[2],"white")) +
  theme_bw() +
  scale_x_reverse() +
  removeGrid(x=T, y=T)
```

``` r
pca1+pca5+pca10+ plot_layout(guides = "collect", axes="collect", axis_titles="collect") & 
  theme(legend.position = 'bottom')
```

![](Readme_files/figure-gfm/1_pca_plot_variants-1.png)<!-- -->

``` r
pca10
```

<img src="Readme_files/figure-gfm/1_pca_plot_final-1.png" width="50%" />

# Figure 2: Differential expressed genes

## -Volcano_function

``` r
getdeg <- function(x){subset(results_list[[x]], padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) %>% data.frame()}

volcano_sk3 <- function(n,
                        col="red",
                        celline="cells",
                        deg=deg) {
xlim <- 12
ylim <- -250
res <- results_list[[n]]
res <- res_shrink_list[[n]] %>% data.frame()

# of deg genes
up <- subset(deg, log2FoldChange > 1) %>% nrow()
down <- subset(deg, log2FoldChange < -1) %>% nrow()
total <- up+down

# points outside the grid
outx <- subset(res, log2FoldChange > xlim | log2FoldChange < -xlim) %>% rownames()
outy <- subset(res, padj < 10^ylim) %>% rownames()

res$outlier <- ifelse(rownames(res) %in% c(outx,outy),"yes","no")
res$deg <- ifelse(rownames(res) %in% rownames(deg),"hypoxic","n.s.") %>% factor()

res <- res %>% arrange(desc(res$deg))

res[outx,"log2FoldChange"] <- ifelse(res[outx,"log2FoldChange"] > xlim,xlim,-xlim)
res[outy,"padj"] <- 10^ylim

volcano_func <- ggplot(res,aes(x=log2FoldChange,y=-log10(padj),color=deg, shape=outlier, fill=deg)) +
  geom_hline(yintercept = 0, linewidth = 0.2) + 
  geom_vline(xintercept = 0, linewidth = 0.2) +
  geom_point(size=1.5, stroke=0.5) +
  scale_shape_manual(values = c(21,3)) + 
  scale_alpha_manual(values = c(0.3,0.6)) + 
  labs(title=paste0("Hypoxic response in ",celline),
       subtitle = paste0("upregulated: ",up,", downregulated: ",down," (total: ",total,")") )+
  theme(plot.title = element_text(size = 1), 
        plot.subtitle = element_text(size = 0.5) )+
  ylab("padj (-log10)") +
  xlab("log2-foldchange") +
  scale_fill_manual(values = alpha(c(lighten(col,0.3),"grey70"),0.5)) + 
  scale_color_manual(values = c(col,"grey40")) + 
  theme_bw() +
  # geom_text_repel(label=res$symbol, color="black") + 
  removeGrid(x=T, y=T)
volcano_func
}
```

## -Volcano_function2

``` r
getdeg <- function(x){subset(results_list[[x]], padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) %>% data.frame()}

volcano_sk4 <- function(n,
                        col="red",
                        col2="blue",
                        celline="cells",
                        deg=deg,
                        deg2=deg2) {
xlim <- 12
ylim <- -250
res <- results_list[[n]]
res <- res_shrink_list[[n]] %>% data.frame()

# of deg genes
hx_up <- subset(deg, log2FoldChange > 1) %>% nrow()
hx_down <- subset(deg, log2FoldChange < -1) %>% nrow()
hx_total <- hx_up+hx_down

# of deg2 genes
up <- subset(deg2, log2FoldChange > 1) %>% nrow()
down <- subset(deg2, log2FoldChange < -1) %>% nrow()
total <- up+down

# points outside the grid
outx <- subset(res, log2FoldChange > xlim | log2FoldChange < -xlim) %>% rownames()
outy <- subset(res, padj < 10^ylim) %>% rownames()

res$outlier <- ifelse(rownames(res) %in% c(outx,outy),"yes","no")
res$deg <- ifelse(rownames(res) %in% rownames(deg2),"different from Kelly",
                  ifelse(rownames(res) %in% rownames(deg),"hypoxic","n.s.")) %>% factor()

res <- res %>% arrange(desc(res$deg))

res[outx,"log2FoldChange"] <- ifelse(res[outx,"log2FoldChange"] > xlim,xlim,-xlim)
res[outy,"padj"] <- 10^ylim

volcano_func <- ggplot(res,aes(x=log2FoldChange,y=-log10(padj),color=deg, shape=outlier, fill=deg)) +
  geom_hline(yintercept = 0, linewidth = 0.2) + 
  geom_vline(xintercept = 0, linewidth = 0.2) +
  geom_point(size=1.5, stroke=0.5) +
  scale_shape_manual(values = c(21,3)) + 
  scale_alpha_manual(values = c(0.3,0.6)) + 
  labs(title=paste0("Hypoxic response in ",celline),
       subtitle = paste0("Hypoxic: up ",hx_up,", down ",hx_down," (total ",hx_total,")","\nDifferent from Kelly: up ",up,", down: ",down," (total: ",total,")")) +
  theme(plot.title = element_text(size = 1), 
        plot.subtitle = element_text(size = 0.5) )+
  ylab("padj (-log10)") +
  xlab("log2-foldchange") +
  scale_fill_manual(values = alpha(c(lighten(c(col2,col),0.3),"grey70"),0.5)) + 
  scale_color_manual(values = c(col2,col,"grey40")) + 
  theme_bw() +
  # geom_text_repel(label=res$symbol, color="black") + 
  removeGrid(x=T, y=T)
volcano_func
}
```

## -Plot Vulcanos

### Simple

``` r
# Simple Volcanos
volcano_Kelly <- volcano_sk3(n="Kelly.Hx.vs.Nx", deg=getdeg("Kelly.Hx.vs.Nx"),col=colors[2], celline="Kelly")
volcano_hif1a <- volcano_sk3(n="Hif1a.Hx.vs.Nx", deg=getdeg("Hif1a.Hx.vs.Nx"),col=colors[4], celline="HIF1A")
volcano_hif2a <- volcano_sk3(n="Hif2a.Hx.vs.Nx", deg=getdeg("Hif2a.Hx.vs.Nx"),col=colors[6], celline="HIF2A")
volcano_hif1b <- volcano_sk3(n="Hif1b.Hx.vs.Nx", deg=getdeg("Hif1b.Hx.vs.Nx"),col=colors[8], celline="HIF1B")

(volcano_Kelly+volcano_hif1b + plot_layout(guides = "collect", axes="collect", axis_titles="collect") ) / 
  (volcano_hif1a+volcano_hif2a + plot_layout(guides = "collect", axes="collect", axis_titles="collect") ) & 
  theme(legend.position = 'right')
```

![](Readme_files/figure-gfm/2_volcanos_plot-1.png)<!-- -->

``` r
n2 <- {}

deg <- subset(results_list[["Kelly.Hx.vs.Nx"]], padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) %>% data.frame()


# Two set volcanos
n <- "Hif1a.Hx.vs.Nx"
deg <- subset(results_list[["Hif1a.Hx.vs.Nx"]], padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) %>% data.frame()
deg2 <- subset(results_list[["Hif1aHxNx.vs.KellyHxNx"]], padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) %>% data.frame()
col <- colors[2]
col2 <- colors[4]
celline="HIF1A"

volcano_Kelly <- volcano_sk3(n="Kelly.Hx.vs.Nx", deg=getdeg("Kelly.Hx.vs.Nx"),col=colors[2], celline="Kelly")
volcano_hif1a <- volcano_sk3(n="Hif1aHxNx.vs.KellyHxNx", deg=getdeg("Hif1aHxNx.vs.KellyHxNx"),col=colors[4], celline="HIF1A") + ggtitle(label="Difference in Effect vs. Kelly")
volcano_hif2a <- volcano_sk3(n="Hif2aHxNx.vs.KellyHxNx", deg=getdeg("Hif2aHxNx.vs.KellyHxNx"),col=colors[6], celline="HIF2A")  + ggtitle(label="Difference in Effect vs. Kelly")
volcano_hif1b <- volcano_sk3(n="Hif1bHxNx.vs.KellyHxNx", deg=getdeg("Hif1bHxNx.vs.KellyHxNx"),col=colors[8], celline="HIF1B")  + ggtitle(label="Difference in Effect vs. Kelly")

(volcano_Kelly+volcano_hif1b + plot_layout(guides = "collect", axes="collect", axis_titles="collect") ) / 
  (volcano_hif1a+volcano_hif2a + plot_layout(guides = "collect", axes="collect", axis_titles="collect") ) & 
  theme(legend.position = 'right')
```

![](Readme_files/figure-gfm/2_volcanos_plot-2.png)<!-- -->

### Advanced

``` r
# Two set volcanos
n <- "Hif1a.Hx.vs.Nx"
deg <- subset(results_list[["Hif1a.Hx.vs.Nx"]], padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) %>% data.frame()
deg2 <- subset(results_list[["Hif1aHxNx.vs.KellyHxNx"]], padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) %>% data.frame()
col <- colors[2]
col2 <- colors[4]
celline="HIF1A"

volcano_Kelly <- volcano_sk3(n="Kelly.Hx.vs.Nx", deg=getdeg("Kelly.Hx.vs.Nx"),col=colors[2], celline="Kelly")
volcano_hif1a <- volcano_sk4(n="Hif1a.Hx.vs.Nx", deg=getdeg("Hif1a.Hx.vs.Nx"), deg2=getdeg("Hif1aHxNx.vs.KellyHxNx"),col=colors[2], col2=colors[4], celline="HIF1A")
volcano_hif2a <- volcano_sk4(n="Hif2a.Hx.vs.Nx", deg=getdeg("Hif2a.Hx.vs.Nx"), deg2=getdeg("Hif2aHxNx.vs.KellyHxNx"),col=colors[2], col2=colors[6], celline="HIF2A")
volcano_hif1b <- volcano_sk4(n="Hif1b.Hx.vs.Nx", deg=getdeg("Hif1b.Hx.vs.Nx"), deg2=getdeg("Hif1bHxNx.vs.KellyHxNx"),col=colors[2], col2=colors[8], celline="HIF1B")

(volcano_Kelly+volcano_hif1b + plot_layout(guides = "collect", axes="collect", axis_titles="collect") ) / 
  (volcano_hif1a+volcano_hif2a + plot_layout(guides = "collect", axes="collect", axis_titles="collect") ) & 
  theme(legend.position = 'right')
```

![](Readme_files/figure-gfm/2_volcanos_plot2-1.png)<!-- -->

## -Venn

``` r
input_list <- main_degs
plt1 <- venn.diagram(
    x = input_list,
    fill = colors[c(2,7,3,5)],
    main.fontface = "bold",
    fontfamily ="Arial",
    category.names = paste(names(input_list),"\n(",input_list %>% summary() %>% .[c(1:length(input_list))],")",sep=""),
    force.unique = TRUE, na = "remove", total.population = TRUE,
    filename = NULL,
    lwd = 2,
    lty = 'blank',
    cat.fontface = "bold",
    cat.fontfamily = "arial")

input_list <- main_degs[c(3,4,1)]
plt2 <- venn.diagram(
    x = input_list,
    fill = colors[c(3,5,2)],
    main.fontface = "bold",
    fontfamily ="Arial",
    category.names = paste(names(input_list),"\n(",input_list %>% summary() %>% .[c(1:length(input_list))],")",sep=""),
    force.unique = TRUE, na = "remove", total.population = TRUE,
    filename = NULL,
    lwd = 2,
    lty = 'blank',
    cat.fontface = "bold",
    cat.fontfamily = "arial")
    
#     main = "Compare Hif KOs",


patchwork::wrap_elements(plt1) / patchwork::wrap_elements(plt2)
```

<img src="Readme_files/figure-gfm/2_venn-1.png" width="50%" />

# Figure 3: Gene Cluster

## Gene Cluster

``` r
# cluster_venn

main_degs %>% names()
```

    ## [1] "Kelly: Hx.vs.Nx" "Hif1b"           "Hif1a"           "Hif2a"

``` r
length(main_degs[[3]])
```

    ## [1] 863

``` r
length(main_degs[[4]])
```

    ## [1] 2856

``` r
venns <- calculate.overlap(main_degs[c(3,4)])
lapply(venns,length)
```

    ## $a1
    ## [1] 863
    ## 
    ## $a2
    ## [1] 2856
    ## 
    ## $a3
    ## [1] 324

``` r
venns %>% unlist() %>% length()
```

    ## [1] 4043

``` r
venns %>% unlist() %>% unique() %>% length()
```

    ## [1] 3395

``` r
venns$a1 <- setdiff(venns$a1,venns$a3)
venns$a2 <- setdiff(venns$a2,venns$a3)
lapply(venns,length)
```

    ## $a1
    ## [1] 539
    ## 
    ## $a2
    ## [1] 2532
    ## 
    ## $a3
    ## [1] 324

``` r
venns %>% unlist() %>% length()
```

    ## [1] 3395

``` r
venns %>% unlist() %>% unique() %>% length()
```

    ## [1] 3395

``` r
res_hif1a_2a$venn <- ifelse(rownames(res_hif1a_2a) %in% venns$a1,"HIF1A",
                      ifelse(rownames(res_hif1a_2a) %in% venns$a2,"HIF2A",
                      ifelse(rownames(res_hif1a_2a) %in% venns$a3,"overlap","interaction")))
res_hif1a_2a$venn %>% table()
```

    ## .
    ##       HIF1A       HIF2A interaction     overlap 
    ##         539        2532        1979         324

``` r
# Cluster Venn
cluster_venn <- ggplot(res_hif1a_2a,aes(x=Hif1aHxNx.vs.KellyHxNx.log2FoldChange, y=Hif2aHxNx.vs.KellyHxNx.log2FoldChange, color=venn, fill=venn)) +
  geom_hline(yintercept = 0, linewidth = 0.1) + 
  geom_vline(xintercept = 0, linewidth = 0.1) +
  geom_point(size=1, stroke=0.5, shape=21) +
  labs(title = "Simple/Venn Cluster") +
  xlab("Hif1a vs. Kelly") +
  ylab("Hif2a vs. Kelly") +
  scale_color_manual(values=c(colors[c(4,6)],"orange",colors[2])) +
  # scale_shape_manual(values = c(21,16)) + 
  scale_fill_manual(values=alpha(c(colors[c(4,6)],"orange",colors[2]),0.2)) +
  theme_bw() +
  removeGrid(x=T, y=T) +
  coord_cartesian(xlim = c(-5, 5),ylim = c(-5,5))

res_hif1a_2a$group %>% table()
```

    ## < table of extent 0 >

``` r
# Manual Cluster

hif_parallel  <- res_hif1a_2a %>% filter(abs(Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange) < 1)

hif2a_up <- res_hif1a_2a %>% filter(Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange > 1 & (Hif2aHxNx.vs.KellyHxNx.log2FoldChange-1 > 2*-Hif1aHxNx.vs.KellyHxNx.log2FoldChange))
hif2a_do <- res_hif1a_2a %>% filter(Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange < -1 & (Hif2aHxNx.vs.KellyHxNx.log2FoldChange < -Hif1aHxNx.vs.KellyHxNx.log2FoldChange+1))

hif1a_up <- res_hif1a_2a %>% filter(Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange > 1 & (-Hif1aHxNx.vs.KellyHxNx.log2FoldChange > (abs(Hif2aHxNx.vs.KellyHxNx.log2FoldChange-1))))

hif1a_do <- res_hif1a_2a %>% filter(Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange < -1 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange-1 > 2*-Hif2aHxNx.vs.KellyHxNx.log2FoldChange))


res_hif1a_2a$group <- ifelse(rownames(res_hif1a_2a) %in% rownames(hif_parallel),"HIF1A_HIF2A",
                      ifelse(rownames(res_hif1a_2a) %in% rownames(hif2a_up),"HIF2A",
                      ifelse(rownames(res_hif1a_2a) %in% rownames(hif2a_do),"HIF2A",
                      ifelse(rownames(res_hif1a_2a) %in% rownames(hif1a_up),"HIF1A",
                      ifelse(rownames(res_hif1a_2a) %in% rownames(hif1a_do),"HIF1A","orange")
                      ))))

cluster <- ggplot(res_hif1a_2a,aes(x=Hif1aHxNx.vs.KellyHxNx.log2FoldChange, y=Hif2aHxNx.vs.KellyHxNx.log2FoldChange, color=group, fill=group)) +
  geom_hline(yintercept = 0, linewidth = 0.1) + 
  geom_vline(xintercept = 0, linewidth = 0.1) +
  geom_abline(intercept=c(1,-1)) +
  geom_abline(slope=c(-1), intercept = 1) +
  annotate("segment", x = c(0,1), y = c(1,0), xend = c(-10,11), yend = c(21,-5),color="black") +
  geom_point(size=1, stroke=0.5, shape=21) +
  labs(title = "Geometric Cluster") +
  xlab("Hif1a vs. Kelly") +
  ylab("Hif2a vs. Kelly") +
  scale_color_manual(values=c(colors[c(4,2,6)],"orange")) +
  # scale_shape_manual(values = c(21,16)) + 
  scale_fill_manual(values=alpha(c(colors[c(4,2,6)],"orange"),0.2)) +
  theme_bw() +
  removeGrid(x=T, y=T) +
  coord_cartesian(xlim = c(-5, 5),ylim = c(-5,5))

cluster_venn + cluster
```

![](Readme_files/figure-gfm/cluster-1.png)<!-- -->

``` r
# only genes with p < 0.05
colnames(res_hif1a_2a)
```

    ##  [1] "ENSEMBL"                              
    ##  [2] "ENTREZ"                               
    ##  [3] "symbol"                               
    ##  [4] "baseMean"                             
    ##  [5] "Kelly.Hx.vs.Nx.log2FoldChange"        
    ##  [6] "Kelly.Hx.vs.Nx.padj"                  
    ##  [7] "Hif1a.Hx.vs.Nx.log2FoldChange"        
    ##  [8] "Hif1a.Hx.vs.Nx.padj"                  
    ##  [9] "Hif2a.Hx.vs.Nx.log2FoldChange"        
    ## [10] "Hif2a.Hx.vs.Nx.padj"                  
    ## [11] "Hif1aHxNx.vs.KellyHxNx.log2FoldChange"
    ## [12] "Hif1aHxNx.vs.KellyHxNx.padj"          
    ## [13] "Hif2aHxNx.vs.KellyHxNx.log2FoldChange"
    ## [14] "Hif2aHxNx.vs.KellyHxNx.padj"          
    ## [15] "Hif1bHxNx.vs.KellyHxNx.log2FoldChange"
    ## [16] "Hif1bHxNx.vs.KellyHxNx.padj"          
    ## [17] "Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange"
    ## [18] "Hif2aHxNx.vs.Hif1aHxNx.padj"          
    ## [19] "venn"                                 
    ## [20] "group"

``` r
res_hif1a_2a_p <- filter(res_hif1a_2a, Hif1aHxNx.vs.KellyHxNx.padj < 0.05 & Hif2aHxNx.vs.KellyHxNx.padj < 0.05)

cluster_p <- ggplot(res_hif1a_2a_p,aes(x=Hif1aHxNx.vs.KellyHxNx.log2FoldChange, y=Hif2aHxNx.vs.KellyHxNx.log2FoldChange, color=group, fill=group)) +
  geom_hline(yintercept = 0, linewidth = 0.1) + 
  geom_vline(xintercept = 0, linewidth = 0.1) +
  geom_abline(intercept=c(1,-1)) +
  geom_abline(slope=c(-1), intercept = 1) +
  annotate("segment", x = c(0,1), y = c(1,0), xend = c(-10,11), yend = c(21,-5),color="black") +
  geom_point(size=1, stroke=0.5, shape=21) +
  labs(title = "Geometric Cluster (all p < 0.05") +
  xlab("Hif1a vs. Kelly") +
  ylab("Hif2a vs. Kelly") +
  scale_color_manual(values=c(colors[c(4,2,6)],"orange")) +
  # scale_shape_manual(values = c(21,16)) + 
  scale_fill_manual(values=alpha(c(colors[c(4,2,6)],"orange"),0.2)) +
  theme_bw() +
  removeGrid(x=T, y=T) +
  coord_cartesian(xlim = c(-5, 5),ylim = c(-5,5))

cluster_p
```

![](Readme_files/figure-gfm/cluster-2.png)<!-- -->

``` r
cluster_venn + cluster
```

![](Readme_files/figure-gfm/cluster-3.png)<!-- -->

``` r
write.xlsx(res_hif1a_2a,"DEG_genes.xlsx")
```

# Table 1: Gene List

``` r
HIF1A_genes <- res_hif1a_2a_p %>% 
  .[order(abs(.$Hif1aHxNx.vs.KellyHxNx.log2FoldChange), decreasing = TRUE),] %>% 
  filter(group =="HIF1A") %>%
  filter(baseMean > 1000)
  
nrow(HIF1A_genes)
```

    ## [1] 380

``` r
HIF1A_genes[1:50,] %>% kable()
```

|  | ENSEMBL | ENTREZ | symbol | baseMean | Kelly.Hx.vs.Nx.log2FoldChange | Kelly.Hx.vs.Nx.padj | Hif1a.Hx.vs.Nx.log2FoldChange | Hif1a.Hx.vs.Nx.padj | Hif2a.Hx.vs.Nx.log2FoldChange | Hif2a.Hx.vs.Nx.padj | Hif1aHxNx.vs.KellyHxNx.log2FoldChange | Hif1aHxNx.vs.KellyHxNx.padj | Hif2aHxNx.vs.KellyHxNx.log2FoldChange | Hif2aHxNx.vs.KellyHxNx.padj | Hif1bHxNx.vs.KellyHxNx.log2FoldChange | Hif1bHxNx.vs.KellyHxNx.padj | Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange | Hif2aHxNx.vs.Hif1aHxNx.padj | venn | group |
|:---|:---|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---|
| ENSG00000114023 | ENSG00000114023 | 26355 | FAM162A | 12808.121 | 2.6749198 | 0.0000000 | -1.0871282 | 0.0000000 | 4.1525578 | 0.0000000 | -3.762048 | 0.0000000 | 1.4776380 | 0.0000000 | 0.4415823 | 0.0460083 | 5.239686 | 0 | overlap | HIF1A |
| ENSG00000176171 | ENSG00000176171 | 664 | BNIP3 | 21434.708 | 3.1965388 | 0.0000000 | -0.3905581 | 0.0021865 | 4.0130817 | 0.0000000 | -3.587097 | 0.0000000 | 0.8165428 | 0.0000005 | -0.0984268 | 0.6234311 | 4.403640 | 0 | HIF1A | HIF1A |
| ENSG00000100314 | ENSG00000100314 | 164633 | CABP7 | 2092.285 | 8.1923082 | 0.0000000 | 4.9574998 | 0.0000000 | 9.3811748 | 0.0000000 | -3.234808 | 0.0000000 | 1.1888666 | 0.0060813 | -1.2450997 | 0.0047254 | 4.423675 | 0 | overlap | HIF1A |
| ENSG00000159208 | ENSG00000159208 | 148523 | CIART | 2103.713 | 2.8415121 | 0.0000000 | -0.0961627 | 0.5848455 | 3.8535705 | 0.0000000 | -2.937675 | 0.0000000 | 1.0120583 | 0.0000005 | -0.0590673 | 0.8220711 | 3.949733 | 0 | overlap | HIF1A |
| ENSG00000074800 | ENSG00000074800 | 2023 | ENO1 | 193555.770 | 1.2878157 | 0.0000000 | -1.6300196 | 0.0000000 | 2.3330043 | 0.0000000 | -2.917835 | 0.0000000 | 1.0451886 | 0.0000209 | 1.0149729 | 0.0000614 | 3.963024 | 0 | overlap | HIF1A |
| ENSG00000079739 | ENSG00000079739 | 5236 | PGM1 | 6602.812 | 1.2982488 | 0.0000000 | -1.1518298 | 0.0000000 | 2.2307331 | 0.0000000 | -2.450079 | 0.0000000 | 0.9324843 | 0.0003620 | 0.4075537 | 0.1536669 | 3.382563 | 0 | HIF1A | HIF1A |
| ENSG00000102144 | ENSG00000102144 | 5230 | PGK1 | 69853.472 | 2.3595632 | 0.0000000 | 0.1341692 | 0.4241579 | 2.9941319 | 0.0000000 | -2.225394 | 0.0000000 | 0.6345687 | 0.0025437 | 0.1551542 | 0.5212432 | 2.859963 | 0 | HIF1A | HIF1A |
| ENSG00000165802 | ENSG00000165802 | 26012 | NSMF | 12663.169 | 0.8783827 | 0.0000000 | -1.2842273 | 0.0000000 | 2.3874643 | 0.0000000 | -2.162610 | 0.0000000 | 1.5090817 | 0.0000000 | -0.3505130 | 0.1519806 | 3.671692 | 0 | overlap | HIF1A |
| ENSG00000291995 | ENSG00000291995 | 23015 | GOLGA8A | 11149.052 | 2.5231385 | 0.0000000 | 0.4428410 | 0.2248031 | 3.6046852 | 0.0000000 | -2.080297 | 0.0000151 | 1.0815466 | 0.0251738 | -0.3801348 | 0.4750074 | 3.161844 | 0 | overlap | HIF1A |
| ENSG00000114480 | ENSG00000114480 | 2632 | GBE1 | 3845.354 | 2.9488275 | 0.0000000 | 0.9458202 | 0.0000000 | 3.8232133 | 0.0000000 | -2.003007 | 0.0000000 | 0.8743858 | 0.0000495 | -1.3279561 | 0.0000000 | 2.877393 | 0 | HIF1A | HIF1A |
| ENSG00000173157 | ENSG00000173157 | 80070 | ADAMTS20 | 1182.692 | 2.5175767 | 0.0000000 | 0.5801152 | 0.0273490 | 3.4119686 | 0.0000000 | -1.937461 | 0.0000000 | 0.8943919 | 0.0111511 | -1.0710456 | 0.0023819 | 2.831853 | 0 | HIF1A | HIF1A |
| ENSG00000291647 | ENSG00000291647 | NA |  | 1641.217 | 2.9399548 | 0.0000000 | 1.0043865 | 0.0000008 | 3.5958931 | 0.0000000 | -1.935568 | 0.0000000 | 0.6559383 | 0.0212520 | -2.2176404 | 0.0000000 | 2.591507 | 0 | HIF1A | HIF1A |
| ENSG00000143590 | ENSG00000143590 | 1944 | EFNA3 | 1307.360 | 1.4326017 | 0.0000000 | -0.3787733 | 0.0541316 | 2.4059099 | 0.0000000 | -1.811375 | 0.0000000 | 0.9733082 | 0.0000953 | -0.1161714 | 0.7046982 | 2.784683 | 0 | HIF1A | HIF1A |
| ENSG00000115657 | ENSG00000115657 | 10058 | ABCB6 | 2501.505 | 0.7838843 | 0.0000000 | -0.9708463 | 0.0000000 | 1.4522205 | 0.0000000 | -1.754730 | 0.0000000 | 0.6683363 | 0.0010309 | -0.3808896 | 0.0797277 | 2.423067 | 0 | HIF1A | HIF1A |
| ENSG00000152256 | ENSG00000152256 | 5163 | PDK1 | 13129.395 | 2.4376427 | 0.0000000 | 0.7481050 | 0.0000005 | 3.9562136 | 0.0000000 | -1.689538 | 0.0000000 | 1.5185709 | 0.0000000 | 1.2052944 | 0.0000000 | 3.208109 | 0 | overlap | HIF1A |
| ENSG00000214063 | ENSG00000214063 | 7106 | TSPAN4 | 1420.589 | 1.7102559 | 0.0000000 | 0.0587888 | 0.7744888 | 2.2767738 | 0.0000000 | -1.651467 | 0.0000000 | 0.5665179 | 0.0252022 | -0.1410303 | 0.6242482 | 2.217985 | 0 | HIF1A | HIF1A |
| ENSG00000169299 | ENSG00000169299 | 55276 | PGM2 | 2050.709 | 0.6660084 | 0.0000000 | -0.9834761 | 0.0000000 | 1.7450178 | 0.0000000 | -1.649485 | 0.0000000 | 1.0790094 | 0.0000000 | -0.5900683 | 0.0008102 | 2.728494 | 0 | overlap | HIF1A |
| ENSG00000177181 | ENSG00000177181 | 284716 | RIMKLA | 2825.763 | 1.4243237 | 0.0000000 | -0.1981704 | 0.0453906 | 1.7606305 | 0.0000000 | -1.622494 | 0.0000000 | 0.3363068 | 0.0105027 | -0.4479529 | 0.0005788 | 1.958801 | 0 | HIF1A | HIF1A |
| ENSG00000083444 | ENSG00000083444 | 5351 | PLOD1 | 5404.734 | 1.5441125 | 0.0000000 | -0.0653180 | 0.6936673 | 2.0711872 | 0.0000000 | -1.609430 | 0.0000000 | 0.5270747 | 0.0099745 | 1.0253415 | 0.0000002 | 2.136505 | 0 | HIF1A | HIF1A |
| ENSG00000054967 | ENSG00000054967 | 84957 | RELT | 1117.456 | -0.1053715 | 0.3010434 | -1.6879889 | 0.0000000 | 0.7857247 | 0.0000000 | -1.582617 | 0.0000000 | 0.8910962 | 0.0000003 | 0.7825556 | 0.0000147 | 2.473714 | 0 | HIF1A | HIF1A |
| ENSG00000130810 | ENSG00000130810 | 56342 | PPAN | 5041.489 | 0.1463942 | 0.1145844 | -1.3920492 | 0.0000000 | 0.6506359 | 0.0000003 | -1.538443 | 0.0000000 | 0.5042417 | 0.0030008 | 0.0157866 | 0.9418571 | 2.042685 | 0 | HIF1A | HIF1A |
| ENSG00000183258 | ENSG00000183258 | 51428 | DDX41 | 7789.756 | 1.1151273 | 0.0000000 | -0.4218511 | 0.0336632 | 1.7056059 | 0.0000000 | -1.536978 | 0.0000000 | 0.5904786 | 0.0295347 | -0.6454242 | 0.0165710 | 2.127457 | 0 | HIF1A | HIF1A |
| ENSG00000087116 | ENSG00000087116 | 9509 | ADAMTS2 | 6978.341 | 1.3396820 | 0.0000000 | -0.1874562 | 0.4112692 | 3.2529619 | 0.0000000 | -1.527138 | 0.0000001 | 1.9132800 | 0.0000000 | 1.7709349 | 0.0000000 | 3.440418 | 0 | overlap | HIF1A |
| ENSG00000291087 | ENSG00000291087 | 1416 | CRYBB2P1 | 3245.348 | 0.5887326 | 0.0000000 | -0.8993326 | 0.0000000 | 1.6770023 | 0.0000000 | -1.488065 | 0.0000000 | 1.0882697 | 0.0000000 | -0.4509091 | 0.0209510 | 2.576335 | 0 | overlap | HIF1A |
| ENSG00000181418 | ENSG00000181418 | 23109 | DDN | 2123.712 | 0.2386186 | 0.0072850 | -1.2248994 | 0.0000000 | 1.1341277 | 0.0000000 | -1.463518 | 0.0000000 | 0.8955091 | 0.0000000 | 0.6749770 | 0.0000410 | 2.359027 | 0 | HIF1A | HIF1A |
| ENSG00000166123 | ENSG00000166123 | 84706 | GPT2 | 3610.724 | 0.3468948 | 0.0169661 | -1.0667671 | 0.0000001 | 1.6030059 | 0.0000000 | -1.413662 | 0.0000002 | 1.2561111 | 0.0000011 | 0.6063258 | 0.0306415 | 2.669773 | 0 | overlap | HIF1A |
| ENSG00000162373 | ENSG00000162373 | 79656 | BEND5 | 1852.010 | 1.8686933 | 0.0000000 | 0.4658387 | 0.0262085 | 2.6733978 | 0.0000000 | -1.402855 | 0.0000005 | 0.8047044 | 0.0035992 | -0.8981269 | 0.0012644 | 2.207559 | 0 | HIF1A | HIF1A |
| ENSG00000171314 | ENSG00000171314 | 5223 | PGAM1 | 34522.854 | 0.5383551 | 0.0000003 | -0.8549299 | 0.0000000 | 1.0966365 | 0.0000000 | -1.393285 | 0.0000000 | 0.5582815 | 0.0055043 | 1.1900214 | 0.0000000 | 1.951566 | 0 | HIF1A | HIF1A |
| ENSG00000165704 | ENSG00000165704 | 3251 | HPRT1 | 1558.059 | -0.8977131 | 0.0000000 | -2.2709809 | 0.0000000 | -0.0651497 | 0.7376740 | -1.373268 | 0.0000000 | 0.8325634 | 0.0000515 | 1.0418365 | 0.0000006 | 2.205831 | 0 | HIF1A | HIF1A |
| ENSG00000069998 | ENSG00000069998 | 27440 | HDHD5 | 8166.301 | 0.3563772 | 0.0000263 | -0.9539115 | 0.0000000 | 1.1183386 | 0.0000000 | -1.310289 | 0.0000000 | 0.7619613 | 0.0000007 | 0.3620188 | 0.0298267 | 2.072250 | 0 | HIF1A | HIF1A |
| ENSG00000165195 | ENSG00000165195 | 5277 | PIGA | 1211.285 | -0.6756653 | 0.0000000 | -1.9621120 | 0.0000000 | -0.0730204 | 0.6483747 | -1.286447 | 0.0000000 | 0.6026448 | 0.0006117 | 0.5437519 | 0.0027528 | 1.889092 | 0 | HIF1A | HIF1A |
| ENSG00000107077 | ENSG00000107077 | 23081 | KDM4C | 2127.336 | 1.0806486 | 0.0000000 | -0.1785491 | 0.0875722 | 1.8062271 | 0.0000000 | -1.259198 | 0.0000000 | 0.7255785 | 0.0000000 | -0.3423431 | 0.0143658 | 1.984776 | 0 | HIF1A | HIF1A |
| ENSG00000104722 | ENSG00000104722 | 4741 | NEFM | 19345.680 | -0.4769474 | 0.0100463 | -1.7302990 | 0.0000000 | 0.3121124 | 0.2828241 | -1.253351 | 0.0007492 | 0.7890599 | 0.0277671 | -0.3631407 | 0.3450017 | 2.042411 | 0 | HIF1A | HIF1A |
| ENSG00000198435 | ENSG00000198435 | 441478 | NRARP | 1158.396 | -0.3450409 | 0.0000329 | -1.5759595 | 0.0000000 | 1.6378908 | 0.0000000 | -1.230919 | 0.0000000 | 1.9829317 | 0.0000000 | 1.0182673 | 0.0000000 | 3.213850 | 0 | overlap | HIF1A |
| ENSG00000175414 | ENSG00000175414 | 285598 | ARL10 | 5805.146 | 1.6448342 | 0.0000000 | 0.4187043 | 0.0000856 | 2.0787476 | 0.0000000 | -1.226130 | 0.0000000 | 0.4339134 | 0.0025046 | -0.8892701 | 0.0000000 | 1.660043 | 0 | HIF1A | HIF1A |
| ENSG00000138336 | ENSG00000138336 | 80312 | TET1 | 1910.720 | 1.3262678 | 0.0000000 | 0.1115671 | 0.5853869 | 1.9718809 | 0.0000000 | -1.214701 | 0.0000021 | 0.6456131 | 0.0118879 | -1.2726230 | 0.0000003 | 1.860314 | 0 | HIF1A | HIF1A |
| ENSG00000157510 | ENSG00000157510 | 134265 | AFAP1L1 | 1016.267 | 2.0789568 | 0.0000000 | 0.8683951 | 0.0000000 | 2.9926345 | 0.0000000 | -1.210562 | 0.0000000 | 0.9136777 | 0.0000021 | -0.7673471 | 0.0001421 | 2.124239 | 0 | HIF1A | HIF1A |
| ENSG00000109107 | ENSG00000109107 | 230 | ALDOC | 7804.931 | 2.3700489 | 0.0000000 | 1.1714388 | 0.0000001 | 3.1096082 | 0.0000000 | -1.198610 | 0.0001162 | 0.7395593 | 0.0147225 | 0.6930218 | 0.0237489 | 1.938169 | 0 | HIF1A | HIF1A |
| ENSG00000067992 | ENSG00000067992 | 5165 | PDK3 | 1826.282 | 1.1202125 | 0.0000000 | -0.0589904 | 0.6340079 | 1.4968928 | 0.0000000 | -1.179203 | 0.0000000 | 0.3766802 | 0.0149725 | -0.5666552 | 0.0002000 | 1.555883 | 0 | HIF1A | HIF1A |
| ENSG00000048028 | ENSG00000048028 | 57646 | USP28 | 2720.259 | 0.2200945 | 0.0033335 | -0.9466727 | 0.0000000 | 1.1181716 | 0.0000000 | -1.166767 | 0.0000000 | 0.8980771 | 0.0000000 | -0.2236986 | 0.1370979 | 2.064844 | 0 | HIF1A | HIF1A |
| ENSG00000196365 | ENSG00000196365 | 9361 | LONP1 | 11682.548 | 0.7728883 | 0.0000000 | -0.3646125 | 0.0264615 | 1.4185936 | 0.0000000 | -1.137501 | 0.0000002 | 0.6457054 | 0.0028933 | 0.2696479 | 0.2565138 | 1.783206 | 0 | HIF1A | HIF1A |
| ENSG00000109519 | ENSG00000109519 | 80273 | GRPEL1 | 1308.634 | 0.0032781 | 0.9757183 | -1.1317726 | 0.0000000 | 0.9223537 | 0.0000000 | -1.135051 | 0.0000000 | 0.9190756 | 0.0000001 | 0.6698318 | 0.0001644 | 2.054126 | 0 | HIF1A | HIF1A |
| ENSG00000131398 | ENSG00000131398 | 3748 | KCNC3 | 1903.200 | 0.6587871 | 0.0000004 | -0.4754457 | 0.0108936 | 1.4764342 | 0.0000000 | -1.134233 | 0.0000082 | 0.8176470 | 0.0008553 | -0.0990200 | 0.7430841 | 1.951880 | 0 | HIF1A | HIF1A |
| ENSG00000004455 | ENSG00000004455 | 204 | AK2 | 12512.800 | 0.1006023 | 0.1577137 | -1.0255227 | 0.0000000 | 1.0699390 | 0.0000000 | -1.126125 | 0.0000000 | 0.9693367 | 0.0000000 | 0.2961433 | 0.0277512 | 2.095462 | 0 | HIF1A | HIF1A |
| ENSG00000160124 | ENSG00000160124 | 131076 | MIX23 | 1710.121 | 0.2588632 | 0.0002513 | -0.8618296 | 0.0000000 | 1.0009414 | 0.0000000 | -1.120693 | 0.0000000 | 0.7420782 | 0.0000000 | 0.1211862 | 0.4188041 | 1.862771 | 0 | HIF1A | HIF1A |
| ENSG00000115457 | ENSG00000115457 | 3485 | IGFBP2 | 14451.065 | -0.3253363 | 0.0087879 | -1.4435181 | 0.0000000 | 0.9856731 | 0.0000000 | -1.118182 | 0.0000021 | 1.3110093 | 0.0000000 | 0.4756764 | 0.0498509 | 2.429191 | 0 | overlap | HIF1A |
| ENSG00000139352 | ENSG00000139352 | 429 | ASCL1 | 4728.131 | -1.5321100 | 0.0000000 | -2.6465436 | 0.0000000 | -0.1767893 | 0.5408597 | -1.114434 | 0.0018100 | 1.3553207 | 0.0000190 | 1.6617711 | 0.0000003 | 2.469754 | 0 | overlap | HIF1A |
| ENSG00000259330 | ENSG00000259330 | 100505573 | INAFM2 | 1305.396 | 0.5907874 | 0.0000000 | -0.5225746 | 0.0000005 | 0.9156838 | 0.0000000 | -1.113362 | 0.0000000 | 0.3248963 | 0.0251662 | -0.2749454 | 0.0627948 | 1.438258 | 0 | HIF1A | HIF1A |
| ENSG00000020633 | ENSG00000020633 | 864 | RUNX3 | 1220.647 | 1.2189673 | 0.0000000 | 0.1081430 | 0.6173585 | 2.0627922 | 0.0000000 | -1.110824 | 0.0000610 | 0.8438250 | 0.0013813 | -0.6277354 | 0.0226062 | 1.954649 | 0 | HIF1A | HIF1A |
| ENSG00000197872 | ENSG00000197872 | 81553 | CYRIA | 2165.458 | 0.9297733 | 0.0000000 | -0.1691507 | 0.4113193 | 2.0053403 | 0.0000000 | -1.098924 | 0.0000438 | 1.0755671 | 0.0000165 | -1.0914347 | 0.0000233 | 2.174491 | 0 | overlap | HIF1A |

``` r
plotCounts_SK(HIF1A_genes[1:9,] %>% rownames()) + labs(title="Hif1A")
```

![](Readme_files/figure-gfm/gene_lists-1.png)<!-- -->

``` r
HIF2A_genes <- res_hif1a_2a_p %>% 
  .[order(abs(.$Hif1aHxNx.vs.KellyHxNx.log2FoldChange), decreasing = TRUE),] %>% 
  filter(group =="HIF2A") %>%
  filter(baseMean > 1000)
nrow(HIF2A_genes)
```

    ## [1] 163

``` r
HIF2A_genes[1:50,] %>% kable()
```

|  | ENSEMBL | ENTREZ | symbol | baseMean | Kelly.Hx.vs.Nx.log2FoldChange | Kelly.Hx.vs.Nx.padj | Hif1a.Hx.vs.Nx.log2FoldChange | Hif1a.Hx.vs.Nx.padj | Hif2a.Hx.vs.Nx.log2FoldChange | Hif2a.Hx.vs.Nx.padj | Hif1aHxNx.vs.KellyHxNx.log2FoldChange | Hif1aHxNx.vs.KellyHxNx.padj | Hif2aHxNx.vs.KellyHxNx.log2FoldChange | Hif2aHxNx.vs.KellyHxNx.padj | Hif1bHxNx.vs.KellyHxNx.log2FoldChange | Hif1bHxNx.vs.KellyHxNx.padj | Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange | Hif2aHxNx.vs.Hif1aHxNx.padj | venn | group |
|:---|:---|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---|
| ENSG00000183691 | ENSG00000183691 | 9241 | NOG | 1680.821 | 12.0596659 | 0.0000000 | 10.2208079 | 0.0000000 | 6.9199550 | 0.0000000 | -1.8388580 | 0.0416464 | -5.1397109 | 0.0000000 | -9.1143915 | 0.0000000 | -3.300853 | 0.0002787 | overlap | HIF2A |
| ENSG00000129757 | ENSG00000129757 | 1028 | CDKN1C | 1908.536 | 7.9525101 | 0.0000000 | 6.2939207 | 0.0000000 | 3.4391792 | 0.0000000 | -1.6585894 | 0.0113011 | -4.5133309 | 0.0000000 | -5.1992419 | 0.0000000 | -2.854741 | 0.0000039 | overlap | HIF2A |
| ENSG00000005108 | ENSG00000005108 | 221981 | THSD7A | 1127.535 | -0.0601233 | 0.7462775 | 1.2139511 | 0.0000002 | -1.0431227 | 0.0000135 | 1.2740744 | 0.0001292 | -0.9829994 | 0.0018408 | -1.7097642 | 0.0000000 | -2.257074 | 0.0000000 | HIF1A | HIF2A |
| ENSG00000211459 | ENSG00000211459 | NA | MT-RNR1 | 80906.395 | 0.9532438 | 0.0000000 | 2.1787628 | 0.0000000 | 0.3686447 | 0.0833038 | 1.2255190 | 0.0000053 | -0.5845992 | 0.0330120 | 1.1419769 | 0.0000126 | -1.810118 | 0.0000000 | HIF1A | HIF2A |
| ENSG00000131378 | ENSG00000131378 | 23180 | RFTN1 | 1155.647 | 1.1364445 | 0.0000000 | 2.3475173 | 0.0000000 | -1.0875627 | 0.0000000 | 1.2110728 | 0.0000000 | -2.2240073 | 0.0000000 | -0.9885731 | 0.0000003 | -3.435080 | 0.0000000 | overlap | HIF2A |
| ENSG00000006468 | ENSG00000006468 | 2115 | ETV1 | 2732.487 | 0.1507106 | 0.4038024 | 1.3533783 | 0.0000000 | -2.4061098 | 0.0000000 | 1.2026677 | 0.0004498 | -2.5568204 | 0.0000000 | -1.6571126 | 0.0000002 | -3.759488 | 0.0000000 | overlap | HIF2A |
| ENSG00000168779 | ENSG00000168779 | 6474 | SHOX2 | 2054.783 | -0.7030353 | 0.0000730 | 0.4984821 | 0.0507013 | -1.4331730 | 0.0000000 | 1.2015174 | 0.0008689 | -0.7301377 | 0.0367874 | -1.1677377 | 0.0005326 | -1.931655 | 0.0000001 | HIF1A | HIF2A |
| ENSG00000116147 | ENSG00000116147 | 7143 | TNR | 1919.320 | 2.3355169 | 0.0000000 | 3.5017676 | 0.0000000 | -1.4168155 | 0.0000000 | 1.1662507 | 0.0017831 | -3.7523324 | 0.0000000 | -3.5931850 | 0.0000000 | -4.918583 | 0.0000000 | overlap | HIF2A |
| ENSG00000189120 | ENSG00000189120 | 80320 | SP6 | 1446.799 | 5.9384507 | 0.0000000 | 7.0889482 | 0.0000000 | 1.3178276 | 0.0002727 | 1.1504975 | 0.0134005 | -4.6206231 | 0.0000000 | -1.0976669 | 0.0091774 | -5.771121 | 0.0000000 | overlap | HIF2A |
| ENSG00000173762 | ENSG00000173762 | 924 | CD7 | 1205.426 | 3.3639202 | 0.0000000 | 4.4931330 | 0.0000000 | -0.3388951 | 0.0919844 | 1.1292128 | 0.0000030 | -3.7028153 | 0.0000000 | -2.3588985 | 0.0000000 | -4.832028 | 0.0000000 | overlap | HIF2A |
| ENSG00000158560 | ENSG00000158560 | 1780 | DYNC1I1 | 1329.032 | 1.5942459 | 0.0000000 | 2.6976504 | 0.0000000 | -0.9576680 | 0.0000000 | 1.1034045 | 0.0000002 | -2.5519139 | 0.0000000 | -1.9637207 | 0.0000000 | -3.655318 | 0.0000000 | overlap | HIF2A |
| ENSG00000213626 | ENSG00000213626 | 81606 | LBH | 2362.312 | 2.3320540 | 0.0000000 | 3.3917924 | 0.0000000 | 0.1770539 | 0.6201698 | 1.0597384 | 0.0225608 | -2.1550001 | 0.0000000 | -1.6580070 | 0.0000297 | -3.214739 | 0.0000000 | overlap | HIF2A |
| ENSG00000144730 | ENSG00000144730 | 54756 | IL17RD | 1065.588 | -0.8998280 | 0.0000000 | 0.1199348 | 0.4843366 | -1.6715395 | 0.0000000 | 1.0197628 | 0.0000029 | -0.7717115 | 0.0002246 | -0.9090940 | 0.0000208 | -1.791474 | 0.0000000 | HIF1A | HIF2A |
| ENSG00000139926 | ENSG00000139926 | 122786 | FRMD6 | 1128.222 | 0.8571116 | 0.0000000 | 1.7967241 | 0.0000000 | -1.5002319 | 0.0000000 | 0.9396125 | 0.0000114 | -2.3573435 | 0.0000000 | -2.2036573 | 0.0000000 | -3.296956 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000138759 | ENSG00000138759 | 80144 | FRAS1 | 3292.396 | 1.2960604 | 0.0000000 | 2.2219719 | 0.0000000 | 0.0562012 | 0.8294481 | 0.9259114 | 0.0019453 | -1.2398592 | 0.0000024 | -1.1107657 | 0.0000498 | -2.165771 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000150961 | ENSG00000150961 | 9871 | SEC24D | 2269.986 | 0.4686477 | 0.0000000 | 1.3944519 | 0.0000000 | -1.3942479 | 0.0000000 | 0.9258042 | 0.0000000 | -1.8628956 | 0.0000000 | -1.0057570 | 0.0000000 | -2.788700 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000175274 | ENSG00000175274 | 9537 | TP53I11 | 10799.099 | 1.3861084 | 0.0000000 | 2.3100250 | 0.0000000 | -0.8990063 | 0.0000000 | 0.9239166 | 0.0000002 | -2.2851147 | 0.0000000 | -1.3625428 | 0.0000000 | -3.209031 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000244405 | ENSG00000244405 | 2119 | ETV5 | 3105.739 | 0.8660303 | 0.0000000 | 1.7890970 | 0.0000000 | -0.8689693 | 0.0000002 | 0.9230667 | 0.0000739 | -1.7349996 | 0.0000000 | -1.0656094 | 0.0000013 | -2.658066 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000197748 | ENSG00000197748 | 80217 | CFAP43 | 1894.875 | 7.9785445 | 0.0000000 | 8.9001244 | 0.0000000 | 4.9838285 | 0.0000000 | 0.9215799 | 0.0216112 | -2.9947160 | 0.0000000 | -5.0617567 | 0.0000000 | -3.916296 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000109762 | ENSG00000109762 | 83891 | SNX25 | 1246.634 | 2.7323672 | 0.0000000 | 3.6413772 | 0.0000000 | 1.8286606 | 0.0000000 | 0.9090101 | 0.0000008 | -0.9037066 | 0.0000002 | -1.4010950 | 0.0000000 | -1.812717 | 0.0000000 | interaction | HIF2A |
| ENSG00000163879 | ENSG00000163879 | 7802 | DNALI1 | 1107.798 | 2.9662626 | 0.0000000 | 2.0860670 | 0.0000000 | 0.9199528 | 0.0000000 | -0.8801956 | 0.0000003 | -2.0463099 | 0.0000000 | -2.7984324 | 0.0000000 | -1.166114 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000163430 | ENSG00000163430 | 11167 | FSTL1 | 2058.247 | 1.2094895 | 0.0000000 | 2.0741318 | 0.0000000 | -0.3088552 | 0.0623844 | 0.8646423 | 0.0000593 | -1.5183447 | 0.0000000 | -1.1851329 | 0.0000000 | -2.382987 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000285438 | ENSG00000285438 | 83595 | SOX7 | 1112.198 | 3.5133003 | 0.0000000 | 4.3762767 | 0.0000000 | 0.5739766 | 0.0001808 | 0.8629764 | 0.0000254 | -2.9393237 | 0.0000000 | -3.1070310 | 0.0000000 | -3.802300 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000153132 | ENSG00000153132 | 1047 | CLGN | 1064.645 | 0.1684716 | 0.1941431 | 1.0099403 | 0.0000000 | -0.4494093 | 0.0150769 | 0.8414687 | 0.0008769 | -0.6178810 | 0.0097392 | -1.0396814 | 0.0000097 | -1.459350 | 0.0000000 | interaction | HIF2A |
| ENSG00000169432 | ENSG00000169432 | 6335 | SCN9A | 1024.486 | 0.5068057 | 0.0013202 | 1.3463143 | 0.0000000 | -1.3659325 | 0.0000000 | 0.8395086 | 0.0129171 | -1.8727383 | 0.0000000 | -1.7414107 | 0.0000000 | -2.712247 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000140577 | ENSG00000140577 | 64784 | CRTC3 | 2120.373 | 0.4621049 | 0.0000000 | 1.2967985 | 0.0000000 | -0.3216799 | 0.0018844 | 0.8346935 | 0.0000000 | -0.7837848 | 0.0000000 | -0.5400293 | 0.0000557 | -1.618478 | 0.0000000 | interaction | HIF2A |
| ENSG00000182022 | ENSG00000182022 | 51363 | CHST15 | 4903.285 | -0.5134735 | 0.0000000 | 0.3183952 | 0.0026743 | -1.4530825 | 0.0000000 | 0.8318686 | 0.0000000 | -0.9396090 | 0.0000000 | -0.6977014 | 0.0000007 | -1.771478 | 0.0000000 | interaction | HIF2A |
| ENSG00000141198 | ENSG00000141198 | 10040 | TOM1L1 | 1755.808 | 0.0164868 | 0.8060105 | 0.8399539 | 0.0000000 | -1.0308364 | 0.0000000 | 0.8234671 | 0.0000000 | -1.0473232 | 0.0000000 | -0.2391372 | 0.0436093 | -1.870790 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000137841 | ENSG00000137841 | 5330 | PLCB2 | 1527.128 | 0.3499819 | 0.0006701 | 1.1510506 | 0.0000000 | -0.2981102 | 0.0537325 | 0.8010687 | 0.0000707 | -0.6480921 | 0.0006882 | -0.4504596 | 0.0242378 | -1.449161 | 0.0000000 | interaction | HIF2A |
| ENSG00000166147 | ENSG00000166147 | 2200 | FBN1 | 2361.478 | 0.3887368 | 0.0025154 | 1.1884482 | 0.0000000 | -0.5931496 | 0.0014144 | 0.7997113 | 0.0024338 | -0.9818865 | 0.0000254 | -0.7598103 | 0.0018540 | -1.781598 | 0.0000000 | interaction | HIF2A |
| ENSG00000231764 | ENSG00000231764 | NA | DLX6-AS1 | 1338.823 | 0.7449452 | 0.0000000 | 1.5388567 | 0.0000000 | -0.0446067 | 0.8428687 | 0.7939116 | 0.0019892 | -0.7895519 | 0.0006938 | -1.1771113 | 0.0000005 | -1.583463 | 0.0000000 | interaction | HIF2A |
| ENSG00000092929 | ENSG00000092929 | 201294 | UNC13D | 3057.350 | 1.2824513 | 0.0000000 | 2.0716364 | 0.0000000 | -0.2590653 | 0.2737492 | 0.7891851 | 0.0143393 | -1.5415166 | 0.0000000 | -1.1809810 | 0.0000249 | -2.330702 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000058668 | ENSG00000058668 | 493 | ATP2B4 | 2587.386 | 0.0529318 | 0.5872882 | 0.8415091 | 0.0000000 | -1.6395347 | 0.0000000 | 0.7885773 | 0.0000056 | -1.6924665 | 0.0000000 | -1.3127029 | 0.0000000 | -2.481044 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000035403 | ENSG00000035403 | 7414 | VCL | 3592.997 | -0.7814983 | 0.0000000 | 0.0045406 | 0.9704556 | -1.2099910 | 0.0000000 | 0.7860388 | 0.0000000 | -0.4284927 | 0.0023122 | 0.5257351 | 0.0001913 | -1.214532 | 0.0000000 | interaction | HIF2A |
| ENSG00000185950 | ENSG00000185950 | 8660 | IRS2 | 2124.727 | 1.3011118 | 0.0000000 | 2.0773448 | 0.0000000 | -0.2522202 | 0.0302621 | 0.7762330 | 0.0000001 | -1.5533320 | 0.0000000 | -1.2604589 | 0.0000000 | -2.329565 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000186469 | ENSG00000186469 | 54331 | GNG2 | 8742.006 | 1.8664939 | 0.0000000 | 2.6328978 | 0.0000000 | -0.3208083 | 0.0153783 | 0.7664039 | 0.0000088 | -2.1873022 | 0.0000000 | -2.4061566 | 0.0000000 | -2.953706 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000114923 | ENSG00000114923 | 6508 | SLC4A3 | 1160.962 | 1.4228971 | 0.0000000 | 2.1881050 | 0.0000000 | 0.5959518 | 0.0001322 | 0.7652079 | 0.0004690 | -0.8269453 | 0.0000305 | -0.6072756 | 0.0036018 | -1.592153 | 0.0000000 | interaction | HIF2A |
| ENSG00000013016 | ENSG00000013016 | 30845 | EHD3 | 1157.292 | 0.0123414 | 0.8760821 | 0.7744821 | 0.0000000 | -0.5913717 | 0.0000000 | 0.7621407 | 0.0000000 | -0.6037131 | 0.0000017 | -0.3741164 | 0.0055186 | -1.365854 | 0.0000000 | interaction | HIF2A |
| ENSG00000104490 | ENSG00000104490 | 83988 | NCALD | 3617.354 | 0.5281467 | 0.0000000 | 1.2882649 | 0.0000000 | -1.0664389 | 0.0000000 | 0.7601182 | 0.0000000 | -1.5945856 | 0.0000000 | -1.4282794 | 0.0000000 | -2.354704 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000147862 | ENSG00000147862 | 4781 | NFIB | 1457.401 | -0.3261211 | 0.0296876 | 0.4294591 | 0.0416447 | -1.1216047 | 0.0000001 | 0.7555803 | 0.0183407 | -0.7954836 | 0.0043157 | -0.5293604 | 0.0708781 | -1.551064 | 0.0000002 | interaction | HIF2A |
| ENSG00000154545 | ENSG00000154545 | 728239 | MAGED4 | 6035.260 | 1.9805949 | 0.0000000 | 2.7340568 | 0.0000000 | 0.8906306 | 0.0000000 | 0.7534619 | 0.0012769 | -1.0899642 | 0.0000001 | -1.3635943 | 0.0000000 | -1.843426 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000229802 | ENSG00000229802 | 7922 | SLC39A7 | 4992.755 | 0.0646786 | 0.4266143 | 0.8140413 | 0.0000000 | -0.7840850 | 0.0000000 | 0.7493626 | 0.0000002 | -0.8487637 | 0.0000000 | -0.2853762 | 0.0595508 | -1.598126 | 0.0000000 | interaction | HIF2A |
| ENSG00000117385 | ENSG00000117385 | 64175 | P3H1 | 1609.489 | 0.6107098 | 0.0000000 | 1.3557441 | 0.0000000 | -0.0031382 | 0.9839532 | 0.7450343 | 0.0000014 | -0.6138480 | 0.0000339 | -0.2048790 | 0.2183609 | -1.358882 | 0.0000000 | interaction | HIF2A |
| ENSG00000149256 | ENSG00000149256 | 26011 | TENM4 | 3274.093 | -0.4345411 | 0.0000691 | 0.2931356 | 0.0625716 | -1.2092143 | 0.0000000 | 0.7276767 | 0.0010945 | -0.7746731 | 0.0001170 | -0.2115988 | 0.3584129 | -1.502350 | 0.0000000 | interaction | HIF2A |
| ENSG00000149809 | ENSG00000149809 | 7108 | TM7SF2 | 1468.200 | 1.0570204 | 0.0000000 | 1.7804095 | 0.0000000 | 0.3276628 | 0.1286652 | 0.7233891 | 0.0173843 | -0.7293576 | 0.0062206 | 0.4047281 | 0.1524848 | -1.452747 | 0.0000003 | interaction | HIF2A |
| ENSG00000182158 | ENSG00000182158 | 64764 | CREB3L2 | 7797.381 | 2.6293443 | 0.0000000 | 3.3455855 | 0.0000000 | 0.1776441 | 0.2549535 | 0.7162412 | 0.0002875 | -2.4517002 | 0.0000000 | -1.5635579 | 0.0000000 | -3.167941 | 0.0000000 | HIF2A | HIF2A |
| ENSG00000182195 | ENSG00000182195 | 23641 | LDOC1 | 2567.935 | 0.8986984 | 0.0000000 | 1.6101640 | 0.0000000 | 0.5241309 | 0.0000192 | 0.7114656 | 0.0000205 | -0.3745675 | 0.0253626 | -0.3743381 | 0.0255375 | -1.086033 | 0.0000000 | interaction | HIF2A |
| ENSG00000100360 | ENSG00000100360 | 11020 | IFT27 | 2189.216 | 0.2765802 | 0.0001990 | 0.9834116 | 0.0000000 | -0.2293783 | 0.0392267 | 0.7068314 | 0.0000005 | -0.5059585 | 0.0002233 | -0.1984178 | 0.1910423 | -1.212790 | 0.0000000 | interaction | HIF2A |
| ENSG00000122507 | ENSG00000122507 | 27241 | BBS9 | 1227.982 | 1.0983279 | 0.0000000 | 1.8029782 | 0.0000000 | 0.4962223 | 0.0009997 | 0.7046503 | 0.0008338 | -0.6021056 | 0.0020723 | -1.1189176 | 0.0000000 | -1.306756 | 0.0000000 | interaction | HIF2A |
| ENSG00000102316 | ENSG00000102316 | 10916 | MAGED2 | 12464.272 | 0.7649346 | 0.0000000 | 1.4665934 | 0.0000000 | 0.0163088 | 0.9367938 | 0.7016588 | 0.0026927 | -0.7486259 | 0.0003533 | -0.3097901 | 0.1781369 | -1.450285 | 0.0000000 | interaction | HIF2A |

``` r
plotCounts_SK(HIF2A_genes[1:9,] %>% rownames()) + labs(title="Hif2A")
```

![](Readme_files/figure-gfm/gene_lists-2.png)<!-- -->

``` r
HIF1A_HIF2A_genes <- res_hif1a_2a_p %>% 
  .[order(abs(.$Hif1aHxNx.vs.KellyHxNx.log2FoldChange), decreasing = TRUE),] %>% 
  filter(group =="HIF1A_HIF2A") %>%
  filter(baseMean > 1000)
nrow(HIF1A_HIF2A_genes)
```

    ## [1] 9

``` r
HIF1A_HIF2A_genes[1:9,] %>% kable()
```

|  | ENSEMBL | ENTREZ | symbol | baseMean | Kelly.Hx.vs.Nx.log2FoldChange | Kelly.Hx.vs.Nx.padj | Hif1a.Hx.vs.Nx.log2FoldChange | Hif1a.Hx.vs.Nx.padj | Hif2a.Hx.vs.Nx.log2FoldChange | Hif2a.Hx.vs.Nx.padj | Hif1aHxNx.vs.KellyHxNx.log2FoldChange | Hif1aHxNx.vs.KellyHxNx.padj | Hif2aHxNx.vs.KellyHxNx.log2FoldChange | Hif2aHxNx.vs.KellyHxNx.padj | Hif1bHxNx.vs.KellyHxNx.log2FoldChange | Hif1bHxNx.vs.KellyHxNx.padj | Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange | Hif2aHxNx.vs.Hif1aHxNx.padj | venn | group |
|:---|:---|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---|
| ENSG00000128285 | ENSG00000128285 | 2847 | MCHR1 | 1441.330 | 6.0238259 | 0.0e+00 | 4.7856609 | 0.0000000 | 4.0754722 | 0.0000000 | -1.2381651 | 0.0013832 | -1.9483537 | 0.0000000 | -2.7343326 | 0.0000000 | -0.7101887 | 0.0971408 | overlap | HIF1A_HIF2A |
| ENSG00000102924 | ENSG00000102924 | 869 | CBLN1 | 1259.609 | -0.3875171 | 1.8e-06 | 0.7434115 | 0.0000000 | 0.0977362 | 0.4597855 | 1.1309286 | 0.0000000 | 0.4852533 | 0.0014055 | 0.8005722 | 0.0000001 | -0.6456753 | 0.0001035 | HIF1A | HIF1A_HIF2A |
| ENSG00000109501 | ENSG00000109501 | 7466 | WFS1 | 1084.739 | 0.8072203 | 0.0e+00 | -0.2938728 | 0.0025945 | -0.4844818 | 0.0000010 | -1.1010931 | 0.0000000 | -1.2917021 | 0.0000000 | -1.3442817 | 0.0000000 | -0.1906090 | 0.2369341 | overlap | HIF1A_HIF2A |
| ENSG00000100033 | ENSG00000100033 | 5625 | PRODH | 1022.754 | 5.2475971 | 0.0e+00 | 4.2839013 | 0.0000000 | 3.8180744 | 0.0000000 | -0.9636958 | 0.0305031 | -1.4295227 | 0.0001199 | -1.3206168 | 0.0007102 | -0.4658269 | 0.3269849 | HIF2A | HIF1A_HIF2A |
| ENSG00000145911 | ENSG00000145911 | 23138 | N4BP3 | 1053.447 | 3.9025379 | 0.0e+00 | 3.0312547 | 0.0000000 | 2.8566776 | 0.0000000 | -0.8712832 | 0.0008396 | -1.0458603 | 0.0000070 | -3.3670552 | 0.0000000 | -0.1745771 | 0.5939696 | HIF2A | HIF1A_HIF2A |
| ENSG00000197381 | ENSG00000197381 | 104 | ADARB1 | 6359.279 | 3.1868643 | 0.0e+00 | 2.3172845 | 0.0000000 | 2.0958326 | 0.0000000 | -0.8695798 | 0.0000000 | -1.0910317 | 0.0000000 | -2.4712562 | 0.0000000 | -0.2214519 | 0.2512870 | HIF2A | HIF1A_HIF2A |
| ENSG00000165973 | ENSG00000165973 | 4745 | NELL1 | 1629.032 | -2.4131601 | 0.0e+00 | -1.6884522 | 0.0000000 | -1.4029088 | 0.0000000 | 0.7247078 | 0.0119803 | 1.0102512 | 0.0000340 | 0.6286834 | 0.0154136 | 0.2855434 | 0.3719025 | HIF2A | HIF1A_HIF2A |
| ENSG00000102996 | ENSG00000102996 | 4324 | MMP15 | 1266.560 | -2.1593740 | 0.0e+00 | -1.4956359 | 0.0000000 | -0.7168216 | 0.0000014 | 0.6637381 | 0.0023198 | 1.4425525 | 0.0000000 | 1.7011626 | 0.0000000 | 0.7788143 | 0.0003390 | HIF2A | HIF1A_HIF2A |
| ENSG00000182771 | ENSG00000182771 | 2894 | GRID1 | 2267.752 | -1.2525679 | 0.0e+00 | -0.8969890 | 0.0000000 | -0.0235387 | 0.8572363 | 0.3555789 | 0.0284194 | 1.2290292 | 0.0000000 | 1.3392019 | 0.0000000 | 0.8734503 | 0.0000000 | HIF2A | HIF1A_HIF2A |

``` r
plotCounts_SK(HIF1A_HIF2A_genes[1:9,] %>% rownames()) + labs(title="HIF1A+HIF2A")
```

![](Readme_files/figure-gfm/gene_lists-3.png)<!-- -->

# Figure 4: Gene Set enrichment

## GO Analysis

``` r
res_hif1a_2a %>% filter(group== "HIF1A") %>% head(n=20) %>% kable()
```

|  | ENSEMBL | ENTREZ | symbol | baseMean | Kelly.Hx.vs.Nx.log2FoldChange | Kelly.Hx.vs.Nx.padj | Hif1a.Hx.vs.Nx.log2FoldChange | Hif1a.Hx.vs.Nx.padj | Hif2a.Hx.vs.Nx.log2FoldChange | Hif2a.Hx.vs.Nx.padj | Hif1aHxNx.vs.KellyHxNx.log2FoldChange | Hif1aHxNx.vs.KellyHxNx.padj | Hif2aHxNx.vs.KellyHxNx.log2FoldChange | Hif2aHxNx.vs.KellyHxNx.padj | Hif1bHxNx.vs.KellyHxNx.log2FoldChange | Hif1bHxNx.vs.KellyHxNx.padj | Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange | Hif2aHxNx.vs.Hif1aHxNx.padj | venn | group |
|:---|:---|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---|
| ENSG00000170525 | ENSG00000170525 | 5209 | PFKFB3 | 8312.9632 | 5.8234937 | 0.0000000 | 0.2920504 | 0.2227293 | 6.4294748 | 0.0000000 | -5.531443 | 0 | 0.6059811 | 0.0560612 | 1.6216771 | 0.0000000 | 6.137424 | 0.0000000 | HIF1A | HIF1A |
| ENSG00000107159 | ENSG00000107159 | 768 | CA9 | 3017.8336 | 10.7042616 | 0.0000000 | 3.9105760 | 0.0000000 | 11.4553409 | 0.0000000 | -6.793686 | 0 | 0.7510793 | 0.3494878 | -1.2066128 | 0.0957458 | 7.544765 | 0.0000000 | HIF1A | HIF1A |
| ENSG00000114023 | ENSG00000114023 | 26355 | FAM162A | 12808.1207 | 2.6749198 | 0.0000000 | -1.0871282 | 0.0000000 | 4.1525578 | 0.0000000 | -3.762048 | 0 | 1.4776380 | 0.0000000 | 0.4415823 | 0.0460083 | 5.239686 | 0.0000000 | overlap | HIF1A |
| ENSG00000176171 | ENSG00000176171 | 664 | BNIP3 | 21434.7075 | 3.1965388 | 0.0000000 | -0.3905581 | 0.0021865 | 4.0130817 | 0.0000000 | -3.587097 | 0 | 0.8165428 | 0.0000005 | -0.0984268 | 0.6234311 | 4.403640 | 0.0000000 | HIF1A | HIF1A |
| ENSG00000228709 | ENSG00000228709 | NA | LINC02575 | 8376.0391 | 9.9795584 | 0.0000000 | 3.7018407 | 0.0000009 | 9.7110362 | 0.0000000 | -6.277718 | 0 | -0.2685222 | 0.7878092 | -2.8390173 | 0.0014617 | 6.009195 | 0.0000000 | HIF1A | HIF1A |
| ENSG00000074800 | ENSG00000074800 | 2023 | ENO1 | 193555.7703 | 1.2878157 | 0.0000000 | -1.6300196 | 0.0000000 | 2.3330043 | 0.0000000 | -2.917835 | 0 | 1.0451886 | 0.0000209 | 1.0149729 | 0.0000614 | 3.963024 | 0.0000000 | overlap | HIF1A |
| ENSG00000159208 | ENSG00000159208 | 148523 | CIART | 2103.7130 | 2.8415121 | 0.0000000 | -0.0961627 | 0.5848455 | 3.8535705 | 0.0000000 | -2.937675 | 0 | 1.0120583 | 0.0000005 | -0.0590673 | 0.8220711 | 3.949733 | 0.0000000 | overlap | HIF1A |
| ENSG00000163536 | ENSG00000163536 | 5274 | SERPINI1 | 881.7979 | 0.4841107 | 0.0208196 | 4.6362397 | 0.0000000 | -0.0677924 | 0.8578828 | 4.152129 | 0 | -0.5519030 | 0.2093907 | -1.6422193 | 0.0000270 | -4.704032 | 0.0000000 | HIF1A | HIF1A |
| ENSG00000101204 | ENSG00000101204 | 1137 | CHRNA4 | 1808.3600 | 11.2833363 | 0.0000000 | 3.8747023 | 0.0000279 | 12.2045196 | 0.0000000 | -7.408634 | 0 | 0.9211833 | 0.5222658 | -2.5145505 | 0.0809627 | 8.329817 | 0.0000000 | HIF1A | HIF1A |
| ENSG00000102144 | ENSG00000102144 | 5230 | PGK1 | 69853.4717 | 2.3595632 | 0.0000000 | 0.1341692 | 0.4241579 | 2.9941319 | 0.0000000 | -2.225394 | 0 | 0.6345687 | 0.0025437 | 0.1551542 | 0.5212432 | 2.859963 | 0.0000000 | HIF1A | HIF1A |
| ENSG00000079739 | ENSG00000079739 | 5236 | PGM1 | 6602.8118 | 1.2982488 | 0.0000000 | -1.1518298 | 0.0000000 | 2.2307331 | 0.0000000 | -2.450079 | 0 | 0.9324843 | 0.0003620 | 0.4075537 | 0.1536669 | 3.382563 | 0.0000000 | HIF1A | HIF1A |
| ENSG00000163516 | ENSG00000163516 | 55139 | ANKZF1 | 4282.9780 | 3.0398887 | 0.0000000 | 0.6596566 | 0.0000669 | 2.8478554 | 0.0000000 | -2.380232 | 0 | -0.1920333 | 0.4795841 | -1.6634039 | 0.0000000 | 2.188199 | 0.0000000 | HIF1A | HIF1A |
| ENSG00000100314 | ENSG00000100314 | 164633 | CABP7 | 2092.2854 | 8.1923082 | 0.0000000 | 4.9574998 | 0.0000000 | 9.3811748 | 0.0000000 | -3.234808 | 0 | 1.1888666 | 0.0060813 | -1.2450997 | 0.0047254 | 4.423675 | 0.0000000 | overlap | HIF1A |
| ENSG00000165802 | ENSG00000165802 | 26012 | NSMF | 12663.1695 | 0.8783827 | 0.0000000 | -1.2842273 | 0.0000000 | 2.3874643 | 0.0000000 | -2.162610 | 0 | 1.5090817 | 0.0000000 | -0.3505130 | 0.1519806 | 3.671692 | 0.0000000 | overlap | HIF1A |
| ENSG00000123095 | ENSG00000123095 | 79365 | BHLHE41 | 591.4573 | 4.7658766 | 0.0000000 | 1.0268082 | 0.0005643 | 7.4469315 | 0.0000000 | -3.739069 | 0 | 2.6810549 | 0.0000000 | -1.8971434 | 0.0000012 | 6.420123 | 0.0000000 | overlap | HIF1A |
| ENSG00000186352 | ENSG00000186352 | 353322 | ANKRD37 | 742.5260 | 2.6555190 | 0.0000000 | 5.6608158 | 0.0000000 | 2.9944297 | 0.0000000 | 3.005297 | 0 | 0.3389107 | 0.2002079 | 0.0346447 | 0.9081232 | -2.666386 | 0.0000000 | HIF1A | HIF1A |
| ENSG00000085563 | ENSG00000085563 | 5243 | ABCB1 | 1046.0696 | -1.9274408 | 0.0000000 | 0.9676360 | 0.0002242 | -2.0925485 | 0.0000000 | 2.895077 | 0 | -0.1651077 | 0.7173453 | 0.3719400 | 0.3326388 | -3.060185 | 0.0000000 | HIF1A | HIF1A |
| ENSG00000261992 | ENSG00000261992 | 57459 | GATAD2B | 453.5022 | 0.2414106 | 0.9269869 | -26.5376503 | 0.0000000 | 20.4094194 | 0.0000000 | -26.779061 | 0 | 20.1680088 | 0.0000026 | -0.1024105 | 0.9860468 | 46.947070 | 0.0000000 | overlap | HIF1A |
| ENSG00000114480 | ENSG00000114480 | 2632 | GBE1 | 3845.3543 | 2.9488275 | 0.0000000 | 0.9458202 | 0.0000000 | 3.8232133 | 0.0000000 | -2.003007 | 0 | 0.8743858 | 0.0000495 | -1.3279561 | 0.0000000 | 2.877393 | 0.0000000 | HIF1A | HIF1A |
| ENSG00000278718 | ENSG00000278718 | 79792 | GSDMD | 364.5132 | 9.9927998 | 0.0000000 | 5.1654055 | 0.0000000 | 7.3256811 | 0.0000000 | -4.827394 | 0 | -2.6671188 | 0.0001387 | -9.0627826 | 0.0000000 | 2.160276 | 0.0049857 | overlap | HIF1A |

``` r
res_hif1a_2a_list_ens <- list("HIF1A" = res_hif1a_2a %>% filter(group == "HIF1A") %>% .[,"ENSEMBL"],
                              "HIF2A" = res_hif1a_2a %>% filter(group == "HIF2A") %>% .[,"ENSEMBL"],
                              "both" = res_hif1a_2a %>% filter(group == "HIF1A_HIF2A") %>% .[,"ENSEMBL"])

res_hif1a_2a_list_ez <- list("HIF1A" = res_hif1a_2a %>% filter(group == "HIF1A") %>% .[,"ENTREZ"],
                              "HIF2A" = res_hif1a_2a %>% filter(group == "HIF2A") %>% .[,"ENTREZ"],
                              "both" = res_hif1a_2a %>% filter(group == "HIF1A_HIF2A") %>% .[,"ENTREZ"])

load(file="GO_analysis/GO_cc_groups.go")

GO_cc_groups_BP <- GO_cc_groups %>% filter(ONTOLOGY=="BP")
dotplot(GO_cc_groups_BP, showCategory=12)
```

<img src="Readme_files/figure-gfm/unnamed-chunk-2-1.png" width="50%" />

## Cluster GO terms

``` r
# Search for clusters
GO_IDs_list <- split(GO_cc_groups_BP@compareClusterResult,f=GO_cc_groups_BP@compareClusterResult$Cluster) %>% lapply('[',,"ID")
names(GO_IDs_list)
```

    ## [1] "HIF1A" "HIF2A" "both"

``` r
# simplifyGOFromMultipleLists(GO_IDs_list[1:2])

simplifyGO(GO_IDs_list[[1]], column_title = paste0("HIF1A (",length(GO_IDs_list[[1]])," GO terms)"))
```

![](Readme_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
simplifyGO(GO_IDs_list[[2]], column_title = paste0("HIF2A (",length(GO_IDs_list[[2]])," GO terms)"))
```

![](Readme_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

## KEGG

``` r
cc_kegg <- compareCluster(geneCluster = res_hif1a_2a_list_ez[1:2],
                          fun = "enrichKEGG",
                          organism     = 'hsa',
                          pvalueCutoff = 0.05)

dotplot(cc_kegg, showCategory=12)
```

<img src="Readme_files/figure-gfm/unnamed-chunk-4-1.png" width="50%" />

``` r
cc_kegg %>% data.frame()
```

    ##    Cluster                             category
    ## 1    HIF1A       Genetic Information Processing
    ## 2    HIF1A       Genetic Information Processing
    ## 3    HIF1A                   Cellular Processes
    ## 4    HIF1A       Genetic Information Processing
    ## 5    HIF1A       Genetic Information Processing
    ## 6    HIF1A       Genetic Information Processing
    ## 7    HIF1A                           Metabolism
    ## 8    HIF1A                           Metabolism
    ## 9    HIF1A                           Metabolism
    ## 10   HIF1A                           Metabolism
    ## 11   HIF1A       Genetic Information Processing
    ## 12   HIF1A                           Metabolism
    ## 13   HIF1A                           Metabolism
    ## 14   HIF1A                           Metabolism
    ## 15   HIF1A                   Cellular Processes
    ## 16   HIF1A                           Metabolism
    ## 17   HIF1A                           Metabolism
    ## 18   HIF1A                           Metabolism
    ## 19   HIF2A                                 <NA>
    ## 20   HIF2A Environmental Information Processing
    ## 21   HIF2A Environmental Information Processing
    ## 22   HIF2A Environmental Information Processing
    ## 23   HIF2A                   Cellular Processes
    ## 24   HIF2A Environmental Information Processing
    ## 25   HIF2A Environmental Information Processing
    ## 26   HIF2A                           Metabolism
    ## 27   HIF2A                   Cellular Processes
    ## 28   HIF2A Environmental Information Processing
    ## 29   HIF2A                   Organismal Systems
    ## 30   HIF2A                       Human Diseases
    ## 31   HIF2A                   Organismal Systems
    ## 32   HIF2A Environmental Information Processing
    ## 33   HIF2A                   Organismal Systems
    ## 34   HIF2A                   Organismal Systems
    ## 35   HIF2A                   Organismal Systems
    ## 36   HIF2A                   Organismal Systems
    ## 37   HIF2A                       Human Diseases
    ## 38   HIF2A Environmental Information Processing
    ##                             subcategory       ID
    ## 1                Replication and repair hsa03030
    ## 2                Replication and repair hsa03440
    ## 3                 Cell growth and death hsa04110
    ## 4                           Translation hsa03008
    ## 5                Replication and repair hsa03460
    ## 6                Replication and repair hsa03410
    ## 7              Global and overview maps hsa01232
    ## 8              Global and overview maps hsa01230
    ## 9                 Nucleotide metabolism hsa00230
    ## 10              Carbohydrate metabolism hsa00030
    ## 11               Replication and repair hsa03430
    ## 12      Metabolism of other amino acids hsa00450
    ## 13             Global and overview maps hsa01200
    ## 14 Metabolism of cofactors and vitamins hsa00670
    ## 15                Cell growth and death hsa04115
    ## 16                Nucleotide metabolism hsa00240
    ## 17              Carbohydrate metabolism hsa00051
    ## 18              Carbohydrate metabolism hsa00010
    ## 19                                 <NA> hsa04820
    ## 20                  Signal transduction hsa04151
    ## 21  Signaling molecules and interaction hsa04512
    ## 22  Signaling molecules and interaction hsa04080
    ## 23             Transport and catabolism hsa04148
    ## 24                  Signal transduction hsa04020
    ## 25  Signaling molecules and interaction hsa04514
    ## 26   Glycan biosynthesis and metabolism hsa00534
    ## 27      Cellular community - eukaryotes hsa04510
    ## 28                  Signal transduction hsa04010
    ## 29                   Circulatory system hsa04270
    ## 30                 Substance dependence hsa05034
    ## 31                       Nervous system hsa04727
    ## 32                  Signal transduction hsa04014
    ## 33                     Digestive system hsa04974
    ## 34         Development and regeneration hsa04360
    ## 35                     Endocrine system hsa04926
    ## 36                       Nervous system hsa04725
    ## 37      Endocrine and metabolic disease hsa04933
    ## 38                  Signal transduction hsa04072
    ##                                                   Description GeneRatio
    ## 1                                             DNA replication    19/612
    ## 2                                    Homologous recombination    15/612
    ## 3                                                  Cell cycle    29/612
    ## 4                           Ribosome biogenesis in eukaryotes    30/612
    ## 5                                      Fanconi anemia pathway    15/612
    ## 6                                        Base excision repair    12/612
    ## 7                                       Nucleotide metabolism    17/612
    ## 8                                 Biosynthesis of amino acids    15/612
    ## 9                                           Purine metabolism    21/612
    ## 10                                  Pentose phosphate pathway     9/612
    ## 11                                            Mismatch repair     7/612
    ## 12                                  Selenocompound metabolism     6/612
    ## 13                                          Carbon metabolism    18/612
    ## 14                                  One carbon pool by folate     9/612
    ## 15                                      p53 signaling pathway    13/612
    ## 16                                      Pyrimidine metabolism    11/612
    ## 17                            Fructose and mannose metabolism     8/612
    ## 18                               Glycolysis / Gluconeogenesis    12/612
    ## 19                               Cytoskeleton in muscle cells    49/830
    ## 20                                 PI3K-Akt signaling pathway    65/830
    ## 21                                   ECM-receptor interaction    24/830
    ## 22                    Neuroactive ligand-receptor interaction    62/830
    ## 23                                              Efferocytosis    33/830
    ## 24                                  Calcium signaling pathway    44/830
    ## 25                                    Cell adhesion molecules    30/830
    ## 26 Glycosaminoglycan biosynthesis - heparan sulfate / heparin     9/830
    ## 27                                             Focal adhesion    35/830
    ## 28                                     MAPK signaling pathway    46/830
    ## 29                         Vascular smooth muscle contraction    25/830
    ## 30                                                 Alcoholism    32/830
    ## 31                                          GABAergic synapse    18/830
    ## 32                                      Ras signaling pathway    37/830
    ## 33                           Protein digestion and absorption    20/830
    ## 34                                              Axon guidance    30/830
    ## 35                                  Relaxin signaling pathway    23/830
    ## 36                                        Cholinergic synapse    21/830
    ## 37       AGE-RAGE signaling pathway in diabetic complications    19/830
    ## 38                          Phospholipase D signaling pathway    25/830
    ##     BgRatio RichFactor FoldEnrichment    zScore       pvalue     p.adjust
    ## 1   36/8865  0.5277778       7.645016 10.878680 1.892361e-13 6.206945e-11
    ## 2   41/8865  0.3658537       5.299498  7.513834 3.823250e-08 6.270131e-06
    ## 3  158/8865  0.1835443       2.658693  5.728554 9.747066e-07 9.520142e-05
    ## 4  168/8865  0.1785714       2.586660  5.653773 1.160993e-06 9.520142e-05
    ## 5   55/8865  0.2727273       3.950535  5.976939 2.890910e-06 1.896437e-04
    ## 6   44/8865  0.2727273       3.950535  5.342602 2.857573e-05 1.562140e-03
    ## 7   85/8865  0.2000000       2.897059  4.785504 5.615534e-05 2.631279e-03
    ## 8   75/8865  0.2000000       2.897059  4.492640 1.521483e-04 5.595842e-03
    ## 9  128/8865  0.1640625       2.376494  4.271522 1.651894e-04 5.595842e-03
    ## 10  31/8865  0.2903226       4.205408  4.868219 1.706049e-04 5.595842e-03
    ## 11  23/8865  0.3043478       4.408568  4.457027 6.668430e-04 1.849149e-02
    ## 12  17/8865  0.3529412       5.112457  4.621553 6.765179e-04 1.849149e-02
    ## 13 116/8865  0.1551724       2.247718  3.683417 9.383859e-04 2.367620e-02
    ## 14  39/8865  0.2307692       3.342760  3.992665 1.090672e-03 2.555290e-02
    ## 15  75/8865  0.1733333       2.510784  3.577860 1.689950e-03 3.435644e-02
    ## 16  58/8865  0.1896552       2.747211  3.635214 1.789520e-03 3.435644e-02
    ## 17  34/8865  0.2352941       3.408304  3.831163 1.794860e-03 3.435644e-02
    ## 18  67/8865  0.1791045       2.594381  3.567157 1.885414e-03 3.435644e-02
    ## 19 232/8865  0.2112069       2.255842  6.229586 3.311671e-08 1.089540e-05
    ## 20 362/8865  0.1795580       1.917809  5.730357 1.496781e-07 2.462204e-05
    ## 21  89/8865  0.2696629       2.880195  5.729419 1.290754e-06 1.415526e-04
    ## 22 368/8865  0.1684783       1.799470  5.034473 2.781391e-06 2.287694e-04
    ## 23 157/8865  0.2101911       2.244993  5.058461 6.577286e-06 4.327854e-04
    ## 24 254/8865  0.1732283       1.850204  4.418493 3.987487e-05 2.186472e-03
    ## 25 158/8865  0.1898734       2.027985  4.190272 1.263612e-04 5.938976e-03
    ## 26  24/8865  0.3750000       4.005271  4.738049 1.896635e-04 7.799910e-03
    ## 27 203/8865  0.1724138       1.841504  3.898128 2.645408e-04 9.670436e-03
    ## 28 300/8865  0.1533333       1.637711  3.611451 5.130278e-04 1.651619e-02
    ## 29 134/8865  0.1865672       1.992672  3.721236 5.924418e-04 1.651619e-02
    ## 30 188/8865  0.1702128       1.817995  3.643392 6.024144e-04 1.651619e-02
    ## 31  89/8865  0.2022472       2.160146  3.535252 1.294367e-03 3.206733e-02
    ## 32 238/8865  0.1554622       1.660449  3.319392 1.364567e-03 3.206733e-02
    ## 33 105/8865  0.1904762       2.034423  3.426904 1.556950e-03 3.414911e-02
    ## 34 184/8865  0.1630435       1.741422  3.266259 1.783445e-03 3.667209e-02
    ## 35 130/8865  0.1769231       1.889666  3.284191 2.018253e-03 3.905912e-02
    ## 36 116/8865  0.1810345       1.933579  3.252834 2.324701e-03 4.057833e-02
    ## 37 101/8865  0.1881188       2.009245  3.278439 2.343429e-03 4.057833e-02
    ## 38 149/8865  0.1677852       1.792068  3.133702 2.802279e-03 4.609748e-02
    ##          qvalue
    ## 1  5.716923e-11
    ## 2  5.775120e-06
    ## 3  8.768551e-05
    ## 4  8.768551e-05
    ## 5  1.746718e-04
    ## 6  1.438813e-03
    ## 7  2.423546e-03
    ## 8  5.154065e-03
    ## 9  5.154065e-03
    ## 10 5.154065e-03
    ## 11 1.703163e-02
    ## 12 1.703163e-02
    ## 13 2.180703e-02
    ## 14 2.353556e-02
    ## 15 3.164409e-02
    ## 16 3.164409e-02
    ## 17 3.164409e-02
    ## 18 3.164409e-02
    ## 19 9.586415e-06
    ## 20 2.166393e-05
    ## 21 1.245464e-04
    ## 22 2.012849e-04
    ## 23 3.807902e-04
    ## 24 1.923787e-03
    ## 25 5.225463e-03
    ## 26 6.862823e-03
    ## 27 8.508622e-03
    ## 28 1.453193e-02
    ## 29 1.453193e-02
    ## 30 1.453193e-02
    ## 31 2.821473e-02
    ## 32 2.821473e-02
    ## 33 3.004641e-02
    ## 34 3.226628e-02
    ## 35 3.436653e-02
    ## 36 3.570322e-02
    ## 37 3.570322e-02
    ## 38 4.055930e-02
    ##                                                                                                                                                                                                                                                                                                                                              geneID
    ## 1                                                                                                                                                                                                                                               5424/4173/4171/79621/4174/4175/2237/23649/5427/5422/4172/10714/5425/5983/54107/4176/5111/1763/10535
    ## 2                                                                                                                                                                                                                                                                       79184/5424/672/641/7516/675/8438/11073/4361/5889/10714/5425/83990/580/79728
    ## 3                                                                                                                                                                                           9088/4173/8318/4171/4174/4175/5591/9319/90381/990/993/1021/64682/10912/4085/23594/4172/1663/7042/4998/113130/11200/1026/157570/4176/5111/81620/545/9134
    ## 4                                                                                                                                                             54913/1736/5822/65083/10799/10171/92856/84135/84916/55651/29107/55226/55127/10885/23160/2091/26354/6949/55272/55341/10528/54552/51096/166378/54433/28987/134430/23195/51602/102157402
    ## 5                                                                                                                                                                                                                                                                     2175/672/641/675/55215/2177/5889/83990/79728/545/2189/29089/116028/2187/57697
    ## 6                                                                                                                                                                                                                                                                                  5424/2237/4913/5427/252969/10714/5425/5983/3146/54107/55775/5111
    ## 7                                                                                                                                                                                                                                                      3251/204/377841/84618/4830/79077/7083/654364/7298/132/159/100/6241/1503/4860/102157402/51727
    ## 8                                                                                                                                                                                                                                                                    2023/5230/5223/84706/230/8277/441531/22934/5634/5091/29968/5832/5214/4144/6120
    ## 9                                                                                                                                                                                                                                     5236/55276/9060/3251/204/377841/84618/114/5198/4830/5634/471/654364/132/159/5141/100/6241/2618/4860/102157402
    ## 10                                                                                                                                                                                                                                                                                                    5236/55276/230/2821/8277/22934/5634/5214/6120
    ## 11                                                                                                                                                                                                                                                                                                              5424/9156/10714/5425/5983/4436/5111
    ## 12                                                                                                                                                                                                                                                                                                                   9060/92935/51540/7296/4141/883
    ## 13                                                                                                                                                                                                                                                    2023/5230/5223/84706/3098/230/2821/8277/441531/22934/5634/5091/29968/2731/5214/2653/2746/6120
    ## 14                                                                                                                                                                                                                                                                                                     4522/2731/471/7298/1719/2653/10797/4144/2618
    ## 15                                                                                                                                                                                                                                                                                5054/55240/5366/1021/10912/6241/11200/1026/55367/545/637/9134/841
    ## 16                                                                                                                                                                                                                                                                                    377841/84618/4830/1723/79077/7083/654364/7298/6241/1503/51727
    ## 17                                                                                                                                                                                                                                                                                                           5209/3098/230/5207/5214/6652/5210/2762
    ## 18                                                                                                                                                                                                                                                                                     2023/5230/5236/55276/5223/3098/230/2821/441531/3939/5214/219
    ## 19                                                                                      1832/25802/3673/4607/22801/6385/1277/3672/3680/27295/1292/1289/1291/8516/7058/2335/1301/6383/22795/3339/22989/4703/288/256076/23500/8736/1293/1290/1281/1298/3696/1634/825/9172/2200/8910/7414/752/633/2318/3728/1462/5318/23345/3693/57644/482/58529/51332
    ## 20 2791/7143/3791/4915/84699/80310/3673/64764/54331/2056/5295/5979/5159/22801/2321/5156/3570/1277/3672/3680/1292/1291/2668/285/8516/3913/7058/5155/9586/54541/2261/6446/1902/2335/3481/3667/56034/94235/118788/256076/8074/1293/2792/2252/4908/1298/3696/2788/57121/5618/2066/1441/26281/7533/5578/842/2064/9170/2323/3693/4515/3918/4804/3718/7448
    ## 21                                                                                                                                                                                                                   7143/3673/22801/6385/1277/3672/3680/1292/1291/8516/3913/7058/22987/2335/3339/80144/158326/256076/1293/1298/3696/3693/3918/7448
    ## 22                         5732/2741/1901/3061/2911/1906/2837/6750/1910/1135/5028/150/4852/147/3360/56923/4923/11255/2834/8862/5031/4922/7425/1902/2692/1392/135/7432/9248/4986/6751/1268/9568/2696/8973/1903/2895/2914/5733/5737/552/2570/59350/116443/5618/1081/2642/116/2900/9127/1140/5025/9170/66004/2832/6870/6915/551/4887/344838/2568/10203
    ## 23                                                                                                                                                                5732/1901/2056/100133941/5627/4953/9844/10461/6566/5468/1374/80824/5031/1846/6446/6376/51761/51454/177/1844/19/5175/89790/6513/343702/8398/5600/9261/10396/10062/3693/121601/4772
    ## 24                                                                                                                      5137/3791/4915/2911/80310/491/1910/5979/5159/2321/5156/493/147/3360/4923/2668/5155/2261/56034/135/8913/8074/778/3707/2252/5733/5737/552/4485/116443/3270/2066/26281/84812/5330/9127/5332/5578/5025/2064/5333/6870/6915/4772
    ## 25                                                                                                                                                                                84189/1002/54413/23705/6385/9369/3680/5789/5818/5797/8516/9378/6383/23114/1364/114798/5010/3683/152404/149461/3696/920/5175/6404/6693/90952/84628/1462/4685/29126
    ## 26                                                                                                                                                                                                                                                                                               9957/9348/2131/266722/26035/222537/64131/2134/9955
    ## 27                                                                                                                                                            330/7143/3791/56924/80310/3673/5295/5159/22801/2321/5156/1277/3672/3680/1292/1291/8516/3913/7058/5155/399694/2335/56034/256076/1293/1298/3696/7414/5829/5578/5602/2064/3693/3918/7448
    ## 28                                                                                                   7048/3791/4915/80310/4208/5979/5159/2321/5156/2668/285/80824/5971/4217/5155/59283/1846/2261/3481/56034/8913/283748/11184/27330/8074/778/11221/2252/4908/27091/1844/4616/2066/26281/8605/5578/5602/1848/5600/2064/9261/3305/2323/4772/4804/5801
    ## 29                                                                                                                                                                                                                        111/3778/59/1906/147/135/5583/4879/283748/2982/778/552/4880/8605/5330/5332/5578/5581/8398/2983/4637/10268/50487/551/10203
    ## 30                                                                                                                                                                          111/2791/4915/84699/64764/54331/4852/8353/8347/9586/399694/1392/94235/135/3017/8329/8357/8344/4129/554313/8365/2792/8343/8339/2788/8360/8367/116443/8348/8350/8358/8336
    ## 31                                                                                                                                                                                                                                                       111/2791/54331/18/3763/92745/94235/6540/6539/9568/778/2792/9001/140679/2788/2570/5578/2568
    ## 32                                                                                                                                             5337/2791/3791/56924/4915/80310/54331/7074/5295/5159/2321/5156/285/23179/5155/399694/2261/3481/56034/94235/4303/283748/2113/8074/2792/2252/4908/2788/64926/26281/8605/5578/5602/8398/2323/50487/4804
    ## 33                                                                                                                                                                                                                                         206358/1308/1277/1292/1289/1291/80781/1301/91522/1310/1303/1295/256076/7512/1293/1290/1281/1298/5547/482
    ## 34                                                                                                                                                                          22885/56924/7474/2051/5362/91584/5295/9037/10512/223117/6092/2043/2046/7224/54437/56896/9423/54361/2044/5578/1949/84628/2049/10507/137970/10154/56963/151449/64221/5365
    ## 35                                                                                                                                                                                                                           7048/111/2791/84699/59/1906/64764/54331/5295/1910/1277/9586/399694/94235/2792/1281/2788/59350/5330/5332/5578/5602/5600
    ## 36                                                                                                                                                                                                                                         111/2791/84699/64764/54331/5295/9132/43/3763/9586/94235/778/8973/2792/2788/1103/5330/5332/5578/3768/6572
    ## 37                                                                                                                                                                                                                                                   7048/2308/1906/5295/1277/2335/7056/1281/177/1958/84812/5330/5332/5578/5581/5602/5600/5333/4772
    ## 38                                                                                                                                                                                                                5337/111/2911/80310/5295/5159/5156/5155/399694/1902/56034/283748/2914/5737/552/57121/8605/5330/5332/5578/26052/9170/551/1606/1609
    ##    Count
    ## 1     19
    ## 2     15
    ## 3     29
    ## 4     30
    ## 5     15
    ## 6     12
    ## 7     17
    ## 8     15
    ## 9     21
    ## 10     9
    ## 11     7
    ## 12     6
    ## 13    18
    ## 14     9
    ## 15    13
    ## 16    11
    ## 17     8
    ## 18    12
    ## 19    49
    ## 20    65
    ## 21    24
    ## 22    62
    ## 23    33
    ## 24    44
    ## 25    30
    ## 26     9
    ## 27    35
    ## 28    46
    ## 29    25
    ## 30    32
    ## 31    18
    ## 32    37
    ## 33    20
    ## 34    30
    ## 35    23
    ## 36    21
    ## 37    19
    ## 38    25

``` r
pathview(gene.data  = res_hif1a_2a_list_ez[[1]],
                     pathway.id = "hsa04115",
                     species    = "hsa")
# ,limit      = list(gene=max(abs(cc_kegg)), cpd=1)
```

# Figure 5: Compare with ChIP-Seq

## Load datasets

``` r
# Own Kelly ChIP-Seq
load("~/S/AG/AG-Scholz-NGS/Daten/Simon/ChIP_Workflow/git_ChIPSeq_Workflow/Kelly_Hx_Annotated.peaks")
# ReMAP
load("~/S/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/git_RNAseq_Kelly_Hx/3A_ChIP-Seq_data/remap_hif1a_peaks_anno_table.peaks")
# SKNBE2
load("~/S/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/git_RNAseq_Kelly_Hx/3A_ChIP-Seq_data/SKNBE2_peaks_anno_table.peaks")
# Schödel
# load()
```

## ChIP Venns

### Own data

``` r
rna_hif1a <- res_hif1a_2a_list_ens[c("HIF1A","both")] %>% unlist() %>% unique()

chip_hif1a_all <- str_detect(names(Kelly_Hx_Annotated),pattern="HIF-1A") %>% Kelly_Hx_Annotated[.] %>% lapply('[',,"geneId") %>% unlist() %>% unique()

# Own data
input_list <- list(RNA = rna_hif1a,
                   ChIP_9929 = Kelly_Hx_Annotated[["9929_HIF-1A"]]$geneId,
                   ChIP_all = chip_hif1a_all)
                   
plt1 <- venn.diagram(
    x = input_list,
   fill = viridis(3),
    main.fontface = "bold",
    fontfamily ="Arial",
    category.names = paste(names(input_list),"\n(",input_list %>% summary() %>% .[c(1:length(input_list))],")",sep=""),
    force.unique = TRUE, na = "remove", total.population = TRUE,
    filename = NULL,
    lwd = 2,
    lty = 'blank',
    cat.fontface = "bold",
    cat.fontfamily = "arial",
   main = "HIF1A")

patchwork::wrap_elements(plt1) 
```

![](Readme_files/figure-gfm/chip_venn-1.png)<!-- -->

``` r
# get gene names
overlap <- calculate.overlap(input_list)
lapply(overlap,length)
```

    ## $a5
    ## [1] 23
    ## 
    ## $a2
    ## [1] 0
    ## 
    ## $a4
    ## [1] 14
    ## 
    ## $a6
    ## [1] 193
    ## 
    ## $a1
    ## [1] 2033
    ## 
    ## $a3
    ## [1] 0
    ## 
    ## $a7
    ## [1] 276

``` r
c(overlap$a4,overlap$a5) %>% unlist() %>% res_hif1a_2a[.,"symbol"] %>% kable()
```

| x         |
|:----------|
| TERT      |
| SLC35F3   |
| PPFIA4    |
| ELFN2     |
| ADARB2    |
| ITGB3     |
| HBA1      |
| TFB1M     |
| SDK1      |
| SP2-AS1   |
| NPIPB4    |
| ZSWIM5    |
| CADM2     |
| YAP1      |
| PGK1      |
| ANKRD37   |
| FOXE3     |
| HTR5A     |
| CRYBB2P1  |
| LNPK      |
| ALDOC     |
| MTFP1     |
| LINC02932 |
| ARRDC2    |
| LINC03047 |
| TMEM45A   |
| LUCAT1    |
| PCSK5     |
| GPI       |
| SFMBT2    |
| PFKFB4    |
| STEAP1B   |
| CROCC2    |
|           |
| HMOX1     |
| TBX3      |
| SCIRT     |

``` r
# HIF2A
lapply(res_hif1a_2a_list_ens,length)
```

    ## $HIF1A
    ## [1] 1736
    ## 
    ## $HIF2A
    ## [1] 3133
    ## 
    ## $both
    ## [1] 334

``` r
rna_hif2a <- res_hif1a_2a_list_ens[c("HIF2A","both")] %>% unlist() %>% unique()
lapply(res_hif1a_2a_list_ens,length)
```

    ## $HIF1A
    ## [1] 1736
    ## 
    ## $HIF2A
    ## [1] 3133
    ## 
    ## $both
    ## [1] 334

``` r
chip_hif2a_all <- str_detect(names(Kelly_Hx_Annotated),pattern="HIF-2A") %>% Kelly_Hx_Annotated[.] %>% lapply('[',,"geneId") %>% unlist() %>% unique()

# Own data
input_list <- list(RNA = rna_hif2a,
                   ChIP_9930 = Kelly_Hx_Annotated[["9930_HIF-2A"]]$geneId,
                   ChIP_all = chip_hif2a_all)
                   
plt1 <- venn.diagram(
    x = input_list,
   fill = viridis(3),
    main.fontface = "bold",
    fontfamily ="Arial",
    category.names = paste(names(input_list),"\n(",input_list %>% summary() %>% .[c(1:length(input_list))],")",sep=""),
    force.unique = TRUE, na = "remove", total.population = TRUE,
    filename = NULL,
    lwd = 2,
    lty = 'blank',
    cat.fontface = "bold",
    cat.fontfamily = "arial",
   main = "HIF2A")

patchwork::wrap_elements(plt1) 
```

![](Readme_files/figure-gfm/chip_venn-2.png)<!-- -->

``` r
# get gene names
overlap <- calculate.overlap(input_list)
lapply(overlap,length)
```

    ## $a5
    ## [1] 106
    ## 
    ## $a2
    ## [1] 0
    ## 
    ## $a4
    ## [1] 24
    ## 
    ## $a6
    ## [1] 352
    ## 
    ## $a1
    ## [1] 3337
    ## 
    ## $a3
    ## [1] 0
    ## 
    ## $a7
    ## [1] 239

``` r
c(overlap$a4,overlap$a5) %>% unlist() %>% res_hif1a_2a[.,"symbol"] %>% kable()
```

| x          |
|:-----------|
| LINC01320  |
| RAG1       |
| EXT1       |
| SMOC2      |
| TNFRSF19   |
| JPH1       |
| INSM1      |
| RASL11B    |
| OTOF       |
| BAHCC1     |
| PDGFC      |
| CNTN4      |
| COL23A1    |
| LINC00624  |
| SNTG2      |
| LINC01341  |
| SOX5       |
| PBXIP1     |
| MYT1L      |
| ZBTB20     |
| LMF1       |
| ZNF350-AS1 |
| ALDH1A2    |
| CALN1      |
| RFTN1      |
| DSP        |
| KCNMA1     |
| RBFOX3     |
| NOG        |
| NTRK2      |
| DSP-AS1    |
| NOD1       |
| PRICKLE2   |
| FGD5       |
| PCDH19     |
| AQP10      |
| SOBP       |
| PLXNA4     |
| KCNK13     |
| POU6F2     |
| GPER1      |
| LARGE1     |
| SOX7       |
| CFAP43     |
| TIAM1      |
| MARCHF3    |
| PIK3R1     |
| MCC        |
| PCDH15     |
| ANKS1A     |
| LRATD2     |
| SNCA       |
| IL6R       |
| SCARB1     |
| GGCT       |
| PRR15      |
| NECAB1     |
| COL5A1     |
| MTUS1      |
|            |
|            |
| KCNJ6      |
| KIF26B     |
| SSUH2      |
| ROBO2      |
| SULF1      |
| PDE10A     |
| RHPN1      |
| NRXN1      |
| SAMD4A     |
| SLC9C2     |
| CFAP46     |
| PTPRQ      |
| CYB5A      |
| LINC02073  |
|            |
|            |
| COL12A1    |
| RMST       |
| NTM        |
| SDK2       |
| TMCC3      |
| CYP2W1     |
| TUBA8      |
| PLEKHG4B   |
|            |
| CCP110     |
| LOXL2-AS1  |
| LINC00967  |
| LINC01250  |
| ANKRD60    |
| MGAT4C     |
| C8orf34    |
| ABCA1      |
| PECAM1     |
| KIF25      |
|            |
|            |
| NRG3       |
| FHAD1      |
| GLI3       |
| LINC01151  |
|            |
| CHST15     |
| CRTC3      |
| SNX25      |
| PLCB4      |
| PRKCE      |
| FLNC       |
| NTNG2      |
| CDK19      |
| ZNF19      |
| TEAD1      |
| MEF2C-AS1  |
| SARDH      |
| SYNE1      |
| ITGB5      |
| SLC2A14    |
| SCAND3     |
| MIR99AHG   |
| EGLN3      |
| MEGF6      |
| LINC03033  |
|            |
|            |
| HMOX1      |
| FOS        |
| TBX3       |
| KRT17      |
| SCIRT      |

``` r
# Compare HIF1A, HIF2A

rna_hif2a <- res_hif1a_2a_list_ens[c("HIF2A","both")] %>% unlist() %>% unique()

chip_hif2a_all <- str_detect(names(Kelly_Hx_Annotated),pattern="HIF-2A") %>% Kelly_Hx_Annotated[.] %>% lapply('[',,"geneId") %>% unlist() %>% unique()

input_list <- c(res_hif1a_2a_list_ens,list(chip_hif1a_all = chip_hif1a_all,
                   chip_hif2a_all = chip_hif2a_all))

plt1 <- venn.diagram(
    x = input_list,
   fill = viridis(5),
    main.fontface = "bold",
    fontfamily ="Arial",
    category.names = paste(names(input_list),"\n(",input_list %>% summary() %>% .[c(1:length(input_list))],")",sep=""),
    force.unique = TRUE, na = "remove", total.population = TRUE,
    filename = NULL,
    lwd = 2,
    lty = 'blank',
    cat.fontface = "bold",
    cat.fontfamily = "arial",
   main = "ChIP all"
   )

input_list <- c(res_hif1a_2a_list_ens,
                list(ChIP_9929_HIF1A = Kelly_Hx_Annotated[["9929_HIF-1A"]]$geneId,
                   ChIP_9930_HIF2A = Kelly_Hx_Annotated[["9930_HIF-2A"]]$geneId))

plt2 <- venn.diagram(
    x = input_list,
   fill = viridis(5),
    main.fontface = "bold",
    fontfamily ="Arial",
    category.names = paste(names(input_list),"\n(",input_list %>% summary() %>% .[c(1:length(input_list))],")",sep=""),
    force.unique = TRUE, na = "remove", total.population = TRUE,
    filename = NULL,
    lwd = 2,
    lty = 'blank',
    cat.fontface = "bold",
    cat.fontfamily = "arial",
   main = "ChIP good samples"
   )

patchwork::wrap_elements(plt1) + patchwork::wrap_elements(plt2) 
```

![](Readme_files/figure-gfm/chip_venn-3.png)<!-- -->

``` r
# get gene names
overlap <- calculate.overlap(input_list)
lapply(overlap,length)
```

    ## $a31
    ## [1] 0
    ## 
    ## $a30
    ## [1] 0
    ## 
    ## $a29
    ## [1] 0
    ## 
    ## $a28
    ## [1] 0
    ## 
    ## $a27
    ## [1] 0
    ## 
    ## $a26
    ## [1] 0
    ## 
    ## $a25
    ## [1] 19
    ## 
    ## $a24
    ## [1] 0
    ## 
    ## $a23
    ## [1] 0
    ## 
    ## $a22
    ## [1] 0
    ## 
    ## $a21
    ## [1] 0
    ## 
    ## $a20
    ## [1] 0
    ## 
    ## $a19
    ## [1] 0
    ## 
    ## $a18
    ## [1] 5
    ## 
    ## $a17
    ## [1] 0
    ## 
    ## $a16
    ## [1] 4
    ## 
    ## $a15
    ## [1] 56
    ## 
    ## $a14
    ## [1] 12
    ## 
    ## $a13
    ## [1] 0
    ## 
    ## $a12
    ## [1] 0
    ## 
    ## $a11
    ## [1] 0
    ## 
    ## $a10
    ## [1] 81
    ## 
    ## $a9
    ## [1] 0
    ## 
    ## $a8
    ## [1] 14
    ## 
    ## $a7
    ## [1] 19
    ## 
    ## $a6
    ## [1] 2
    ## 
    ## $a5
    ## [1] 286
    ## 
    ## $a4
    ## [1] 110
    ## 
    ## $a3
    ## [1] 328
    ## 
    ## $a2
    ## [1] 3021
    ## 
    ## $a1
    ## [1] 1698

``` r
# c(overlap$a4,overlap$a5) %>% unlist() %>% res_hif1a_2a[.,"symbol"] %>% kable()
```

### SKNBE2

``` r
# Generate SKNBE2 min list
input_list <- SKNBE2_genes_list <- SKNBE2_peaks_anno_table %>% lapply('[',,"geneId")

plt1 <- venn.diagram(
    x = input_list,
   fill = viridis(3),
    main.fontface = "bold",
    fontfamily ="Arial",
    category.names = paste(names(input_list),"\n(",input_list %>% summary() %>% .[c(1:length(input_list))],")",sep=""),
    force.unique = TRUE, na = "remove", total.population = TRUE,
    filename = NULL,
    lwd = 2,
    lty = 'blank',
    cat.fontface = "bold",
    cat.fontfamily = "arial",
   main = "SKNBE2")
patchwork::wrap_elements(plt1) 
```

![](Readme_files/figure-gfm/chip_venn_sknbe2-1.png)<!-- -->

``` r
overlap <- calculate.overlap(input_list)
lapply(overlap,length)
```

    ## $a5
    ## [1] 686
    ## 
    ## $a2
    ## [1] 24
    ## 
    ## $a4
    ## [1] 23
    ## 
    ## $a6
    ## [1] 1643
    ## 
    ## $a1
    ## [1] 26
    ## 
    ## $a3
    ## [1] 625
    ## 
    ## $a7
    ## [1] 1339

``` r
SKNBE2_genes <- c(overlap$a1,overlap$a2,overlap$a3,overlap$a4) %>% unlist() %>% unique()



# Compare with RNA-Seq
input_list <- c(res_hif1a_2a_list_ens[1],list(SKNBE2_chip = SKNBE2_genes),res_hif1a_2a_list_ens[2:3])

plt1 <- venn.diagram(
    x = input_list,
   fill = viridis(10)[c(3,10,4,5)],
    main.fontface = "bold",
    fontfamily ="Arial",
    category.names = paste(names(input_list),"\n(",input_list %>% summary() %>% .[c(1:length(input_list))],")",sep=""),
    force.unique = TRUE, na = "remove", total.population = TRUE,
    filename = NULL,
    lwd = 2,
    lty = 'blank',
    cat.fontface = "bold",
    cat.fontfamily = "arial",
   main = "RNA vs- ChIP (SKNBE2)")
patchwork::wrap_elements(plt1) 
```

![](Readme_files/figure-gfm/chip_venn_sknbe2-2.png)<!-- -->

``` r
# get gene names
overlap <- calculate.overlap(input_list)
lapply(overlap,length)
```

    ## $a6
    ## [1] 0
    ## 
    ## $a12
    ## [1] 0
    ## 
    ## $a11
    ## [1] 0
    ## 
    ## $a5
    ## [1] 0
    ## 
    ## $a7
    ## [1] 0
    ## 
    ## $a15
    ## [1] 36
    ## 
    ## $a4
    ## [1] 0
    ## 
    ## $a10
    ## [1] 0
    ## 
    ## $a13
    ## [1] 75
    ## 
    ## $a8
    ## [1] 8
    ## 
    ## $a2
    ## [1] 0
    ## 
    ## $a9
    ## [1] 1700
    ## 
    ## $a14
    ## [1] 537
    ## 
    ## $a1
    ## [1] 3058
    ## 
    ## $a3
    ## [1] 326

``` r
c(overlap$a4,overlap$a5) %>% unlist() %>% res_hif1a_2a[.,"symbol"] %>% kable()
```

| x   |
|:----|

### ReMAP

``` r
# Generate SKNBE2 min list
input_list <- remap_genes_list <- remap_hif1a_peaks_anno_table %>% lapply('[',,"geneId")
names(remap_genes_list)
```

    ## [1] "ReMap_hif1a_RCC10"      "ReMap_hif1a_K-562"      "ReMap_hif1a_HUVEC-C"   
    ## [4] "ReMap_hif1a_macrophage" "ReMap_hif1a_501-mel"

``` r
plt1 <- venn.diagram(
    x = input_list,
   fill = viridis(length(input_list)),
    main.fontface = "bold",
    fontfamily ="Arial",
    category.names = paste(names(input_list),"\n(",input_list %>% summary() %>% .[c(1:length(input_list))],")",sep=""),
    force.unique = TRUE, na = "remove", total.population = TRUE,
    filename = NULL,
    lwd = 2,
    lty = 'blank',
    cat.fontface = "bold",
    cat.fontfamily = "arial",
   main = "SKNBE2")
patchwork::wrap_elements(plt1) 
```

![](Readme_files/figure-gfm/chip_venn_remap-1.png)<!-- -->

``` r
overlap <- calculate.overlap(input_list)
names(overlap)
```

    ##  [1] "a31" "a30" "a29" "a28" "a27" "a26" "a25" "a24" "a23" "a22" "a21" "a20"
    ## [13] "a19" "a18" "a17" "a16" "a15" "a14" "a13" "a12" "a11" "a10" "a9"  "a8" 
    ## [25] "a7"  "a6"  "a5"  "a4"  "a3"  "a2"  "a1"

``` r
# lapply(overlap,length)
# SKNBE2_genes <- overlap %>% unlist() %>% unique()


# Compare with RNA-Seq
input_list <- c(res_hif1a_2a_list_ens[1],remap_genes_list[c(1,3:5)])

plt1 <- venn.diagram(
    x = input_list,
   fill = viridis(10)[c(10,3:6)],
    main.fontface = "bold",
    fontfamily ="Arial",
    category.names = paste(names(input_list),"\n(",input_list %>% summary() %>% .[c(1:length(input_list))],")",sep=""),
    force.unique = TRUE, na = "remove", total.population = TRUE,
    filename = NULL,
    lwd = 2,
    lty = 'blank',
    cat.fontface = "bold",
    cat.fontfamily = "arial",
   main = "RNA vs- ChIP (SKNBE2)")
patchwork::wrap_elements(plt1) 
```

![](Readme_files/figure-gfm/chip_venn_remap-2.png)<!-- -->

``` r
# get gene names
overlap <- calculate.overlap(input_list)
lapply(overlap,length)
```

    ## $a31
    ## [1] 79
    ## 
    ## $a30
    ## [1] 28
    ## 
    ## $a29
    ## [1] 15
    ## 
    ## $a28
    ## [1] 29
    ## 
    ## $a27
    ## [1] 9
    ## 
    ## $a26
    ## [1] 309
    ## 
    ## $a25
    ## [1] 221
    ## 
    ## $a24
    ## [1] 280
    ## 
    ## $a23
    ## [1] 31
    ## 
    ## $a22
    ## [1] 21
    ## 
    ## $a21
    ## [1] 169
    ## 
    ## $a20
    ## [1] 30
    ## 
    ## $a19
    ## [1] 31
    ## 
    ## $a18
    ## [1] 36
    ## 
    ## $a17
    ## [1] 12
    ## 
    ## $a16
    ## [1] 146
    ## 
    ## $a15
    ## [1] 587
    ## 
    ## $a14
    ## [1] 390
    ## 
    ## $a13
    ## [1] 240
    ## 
    ## $a12
    ## [1] 29
    ## 
    ## $a11
    ## [1] 254
    ## 
    ## $a10
    ## [1] 357
    ## 
    ## $a9
    ## [1] 60
    ## 
    ## $a8
    ## [1] 116
    ## 
    ## $a7
    ## [1] 170
    ## 
    ## $a6
    ## [1] 185
    ## 
    ## $a5
    ## [1] 2742
    ## 
    ## $a4
    ## [1] 2854
    ## 
    ## $a3
    ## [1] 625
    ## 
    ## $a2
    ## [1] 1220
    ## 
    ## $a1
    ## [1] 1040

``` r
c(overlap$a31) %>% unlist() %>% res_hif1a_2a[.,"symbol"] %>% kable()
```

| x         |
|:----------|
| FAM162A   |
| PGK1      |
| ANKZF1    |
| ANKRD37   |
| GBE1      |
| ABCB6     |
| PPAN      |
| EFNA3     |
| ADORA2B   |
| DDX41     |
| KDM4C     |
| LNPK      |
| PIGA      |
| HERC3     |
| AK2       |
| USP28     |
| TRIM9     |
| PPFIA4    |
| PRELID2   |
| GRPEL1    |
| DIPK2A    |
| PPP1R3E   |
|           |
| BNIP3L    |
| ALDOC     |
| ZNF160    |
| BICDL2    |
| MTFP1     |
| ARRDC2    |
| SLC25A36  |
| LINC03047 |
| NUDT18    |
| CFAP96    |
| LUCAT1    |
| ANGPTL4   |
| CPNE5     |
| PPM1H     |
| CACNA1C   |
| FOSL2     |
| NME1      |
| EIF4EBP1  |
|           |
| VDAC1     |
| TSR1      |
| NAT10     |
| PFKP      |
| NDC1      |
| NEMP1     |
| HSPA9     |
| WDR43     |
| SLC7A5    |
| ZNF770    |
| CHSY1     |
| NAV1      |
| IMP3      |
| TMEM248   |
| OPA3      |
| DLGAP1    |
| TRMU      |
| SAP30     |
| DDX50     |
| DENND1A   |
| FEM1C     |
| DPCD      |
| FAM117B   |
| PDE12     |
| TOMM20    |
| PPFIA3    |
| SLC25A26  |
| MRPL17    |
| SAP30-DT  |
| SERGEF    |
| METTL26   |
| DUS3L     |
| PFKFB4    |
| KLHL3     |
| SLC8A2    |
| THAP8     |
| EPHB1     |

# 

# \###Old code

## Enhanced volcano

## Volcanos

``` r
n <- "Kelly.Hx.vs.Nx"
xlim <- 12
ylim <- -250
res <- results_list[[n]]
res <- res_shrink_list[[n]] %>% data.frame()

# of deg genes
up <- subset(results_list[[n]], padj< 0.05 & log2FoldChange > 1) %>% nrow()
down <- subset(results_list[[n]], padj< 0.05 & log2FoldChange < -1) %>% nrow()
total <- up+down
deg <- subset(results_list[[n]], padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) %>% data.frame()

# points outside the grid
outx <- subset(res, log2FoldChange > xlim | log2FoldChange < -xlim) %>% rownames()
outy <- subset(res, padj < 10^ylim) %>% rownames()

res$outlier <- ifelse(rownames(res) %in% c(outx,outy),"yes","no")
res$deg <- ifelse(rownames(res) %in% rownames(deg),"yes","no")

res[outx,"log2FoldChange"] <- ifelse(res[outx,"log2FoldChange"] > xlim,xlim,-xlim)
res[outy,"padj"] <- 10^ylim

volcano_kelly <- ggplot(res,aes(x=log2FoldChange,y=-log10(padj),color=deg, shape=outlier)) +
  geom_hline(yintercept = 0, linewidth = 0.2) + 
  geom_vline(xintercept = 0, linewidth = 0.2) +
  geom_point(size=2, stroke=1) +
  labs(title = "Hypoxic response in Kelly (l2FC > 1, p<0.05") +
  ylab("padj (-log10)") +
  xlab("log2-foldchange") +
  scale_shape_manual(values = c(16,21)) + 
  scale_color_manual(values = c("grey80","orchid4")) + 
  scale_fill_manual(values=c(colors[1],"white",colors[2],"white")) +
  theme_bw() +
  # geom_text_repel(label=res$symbol, color="black") + 
  removeGrid(x=T, y=T)

volcano_kelly
```

![](Readme_files/figure-gfm/2_volcanos-1.png)<!-- -->

``` r
# Hif1a
n <- "Hif1a.Hx.vs.Nx"
n2 <- "Hif1aHxNx.vs.KellyHxNx"
col <- colors[4]
```

``` r
vol_kelly <- EnhancedVolcano(res,
    colAlpha=0.2,
    lab = res[["symbol"]],
    x = 'log2FoldChange',
    y = 'padj',
    col=c("grey","grey","grey","orchid4"),
    title = "Hypoxic response in kelly cells",
    titleLabSize = 12,
    subtitle = paste0("upregulated: ",up,", downregulated: ",down,"\n(total: ",total,")"),
    subtitleLabSize = 10,
    caption = NULL,
    # xlim = c(-11,10),
    # ylim = c(0,50),
    pCutoff = 0.05,
    FCcutoff = 0.5,
    maxoverlapsConnectors = 30,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = "grey70",
    legendLabels=c('ns','ns','ns',
      'padj < 0.05 & Log2FC > 0.5'),
    labSize = 4,
    axisLabSize = 12,
    legendLabSize = 12,
    legendIconSize = 4,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    pointSize = 2
)
vol_kelly
```

## Cluster genes

``` r
# Hypoxia up:
## Volcano

deg <- res_hif1a_2a %>% filter(Kelly.Hx.vs.Nx.log2FoldChange > 1) %>% mutate(log2FoldChange = Kelly.Hx.vs.Nx.log2FoldChange)
volcano_kelly_up <- volcano_sk3(n="Kelly.Hx.vs.Nx",n2="Kelly.Hx.vs.Nx", deg=deg,col=colors_vul[5], celline="Kelly") + labs(title=paste0("DEG genes Hx upregulated"), subtitle = paste0("upregulated: ",nrow(deg)))


# Hif1a ~ Kelly, color: 
color3 <- cut(deg$Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange,c(-Inf,seq(-3,3,by=1),+Inf))
cc <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=8))
cc <- viridis(8)
cc[4:5] <- "grey40"
cluster_h1vK <- ggplot(deg,aes(x=Kelly.Hx.vs.Nx.log2FoldChange, y=Hif1a.Hx.vs.Nx.log2FoldChange, color=color3)) +
  geom_hline(yintercept=c(0)) +
  geom_vline(xintercept=c(0)) +
  geom_abline(intercept=c(1,-1)) +
  geom_point(alpha=0.5) +
  # scale_color_viridis_d(option = 'C') +
  scale_colour_manual(name="log2-FC",values=cc) +
  ggtitle(label="Kelly: Hx vs. Nx") + 
  coord_cartesian(xlim = c(-5, 12),ylim = c(-5,12))
cluster_h1vK

# Hif1a ~ Hif2a, color: 
color3 <- cut(deg$Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange,c(-Inf,seq(-3,3,by=1),+Inf))
cc <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=8))
cc <- viridis(8)
cc[4:5] <- "grey40"
cluster_h2v1 <- ggplot(deg,aes(x=Hif1aHxNx.vs.KellyHxNx.log2FoldChange, y=Hif2aHxNx.vs.KellyHxNx.log2FoldChange, color=color3)) +
  geom_hline(yintercept=c(0)) +
  geom_vline(xintercept=c(0)) +
  geom_abline(intercept=c(1,-1)) +
  geom_point(alpha=0.5) +
  # scale_color_viridis_d(option = 'C') +
  scale_colour_manual(name="log2-FC",values=cc) +
  ggtitle(label="Kelly: Hx vs. Nx") + 
  coord_cartesian(xlim = c(-10, 10),ylim = c(-10,10))
cluster_h2v1


# hif1a = 0, Hiff2A > 1
subset(deg, Hif2aHxNx.vs.KellyHxNx.log2FoldChange >1 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange > -1 & Hif1aHxNx.vs.KellyHxNx.log2FoldChange < 1))
goi <- "ENSG00000135636"
plotCounts_anno(goi) + labs(title="Hif2A target")
point_xy <- data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF2A",gene=goi)

# Hif1a > 1, Hif2A > 1 Hif2:Hif1 > 1
subset(deg, Hif2aHxNx.vs.KellyHxNx.log2FoldChange >1 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange > 1 & Hif1aHxNx.vs.KellyHxNx.log2FoldChange > 1) & Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange > 1)
goi <- "ENSG00000240801"
plotCounts_anno(goi) + labs(title="Hif2A target")
point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF2A",gene=goi))

# Hif1a > 1, Hif2A > 1 Hif2:Hif1 < 1
subset(deg, Hif2aHxNx.vs.KellyHxNx.log2FoldChange >1 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange > 1 & Hif1aHxNx.vs.KellyHxNx.log2FoldChange > 1) & Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange < 1)
goi <- "ENSG00000124205"
plotCounts_anno(goi) + labs(title="Hif1A & Hif2A target")
point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF1A & HIF2A",gene=goi))

# Hif1a > 1, Hif2A > 1, Hif2:Hif1 < -1
subset(deg, Hif2aHxNx.vs.KellyHxNx.log2FoldChange >1 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange > 1 & Hif1aHxNx.vs.KellyHxNx.log2FoldChange > 1) & Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange < -1)
goi <- "ENSG00000249815"
plotCounts_anno(goi) + labs(title="Hif1A")
point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF1A",gene=goi))

# Hif1a > 1, Hif2A < 1, Hif1a > Hif2A
subset(deg, abs(Hif2aHxNx.vs.KellyHxNx.log2FoldChange) < 1 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange > 1 & Hif1aHxNx.vs.KellyHxNx.log2FoldChange > 1) & Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange < -1)
goi <- "ENSG00000186352"
plotCounts_anno(goi) + labs(title="Hif1A")
point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF1A",gene=goi))
goi <- "ENSG00000099994"
plotCounts_anno(goi) + labs(title="Hif1A")

point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF1A",gene=goi))

# Hif1a > 1, Hif2A < 1, Hif1a = Hif2A
subset(deg, Hif2aHxNx.vs.KellyHxNx.log2FoldChange < 0 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange > 1 & Hif1aHxNx.vs.KellyHxNx.log2FoldChange > 1) & ((-Hif2aHxNx.vs.KellyHxNx.log2FoldChange - Hif1aHxNx.vs.KellyHxNx.log2FoldChange) < 1))
goi <- "ENSG00000125384"
plotCounts_anno(goi) + labs(title="Hif2A")
point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF2A",gene=goi))
goi <- "ENSG00000099994"
point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF1A",gene=goi))

# Hif1a < 1, Hif2A < 0
subset(deg, Hif2aHxNx.vs.KellyHxNx.log2FoldChange < 0 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange < 1 & Hif1aHxNx.vs.KellyHxNx.log2FoldChange > -1))
goi <- "ENSG00000129675"
plotCounts_anno(goi) + labs(title="Hif2A")
point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF2A",gene=goi))

# Hif1a < -1, Hif2A < 0
subset(deg, Hif2aHxNx.vs.KellyHxNx.log2FoldChange < 0 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange < -1 & Hif1aHxNx.vs.KellyHxNx.log2FoldChange < -1) & Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange < -1)
goi <- "ENSG00000023445"
plotCounts_anno(goi) + labs(title="Hif2A")
point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF2A",gene=goi))

# Hif1a < -1, Hif2A < -1, Hif2a:Hif1a < 1
subset(deg, Hif2aHxNx.vs.KellyHxNx.log2FoldChange < -1 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange < -1 & Hif1aHxNx.vs.KellyHxNx.log2FoldChange < -1) & abs(Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange) < 1)
goi <- "ENSG00000136235"
plotCounts_anno(goi) + labs(title="Hif1A & Hif2A target")
point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF1A & HIF2A",gene=goi))


# Hif1a < -1, Hif2A < -1, Hif2a:Hif1a < 1
subset(deg, Hif2aHxNx.vs.KellyHxNx.log2FoldChange < -1 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange < -1 & Hif1aHxNx.vs.KellyHxNx.log2FoldChange < -1) & Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange> 1)
goi <- "ENSG00000278718"
plotCounts_anno(goi) + labs(title="Hif1A")
point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF1A",gene=goi))

# Hif1a < -1, Hif2A < -1, Hif2a:Hif1a > 1
subset(deg, abs(Hif2aHxNx.vs.KellyHxNx.log2FoldChange) < 1 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange < -1 & Hif1aHxNx.vs.KellyHxNx.log2FoldChange < -1) & Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange> 1)
goi <- "ENSG00000170525"
plotCounts_anno(goi) + labs(title="Hif1A")
point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF1A",gene=goi))

# Hif1a < -1, Hif2A > 1
subset(deg, Hif2aHxNx.vs.KellyHxNx.log2FoldChange > 1 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange < -1 & Hif1aHxNx.vs.KellyHxNx.log2FoldChange < -1) & (-Hif2aHxNx.vs.KellyHxNx.log2FoldChange - Hif1aHxNx.vs.KellyHxNx.log2FoldChange) < 1)
goi <- "ENSG00000137819"
plotCounts_anno(goi) + labs(title="Hif1A")
point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF1A",gene=goi))

# Hif1a < -1, Hif2A > 1
subset(deg, Hif2aHxNx.vs.KellyHxNx.log2FoldChange > 1 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange < -1 & Hif1aHxNx.vs.KellyHxNx.log2FoldChange < -1) & (-Hif2aHxNx.vs.KellyHxNx.log2FoldChange - Hif1aHxNx.vs.KellyHxNx.log2FoldChange) < 1) %>% arrange(Hif1aHxNx.vs.KellyHxNx.log2FoldChange)
goi <- "ENSG00000285321"
plotCounts_anno(goi) + labs(title="Hif1A")
point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="HIF1A",gene=goi))

# Hif1a -> Hif2A
goi <- subset(deg, Hif2aHxNx.vs.KellyHxNx.log2FoldChange > 5 & (Hif1aHxNx.vs.KellyHxNx.log2FoldChange < 0 & Hif1aHxNx.vs.KellyHxNx.log2FoldChange < 0) & (-Hif2aHxNx.vs.KellyHxNx.log2FoldChange - Hif1aHxNx.vs.KellyHxNx.log2FoldChange) < 1) %>% arrange(Hif1aHxNx.vs.KellyHxNx.log2FoldChange) 
plotCounts_SK(goi %>% rownames()) + labs(title="Hif1A")

# Hif1a -> Hif2A
goi <- subset(deg, Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange > 1 & (Hif2aHxNx.vs.KellyHxNx.log2FoldChange + Hif1aHxNx.vs.KellyHxNx.log2FoldChange) > 1 & (Hif2aHxNx.vs.KellyHxNx.log2FoldChange + Hif1aHxNx.vs.KellyHxNx.log2FoldChange) < 2) %>% arrange(Hif1aHxNx.vs.KellyHxNx.log2FoldChange) 

plotCounts_SK(goi[1:9,] %>% rownames()) + labs(title="Hif1A")

goi <- "ENSG00000127824"
plotCounts_anno(goi=goi)
point_xy <- rbind(point_xy,data.frame(x=deg[goi,"Hif1aHxNx.vs.KellyHxNx.log2FoldChange"],y=deg[goi,"Hif2aHxNx.vs.KellyHxNx.log2FoldChange"], target="opposing",gene=goi))

cc <- c(viridis(8,option = 'A'),colors[c(4,2,6)])
cc <- c(rep("grey40",8),colors[c(4,2,6)], "orange")

point_xy$symbol <- mcols(dds)[point_xy$gene,"symbol"]

cluster_h2v1 + geom_point(data = point_xy, aes(x=x, y=y, color=target), size=5) + scale_colour_manual(name="log2-FC",values=cc) + geom_text_repel(data = point_xy, aes(x=x, y=y), label=point_xy$symbol, colour="black")
```

### Other cluster

``` r
# Hif1a ~ Hif2a, color: Kelly_Hx
color1 <- cut(res_hif1a_2a$Kelly.Hx.vs.Nx.log2FoldChange,c(-Inf,seq(-3,3,by=1),+Inf))
cc <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=8))
cc <- viridis(8)
cc[4:5] <- "grey40"
cluster_hx <- ggplot(res_hif1a_2a,aes(x=Hif1a.Hx.vs.Nx.log2FoldChange, y=Hif2a.Hx.vs.Nx.log2FoldChange, color=color1)) +
  geom_hline(yintercept=c(0)) +
  geom_vline(xintercept=c(0)) +
  geom_abline(intercept=c(1,-1)) +
  geom_point(alpha=0.5) +
  # scale_color_viridis_d(option = 'C') +
  scale_colour_manual(name="log2-FC",values=cc) +
  ggtitle(label="Kelly: Hx vs. Nx") + 
  coord_cartesian(xlim = c(-5, 12),ylim = c(-5,12))

# Hif1a ~ Hif2a, color: Hif2a vs. Hif1a
color2 <- cut(res_hif1a_2a$Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange,c(-Inf,seq(-3,3,by=1),+Inf))
clusterh1v2 <- ggplot(res_hif1a_2a,aes(x=Hif1a.Hx.vs.Nx.log2FoldChange, y=Hif2a.Hx.vs.Nx.log2FoldChange, color=color2)) +
  geom_hline(yintercept=c(0)) +
  geom_vline(xintercept=c(0)) +
  geom_abline(intercept=c(1,-1)) +
  geom_point(alpha=0.5) +
  scale_colour_manual(name="log2-FC",values=cc) +
  # scale_color_viridis_d(option = 'C') +
  ggtitle(label="Hif1a vs. Hif2a") + 
  coord_cartesian(xlim = c(-5, 12),ylim = c(-5,12))

cluster_hx + clusterh1v2 + plot_layout(guides = "collect", axes="collect", axis_titles="collect") & 
  theme(legend.position = 'right')


# Hif1a ~ Kelly, color: 
color3 <- cut(res_hif1a_2a$Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange,c(-Inf,seq(-3,3,by=1),+Inf))
cc <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=8))
cc <- viridis(8)
cc[4:5] <- "grey40"
cluster_h1vK <- ggplot(res_hif1a_2a,aes(x=Kelly.Hx.vs.Nx.log2FoldChange, y=Hif1a.Hx.vs.Nx.log2FoldChange, color=color3)) +
  geom_hline(yintercept=c(0)) +
  geom_vline(xintercept=c(0)) +
  geom_abline(intercept=c(1,-1)) +
  geom_point(alpha=0.5) +
  # scale_color_viridis_d(option = 'C') +
  scale_colour_manual(name="log2-FC",values=cc) +
  ggtitle(label="Kelly: Hx vs. Nx") + 
  coord_cartesian(xlim = c(-5, 12),ylim = c(-5,12))
cluster_h1vK


color3 <- cut(res_hif1a_2a$Hif2aHxNx.vs.Hif1aHxNx.log2FoldChange,c(-Inf,seq(-3,3,by=1),+Inf))
clusterh1vK <- ggplot(res_hif1a_2a,aes(x=Hif1a.Hx.vs.Nx.log2FoldChange, y=Hif2a.Hx.vs.Nx.log2FoldChange, color=color3)) +
  geom_hline(yintercept=c(0)) +
  geom_vline(xintercept=c(0)) +
  geom_abline(intercept=c(1,-1)) +
  geom_point(alpha=0.5) +
  scale_colour_manual(name="log2-FC",values=cc) +
  # scale_color_viridis_d(option = 'C') +
  ggtitle(label="Hif1a vs. Hif2a") + 
  coord_cartesian(xlim = c(-5, 12),ylim = c(-5,12))
clusterh1vK
```
