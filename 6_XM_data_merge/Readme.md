Master_Table
================

- [Input data](#input-data)
- [Kelly: fresh contrasts from the new
  dds](#kelly-fresh-contrasts-from-the-new-dds)
- [SK-N-AS / SH-SY5Y: independent per-line Hx vs
  Nx](#sk-n-as--sh-sy5y-independent-per-line-hx-vs-nx)
- [Combine into master table](#combine-into-master-table)
  - [Overlap of the Ensembl IDs](#overlap-of-the-ensembl-ids)
- [Target classification (Kelly
  only)](#target-classification-kelly-only)
- [Hypoxia baseline (SK-N-AS /
  SH-SY5Y)](#hypoxia-baseline-sk-n-as--sh-sy5y)
  - [Target vs. baseline](#target-vs-baseline)
  - [Venn: Target vs. Hx_Baseline](#venn-target-vs-hx_baseline)
- [Literature coverage](#literature-coverage)
  - [Candidates: strong, Kelly-specific, undescribed in
    hypoxia](#candidates-strong-kelly-specific-undescribed-in-hypoxia)
- [Master counts table](#master-counts-table)
- [Export](#export)

This document merges the three RNA-seq / microarray data sets of the
project into one **master table**, classifies every gene by its HIF
dependency, and exports the result as Excel. Code blocks are collapsed -
click *Code* to expand.

<details>

<summary>

Code
</summary>

``` r
library(dplyr)
library(DESeq2)
library(tidyverse)
library(knitr)
library(kableExtra)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(VennDiagram)
library(patchwork)

# Suppress the .log files the VennDiagram package writes on every call.
invisible(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))

colors <- c("lavenderblush3", "lavenderblush4", "#90caf9", "#1976d2",
            "#82e0aa", "#239b56", "#f8c471", "#b9770e")
```

</details>

# Input data

Three sources go into the master table: the Kelly line with its three
HIF knockouts, the two MYCN-non-amplified reference lines, and the
high-risk microarray gene list.

<details>

<summary>

Code
</summary>

``` r
dds_kelly <- readRDS(
  file = "/Users/simonkelterborn/Library/CloudStorage/GoogleDrive-simon.kelterborn@gmail.com/Meine Ablage/! OneDrive Charité/dds_filtered_with_metadata.rds"
)

load(
  file = "/Users/simonkelterborn/Library/CloudStorage/GoogleDrive-simon.kelterborn@gmail.com/Meine Ablage/! OneDrive Charité/SKNAS_SHSY5Y.dds"
)
dds_SKNAS <- dds
rm(dds)

HR_table <- readxl::read_xlsx("/Users/simonkelterborn/Library/CloudStorage/GoogleDrive-simon.kelterborn@gmail.com/Meine Ablage/! OneDrive Charité/HR_genes_Microarray.xlsx")

tab(data.frame(
  Data_set = c("Kelly (WT + HIF1a/HIF1b/HIF2a KO)",
               "SK-N-AS / SH-SY5Y",
               "High-risk microarray"),
  File     = c("dds_filtered_with_metadata.rds", "SKNAS_SHSY5Y.dds",
               "HR_genes_Microarray.xlsx"),
  Genes    = c(nrow(dds_kelly), nrow(dds_SKNAS), nrow(HR_table)),
  Samples  = c(ncol(dds_kelly), ncol(dds_SKNAS), NA)
), "Input data sets")
```

</details>

**Input data sets**

| Data_set | File | Genes | Samples |
|:---|:---|---:|---:|
| Kelly (WT + HIF1a/HIF1b/HIF2a KO) | dds_filtered_with_metadata.rds | 64537 | 88 |
| SK-N-AS / SH-SY5Y | SKNAS_SHSY5Y.dds | 76045 | 12 |
| High-risk microarray | HR_genes_Microarray.xlsx | 27855 | NA |

# Kelly: fresh contrasts from the new dds

DESeq2 is run once with the full interaction design
`~genotype + condition + genotype:condition`. All contrasts of interest
are then derived from that single fit as weighted coefficient vectors,
so every comparison shares one dispersion estimate. The coefficient
names are resolved by pattern instead of being hardcoded, because DESeq2
rewrites factor level separators (`HIF1a-KO` becomes `HIF1a.KO`).

<details>

<summary>

Code
</summary>

``` r
colData(dds_kelly)$genotype <- relevel(factor(colData(dds_kelly)$genotype), "WT")
colData(dds_kelly)$condition <- relevel(factor(colData(dds_kelly)$condition), "Nx")

design(dds_kelly) <- ~genotype + condition + genotype:condition
dds_kelly <- DESeq(dds_kelly)

# entrezid is a 1-to-many list-column (some genes map to multiple Entrez IDs) -
# flatten to a single string so results/data.frame conversions don't explode rows.
entrez_flat <- vapply(mcols(dds_kelly)$entrezid, function(x) {
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}, character(1))

rn <- resultsNames(dds_kelly)

# Resolve the actual (possibly sanitized, e.g. "HIF1a-KO" -> "HIF1a.KO") coefficient
# names by pattern rather than hardcoding them, since DESeq2 rewrites factor level
# separators and level order can change if the input dds changes.
find_coef <- function(pattern) {
  hits <- grep(pattern, rn, value = TRUE)
  if (length(hits) != 1) {
    stop("Expected exactly 1 match for pattern '", pattern, "', found: ",
         paste(hits, collapse = ", "))
  }
  hits
}

coef_geno <- list(
  Hif1a = find_coef("^genotype.*HIF1a.*_vs_WT$"),
  Hif1b = find_coef("^genotype.*HIF1b.*_vs_WT$"),
  Hif2a = find_coef("^genotype.*HIF2a.*_vs_WT$")
)
coef_cond <- find_coef("^condition_Hx_vs_Nx$")
coef_int <- list(
  Hif1a = find_coef("^genotype.*HIF1a.*\\.condition"),
  Hif1b = find_coef("^genotype.*HIF1b.*\\.condition"),
  Hif2a = find_coef("^genotype.*HIF2a.*\\.condition")
)

# Named weight vectors, resolved to the real coefficient names above.
# Mirrors the full contrast set from 2B_DGE/2B_DGE.Rmd, computed fresh here.
contrast_defs <- list(
  Kelly.Hx.vs.Nx = c(setNames(1, coef_cond)),
  Hif1a.Hx.vs.Nx = c(setNames(1, coef_cond), setNames(1, coef_int$Hif1a)),
  Hif1b.Hx.vs.Nx = c(setNames(1, coef_cond), setNames(1, coef_int$Hif1b)),
  Hif2a.Hx.vs.Nx = c(setNames(1, coef_cond), setNames(1, coef_int$Hif2a)),

  Nx.Hif1a.vs.Kelly = c(setNames(1, coef_geno$Hif1a)),
  Nx.Hif1b.vs.Kelly = c(setNames(1, coef_geno$Hif1b)),
  Nx.Hif2a.vs.Kelly = c(setNames(1, coef_geno$Hif2a)),

  Hx.Hif1a.vs.Kelly = c(setNames(1, coef_geno$Hif1a), setNames(1, coef_int$Hif1a)),
  Hx.Hif1b.vs.Kelly = c(setNames(1, coef_geno$Hif1b), setNames(1, coef_int$Hif1b)),
  Hx.Hif2a.vs.Kelly = c(setNames(1, coef_geno$Hif2a), setNames(1, coef_int$Hif2a)),

  Hx.Hif2a.vs.Hif1a = c(setNames(-1, coef_geno$Hif1a), setNames(-1, coef_int$Hif1a),
                        setNames(1, coef_geno$Hif2a), setNames(1, coef_int$Hif2a)),
  Hx.Hif1b.vs.Hif1a = c(setNames(-1, coef_geno$Hif1a), setNames(-1, coef_int$Hif1a),
                        setNames(1, coef_geno$Hif1b), setNames(1, coef_int$Hif1b)),
  Hx.Hif1b.vs.Hif2a = c(setNames(-1, coef_geno$Hif2a), setNames(-1, coef_int$Hif2a),
                        setNames(1, coef_geno$Hif1b), setNames(1, coef_int$Hif1b)),

  Hif1aHxNx.vs.KellyHxNx = c(setNames(1, coef_int$Hif1a)),
  Hif2aHxNx.vs.KellyHxNx = c(setNames(1, coef_int$Hif2a)),
  Hif1bHxNx.vs.KellyHxNx = c(setNames(1, coef_int$Hif1b)),
  Hif2aHxNx.vs.Hif1aHxNx = c(setNames(-1, coef_int$Hif1a), setNames(1, coef_int$Hif2a)),

  Hx.Hif1b.vs.Hif12a = c(setNames(-0.5, coef_geno$Hif1a), setNames(-0.5, coef_geno$Hif2a),
                         setNames(1, coef_geno$Hif1b),
                         setNames(-0.5, coef_int$Hif1a), setNames(-0.5, coef_int$Hif2a),
                         setNames(1, coef_int$Hif1b)),

  Hx.Kelly.vs.allHIFs = c(setNames(1/3, coef_geno$Hif1a), setNames(1/3, coef_geno$Hif2a),
                          setNames(1/3, coef_geno$Hif1b),
                          setNames(1/3, coef_int$Hif1a), setNames(1/3, coef_int$Hif2a),
                          setNames(1/3, coef_int$Hif1b)),

  # NOTE: kept "NOT confirmed" caveat from the original 2B_DGE.Rmd analysis
  Hx.vs.Nx = c(setNames(1, coef_cond),
               setNames(1/4, coef_int$Hif1a), setNames(1/4, coef_int$Hif2a),
               setNames(1/4, coef_int$Hif1b))
)

build_contrast_vec <- function(weights) {
  v <- setNames(rep(0, length(rn)), rn)
  v[names(weights)] <- unname(weights)
  unname(v)
}

results_list <- lapply(contrast_defs, function(w) {
  r1 <- results(dds_kelly, contrast = build_contrast_vec(w))
  r1$symbol  <- mcols(dds_kelly)$symbol
  r1$ENSEMBL <- mcols(dds_kelly)$gene_id
  r1$ENTREZ  <- entrez_flat
  r1
})

tab(data.frame(
  Contrast = names(results_list),
  DEGs     = sapply(results_list, function(x) sum(x$padj < 0.05, na.rm = TRUE)),
  row.names = NULL
), "Significant genes per contrast (padj < 0.05)")
```

</details>

**Significant genes per contrast (padj \< 0.05)**

| Contrast               |  DEGs |
|:-----------------------|------:|
| Kelly.Hx.vs.Nx         | 16675 |
| Hif1a.Hx.vs.Nx         | 16883 |
| Hif1b.Hx.vs.Nx         |  7465 |
| Hif2a.Hx.vs.Nx         | 11905 |
| Nx.Hif1a.vs.Kelly      |   623 |
| Nx.Hif1b.vs.Kelly      |  5717 |
| Nx.Hif2a.vs.Kelly      |  2424 |
| Hx.Hif1a.vs.Kelly      |  6956 |
| Hx.Hif1b.vs.Kelly      | 15787 |
| Hx.Hif2a.vs.Kelly      | 11524 |
| Hx.Hif2a.vs.Hif1a      | 14604 |
| Hx.Hif1b.vs.Hif1a      | 16341 |
| Hx.Hif1b.vs.Hif2a      | 13268 |
| Hif1aHxNx.vs.KellyHxNx |  3050 |
| Hif2aHxNx.vs.KellyHxNx |  7976 |
| Hif1bHxNx.vs.KellyHxNx |  9584 |
| Hif2aHxNx.vs.Hif1aHxNx |  9882 |
| Hx.Hif1b.vs.Hif12a     | 15137 |
| Hx.Kelly.vs.allHIFs    | 12169 |
| Hx.vs.Nx               | 18937 |

# SK-N-AS / SH-SY5Y: independent per-line Hx vs Nx

These two lines are analysed separately from Kelly - each gets its own
DESeq2 fit on the subset of its own samples with the simple design
`~Condition`, so their hypoxia response is estimated independently of
the Kelly model.

<details>

<summary>

Code
</summary>

``` r
get_line_hxvsnx <- function(dds_all, celltype_value, label) {
  dds_sub <- dds_all[, colData(dds_all)$Celltype == celltype_value]
  colData(dds_sub)$Condition <- relevel(factor(colData(dds_sub)$Condition), "Normoxia")
  design(dds_sub) <- ~Condition
  dds_sub <- DESeq(dds_sub)
  res <- results(dds_sub, name = "Condition_Hypoxia_vs_Normoxia")

  out <- data.frame(
    Ensembl = rownames(res),
    baseMean = res$baseMean,
    log2FC   = res$log2FoldChange,
    padj     = res$padj
  )
  colnames(out)[2:4] <- paste0(c("baseMean_", "log2FC_", "padj_"), label)
  out
}

sknas_df  <- get_line_hxvsnx(dds_SKNAS, "SK-N-AS", "SKNAS")
shsy5y_df <- get_line_hxvsnx(dds_SKNAS, "SH-SY5Y", "SHSY5Y")
```

</details>

# Combine into master table

All three sources are joined on Ensembl gene ID, with Kelly as the
reference gene space. Every contrast contributes a `log2FoldChange` and
a `padj` column.

<details>

<summary>

Code
</summary>

``` r
res_table <- lapply(results_list, function(x) data.frame(x)[, c("log2FoldChange", "padj")])
res_table <- do.call('cbind', res_table)

master_table <- data.frame(
  Ensembl        = mcols(dds_kelly)$gene_id,
  ENTREZ         = entrez_flat,
  symbol         = mcols(dds_kelly)$symbol,
  baseMean_Kelly = results_list[["Kelly.Hx.vs.Nx"]]$baseMean
)

master_table <- cbind(master_table, res_table)

master_table <- master_table %>%
  left_join(sknas_df, by = "Ensembl") %>%
  left_join(shsy5y_df, by = "Ensembl") %>%
  left_join(HR_table %>% dplyr::rename(p_value_HR = p_value), by = "Ensembl")
```

</details>

## Overlap of the Ensembl IDs

How well the three sources cover the Kelly gene space:

<details>

<summary>

Code
</summary>

``` r
ens_kelly <- master_table$Ensembl

tab(data.frame(
  Data_set          = c("Kelly (reference)", "SK-N-AS / SH-SY5Y", "High-risk microarray"),
  Genes             = c(length(ens_kelly), nrow(sknas_df), nrow(HR_table)),
  Shared_with_Kelly = c(NA, sum(sknas_df$Ensembl %in% ens_kelly),
                        sum(HR_table$Ensembl %in% ens_kelly)),
  Not_in_Kelly      = c(NA, sum(!sknas_df$Ensembl %in% ens_kelly),
                        sum(!HR_table$Ensembl %in% ens_kelly))
), "Ensembl ID overlap with the Kelly gene space")
```

</details>

**Ensembl ID overlap with the Kelly gene space**

| Data_set             | Genes | Shared_with_Kelly | Not_in_Kelly |
|:---------------------|------:|------------------:|-------------:|
| Kelly (reference)    | 64537 |                NA |           NA |
| SK-N-AS / SH-SY5Y    | 76045 |             57785 |        18260 |
| High-risk microarray | 27855 |             25398 |         2457 |

Being present is not the same as being testable: DESeq2’s independent
filtering drops low-count genes, which leaves their `padj` as `NA`. That
distinction matters for the baseline classification below, so it is
counted separately here.

<details>

<summary>

Code
</summary>

``` r
coverage <- function(label) {
  bm <- master_table[[paste0("baseMean_", label)]]
  pa <- master_table[[paste0("padj_", label)]]
  data.frame(Line                     = label,
             Not_in_data_set          = sum(is.na(bm)),
             Present_but_not_testable = sum(!is.na(bm) & is.na(pa)),
             Testable                 = sum(!is.na(pa)))
}

tab(rbind(coverage("SKNAS"), coverage("SHSY5Y")),
    "Gene coverage of the two baseline lines")
```

</details>

**Gene coverage of the two baseline lines**

| Line   | Not_in_data_set | Present_but_not_testable | Testable |
|:-------|----------------:|-------------------------:|---------:|
| SKNAS  |            6752 |                    36577 |    21208 |
| SHSY5Y |            6752 |                    37072 |    20713 |

<details>

<summary>

Code
</summary>

``` r
master_table %>%
  dplyr::select(Ensembl, symbol, baseMean_Kelly,
                Kelly.Hx.vs.Nx.log2FoldChange, Kelly.Hx.vs.Nx.padj,
                log2FC_SKNAS, log2FC_SHSY5Y) %>%
  head(5) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~signif(.x, 3))) %>%
  tab(sprintf("master_table preview - first 5 of %s genes, 7 of %d columns",
              format(nrow(master_table), big.mark = ","), ncol(master_table)))
```

</details>

**master_table preview - first 5 of 64,537 genes, 7 of 52 columns**

| Ensembl | symbol | baseMean_Kelly | Kelly.Hx.vs.Nx.log2FoldChange | Kelly.Hx.vs.Nx.padj | log2FC_SKNAS | log2FC_SHSY5Y |
|:---|:---|---:|---:|---:|---:|---:|
| ENSG00000000003 | TSPAN6 | 1.21e+03 | 0.111 | 0.7820 | 0.0795 | 0.0862 |
| ENSG00000000005 | TNMD | 2.01e-02 | 0.227 | 0.9930 | 0.9080 | NA |
| ENSG00000000419 | DPM1 | 1.48e+03 | -0.271 | 0.0128 | -0.1050 | -0.1570 |
| ENSG00000000457 | SCYL3 | 7.48e+02 | 0.390 | 0.0000 | 0.1250 | 0.0748 |
| ENSG00000000460 | C1orf112 | 1.15e+03 | -0.989 | 0.0000 | -0.8430 | -0.6820 |

# Target classification (Kelly only)

A gene counts as a **HIF target** if it is induced by hypoxia in WT
Kelly (`padj < 0.05` and `log2FC > 1`) *and* that induction is lost in
the respective knockout, i.e. the interaction term is significantly
negative (`padj < 0.05` and `log2FC < -1`). Genes meeting both criteria
for HIF1A and HIF2A are labelled `Hif1a_Hif2a`. Only Kelly data is used
here.

<details>

<summary>

Code
</summary>

``` r
padj_cut <- 0.05
lfc_cut  <- 1

sig_up   <- function(lfc, padj) !is.na(padj) & padj < padj_cut & !is.na(lfc) & lfc >  lfc_cut
sig_down <- function(lfc, padj) !is.na(padj) & padj < padj_cut & !is.na(lfc) & lfc < -lfc_cut

hx_up <- sig_up(master_table$Kelly.Hx.vs.Nx.log2FoldChange,
                master_table$Kelly.Hx.vs.Nx.padj)

is_hif1a <- hx_up & sig_down(master_table$Hif1aHxNx.vs.KellyHxNx.log2FoldChange,
                             master_table$Hif1aHxNx.vs.KellyHxNx.padj)
is_hif2a <- hx_up & sig_down(master_table$Hif2aHxNx.vs.KellyHxNx.log2FoldChange,
                             master_table$Hif2aHxNx.vs.KellyHxNx.padj)

master_table$Target <- dplyr::case_when(
  is_hif1a & is_hif2a ~ "Hif1a_Hif2a",
  is_hif1a            ~ "Hif1a",
  is_hif2a            ~ "Hif2a",
  TRUE                ~ NA_character_
)

master_table <- master_table %>% dplyr::relocate(Target, .after = symbol)

# Explicit labels + level order, so NA does not surface as the string "NA." and
# the categories read in a sensible order rather than alphabetically.
target_f <- factor(ifelse(is.na(master_table$Target), "no target", master_table$Target),
                   levels = c("Hif1a", "Hif2a", "Hif1a_Hif2a", "no target"))

tab(as.data.frame(table(Target = target_f)),
    sprintf("Target classification (%s genes Hx-induced in WT Kelly)",
            format(sum(hx_up), big.mark = ",")))
```

</details>

**Target classification (5,882 genes Hx-induced in WT Kelly)**

| Target      |  Freq |
|:------------|------:|
| Hif1a       |   322 |
| Hif2a       |  1472 |
| Hif1a_Hif2a |    43 |
| no target   | 62700 |

# Hypoxia baseline (SK-N-AS / SH-SY5Y)

SK-N-AS and SH-SY5Y are both **MYCN-non-amplified** neuroblastoma lines
and serve as the baseline hypoxia response; Kelly is MYCN-amplified. The
same up-criterion as for `Target` is applied, so both columns are
directly comparable.

| Value | Meaning |
|----|----|
| `SKNAS_SHSY5Y` | Hx-induced in both baseline lines |
| `SKNAS` / `SHSY5Y` | Hx-induced in that line only |
| `none` | testable in at least one line, but not Hx-induced - the Kelly-specific candidates |
| `NA` | testable in neither line - unknown, **not** “no response” |

Keeping `NA` separate from `none` matters: without it, genes that were
never measured in the reference lines would masquerade as
Kelly-specific.

<details>

<summary>

Code
</summary>

``` r
up_sknas  <- sig_up(master_table$log2FC_SKNAS,  master_table$padj_SKNAS)
up_shsy5y <- sig_up(master_table$log2FC_SHSY5Y, master_table$padj_SHSY5Y)

testable <- !is.na(master_table$padj_SKNAS) | !is.na(master_table$padj_SHSY5Y)

master_table$Hx_Baseline <- dplyr::case_when(
  up_sknas & up_shsy5y ~ "SKNAS_SHSY5Y",
  up_sknas             ~ "SKNAS",
  up_shsy5y            ~ "SHSY5Y",
  testable             ~ "none",
  TRUE                 ~ NA_character_
)

master_table <- master_table %>% dplyr::relocate(Hx_Baseline, .after = Target)

baseline_f <- factor(ifelse(is.na(master_table$Hx_Baseline), "no data",
                            master_table$Hx_Baseline),
                     levels = c("none", "SKNAS", "SHSY5Y", "SKNAS_SHSY5Y", "no data"))

tab(as.data.frame(table(Hx_Baseline = baseline_f)), "Baseline hypoxia response")
```

</details>

**Baseline hypoxia response**

| Hx_Baseline  |  Freq |
|:-------------|------:|
| none         | 20857 |
| SKNAS        |  1143 |
| SHSY5Y       |   565 |
| SKNAS_SHSY5Y |   437 |
| no data      | 41535 |

## Target vs. baseline

Rows are the HIF targets from Kelly, columns their behaviour in the
MYCN-non-amplified lines. The `none` column holds the Kelly-specific
targets, the `NA` column the ones with no baseline data.

<details>

<summary>

Code
</summary>

``` r
cross <- table(Target = target_f, Baseline = baseline_f)

tab(as.data.frame.matrix(cross),
    sprintf("Target x Hx_Baseline - %s HIF targets are Kelly-specific",
            format(sum(!is.na(master_table$Target) &
                       master_table$Hx_Baseline == "none", na.rm = TRUE),
                   big.mark = ",")))
```

</details>

**Target x Hx_Baseline - 922 HIF targets are Kelly-specific**

|             |  none | SKNAS | SHSY5Y | SKNAS_SHSY5Y | no data |
|:------------|------:|------:|-------:|-------------:|--------:|
| Hif1a       |    82 |    33 |     63 |           81 |      63 |
| Hif2a       |   823 |   141 |    107 |           62 |     339 |
| Hif1a_Hif2a |    17 |     6 |      9 |            3 |       8 |
| no target   | 19935 |   963 |    386 |          291 |   41125 |

## Venn: Target vs. Hx_Baseline

Left: both reference lines pooled into one set. Right: the same
comparison with SK-N-AS and SH-SY5Y resolved separately.

<details>

<summary>

Code
</summary>

``` r
# Restrict to genes that are testable in at least one baseline line. Otherwise
# genes with no baseline data would sit outside the SKNAS/SHSY5Y circles and
# would be read as "Kelly-specific" when they are simply unmeasured.
venn_universe <- master_table %>% dplyr::filter(!is.na(Hx_Baseline))

ens_where <- function(keep) venn_universe$Ensembl[keep]

venn_sets <- list(
  Hif1a  = ens_where(venn_universe$Target %in% c("Hif1a", "Hif1a_Hif2a")),
  Hif2a  = ens_where(venn_universe$Target %in% c("Hif2a", "Hif1a_Hif2a")),
  SKNAS  = ens_where(venn_universe$Hx_Baseline %in% c("SKNAS", "SKNAS_SHSY5Y")),
  SHSY5Y = ens_where(venn_universe$Hx_Baseline %in% c("SHSY5Y", "SKNAS_SHSY5Y"))
)

venn_plot <- function(input_list, main, fill) {
  venn.diagram(
    x = input_list,
    category.names = paste0(names(input_list), "\n(", lengths(input_list), ")"),
    force.unique = TRUE, na = "remove",
    filename = NULL,
    main = main, main.fontface = "bold",
    lwd = 2,
    lty = "blank",
    fill = fill,
    cat.fontface = "bold",
    cex = 0.8, cat.cex = 0.8
  )
}

# Simplified view: both baseline lines pooled into one "Hx in baseline" set.
venn_simple <- list(
  Hif1a = venn_sets$Hif1a,
  Hif2a = venn_sets$Hif2a,
  Baseline_Hx = union(venn_sets$SKNAS, venn_sets$SHSY5Y)
)

plt_simple <- venn_plot(venn_simple, "HIF targets vs. baseline Hx response",
                        colors[c(4, 6, 1)])
plt_full   <- venn_plot(venn_sets, "HIF targets vs. SK-N-AS / SH-SY5Y",
                        colors[c(4, 6, 7, 3)])

wrap_elements(plt_simple) + wrap_elements(plt_full)
```

</details>

![](Readme_files/figure-gfm/target_baseline_venn-1.png)<!-- -->

*Venn universe: 23,002 genes testable in at least one baseline line,
41,535 genes excluded for lack of baseline data.*

# Literature coverage

Two columns describing how well a gene is already studied. Both come
from NCBI’s curated `gene2pubmed` mapping (Entrez GeneID -\> PubMed ID)
rather than from per-gene PubMed text searches: it keys on the `ENTREZ`
column we already have, so ambiguous symbols (`CS`, `AR`, `MARCH1`,
`SEPT9`) cannot mismatch, and it needs one download instead of 64,537
API calls.

- `pubmed_total` - curated publications for the gene
- `pubmed_hypoxia` - of those, publications in the hypoxia / HIF
  literature

The second column is the useful one here. `pubmed_total` mostly measures
how long a gene has been known; `pubmed_hypoxia` answers whether *this*
biology is already described. A strongly induced, Kelly-specific target
with many publications but none in hypoxia is the interesting case.

Genes without an Entrez ID stay `NA` - unknown coverage, **not** zero.
Only genes that have an Entrez ID and no linked publication get a real
`0`.

<details>

<summary>

Code
</summary>

``` r
cache_dir <- here::here("6_XM_data_merge", "cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

g2p_file <- file.path(cache_dir, "gene2pubmed.gz")
pmid_file <- file.path(cache_dir, "hypoxia_pmids.rds")

if (!file.exists(g2p_file)) {
  download.file("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz",
                g2p_file, mode = "wb", quiet = TRUE)
}

# Human subset only; the full file is ~65M rows, so filter while decompressing.
g2p <- data.table::fread(
  cmd = sprintf("gzcat %s | awk -F'\t' '$1==9606'", shQuote(g2p_file)),
  header = FALSE, col.names = c("tax_id", "GeneID", "PubMed_ID"))

hypoxia_query <- paste0(
  'hypoxia[Title/Abstract] OR "hypoxia-inducible factor"[Title/Abstract] ',
  'OR HIF-1[Title/Abstract] OR HIF-2[Title/Abstract] OR hypoxia[MeSH Terms]')

if (file.exists(pmid_file)) {
  hyp_ids <- readRDS(pmid_file)
} else {
  # esearch and efetch both cap at 10,000 records per request, so the result set
  # is split by publication year, and years above the cap by month.
  qenc <- URLencode(paste0("(", hypoxia_query, ")"), reserved = TRUE)
  esearch <- function(mind, maxd, retmax) {
    u <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                "?db=pubmed&retmode=json&retmax=", retmax, "&term=", qenc,
                "&datetype=pdat&mindate=", mind, "&maxdate=", maxd)
    for (att in 1:4) {
      j <- tryCatch(jsonlite::fromJSON(u), error = function(e) NULL)
      if (!is.null(j) && !is.null(j$esearchresult$count)) return(j$esearchresult)
      Sys.sleep(1.5)
    }
    stop("esearch failed for ", mind, "-", maxd)
  }
  window_ids <- function(mind, maxd) {
    n <- as.integer(esearch(mind, maxd, 0)$count)
    if (n == 0) return(character(0))
    if (n > 9500) return(NULL)
    Sys.sleep(0.36)
    as.character(esearch(mind, maxd, 10000)$idlist)
  }
  hyp_ids <- character(0)
  for (y in 1945:as.integer(format(Sys.Date(), "%Y"))) {
    Sys.sleep(0.36)
    got <- window_ids(y, y)
    if (is.null(got)) {
      for (m in 1:12) {
        Sys.sleep(0.36)
        md <- sprintf("%d/%02d", y, m)
        got_m <- window_ids(md, md)
        if (is.null(got_m)) stop("month window still above the cap: ", md)
        hyp_ids <- c(hyp_ids, got_m)
      }
    } else {
      hyp_ids <- c(hyp_ids, got)
    }
  }
  hyp_ids <- unique(as.integer(hyp_ids))
  saveRDS(hyp_ids, pmid_file)
}

# openxlsx writes NA_character_ as the literal string "NA", so a table that has
# been through the Excel round trip carries text where a missing value belongs.
# Normalise before anything keys off is.na().
entrez <- master_table$ENTREZ
entrez[entrez %in% c("NA", "")] <- NA_character_
has_entrez <- !is.na(entrez)

# ENTREZ is ";"-joined where a gene maps to several Entrez IDs - expand to one row
# per (gene, GeneID) and take the UNION of their PMIDs, so nothing is counted twice.
map <- data.table::data.table(row_id = seq_len(nrow(master_table)), ENTREZ = entrez)
map <- map[!is.na(ENTREZ)]
map <- map[, .(GeneID = as.integer(unlist(strsplit(ENTREZ, ";")))), by = row_id]
map <- map[!is.na(GeneID)]

hit <- merge(map, g2p[, .(GeneID, PubMed_ID)], by = "GeneID", allow.cartesian = TRUE)
hit <- unique(hit[, .(row_id, PubMed_ID)])
hit[, is_hyp := PubMed_ID %in% hyp_ids]
agg <- hit[, .(total = .N, hypoxia = sum(is_hyp)), by = row_id]

master_table$pubmed_total   <- ifelse(has_entrez, 0L, NA_integer_)
master_table$pubmed_hypoxia <- ifelse(has_entrez, 0L, NA_integer_)
master_table$pubmed_total[agg$row_id]   <- agg$total
master_table$pubmed_hypoxia[agg$row_id] <- agg$hypoxia

master_table <- master_table %>% dplyr::relocate(pubmed_total, pubmed_hypoxia,
                                                 .after = Hx_Baseline)

tab(data.frame(
  Group = c("no Entrez ID (NA - unknown)", "Entrez ID, no publication",
            "at least 1 publication", "at least 1 hypoxia publication"),
  Genes = c(sum(!has_entrez),
            sum(master_table$pubmed_total == 0, na.rm = TRUE),
            sum(master_table$pubmed_total > 0, na.rm = TRUE),
            sum(master_table$pubmed_hypoxia > 0, na.rm = TRUE))
), sprintf("Literature coverage (%s hypoxia PMIDs, %s curated human gene-paper links)",
           format(length(hyp_ids), big.mark = ","), format(nrow(g2p), big.mark = ",")))
```

</details>

**Literature coverage (210,214 hypoxia PMIDs, 3,322,637 curated human
gene-paper links)**

| Group                          | Genes |
|:-------------------------------|------:|
| no Entrez ID (NA - unknown)    | 35994 |
| Entrez ID, no publication      |  1993 |
| at least 1 publication         | 26550 |
| at least 1 hypoxia publication |  7840 |

Sanity check - the genes with the most hypoxia literature are the ones
you would expect, which is what makes the column trustworthy:

**Best covered genes in the hypoxia literature**

| symbol | pubmed_total | pubmed_hypoxia |
|:-------|-------------:|---------------:|
| HIF1A  |         9767 |           8254 |
| VEGFA  |        16686 |           1675 |
| VHL    |         2236 |            827 |
| CA9    |         1333 |            605 |
| EPAS1  |          638 |            532 |
| MTOR   |        13389 |            518 |

## Candidates: strong, Kelly-specific, undescribed in hypoxia

HIF targets that respond only in Kelly, are strongly induced, and are
established genes (at least 20 publications) with **no** hypoxia
literature at all.

**Strongly induced Kelly-specific targets with no hypoxia literature**

| symbol | Target | pubmed_total | pubmed_hypoxia | log2FC |
|:-------|:-------|-------------:|---------------:|-------:|
| ATP8A2 | Hif2a  |           45 |              0 |   8.69 |
| LRAT   | Hif1a  |           56 |              0 |   8.33 |
| FOXC2  | Hif2a  |          170 |              0 |   7.29 |
| HTR4   | Hif2a  |          107 |              0 |   7.23 |
| SASH1  | Hif2a  |          112 |              0 |   7.13 |
| TSPAN8 | Hif2a  |          141 |              0 |   7.11 |
| NFE2L3 | Hif2a  |           76 |              0 |   7.09 |
| MGAT4C | Hif2a  |           25 |              0 |   6.42 |
| GLI3   | Hif2a  |          457 |              0 |   6.29 |
| C4B    | Hif2a  |          216 |              0 |   6.28 |
| AKR1D1 | Hif2a  |           48 |              0 |   6.27 |
| HERC6  | Hif2a  |           29 |              0 |   6.17 |

# Master counts table

The same table plus the per-sample normalized counts of the Kelly
experiment. Sample columns are relabelled
`counts_genotype_condition_replicate`, so the `counts_` prefix separates
them from the statistics columns in one selection.

<details>

<summary>

Code
</summary>

``` r
sample_labels <- colData(dds_kelly) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::group_by(genotype, condition) %>%
  dplyr::mutate(replicate = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(label = paste("counts", genotype, condition, replicate, sep = "_"))

norm_counts <- counts(dds_kelly, normalized = TRUE)
colnames(norm_counts) <- sample_labels$label[match(colnames(norm_counts), sample_labels$sample_id)]

counts_df <- norm_counts %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Ensembl")

master_counts_table <- master_table %>% left_join(counts_df, by = "Ensembl")

count_cols <- grep("^counts_", colnames(master_counts_table), value = TRUE)

master_counts_table %>%
  dplyr::select(symbol, Target, Hx_Baseline, dplyr::all_of(head(count_cols, 4))) %>%
  head(5) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 1))) %>%
  tab(sprintf("master_counts_table preview - %d columns total, %d of them counts",
              ncol(master_counts_table), length(count_cols)))
```

</details>

**master_counts_table preview - 144 columns total, 88 of them counts**

| symbol | Target | Hx_Baseline | counts_WT_Nx_1 | counts_WT_Hx_1 | counts_HIF1a-KO_Hx_1 | counts_HIF2a-KO_Hx_1 |
|:---|:---|:---|---:|---:|---:|---:|
| TSPAN6 | NA | none | 1485.6 | 1473.7 | 1382.9 | 1313.8 |
| TNMD | NA | NA | 0.0 | 0.0 | 0.0 | 0.0 |
| DPM1 | NA | none | 1694.3 | 1577.5 | 1491.4 | 1361.1 |
| SCYL3 | NA | none | 821.4 | 847.9 | 1064.3 | 692.4 |
| C1orf112 | NA | none | 1454.2 | 715.2 | 727.6 | 874.8 |

# Export

<details>

<summary>

Code
</summary>

``` r
# Absolute paths via here(), so the exports always land in 6_XM_data_merge/ - knit
# sets the wd to this Rmd's folder, but running the chunks from the console uses
# the project root, which previously produced a second copy of both files there.
out_dir <- here::here("6_XM_data_merge")

openxlsx::write.xlsx(master_table,
                     file.path(out_dir, "master_table.xlsx"), rowNames = FALSE)
openxlsx::write.xlsx(master_counts_table,
                     file.path(out_dir, "master_counts_table.xlsx"), rowNames = FALSE)
```

</details>

- [master_table.xlsx](master_table.xlsx) - 64,537 genes x 56 columns,
  27.4 MB (statistics only)
- [master_counts_table.xlsx](master_counts_table.xlsx) - 64,537 genes x
  144 columns, 64.7 MB (statistics + per-sample normalized counts)
