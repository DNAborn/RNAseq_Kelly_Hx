Master_Table_Dive
================

- [Input](#input)
- [Reduction](#reduction)
- [Reduced table](#reduced-table)

Follow-up analysis on the annotated master table built in
[Master_Table.Rmd](Master_Table.Rmd) (see [its rendered
output](Readme.md) for how the annotation columns are defined). This
document starts by reducing the table to the genes that are actually
worth looking at. Code blocks are collapsed - click *Code* to expand.

<details>

<summary>

Code
</summary>

``` r
library(dplyr)
library(tidyverse)
library(knitr)
library(ggplot2)
library(openxlsx)
```

</details>

# Input

`master_counts_table.xlsx` carries everything: the DESeq2 contrasts, the
annotation columns, and the per-sample normalized counts. Reading 65 MB
of Excel on every knit is slow, so it is cached as an RDS and only
re-read when the xlsx is newer.

<details>

<summary>

Code
</summary>

``` r
in_dir  <- here::here("6_XM_data_merge")
xlsx_in <- file.path(in_dir, "master_counts_table.xlsx")
rds_cache <- file.path(in_dir, "cache", "master_counts_table.rds")

dir.create(dirname(rds_cache), showWarnings = FALSE, recursive = TRUE)
if (file.exists(rds_cache) &&
    file.mtime(rds_cache) > file.mtime(xlsx_in)) {
  master <- readRDS(rds_cache)
} else {
  master <- as.data.frame(readxl::read_xlsx(xlsx_in))
  saveRDS(master, rds_cache)
}

count_cols <- grep("^counts_", colnames(master), value = TRUE)

tab(data.frame(
  Item = c("genes", "columns", "of them per-sample counts", "annotation columns"),
  n    = c(nrow(master), ncol(master), length(count_cols),
           sum(colnames(master) %in% c("Target", "Hx_Baseline", "pubmed_total",
                                       "pubmed_hypoxia", "pubmed_source", "HR_risk")))
), "master_counts_table as loaded")
```

</details>

**master_counts_table as loaded**

| Item                      |     n |
|:--------------------------|------:|
| genes                     | 64537 |
| columns                   |   146 |
| of them per-sample counts |    88 |
| annotation columns        |     6 |

# Reduction

Most of the table is not expressed: 61% of all genes have a
`baseMean_Kelly` of exactly zero. Those rows carry no information about
Kelly biology no matter how well studied the gene is, so expression is
the only filter applied here.

Literature coverage is deliberately **not** used as a cutoff.
`pubmed_total` and `pubmed_hypoxia` stay in the table as a ranking aid -
a strongly expressed but unstudied gene is a candidate, not noise, and
filtering on publication counts would preferentially discard exactly
those.

<details>

<summary>

Code
</summary>

``` r
min_basemean <- 10

expressed <- !is.na(master$baseMean_Kelly) & master$baseMean_Kelly >= min_basemean
dive <- master[expressed, ]

tab(data.frame(
  Step = c("master_counts_table", sprintf("baseMean_Kelly >= %d", min_basemean)),
  Genes = c(nrow(master), nrow(dive)),
  Share = c("100%", sprintf("%.0f%%", 100 * nrow(dive) / nrow(master)))
), "Reduction")
```

</details>

**Reduction**

| Step                  | Genes | Share |
|:----------------------|------:|:------|
| master_counts_table   | 64537 | 100%  |
| baseMean_Kelly \>= 10 | 18471 | 29%   |

What survives, by annotation. The reduction removes three quarters of
the rows but keeps the large majority of the annotated candidates, which
is the point:

**Effect of the reduction per annotation set**

| Set                     | before | after | kept |
|:------------------------|-------:|------:|:-----|
| all genes               |  64537 | 18471 | 29%  |
| HIF targets             |   1837 |  1311 | 71%  |
| Kelly-specific targets  |    922 |   685 | 74%  |
| HR adverse              |   6081 |  4598 | 76%  |
| with hypoxia literature |   7945 |  5911 | 74%  |

The genes dropped despite being annotated are the ones the filter is
meant to catch - annotated but silent in this cell line:

**Best-studied HIF targets removed for lack of expression**

| symbol | Target      | baseMean | pubmed_total | pubmed_hypoxia |
|:-------|:------------|---------:|-------------:|---------------:|
| TLR2   | Hif2a       |     2.09 |         3812 |             22 |
| HFE    | Hif2a       |     1.92 |         2259 |             10 |
| NOD2   | Hif2a       |     2.97 |         1997 |              2 |
| POSTN  | Hif2a       |     2.79 |         1016 |             21 |
| SOST   | Hif2a       |     3.61 |          961 |              8 |
| DCN    | Hif2a       |     5.09 |          948 |              6 |
| ABCA4  | Hif1a_Hif2a |     6.61 |          754 |              1 |
| UGT1A1 | Hif2a       |     1.08 |          753 |              0 |

# Reduced table

**Top Kelly-specific targets in the reduced table**

| symbol      | Target | baseMean | log2FC | HR_risk | pubmed_total | pubmed_hypoxia |
|:------------|:-------|---------:|-------:|:--------|-------------:|---------------:|
| CALHM4      | Hif2a  |      665 |  13.14 | ns      |           13 |              0 |
| FLT1        | Hif2a  |      275 |  12.07 | ns      |         2122 |            151 |
| CHRNA4      | Hif1a  |     1622 |  11.40 | adverse |          218 |              2 |
| IL6R        | Hif2a  |      179 |  10.84 | ns      |          893 |              8 |
| SLC9C2      | Hif2a  |      133 |  10.15 | ns      |            5 |              0 |
| STEAP1B-AS1 | Hif2a  |       26 |   9.06 | NA      |            3 |              0 |
| SLCO2A1     | Hif2a  |      187 |   8.80 | ns      |          211 |              2 |
| FRG2B       | Hif2a  |       46 |   8.80 | NA      |            3 |              0 |
| ATP8A2      | Hif2a  |       51 |   8.69 | ns      |           45 |              0 |
| LINGO1-AS1  | Hif2a  |      115 |   8.68 | NA      |            2 |              0 |

<details>

<summary>

Code
</summary>

``` r
out_file <- here::here("6_XM_data_merge", "master_table_dive.xlsx")
openxlsx::write.xlsx(dive, out_file, rowNames = FALSE)
```

</details>

- [master_table_dive.xlsx](master_table_dive.xlsx) - 18,471 genes x 146
  columns, 31.5 MB
