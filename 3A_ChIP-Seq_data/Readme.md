RNA-Seq Kelly Hx 3A ChIP-Seq
================
Kelterborn
2024-06-19

- [0. Load](#0-load)
  - [- Load R librarys](#--load-r-librarys)
  - [- Load dds](#--load-dds)
  - [- functions](#--functions)
- [ChIP-Seq datasets](#chip-seq-datasets)

# 0. Load

## - Load R librarys

## - Load dds

## - functions

# ChIP-Seq datasets

Literature

<table style="width:99%;">
<colgroup>
<col style="width: 12%" />
<col style="width: 20%" />
<col style="width: 9%" />
<col style="width: 55%" />
</colgroup>
<thead>
<tr class="header">
<th>Author</th>
<th>cells</th>
<th>Ab</th>
<th>link</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><strong>Schödel et al., 2011</strong></td>
<td>MCF-7 (breast)</td>
<td>Hif1A, Hif2A, Hif1B</td>
<td><a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3374576/"
class="uri">https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3374576/</a></td>
</tr>
<tr class="even">
<td><strong>Andrysik et al., 2021</strong></td>
<td>HCT116 (colon), RKO (colon), A549, and H460</td>
<td>HIF1A</td>
<td><a href="https://www.nature.com/articles/s41467-021-21687-2#Sec11"
class="uri">https://www.nature.com/articles/s41467-021-21687-2#Sec11</a>
<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145157"
class="uri">https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145157</a></td>
</tr>
<tr class="odd">
<td><strong>James A Smythies</strong></td>
<td>HKC-8, RCC4, HepG2</td>
<td>HIF1A, HIF2A, HIF1B</td>
<td><p><a
href="https://www.embopress.org/doi/pdf/10.15252/embr.201846401"
class="uri">https://www.embopress.org/doi/pdf/10.15252/embr.201846401</a>
GSE120885, GSE120886 and GSE120887</p>
<p><a
href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120887"
class="uri">https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120887</a></p></td>
</tr>
</tbody>
</table>

**ReMAP:**

Hif1a: <https://remap.univ-amu.fr/target_page/HIF1A:9606>

Hif2a: <https://remap.univ-amu.fr/target_page/EPAS1:9606>

Hif1b: <https://remap.univ-amu.fr/target_page/ARNT:9606>

``` r
list.files(pdir)
```

    ## [1] "3A_ChIP-Seq_data_bu.Rmd"      "3A_ChIP-Seq_data.Rmd"        
    ## [3] "Readme.md"                    "ReMAP_ChIP_Hif1a.xlsx"       
    ## [5] "ReMAP_ChIP_Hif1b.xlsx"        "ReMAP_ChIP_Hif2a.xlsx"       
    ## [7] "ReMAP_ChIP-Seq_datasets.xlsx"

``` r
hif1a_datasets <- read_xlsx(paste(pdir,"ReMAP_ChIP_Hif1a.xlsx",sep="/"), )
hif2a_datasets <- read_xlsx(paste(pdir,"ReMAP_ChIP_Hif2a.xlsx",sep="/"), )
hif1b_datasets <- read_xlsx(paste(pdir,"ReMAP_ChIP_Hif1b.xlsx",sep="/"), )

celllines <- data.frame("Hif1A" = (hif1a_datasets$Biotype %>% factor() %>% levels() %>% paste(collapse = " ")),
                       "Hif2A" = (hif2a_datasets$Biotype %>% factor() %>% levels() %>% paste(collapse = " ")),
                       "Hif1B" = (hif1b_datasets$Biotype %>% factor() %>% levels() %>% paste(collapse = " "))) %>% t()

celllines %>% kable()
```

|       |                                                                                                 |
|:------|:------------------------------------------------------------------------------------------------|
| Hif1A | 501-mel 786-O BEAS-2B ccRCC HUVEC-C K-562 LNCaP macrophage MDA-MB-231 NCI-H1299 PC-3 RCC10 U2OS |
| Hif2A | 501-mel 786-O ccRCC HUVEC-C K-562 macrophage PC-3                                               |
| Hif1B | 501-mel A-549 MCF-7 PC-3 RCC10 RCC4 SK-MEL-28 T-47D                                             |
