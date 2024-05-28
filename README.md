# Kelly Hypoxia RNA-Seq

RNA-Seq of Hif1a, Hif2a & Hif1b gene knock-outs

1.) [Data processing](1_data_processing)

### combine sample list

::: {style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:400px; "}
|     | experiment | RNAs | conditions                                            | date                                        | seq_id | Seq_runs |
|:----|:-----------|-----:|:------------------------------------------------------|:--------------------------------------------|:-------|---------:|
| 3   | Katharina  |   16 | Kelly_Nx Kelly_Hx HIF1A_Hx HIF2A_Hx                   | 2018-09-13 2018-09-14                       | P557   |       16 |
| 1   | Simon      |   22 | Kelly_Nx Kelly_Hx HIF1A_Nx HIF1A_Hx HIF1B_Nx HIF1B_Hx | 2017-05-04 2021-06-16 2021-08-25 2021-08-27 | P2041  |       22 |
| 2   | Ulrike     |   50 | Kelly_Nx Kelly_Hx HIF1A_Nx HIF1A_Hx HIF2A_Nx HIF2A_Hx | 2023-06-02 2023-06-08 2023-06-15 2023-06-28 | P3302  |      150 |
:::

### DESeq2 Design

design = \~experiment+genotype+treatment+genotype:treatment

MA plot & Dispersion

<img src="1_data_processing/Readme_files/figure-gfm/dds_design-1.png" width="50%"/><img src="1_data_processing/Readme_files/figure-gfm/dds_design-2.png" width="50%"/>

Transformations

<img src="1_data_processing/Readme_files/figure-gfm/pre_trans_fig, figures-side-1.png" width="33%"/><img src="1_data_processing/Readme_files/figure-gfm/pre_trans_fig, figures-side-2.png" width="33%"/><img src="1_data_processing/Readme_files/figure-gfm/pre_trans_fig, figures-side-3.png" width="33%"/>

Sample distance

<img src="1_data_processing/Readme_files/figure-gfm/pre_sample_dist-1.png" width="100%"/>

Principal component analysis

<img src="1_data_processing/Readme_files/figure-gfm/pca-1.png" width="80%"/>

Plot example counts

<img src="1_data_processing/Readme_files/figure-gfm/example_counts-1.png" width="50%"/><img src="1_data_processing/Readme_files/figure-gfm/example_counts-2.png" width="50%"/>

2.  

    A)  [network analysis](2A_WGCNA)

    B)  [Differential gene expression](2B_DGE)
