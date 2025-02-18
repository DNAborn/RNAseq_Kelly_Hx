# use all genes in dataset
genes_univ_ens <- results_list[[1]]$ENSEMBL
genes_univ_ens %>% length()

genes_univ <- results_list[[1]]$ENTREZ
genes_univ %>% length()
genes_univ %>% is.na() %>% summary()

# Folder to save
go_dir <- "/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/git_RNAseq_Kelly_Hx/4_manuscript/GO_analysis"

# Enrich GO terms (SK-groups)
cc_go <- compareCluster(geneCluster = res_hif1a_2a_list_ens,
                     fun = "enrichGO",
                     universe = genes_univ_ens,
                     OrgDb = "org.Hs.eg.db",
                     keyType = "ENSEMBL",
                     ont = "ALL",
                     # pvalueCutoff = 0.05,
                     # minGSSize = 10,
                     maxGSSize = 1000,
                     pAdjustMethod = "fdr")
cc_go <- simplify(cc_go)
# cc_go <- setReadable(cc_go, OrgDb = org.Hs.eg.db)
GO_cc_groups <- mutate(cc_go, l2FoldEnrichment = -log2(parse_ratio(GeneRatio) / parse_ratio(BgRatio))) %>%
  arrange(desc(l2FoldEnrichment))

go_dir <- "/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/git_RNAseq_Kelly_Hx/4_manuscript/GO_analysis"
save(GO_cc_groups, file = paste0(go_dir,"/GO_cc_groups.go"))



# Enrich GO terms (HS grous)
cc_go <- compareCluster(geneCluster = genes_holger_list[c("HIF1a","HIF2a","HIF1a_HIF2a")],
                        fun = "enrichGO",
                        universe = genes_univ_ens,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "ENSEMBL",
                        ont = "ALL",
                        # pvalueCutoff = 0.05,
                        # minGSSize = 10,
                        maxGSSize = 1000,
                        pAdjustMethod = "fdr")
cc_go <- simplify(cc_go)
# cc_go <- setReadable(cc_go, OrgDb = org.Hs.eg.db)
GO_cc_groups <- mutate(cc_go, l2FoldEnrichment = -log2(parse_ratio(GeneRatio) / parse_ratio(BgRatio))) %>%
  arrange(desc(l2FoldEnrichment))

go_dir <- "/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/git_RNAseq_Kelly_Hx/4_manuscript/GO_analysis"
save(GO_cc_groups, file = paste0(go_dir,"/GO_cc_groups_hs.go"))


# Enrich GO terms (HS grous)
cc_go <- compareCluster(geneCluster = c(res_hif1a_2a_list_ens[1:2],genes_holger_list[c("HIF1a","HIF2a")]),
                        fun = "enrichGO",
                        universe = genes_univ_ens,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "ENSEMBL",
                        ont = "ALL",
                        # pvalueCutoff = 0.05,
                        # minGSSize = 10,
                        maxGSSize = 1000,
                        pAdjustMethod = "fdr")
cc_go <- simplify(cc_go)
# cc_go <- setReadable(cc_go, OrgDb = org.Hs.eg.db)
GO_cc_groups <- mutate(cc_go, l2FoldEnrichment = -log2(parse_ratio(GeneRatio) / parse_ratio(BgRatio))) %>%
  arrange(desc(l2FoldEnrichment))

go_dir <- "/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/git_RNAseq_Kelly_Hx/4_manuscript/GO_analysis"
save(GO_cc_groups, file = paste0(go_dir,"/GO_cc_groups_both.go"))




# KEGG

search_kegg_organism('sapiens', by='scientific_name')



cc_kegg <- compareCluster(geneCluster = res_hif1a_2a_list_ez[1:2],
                          fun = "enrichKEGG",
                          organism     = 'hsa',
                          pvalueCutoff = 0.05)



