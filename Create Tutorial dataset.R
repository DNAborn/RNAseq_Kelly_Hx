load(file=paste(s,"AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/data/tximeta.txm", sep="/"))

colData(gse)
colData(gse) <- colData(gse)[,c("treatment", "genotype","condition","experiment")]
gse <- gse[,c("P557_Kelly_Nx_25","P2041_Kelly_Nx_33","P3302_Kelly_Nx_119",
  "P557_Kelly_Hx_30","P2041_Kelly_Hx_39","P3302_Kelly_Hx_128",
  "P3302_HIF1A_Nx_64","P3302_HIF1A_Nx_98","P3302_HIF1A_Nx_150",
  "P557_HIF1A_Hx_27","P2041_HIF1A_Hx_41","P3302_HIF1A_Hx_133",
  "P3302_HIF2A_Nx_67","P3302_HIF2A_Nx_151","P3302_HIF2A_Nx_171",
  "P557_HIF2A_Hx_28","P3302_HIF2A_Hx_89","P3302_HIF2A_Hx_203",
  "P2041_HIF1B_Nx_34","P2041_HIF1B_Nx_35","P2041_HIF1B_Nx_45",
  "P2041_HIF1B_Hx_50","P2041_HIF1B_Hx_43","P2041_HIF1B_Hx_48" )]

save(gse, file=paste(s,"AG/AG-Scholz-NGS/Daten/Simon/RNA-Seq_Kelly_all/data/tximeta_tutorial.txm", sep="/"))

