
# libraries
library(tidyverse)
library(tximeta)
library(tximport)
library(AnnotationHub)
library(DESeq2)

# ENSEMBL 111 (TS)
dirsalmon <- "/mnt/s/AG/AG-Scholz-NGS/Daten/Tobias/Salmon_TS/human"
indexdir <- "/mnt/s/AG/AG-Scholz-NGS/Daten/Tobias/Salmon_TS/human/index/human_ensh38_index"
fastaPath <- file.path(dirsalmon, "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
gtfPath <- file.path(dirsalmon,"Homo_sapiens.GRCh38.111.gff3.gz")

file.exists(indexdir, fastaPath, gtfPath)

# Load Sample table
dirtest <- "/mnt/s/AG/AG-Scholz-NGS/Daten/Tobias/Simon_test/"
load(file=paste(dirtest,"sample_table_Kelly.RDS", sep="/"))
coldata

# run tximeta
se <- tximeta(coldata[1:10,])


#################################################################
## compare different linkedTxomes

# clear the entire linkedTxome table 
bfcloc <- getTximetaBFC()
bfc <- BiocFileCache(bfcloc)
show(bfc)
bfcinfo()$fpath
bfcquery(bfc, "Homo")
bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
bfcremove(bfc, bfcquery(bfc, "Homo")$rid)

# Use Hub (online)
se_hub <- tximeta(coldata[1:10,], useHub=T)
se_hub
tx <- dim(se_hub) %>% as.data.frame()
colnames(tx) <- "hub"

# Local use gff3
makeLinkedTxome(indexDir=indexdir,
                source="LocalEnsembl",
                organism="Homo sapiens",
                release="111",
                genome="GRCh38",
                fasta=fastaPath,
                gtf=gtfPath,
                write=FALSE)
bfcquery(bfc, "Homo")
se_local_gff3 <- tximeta(coldata[1:10,], useHub=F)
tx <- bind_cols(tx,dim(se_local_gff3))
colnames(tx)[2] <- "local_gff3"

# Local use gtf
gtfPath <- file.path(dirsalmon,"Homo_sapiens.GRCh38.111.gtf.gz")
file.exists(indexdir, fastaPath, gtfPath)

bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
bfcremove(bfc, bfcquery(bfc, "Homo")$rid)

makeLinkedTxome(indexDir=indexdir,
                source="LocalEnsembl",
                organism="Homo sapiens",
                release="111",
                genome="GRCh38",
                fasta=fastaPath,
                gtf=gtfPath,
                write=FALSE)

se_local_gtf <- tximeta(coldata[1:10,], useHub=F)
bfcquery(bfc, "Homo")
tx <- bind_cols(tx,dim(se_local_gtf))
colnames(tx)[3] <- "local_gtf"
tx
