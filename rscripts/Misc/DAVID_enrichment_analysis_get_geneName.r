
load("../data//ENSEMBL_GENE_ID.RData")
head(ENSEMBL_GENE_ID)

Ensembl <- unlist(ENSEMBL_GENE_ID)

Ensembl_df  <- data.frame(ENSEMBL_GENE_ID)

toptableData <- read.table("../data/limma//normal_vs_adenoma_santina.tsv")
head(toptable)



#how to get ENSEMBL gene ids from toptable coordinates (chr:start-end)? 

coordDMRprobe <-
  droplevels(na.omit(coord[cginame[cginame$cginame %in%
                                     rownames(DMR),]$Probe_ID,]))