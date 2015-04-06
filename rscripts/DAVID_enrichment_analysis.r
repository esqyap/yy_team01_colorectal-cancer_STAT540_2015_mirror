library(IlluminaHumanMethylation450k.db)
load("../data/chrlen.RData")
load("../data/coord.RData")

#test data
load("../data//limma/normH_normC_dma.Rdata")
DMR1.1  <- subset(normH_normC_dma, normH_normC_dma$P.Value < 1e-5)
str(DMR1.1)
#test data end


#Get ENSEMBL_GENE_ID from 'IlluminaHumanMethylation450kENTREZID'

x <- IlluminaHumanMethylation450kENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(x)
# Convert to a list
ENSEMBL_GENE_ID <- as.list(x[mapped_genes])
head(ENSEMBL_GENE_ID)
str(ENSEMBL_GENE_ID)

#subset  <- subset(ENSEMBL_GENE_ID, nrow= 100)
#df  <- data.frame(ENSEMBL_GENE_ID)
df <- data.frame(attr(ENSEMBL_GENE_ID, "row.names"), check.rows= T)

#ENSEMBL_GENE_ID  <- as.character(ENSEMBL_GENE_ID)

#ENSEMBL_GENE_ID  <- as.data.frame(ENSEMBL_GENE_ID)

# save chr in .RData
save(ENSEMBL_GENE_ID, file = "../data/ENSEMBL_GENE_ID.RData")

