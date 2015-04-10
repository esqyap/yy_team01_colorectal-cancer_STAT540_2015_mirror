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

#df  <- data.frame(ENSEMBL_GENE_ID)
df <- data.frame(attr(ENSEMBL_GENE_ID, "row.names"), check.rows= T)

# save chr in .RData
save(ENSEMBL_GENE_ID, file = "../data/ENSEMBL_GENE_ID.RData")

Ensembl <- unlist(ENSEMBL_GENE_ID)

write(ENSEMBL_GENE_ID, file='../data/ENSEMBL_GENE_ID.tsv', sep='\t')
