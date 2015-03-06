library(GEOquery)
library(wateRmelon)
library(IlluminaHumanMethylation450k.db)

GSE48684 <- getGEO("GSE48684", destdir="./") # get data save to home dir
show(GSE48684)


#get the methylation beta value # only the first 6 rows
methyl_data <- head(exprs(GSE48684))

#get the metadata
methyl_metadata <- pData(phenoData(GSE48684))


#save Rdata
save(methyl_data, methyl_metadata, file = "GSE48684_head.Rdata")

