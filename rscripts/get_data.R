library(GEOquery)
library(wateRmelon)
library(IlluminaHumanMethylation450k.db)

GSE48684 <- getGEO("GSE48684", destdir="../data") # get data save to "data" dir
show(GSE48684)


#get the methylation beta value # only the first 6 rows
methyl_data <- head(exprs(GSE48684))

#get the full data (commented out)
#methyl_data <- exprs(GSE48684)

#get the metadata
methyl_metadata <- pData(phenoData(GSE48684))


#save Rdata
save(methyl_data, methyl_metadata, file = "../data/GSE48684_head.Rdata")

