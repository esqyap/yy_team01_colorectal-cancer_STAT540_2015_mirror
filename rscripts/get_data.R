#####################################################
# Beryl Zhuang
#
# Download GSE48684 dataset and save the data
# and metadata as .Rdata
#
# Output: GSE48684_full.Rdata
#         methyl_data: contain beta value for all probes (no filter)
#         methyl_metadata: metadata for all samples
#####################################################

library(GEOquery)

GSE48684 <- getGEO("GSE48684", destdir="../data") # get data save to "data" dir # don't commit this data, it's large
show(GSE48684)


#get the methylation beta value for the full data
methyl_data <- exprs(GSE48684)

#get the metadata
methyl_metadata <- pData(phenoData(GSE48684))


#save Rdata
save(methyl_data, methyl_metadata, file = "../data/GSE48684_full.Rdata")


#####################################################
# End of script
#####################################################
