#####################################################
# Beryl Zhuang
#
# Download GSE48684 dataset and save the data
# and metadata as .Rdata
#
# Output: "../data/GSE48684_raw.Rdata", and "../data/GSE48684_metadata.Rdata"
#         methyl_data: contain beta value for all probes (no filter)
#         methyl_metadata: metadata for all samples
#####################################################

library(GEOquery)


filename <- "../data/GSE48684_series_matrix.txt.gz"


if(file.exists(filename)){
  GSE48684 <- getGEO(filename = filename)  # load previousely downloaded data
  show(GSE48684)
} else {
  # get data save to "data" dir # don't commit this data, it's large
  # it should download GSE48684_series_matrix.txt.gz, and GPL13534.soft
  # loading takes some time
  GSE48684 <- getGEO("GSE48684", destdir="../data")
  show(GSE48684)
}


#get the methylation beta value for the full data
methyl_data_raw <- exprs(GSE48684)

#get the metadata
methyl_metadata <- pData(phenoData(GSE48684))


#save Rdata
save(methyl_data_raw, file = "../data/GSE48684_raw.Rdata")  # large file, dont push
save(methyl_metadata, file = "../data/GSE48684_metadata.Rdata")

#####################################################
# End of script
#####################################################
