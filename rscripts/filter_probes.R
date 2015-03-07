#####################################################
# Ka Ming Nip and Beryl Zhuang
#
#
# Output: raw_data_filter ("../data/GSE48684_raw_filtered.Rdata")  
#         - contain raw beta value for cgi and non chrx probes
#####################################################

library(GEOquery)
library(IlluminaHumanMethylation450k.db)



#####################################################
# Ka Ming
#
# Extracting probe ID to CpG islands association
#####################################################

# NOTE: I think this includes shelf, shore, and island
# @TODO: check how to keep island section only
cginame <- as.data.frame(IlluminaHumanMethylation450kCPGINAME)
names(cginame) <- c('Probe_ID', 'cginame')
rownames(cginame) <- cginame$Probe_ID
length(levels(factor(cginame$cginame)))   # No. of CGIs 27176

# Exclude probe ID in chrX regions
chrx <- grep("^chrX:", cginame$cginame, ignore.case=TRUE)
cginame <- cginame[-chrx, ]
length(levels(factor(cginame$cginame)))   # No. of CGIs 26403

#save(cginame, file = "../data/cgi_non_chrx_probes.Rdata")

#####################################################
# Beryl Zhuang
#
# load the raw methyl data "../data/GSE48684_raw.Rdata" 
# (output from get_data.R) 
# 
# save the "../data/GSE48684_raw_filtered.Rdata" for
# probes of cgi and non chrx raw methyl data
#####################################################
filename  <- "../data/GSE48684_raw.Rdata"

if(!file.exists(filename)){
  print("you haven't downloaded the raw data yet! run get_data.R")
} else {
  load(filename)
  raw_data <- methyl_data_raw
  # filter the raw data by cgi and non chrX
  raw_data_filter <- raw_data[cginame$Probe_ID, ]
  # save the raw filtered methyl data for normalization, clustering
  save(raw_data_filter, file = "../data/GSE48684_raw_filtered.Rdata")
}


#####################################################
# End of script
#####################################################



