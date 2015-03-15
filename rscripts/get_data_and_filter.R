#####################################################
# Beryl Zhuang
#
# Download GSE48684 dataset and save the data
# and metadata as .Rdata
#
# Output: methyl_data_raw, file = "../data/GSE48684_raw.Rdata
#					-contain beta value for all probes (no filter)
#					raw_data_filter, file = "../data/GSE48684_raw_filtered.Rdata"
#					-contain beta value for filtered probes (CGI and non chrx)
#         methyl_metadata, file = "../data/GSE48684_metadata.Rdata"
#					-metadata for all samples (RAW)
#					cginame, file = "../data/cgi_non_chrx_probes.Rdata"
#					-probe info of now chrx CGI probes
#					metadata, file = "../data/metadata.Rdata"
#					-processed metadata, for limma and other
#					
#####################################################

library(GEOquery)
if(file.exists("../data/GSE48684_raw.Rdata")){
	Print("raw data downloaded, proceed to the next step.")
} else {
	filename <- "../data/GSE48684_series_matrix.txt.gz"
	if(file.exists(filename)){
	  GSE48684 <- getGEO(filename = filename)  # load previousely downloaded data
	} else {
	  # get data save to "data" dir # don't commit this data, it's large
	  # it should download GSE48684_series_matrix.txt.gz, and GPL13534.soft
	  # loading takes some time
		GSE48684 <- getGEO("GSE48684", destdir="../data")
	}
		#get the methylation beta value for the full data
		methyl_data_raw <- exprs(GSE48684)
		#get the metadata
		methyl_metadata <- pData(phenoData(GSE48684))
		#save Rdata
		save(methyl_data_raw, file = "../data/GSE48684_raw.Rdata")  # large file, don't push
		save(methyl_metadata, file = "../data/GSE48684_metadata.Rdata")
}


#####################################################
# Ka Ming Nip and Beryl Zhuang
#
#
# load the raw methyl data "../data/GSE48684_raw.Rdata" 
# 
# save the "../data/GSE48684_raw_filtered.Rdata" for
# probes of cgi and non chrx raw methyl data
#####################################################
library(IlluminaHumanMethylation450k.db)

filename  <- "../data/GSE48684_raw.Rdata"


if(file.exists("../data/GSE48684_raw_filtered.Rdata")){
	print("filtered raw methylation data saved")
} else {
	
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
	
	save(cginame, file = "../data/cgi_non_chrx_probes.Rdata")
	
	load(filename)
	raw_data <- methyl_data_raw
	# filter the raw data by cgi and non chrX
	raw_data_filter <- raw_data[cginame$Probe_ID, ]
	# save the raw filtered methyl data for normalization, clustering
	save(raw_data_filter, file = "../data/GSE48684_raw_filtered.Rdata")
	print("filtered raw methylation data saved")
}




#####################################################
# Beryl Zhuang
#
# process the raw metadata, select relevant columns
# Output: "../data/metadata.Rdata": metadata
#
#####################################################




if(file.exists("../data/metadata.Rdata")){
	print("metadata processed")
} else {
	library(dplyr)
	filename <- "../data/GSE48684_metadata_raw.Rdata"
	load(filename)
	
	# select the relavent columns form the raw metadata
	metadata <- methyl_metadata_raw[, c("description", "title", 
																			"geo_accession", "source_name_ch1", 
																			"characteristics_ch1.1", "characteristics_ch1.3",
																			"characteristics_ch1.4")]
	
	# change the column names
	colnames(metadata) <- c("group", "title", "geo_accession", "tissue",
													"colon_region", "gender", "stage")
	
	# data cleaning
	metadata <- metadata %>%
		mutate(gender = gsub("gender: ", "", gender)) %>%
		mutate(colon_region = gsub("colon region: ", "", colon_region)) %>%
		mutate(stage = gsub("Stage: ", "", stage))
	
	save(metadata, file = "../data/metadata.Rdata")
	print("metadata is processed")
}

#####################################################
# End of script
#####################################################
