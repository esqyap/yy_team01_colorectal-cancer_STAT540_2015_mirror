#####################################################
# Beryl Zhuang
#
# process the raw metadata, select relevant columns
# Output: "../data/metadata.Rdata": metadata
#
#####################################################

library(dplyr)
filename <- "../data/GSE48684_metadata_raw.Rdata"

if(!file.exists(filename)){
  print("you haven't downloaded the raw metadata yet! pull Git or run get_data.R")
} else {
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
}

#####################################################
# End of script
#####################################################