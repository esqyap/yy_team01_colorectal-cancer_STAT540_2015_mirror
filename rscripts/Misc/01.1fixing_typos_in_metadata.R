#####################################################
# Beryl Zhuang
#
# process the raw metadata, select relevant columns
#
#
# Rashedul Islam
#
# Fixing typos: 'right' to 'Right', 'left' to 'Left', 'colon' to 'Colon', 'male' to 'Male', 'feMale' to 'Female'.
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
  
  #Fixing typos in metadata
  #change 'right' to 'Right' 
  metadata$colon_region <- as.character(gsub("right","Right", as.character(metadata$colon_region)))
  #change 'left' to 'Left'
  metadata$colon_region <- as.character(gsub("left","Left", as.character(metadata$colon_region)))
  #change 'colon' to 'Colon'
  metadata$colon_region <- as.character(gsub("colon","Colon", as.character(metadata$colon_region)))
  #change 'male' to 'Male'
  metadata$gender <- as.character(gsub("male","Male", as.character(metadata$gender)))
  #change 'feMale' to 'Female
  metadata$gender <- as.character(gsub("feMale","Female", as.character(metadata$gender)))
  
  save(metadata, file = "../data/metadata.Rdata")
  print("metadata is processed")
}
#####################################################
# End of script
#####################################################
  
  