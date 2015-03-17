#####################################################
# Ka Ming Nip
# 
# Script to:
#   1. normalize filtered beta-values
#   2. convert beta-values to M-values
#   3. aggregate both beta-values and M-values to CpG islands
# 
# Input Files: (1)
#   raw_data_filter : "../data/GSE48684_raw_filtered.Rdata"
# 
# Output Files: (4)
#   beta.norm       : "../data/GSE48684_raw_filtered.beta.norm.Rdata"
#   beta.norm.CGI   : "../data/GSE48684_raw_filtered.beta.norm.cgi.Rdata"
#   M.norm          : "../data/GSE48684_raw_filtered.m.norm.Rdata"
#   M.norm.CGI      : "../data/GSE48684_raw_filtered.m.norm.cgi.Rdata"
#   
#####################################################

library(wateRmelon)
library(IlluminaHumanMethylation450k.db)

myGseId <- "GSE48684"
dataDir <- "../data" # our working dir is `./rscripts`

# input file
filteredFilePath <- paste(dataDir, "GSE48684_raw_filtered.Rdata", sep="/")

# output files
betaNormFilePath <- paste(dataDir, "GSE48684_raw_filtered.beta.norm.Rdata", sep="/")
betaNormCgiFilePath <- paste(dataDir, "GSE48684_raw_filtered.beta.norm.cgi.Rdata", sep="/")
mNormFilePath <- paste(dataDir, "GSE48684_raw_filtered.m.norm.Rdata", sep="/")
mNormCgiFilePath <- paste(dataDir, "GSE48684_raw_filtered.m.norm.cgi.Rdata", sep="/")

if (!file.exists(filteredFilePath)) {
  # can't find input file, quit with error message
  stop(paste0("ERROR: Cannot find: ", filteredFilePath))  
}

# load input file
load(filteredFilePath)

# assert more than one rows are loaded
stopifnot(nrow(raw_data_filter) > 0)

# to flag whether output files are stale or not
areOutputFilesStale <- FALSE

# If a file is designated as "stale", all dependent output files are also stale!
# From now on, re-create any stale output file(s).

# quantile-normalization on the beta values
if (!areOutputFilesStale && file.exists(betaNormFilePath)) {
  load(betaNormFilePath)
} else {
  areOutputFilesStale <- TRUE
  
  beta.norm <- betaqn(as.matrix(raw_data_filter)[cginame$Probe_ID, ])
  str(beta.norm, max.level = 0)
  save(beta.norm, file = betaNormFilePath)
}

# convert beta-values to M-values
if (!areOutputFilesStale && file.exists(mNormFilePath)) {
  load(mNormFilePath)
} else {
  areOutputFilesStale <- TRUE
    
  M.norm <- beta2m(beta.norm)
  save(M.norm, file = mNormFilePath)
}

# aggregate per-probe beta-values to CGIs
if (!areOutputFilesStale && file.exists(betaNormCgiFilePath)){
  load(betaNormCgiFilePath)
} else {
  areOutputFilesStale <- TRUE
    
  beta.norm.CGI <- aggregate(beta.norm, by = list(cginame$cginame), mean, na.rm = T)
  rownames(beta.norm.CGI) <- beta.norm.CGI[, "Group.1"]
  beta.norm.CGI <- subset(beta.norm.CGI, select = - Group.1)
  str(beta.norm.CGI, max.level = 0)
  save(beta.norm.CGI, file = betaNormCgiFilePath)
}

# aggregate per-probe M-values to CGIs
if (!areOutputFilesStale && file.exists(mNormCgiFilePath)){
  load(mNormCgiFilePath)
} else {
  areOutputFilesStale <- TRUE
  
  M.norm.CGI <- aggregate(M.norm, by = list(cginame$cginame), mean, na.rm = T)
  rownames(M.norm.CGI) <- M.norm.CGI[, "Group.1"]
  M.norm.CGI <- subset(M.norm.CGI, select = - Group.1)
  str(M.norm.CGI, max.level = 0)
  save(M.norm.CGI, file = mNormCgiFilePath)
}

# assert sizes are the same
stopifnot(nrow(beta.norm) == nrow(M.norm))
stopifnot(nrow(beta.norm.CGI) == nrow(M.norm.CGI))

#####################################################
# End of script
#####################################################