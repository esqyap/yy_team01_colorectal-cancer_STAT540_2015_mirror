# Ka Ming Nip (kmnip)

# Script to:
# 1. download dataset from GEO
# 2. downsize dataset to CpG-associated, non-chrX probes only
# 3. normalize data
# 4. aggregate data to group together probes on the same CpG island

library(GEOquery)
library(wateRmelon)
library(IlluminaHumanMethylation450k.db)
library(ggplot2)
library(dplyr)

myGseId <- 'GSE48684'
rawDumpFileName <- "raw.Rdata"
mBetaNormCgiFileName <- "norm.Rdata"


if(file.exists(rawDumpFileName)){ # if previously downloaded
  load(rawDumpFileName)
} else { # if downloading for the first time
  
  # I had to add this line to bypass a ftp issue. FYI: 
  # https://support.rstudio.com/hc/communities/public/questions/201327633-Bug-report-RStudio-causes-download-file-to-fail-with-FTP-downloads
  setInternet2(use=FALSE)
  
  system.time( myGse <- getGEO(myGseId, destdir="./") )
  system.time( show(myGse) )
  
  # Extract expression matrices (turn into data frames at once) 
  system.time( raw.dat <- as.data.frame(exprs(myGse[[1]])) )
  
  # Obtain the meta-data for the samples
  system.time( raw.meta <- pData(phenoData(myGse[[1]])) )
  
  # save the data to avoid future re-downloading
  save(raw.dat, raw.meta, file = rawDumpFileName)
}

# Extracting probe ID to CpG islands association
# NOTE: I think this includes shelf, shore, and island
# @TODO: check how to keep island section only
cginame <- as.data.frame(IlluminaHumanMethylation450kCPGINAME)
names(cginame) <- c('Probe_ID', 'cginame')
rownames(cginame) <- cginame$Probe_ID
length(levels(factor(cginame$cginame)))   # No. of CGIs

# Exclude probe ID in chrX regions
cginame <- filter(cginame, !grepl("^chrX:", cginame, ignore.case=TRUE))
length(levels(factor(cginame$cginame)))   # No. of CGIs  

# filter probes, then perform quantile normalization
system.time(beta.norm <- betaqn(as.matrix(raw.dat)[cginame$Probe_ID, ]))
str(beta.norm, max.level = 0)

# convert to M values
M.norm <- beta2m(beta.norm)

# aggregate probes to CGIs
beta.CGI <- aggregate(beta.norm, by = list(cginame$cginame), mean, na.rm = T)
rownames(beta.CGI) <- beta.CGI[, "Group.1"]
beta.CGI <- subset(beta.CGI, select = - Group.1)
str(beta.CGI, max.level = 0)

M.CGI <- aggregate(M.norm, by = list(cginame$cginame), mean, na.rm = T)
rownames(M.CGI) <- M.CGI[, "Group.1"]
M.CGI <- subset(M.CGI, select = - Group.1)
str(M.CGI, max.level = 0)

# save the raw meta data, normalized/aggregated beta and M values to file
save(raw.meta, beta.CGI, M.CGI, file = mBetaNormCgiFileName)

# EOF