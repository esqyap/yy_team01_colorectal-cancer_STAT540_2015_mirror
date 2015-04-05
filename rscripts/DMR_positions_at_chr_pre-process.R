#####################################################
# Rashedul Islam
#
# Input: IlluminaHumanMethylation450k.db
#
# Output: chrlen, file = "../data/chrlen.RData"
#           - contain length of each chromosome
#         chr, file = "../data/chr.RData
#           - get the probe identifiers that are mapped to chromosome
#         coord, file = "../data/coord.RData
#					  - get chromosome coordinate of each probe
#         cginame, file = "../data/cginame.RData
#           -Extract probe ID to CpG islands association
#
# Note: codes borrowed from Seminar08
#####################################################

source('http://bioconductor.org/biocLite.R')
biocLite('GEOquery')
biocLite('wateRmelon')
biocLite("IlluminaHumanMethylation450k.db")

library(GEOquery)
library(wateRmelon)
library(IlluminaHumanMethylation450k.db)

if(file.exists("../data/chrlen.Rdata")){
  print("chrlen.RData data saved")
} else {  
  # get the length of chromosome 1-22 and Y. In our study we ignored X chromosome.
  chrlen <-
  unlist(as.list(IlluminaHumanMethylation450kCHRLENGTHS)[c(as.character(1:22),
                                                           "Y")])   
  chrlen <- data.frame(chr = factor(names(chrlen)), length = chrlen)
  # save chr in .RData
  save(chrlen, file = "../data/chrlen.RData")
  print("chrlen.RData saved")
}

if(file.exists("../data/chr.Rdata")){
  print("chr.Rdata saved")
} else {  
  chr <- IlluminaHumanMethylation450kCHR        # get the chromosome of each probe
  # get the probe identifiers that are mapped to chromosome
  chr <- unlist(as.list(chr[mappedkeys(chr)]))
  # save chr in .RData
  save(chr, file = "../data/chr.RData")
  print("chr.Rdata saved")
} 

if(file.exists("../data/coord.RData")){
  print("coord.RData saved")
} else {  
  # get chromosome coordinate of each probe
  coord <- IlluminaHumanMethylation450kCPGCOORDINATE
  # get the probe identifiers that are mapped to coordinate
  coord <- unlist(as.list(coord[mappedkeys(coord)]))      
  coord <- data.frame(chr = chr[intersect(names(chr), names(coord))],
                    coord = coord[intersect(names(chr), names(coord))])
  # save coord in .RData
  save(coord, file = "../data/coord.RData")
  print("coord.RData saved")
}

if(file.exists("../data/cginame.RData")){
  print("cginame.RData saved")
} else {
  # Extracting probe ID to CpG islands association
  cginame <- as.data.frame(IlluminaHumanMethylation450kCPGINAME)
  names(cginame) <- c('Probe_ID', 'cginame')
  rownames(cginame) <- cginame$Probe_ID
  length(levels(factor(cginame$cginame)))   # No. of CGIs
  # save cginame as .RData
  save(cginame, file = "../data/cginame.RData")
  print("cginame.RData saved")
} 

#####################################################
# End of script
#####################################################
