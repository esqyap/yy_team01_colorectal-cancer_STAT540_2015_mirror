#####################################################
# Ka Ming Nip
# 

library(dplyr)
library(limma)

# loading the norm. CGI.
M.norm.CGI.path <- "../data/GSE48684_raw_filtered.m.norm.cgi.Rdata"
load(M.norm.CGI.path)

metadata.path <- "../data/metadata.Rdata"
load(metadata.path)

# sanity
head(M.norm.CGI)
head(metadata)

# function to perform limma to generate top table
limmaTopTable <- function(dat, des){
  myFit <- lmFit(dat, des)
  myEbFit <- eBayes(myFit)
  myTopTable <- topTable(myEbFit, number=nrow(dat), colnames(coef(myEbFit)))
  return(myTopTable)
}

# design matrix for limma on group
des <- metadata %>% 
  select(group, colon_region, gender, stage)
rownames(des) <- metadata$geo_accession


desMat.group <- model.matrix(~group, des)

# remove NA from data
nums <- sapply(M.norm.CGI, is.numeric)
dat <- na.omit(M.norm.CGI[ , nums])

topTable.group <- limmaTopTable(dat, desMat.group)


# design matrix for group and gender
des.group.gender <- des %>%
  select(group, gender)
desMat.group.gender <- model.matrix(~group+gender, des.group.gender)

dat.group.gender <- dat[, rownames(des.group.gender)]

topTable.group.gender <- limmaTopTable(dat.group.gender, desMat.group.gender)


# design matrix for cancer and stage
des.cancer.stage <- droplevels(subset(des, group == "cancer"))
desMat.cancer.stage <- model.matrix(~stage, des.cancer.stage)

dat.cancer.stage <- dat[, rownames(des.cancer.stage)]

#myFit <- lmFit(dat.cancer.stage, desMat.cancer.stage)
#myEbFit <- eBayes(myFit)
#myTopTable <- topTable(myEbFit, number=nrow(dat), colnames(coef(myEbFit)))

topTable.cancer.stage <- limmaTopTable(dat.cancer.stage, desMat.cancer.stage)

# WORKS UP TO HERE

# design matrix for limma on group and colon_region
des.no.unknown.region <- droplevels(subset(des, colon_region != 'unknown'))
desMat.group.plus.region <- model.matrix(~group+colon_region, des.no.unknown.region)
desMat.group.star.region <- model.matrix(~group*colon_region, des.no.unknown.region)

dat.no.unknown.region <- dat[, rownames(des.no.unknown.region)]

topTable.group.plus.region <- limmaTopTable(dat.no.unknown.region, desMat.group.plus.region)
topTable.group.star.region <- limmaTopTable(dat.no.unknown.region, desMat.group.star.region)

#####################################################
# End of script
#####################################################