#####################################################
# Ka Ming Nip
# 

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
  select(group, geo_accession, colon_region, gender, stage)
desMat.group <- model.matrix(~group, des)

# remove NA from data
dat <- na.omit(M.norm.CGI)

topTable.group <- limmaTopTable(dat, desMat.group)

# WORKS UP TO HERE

# design matrix for limma on group and colon_region
des.no.unknown.region <- des %>%
  filter(colon_region != 'unknown') %>%
  droplevels()
desMat.group.plus.region <- model.matrix(~group+colon_region, des.no.unknown.region)
desMat.group.star.region <- model.matrix(~group*colon_region, des.no.unknown.region)

dat.no.unknown.region <- dat[, as.character(des.no.unknown.region$geo_accession)]

topTable.group.plus.region <- limmaTopTable(dat.no.unknown.region, desMat.group.plus.region)
topTable.group.star.region <- limmaTopTable(dat.no.unknown.region, desMat.group.star.region)



# design matrix for cancer and stage
des.cancer.stage <- des %>%
  filter(group == "cancer") %>%
  select(geo_accession, stage)
desMat.cancer.stage <- model.matrix(~stage, des.cancer.stage)

dat.cancer.stage <- dat[, as.character(des.cancer.stage$geo_accession)]

topTable.cancer.stage <- limmaTopTable(dat.cancer.stage, desMat.cancer.stage)



# design matrix for group and gender
des.group.gender <- des %>%
  select(group, geo_accession, gender)
desMat.group.gender <- model.matrix(~group+gender, des.group.gender)

dat.group.gender <- dat[, as.character(des.group.gender$geo_accession)]

topTable.group.gender <- limmaTopTable(dat.group.gender, desMat.group.gender)



#####################################################
# End of script
#####################################################