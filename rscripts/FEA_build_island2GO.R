#' Beryl Zhuang
#' April 3, 2015
#' Disclaimer:
#' The following script is taken from
#' https://github.com/sjackman/stat540-project/blob/master/topGO.R
#' island2GO contains GO identifiers, which is used to build the topGOdata object
#' for the functional enrichment analysis


library(IlluminaHumanMethylation450k.db)
library(topGO)

#### ALL ISLANDS 
#GO annotations for all Islands on the 450K
GO <- as.list(IlluminaHumanMethylation450kGO2PROBE)
head(GO) #GO<-GO2probe, want probe to GO then Island to GO
probe2GO<-inverseList(GO)

#Summaize GO groups associated with each island (object for topGO function)
Island<-as.data.frame(IlluminaHumanMethylation450kCPGINAME)
#lookup all the GO of each Island and store as list
islGO<-function(x) probe2GO[[Island[x,1]]]
isl<-as.list(1:nrow(Island))
#Island GO data (from probe GO data)
Island2GO<-lapply(isl,islGO)


names(Island2GO)<-Island$cpgiview.ucscname
save(Island2GO, file="../data/island2go.Rdata")

# load(file="island2go.Rdata")