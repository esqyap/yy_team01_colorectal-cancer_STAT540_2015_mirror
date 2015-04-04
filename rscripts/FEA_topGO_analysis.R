#beryl still working on this
library(IlluminaHumanMethylation450k.db)
library(topGO)


load("../data/FEA_island2go.Rdata")

# get all the cgi island names
island<-as.data.frame(IlluminaHumanMethylation450kCPGINAME)
island_names<-island$cpgiview.ucscname
head(island_names)



####################################
# functions
####################################

#' functions for a list of predefined genes
#' e.g. genes found by multiple DMA
 
predefinedGenes <- function(interesting_genes){
# take a list of interesting genes e.g. "chrY:9363680-9363943"
# to generate a geneList object to build a topGOdata
  gene_list <- factor(as.integer(island_names %in% interesting_genes))
  names(gene_list) <- island_names
  return(gene_list)
}

makeTopGODataPreDefined <- function(predefined_gene_list){
  # take a named list of values (names are the 
  # chr coordinates) to generate a 
  # topGOdata object
  new("topGOdata", 
      description = "CRC Differential methylation analysis",
      ontology = "MF", 
      allGenes = predefined_gene_list, 
      annot = annFUN.gene2GO, 
      gene2GO = Island2GO)
}


#e.g.
x <- predefinedGenes("chrY:9363680-9363943")
makeTopGODataPreDefined(x)

###############################################
#' for results straight from the toptable (no filter)
#' filter is applied when building the topGOdata object

thresholdFun <- function (value) {
  return(value < 0.01)
}

makeTopGOData <- function(gene_list){
  # take a named numeric vector (names are the 
  # chr coordinates) to generate a 
  # topGOdata object
  new("topGOdata", 
      description = "CRC Differential methylation analysis",
      ontology = "MF", 
      allGenes = gene_list, 
      geneSel = thresholdFun,
#      nodeSize = 5,
      annot = annFUN.gene2GO, 
      gene2GO = Island2GO)
}

#'  Once the genes are annotated to the each GO term and the true
#'  path rule is applied the nodes with less than nodeSize annotated genes are removed from the GO hierarchy.
#'  We found that values between 5 and 10 for the nodeSize parameter yield more stable results. The default
#'  value for the nodeSize parameter is 1, meaning that no pruning is performed.

####################################
### example
####################################
load("../data/GSE48684_raw_filtered.m.norm.cgi.Rdata")

# make a random list
x <- row.names(M.norm.CGI)[1:1000]
lm_list <- runif(1000, 0, 1)

names(lm_list) <- x

test_GO <- makeTopGOData(lm_list)


### TESTS
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
###Summary of top GO groups
allRes_APL <- GenTable(GOdata, classicFisher = resultFisher,
                       classicKS = resultKS, elimKS = resultKS.elim,orderBy = "elimKS", 
                       ranksOf = "classicFisher", topNodes = 10)
save(allRes_ALL, file="allRes_ALL.R") 

