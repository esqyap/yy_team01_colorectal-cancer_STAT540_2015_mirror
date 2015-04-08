#########################################################
# Beryl
#' the test codes are taken from 
#' https://github.com/sjackman/stat540-project/blob/master/topGO.R
#' with modification
###########################################################
#' functional enrichment analysis
#' with output tables and intermediate Rdata
#' compare differentially methylated CGI among normal, adenoma, and cancer.
#' output genes associtated with the CGIs
#' output enriched GO terms
#' 
library(IlluminaHumanMethylation450k.db)
library(topGO)

#####################inputs
# load island to GO
island_GO_file <- "../data/FEA_island2go.Rdata"

if(file.exists(island_GO_file)){
	load(island_GO_file)
	} else {
	source("FEA_build_island2GO.R")
	load(island_GO_file)
}

# load the toptables
normal_adenoma <- read.delim("../data/topTables/normal_vs_adenoma_santina.tsv")
normal_cancer <- read.delim("../data/topTables/normal_vs_cancer_santina.tsv")
adenoma_cancer <- read.delim("../data/topTables/adenoma_vs_cancer_santina.tsv")

# set FDR cutoff value
cutoff <- 1e-4

## save the files
file_dir <- paste0("../data/FEA/", as.character(cutoff), "/")
dir.create(path = file_dir, showWarnings = F)


#######################################



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


###############################################
###############################################


getChr <- function(tb){
	candidate_list <- as.character(rownames(tb)[which(tb$q.value < cutoff)])
	return(candidate_list)
}

NA_NC <- intersect(getChr(normal_adenoma), getChr(normal_cancer))

NA_AC <- intersect(getChr(adenoma_cancer), getChr(normal_adenoma))

NC_AC <- intersect(getChr(adenoma_cancer), getChr(normal_cancer))

all_groups <- intersect(NA_AC, NA_NC)




# make GOdata
x <- predefinedGenes(all_groups)
GOdata <- makeTopGODataPreDefined(x)

### TESTS
fisher_file <- paste0(file_dir, "/result_fisher.Rdata")
if(file.exists(fisher_file)){
	load(file = fisher_file)
} else {
	result_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
	save(result_fisher, file = paste0(file_dir, "/result_fisher.Rdata"))
}

KS_file <- paste0(file_dir, "/result_KS.Rdata")
if(file.exists(KS_file)){
	load(file = KS_file)
	} else {
	result_KS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
	save(result_KS, file = KS_file)
}


KS.elim_file <- paste0(file_dir, "/result_KS_elim.Rdata")
if(file.exists(KS.elim_file)){
	load(file = KS.elim_file)
} else {
	result_KS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
	save(result_KS.elim, file = KS.elim_file)
}
###Summary of top GO groups
all_tests_result <- GenTable(GOdata, classicFisher = result_fisher,
											 classicKS = result_KS, elimKS = result_KS.elim, orderBy = "elimKS", 
											 ranksOf = "classicFisher", topNodes = 20)


##### save all the files
save(all_groups, file = paste0(file_dir, "/candidate_chr.Rdata"))
save(all_tests_result, file = paste0(file_dir, "/all_tests_result.Rdata"))

head(all_tests_result)

write.table(all_tests_result, 
						file = paste0(file_dir, "/enrichment_table.tsv"), 
						row.names = F, col.names = T)
write.table(as.data.frame(all_groups), 
						file = paste0(file_dir, "/candidates.tsv"), 
						row.names = F, col.names = F)


###############################################
# from chr coordinates to gene symbol
load(paste0(file_dir, "/candidate_chr.Rdata"))
sym <- as.data.frame(IlluminaHumanMethylation450kSYMBOL)

# get probe.id
probe_id <- island$cpgiview.Probe_ID[which(island$cpgiview.ucscname %in% all_groups)]

# probe id to gene symbol
gene_symbol <- sort(unique(sym$symbol[which(sym$probe_id %in% probe_id)]))

write.table(as.data.frame(gene_symbol), 
						file = paste0(file_dir, "/genes.txt"), 
						row.names = F, col.names = F)


################## making plots
#' Rectangles indicate the top most significant terms. Rectangle color represents the relative significance,
#' ranging from dark red (most significant) to bright yellow (least significant). For each node, some basic information
#' is displayed. The first two lines show the GO identifier and a trimmed GO name. In the third line the raw p-value
#' is shown. The forth line is showing the number of significant genes and the total number of genes annotated to
#' the respective GO term.

printGraph(GOdata, result_KS, firstSigNodes = 10, 
					 result_fisher, fn.prefix = paste0(file_dir, "KS_fisher"), useInfo = "all", 
					 pdfSW = TRUE)

# save both pdf and png
nodes <- 5
pdf(file = paste0(file_dir, "/Fisher_top", nodes, "nodes.pdf"))
showSigOfNodes(GOdata, score(result_fisher), 
							 firstSigNodes = nodes, 
							 useInfo = "all")
dev.off()

png(filename = paste0(file_dir, "/Fisher_top", nodes, "nodes.png"),
		width = 1400, height = 1400, res = 300)
showSigOfNodes(GOdata, score(result_fisher), 
							 firstSigNodes = nodes, 
							 useInfo = "all")
dev.off()


nodes <- 15
pdf(file = paste0(file_dir, "/KS_top", nodes, "nodes.pdf"))
showSigOfNodes(GOdata, score(result_KS), 
							 firstSigNodes = nodes, 
							 useInfo = "all")
dev.off()

png(filename = paste0(file_dir, "/KS_top", nodes, "nodes.png"),
		width = 1400, height = 1400, res = 300)
showSigOfNodes(GOdata, score(result_KS), 
							 firstSigNodes = nodes, 
							 useInfo = "all")
dev.off()

# ###select some genes?
# x <- all_tests_result[c(1,9, 11, 14), ]
# write.table(x, sep="\t", file = paste0(file_dir, "/seleted_top_GO.tsv", 
# 						col.names = TRUE, row.names = F))
# 
# ### 

						

# part 2
###############################################
#' for results straight from the toptable (no filter)
#' filter is applied when building the topGOdata object

cutoff <- 1e-4
thresholdFun <- function (value) {
  return(value < cutoff)
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

#'  from topGO manual:
#'  Once the genes are annotated to the each GO term and the true
#'  path rule is applied the nodes with less than nodeSize annotated genes are removed from the GO hierarchy.
#'  We found that values between 5 and 10 for the nodeSize parameter yield more stable results. The default
#'  value for the nodeSize parameter is 1, meaning that no pruning is performed.

####################################
### example
####################################
# load("../data/GSE48684_raw_filtered.m.norm.cgi.Rdata")
# 
# # make a random list
# x <- row.names(M.norm.CGI)[1:1000]
# lm_list <- runif(1000, 0, 1)
# 
# names(lm_list) <- x
# 
# test_GO <- makeTopGOData(lm_list)
# 
# GOdata <- makeTopGOData(lm_list)
# ### TESTS
# result_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
# result_KS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
# result_KS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
# ###Summary of top GO groups
# all_groups <- GenTable(GOdata, classicFisher = result_fisher,
# 											 classicKS = result_KS, elimKS = result_KS.elim,orderBy = "elimKS", 
# 											 ranksOf = "classicFisher", topNodes = 10)
# 
# head(all_groups)


###################################
