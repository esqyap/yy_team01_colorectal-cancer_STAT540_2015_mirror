#########################################################
# Beryl
#' the test codes are taken from 
#' https://github.com/sjackman/stat540-project/blob/master/topGO.R
#' with modification
###########################################################
#' functional enrichment analysis
#' with output tables and intermediate Rdata
#' compare differentially methylated CGI among normal H and normal C groups
#' output genes associtated with the CGIs
#' output enriched GO terms


library(IlluminaHumanMethylation450k.db)
library(topGO)

######### inputs
# load island to GO
island_GO_file <- "../data/FEA_island2go.Rdata"

if(file.exists(island_GO_file)){
	load(island_GO_file)
} else {
	source("FEA_build_island2GO.R")
	load(island_GO_file)
}

# load the toptables
normal_HC <- read.delim("../data/topTables/normalC_vs_normalH_santina.tsv")
# set FDR cutoff value
cutoff <- 1e-4

## save the files
file_dir <- paste0("../data/FEA/normal_HC_", as.character(cutoff), "/")
dir.create(path = file_dir, showWarnings = F)

############

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

# get the candidate genes
all_groups <- getChr(normal_HC)


# make GOdata
x <- predefinedGenes(all_groups)
GOdata <- makeTopGODataPreDefined(x)
# 
### TESTS
result_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
result_KS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
result_KS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
###Summary of top GO groups
all_tests_result <- GenTable(GOdata, classicFisher = result_fisher,
														 classicKS = result_KS, elimKS = result_KS.elim, orderBy = "elimKS", 
														 ranksOf = "classicFisher", topNodes = 20)


##### save all the files
save(all_groups, file = paste0(file_dir, "/candidate_chr.Rdata"))
save(result_fisher, file = paste0(file_dir, "/result_fisher.Rdata"))
save(result_KS, file = paste0(file_dir, "/result_KS.Rdata"))
save(result_KS.elim, file = paste0(file_dir, "/result_KS_elim.Rdata"))
save(all_tests_result, file = paste0(file_dir, "/all_tests_result.Rdata"))

head(all_tests_result)

write.table(as.data.frame(all_groups), 
						file = paste0(file_dir, "/candidates.tsv"), 
						row.names = F, col.names = F)

write.table(all_tests_result, 
						file = paste0(file_dir, "/enrichment_table.tsv"), 
						row.names = F, col.names = T)

#####################################################################
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
