#########################################################
# Beryl
#' the codes are taken from 
#' https://github.com/sjackman/stat540-project/blob/master/topGO.R
#' with modification
###########################################################
#' functional enrichment analysis
#' with output tables and intermediate Rdata
#' compare differentially methylated CGI among normal H and normal C groups
#' output enriched GO terms and genes





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


###############################################
###############################################
# load the toptables
normal_HC <- read.delim("../data/topTables/normalC_vs_normalH_santina.tsv")
cutoff <- 1e-4

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
file_dir <- paste0("../data/FEA/normal_HC_", as.character(cutoff), "/")
dir.create(path = file_dir, showWarnings = F)
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


###################################
# get genes
load(paste0(file_dir, "/candidate_chr.Rdata"))
myInterestingIslands <- all_groups

## Genes in top Islands
x <- IlluminaHumanMethylation450kSYMBOL
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
xx <- as.data.frame(x[mapped_probes])
gen.isl<-merge(island, xx, by.x="cpgiview.Probe_ID", by.y="probe_id")

gen.isl[1]<-NULL
gen.isl<-unique(gen.isl) #21263 Islands associated with 14770 genes

# function to pull out genes associated with top islands

int.genes<-gen.isl[gen.isl$cpgiview.ucscname %in% myInterestingIslands, 2]

lapply(int.genes, write, paste0(file_dir, "genes.txt"), append=TRUE)