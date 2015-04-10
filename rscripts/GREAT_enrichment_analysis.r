#####################################################
# Rashedul Islam
#
# Input: TtopTable data of four differential methylation analysis.
#
# Output: Genomic coordinates at FDR 1e-4 in TSV file. 
#####################################################

adenoma_vs_cancer  <- 
  
p1  <- read.table("../data//limma/adenoma_vs_cancer_santina.tsv")
p1_subset  <- subset(p1, p1$q.value < 1e-4)
#1011 CGIs
p2  <- read.table("../data//limma/normalC_vs_normalH_santina.tsv")
p2_subset  <- subset(p2, p2$q.value < 1e-4)
#8 CGIs
p3  <- read.table("../data//limma/normal_vs_adenoma_santina.tsv")
p3_subset  <- subset(p3, p3$q.value < 1e-4)
#13276 CGIs
p4  <- read.table("../data//limma/normal_vs_cancer_santina.tsv")
p4_subset  <- subset(p4, p4$q.value < 1e-4)
#8247 CGIs

#look at the number of rows or coordinates
str(p1_subset)
str(p2_subset)
str(p3_subset)
str(p4_subset)

#take row names for adenoma_vs_cancer
takeCoord_1  <- row.names(p1_subset)
str(takeCoord_1)

#save row names as tsv file
takeCoordTabSep_1  <- write.table(takeCoord_1, file= "../data/adenoma_vs_cancer_Coord.tsv", quote=FALSE, sep='\t', col.names = NA)

#2take row names for normalC_vs_normalH
takeCoord_2  <- row.names(p2_subset)
str(takeCoord_2)

#save row names as tsv file
takeCoordTabSep_2  <- write.table(takeCoord_2, file= "../data/normalC_vs_normalH_Coord.tsv", quote=FALSE, sep='\t', col.names = NA)

#3take row names for normal_vs_adenoma
takeCoord_3  <- row.names(p3_subset)
str(takeCoord_3)

#save row names as tsv file
takeCoordTabSep_3  <- write.table(takeCoord_3, file= "../data/normal_vs_adenoma_Coord.tsv", quote=FALSE, sep='\t', col.names = NA)

#4take row names for normal_vs_cancer
takeCoord_4  <- row.names(p4_subset)
str(takeCoord_4)

#save row names as tsv file
takeCoordTabSep_4  <- write.table(takeCoord_4, file= "../data/normal_vs_cancer.tsv", quote=FALSE, sep='\t', col.names = NA)


###################################################################
#End of script
###################################################################