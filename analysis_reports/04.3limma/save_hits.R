# how many DMR are there at FDR < 1e-05?
norm_cancer_dmr <- subset(norm_cancer_dma, adj.P.Val < 1e-04)
nrow(norm_cancer_dmr)

norm_adenoma_dmr <- subset(norm_adenoma_dma, adj.P.Val < 1e-04)
nrow(norm_adenoma_dmr)

adenoma_cancer_dmr <- subset(adenoma_cancer_dma, adj.P.Val < 1e-04)
nrow(adenoma_cancer_dmr)

# save the row.names in data.frame
norm_cancer_dmr_list <- data.frame(row.names(norm_cancer_dmr))
write.table(norm_cancer_dmr_list, "norm_cancer_dmr_list.tsv",
            row.names=F, col.names=F)

norm_adenoma_dmr_list <- data.frame(row.names(norm_adenoma_dmr))
write.table(norm_adenoma_dmr_list, "norm_adenoma_dmr_list.tsv",
            row.names=F, col.names=F)

adenoma_cancer_dmr_list <- data.frame(row.names(adenoma_cancer_dmr))
write.table(adenoma_cancer_dmr_list, "adenoma_cancer_dmr_list.tsv",
            row.names=F, col.names=F)