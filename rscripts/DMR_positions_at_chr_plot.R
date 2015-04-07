#####################################################
# Rashedul Islam
#
# Input: DMRs obtained from four different DMA experiments
# 
# Output: DMR plots
# 
#####################################################
library(ggplot2)

#laod pre-processed information
load("../data//cginame.RData")
load("../data//chr.RData")
load("../data//chrlen.RData")
load("../data//coord.RData")


#subset test data 
#DMR1.1  <- subset(adenoma_cancer_dma, adenoma_cancer_dma$P.Value < 1e-3)
#str(DMR1.1)

p1  <- read.table("../data//limma/adenoma_vs_cancer_santina.tsv")
p1_subset  <- subset(p1, p1$q.value < 1e-4)
p2  <- read.table("../data//limma/normalC_vs_normalH_santina.tsv")
p2_subset  <- subset(p2, p2$q.value < 1e-4)
p3  <- read.table("../data//limma/normal_vs_adenoma_santina.tsv")
p3_subset  <- subset(p3, p3$q.value < 1e-4)
p4  <- read.table("../data//limma/normal_vs_cancer_santina.tsv")
p4_subset  <- subset(p4, p4$q.value < 1e-4)

drawplot <- function(DMR){
  # coordinates of probes in DM CGIs
  coordDMRprobe <-
    droplevels(na.omit(coord[cginame[cginame$cginame %in%
                                       rownames(DMR),]$Probe_ID,]))
#  head(coordDMRprobe)
  (coord.plot <- ggplot(data = coordDMRprobe) + 
     geom_linerange(aes(factor(chr, levels = c("Y", as.character(22:1))),
                        ymin = 0, ymax = length), data = chrlen, alpha = 0.5) + 
     geom_point(aes(x = factor(chr,
                               levels = c("Y", as.character(22:1))), y = coord),
                position = position_jitter(width = 0.03), na.rm = T) + 
     ggtitle("DMR positions on chromosomes") + 
     ylab("Position of DMRs") +
     xlab("chr") +
     coord_flip() + 
     theme_bw())
  return(coord.plot)
}

# now plot coordinates of MDR for probes
Plot_adenoma_vs_cancer  <- drawplot(p1_subset)
Plot_adenoma_vs_cancer
Plot_normalC_vs_normalH <- drawplot(p2_subset)
Plot_normalC_vs_normalH
Plot_normal_vs_adenoma  <- drawplot(p3_subset)
Plot_normal_vs_adenoma
Plot_normal_vs_cancer <- drawplot(p4_subset)
Plot_normal_vs_cancer


ggsave(plot = Plot_adenoma_vs_cancer, filename = "../analysis_reports/06positions_at_DMR/Plot_adenoma_vs_cancer.png", width = 16, height = 10.67, units = "in")
ggsave(plot = Plot_normalC_vs_normalH, filename = "../analysis_reports/06positions_at_DMR/Plot_normalC_vs_normalH.png", width = 16, height = 10.67, units = "in")
ggsave(plot = Plot_normal_vs_adenoma, filename = "../analysis_reports/06positions_at_DMR/Plot_normal_vs_adenoma.png", width = 16, height = 10.67, units = "in")
ggsave(plot = Plot_normal_vs_cancer, filename = "../analysis_reports/06positions_at_DMR/Plot_normal_vs_cancer.png", width = 16, height = 10.67, units = "in")

#####################################################
# End of script
#####################################################
