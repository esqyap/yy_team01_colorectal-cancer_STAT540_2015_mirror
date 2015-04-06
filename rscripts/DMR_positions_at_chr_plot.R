#####################################################
# Rashedul Islam
#
#
# load the topTable data "../data/topTable.*.RData" 
# 
# save the plots
# 
#####################################################
library(ggplot2)

#laod pre-processed information
load("../data//cginame.RData")
load("../data//chr.RData")
load("../data//chrlen.RData")
load("../data//coord.RData")


#subset data
#DMR1.1  <- subset(adenoma_cancer_dma, adenoma_cancer_dma$P.Value < 1e-3)
#str(DMR1.1)

#####test data end

p1  <- read.table("../data//limma/adenoma_vs_cancer_santina.tsv")
p2  <- read.table("../data//limma/normalC_vs_normalH_santina.tsv")
p3  <- read.table("../data//limma/normal_vs_adenoma_santina.tsv")
p4  <- read.table("../data//limma/normal_vs_cancer_santina.tsv")

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

#str(p4)
#plot1  <- drawplot("../data//limma/")
#plot
# now plot coordinates of MDR for probes
Plot_adenoma_vs_cancer  <- drawplot(p1)
Plot_adenoma_vs_cancer
Plot_normalC_vs_normalH <- drawplot(p2)
Plot_normalC_vs_normalH
Plot_normal_vs_adenoma  <- drawplot(p3)
Plot_normal_vs_adenoma
Plot_normal_vs_cancer <- drawplot(p4)
Plot_normal_vs_cancer

#ggsave(plot = coord.plot, filename = "../figures/coord.plot.png", width = 16, height = 10.67, units = "in")
ggsave(plot = Plot_adenoma_vs_cancer, filename = "../figures/Plot_adenoma_vs_cancer.png", width = 16, height = 10.67, units = "in")
ggsave(plot = Plot_normalC_vs_normalH, filename = "../figures/Plot_normalC_vs_normalH.png", width = 16, height = 10.67, units = "in")
ggsave(plot = Plot_normal_vs_adenoma, filename = "../figures/Plot_normal_vs_adenoma.png", width = 16, height = 10.67, units = "in")
ggsave(plot = Plot_normal_vs_cancer, filename = "../figures/Plot_normal_vs_cancer.png", width = 16, height = 10.67, units = "in")

#####################################################
# End of script
#####################################################
