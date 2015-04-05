################################
# Eva's script
# Plot boxplot without reordering samples

################################

library(reshape2)
library(ggplot2)

# load m values
load("../data/GSE48684_raw_filtered.m.norm.cgi.Rdata")
str(M.norm.CGI, max.level=0)

# load metadata
load("../data/metadata.Rdata")
str(metadata, max.level=0)

# rename sample labels for the dataset
colnames(M.norm.CGI) <- paste(metadata$group, 
                              gsub("GSM", "", colnames(M.norm.CGI)), sep="_")

# remove probes with NA
M.norm.CGI.rmna <- M.norm.CGI[complete.cases(M.norm.CGI), ]

# make data frame tall and skinny
m.cgi.tall <- melt(M.norm.CGI.rmna, variable.name="Samples", value.name="MValues")
m.cgi.tall <- data.frame(m.cgi.tall, Group=gsub("_\\d+", "", m.cgi.tall$Samples))

# order the samples by group
m.cgi.tall <- m.cgi.tall[order(m.cgi.tall$Group), ]

# make ggplot stop reordering my genes
m.cgi.tall$Samples <- as.character(m.cgi.tall$Samples)
m.cgi.tall$Samples <- factor(m.cgi.tall$Samples, levels=unique(m.cgi.tall$Samples))

# plot distribution of CGI M values
(p <- ggplot(m.cgi.tall, aes(x=Samples, y=MValues)) +
  geom_boxplot(aes(fill=factor(Group))) + theme_bw() +
  theme(axis.text.x = element_blank()) + 
  xlab("Samples") + ylab("M values") + 
  ggtitle("Distribution of CGI M values"))

# add code to save the boxplot
p <- p +
	theme(text = element_text(size=28))
ggsave(plot = p, filename = "../figures/dataQC_boxplot.png", width = 16, height = 10.67, units = "in")


################################
# End of script
################################