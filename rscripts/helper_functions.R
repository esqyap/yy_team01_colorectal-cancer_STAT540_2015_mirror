########################################################
# Beryl Zhuang
#
# helper functions:
# 1. prepareData
# 2. subsetData
# 3. plotHeatmap
# 4. plotScatter
########################################################

load("../data/metadata.Rdata")
########################################################
#####prepareData: data aggregation
########################################################
prepareData <- function(df, x, design, col_names = c("chr", "value")){
	# takes a list of  (length can be variable) chr/ probes
	# and extract the value from M.norm.CGI, and 
	# metadata from design (metadata), return a tall table
	#
	# args
	#   df: e.g. norm.M.CGI that contains methylation value
	#   x is a list of chr coordinates/probes or other from differential analysis
	# return
	#   a tall table with metadata and value of the M CGI value
	# load packages
	#data aggregation
	temp <- df %>%
		filter(rownames(df) %in% x)
	rownames(temp) <- x
	## transpose as a tall table
	temp2 <- as.data.frame(t(temp))
	## add the rownames name as the new col sample
	temp2$group <- row.names(temp2)
	## melt the table by sample and update colnames
	temp3 <- reshape2::melt(temp2, id=c("group"))
	colnames(temp3) <- c("group", col_names)
	## convert row.names in design as sample name
	temp4 <- design %>%
		mutate(group = rownames(design))
	## dplyr::left_join to add the metadata
	return(suppressMessages(dplyr::left_join(temp3, temp4)))
}


########################################################
####### subsetData
########################################################
# subset the data matrix by a candidate list of chr coordinates, probes etc
# return an ordered matrix (for heatmap)
subsetData <- function(data, candidate_list, design = metadata){
	# take a list of chr coordinates/ probes and get a subset
	## data: the data matrix for methylation values (e.g. normalized M CGI)
	## candidate list: the list of candidates of DE genes, chr coordinates etc
	## design: the metadata df
	
	df <- data[which(row.names(data) %in% candidate_list), ]
	#order the columns by group and colon region
	order_by_group <- design$geo_accession[order(design$group, design$colon_region)]
	return(df[, order_by_group])
}



########################################################
##### plotHeatmap
########################################################
# take a ordered matrix and produce a heatmap
plotHeatmap <- function(x, title = "", legend = "group", 
												size = 2, names = F, 
												row_dendrogram = F, col_dendrogram = F){
	## x is an ordered matrix
	## title is the plot title
	## size is the font size for the rows
	## legend: for group legends, colnames in the design, e.g. c("group", "colon_region")
	## size: the x y axis font size
	## names: turn on(T) or off(F) the row and col names
	## row_dendrogram : turn on(T) or off(F) for the row clustering 
	## col_dendrogram : turn on(T) or off(F) for the column clustering 
	# load packages
	# load plyr, dplyr if it's not loaded
	if("package:RColorBrewer" %in% search() == FALSE)
		library(RColorBrewer)
	if("package:pheatmap" %in% search() == FALSE)
		library(pheatmap)
	#get color palette
	# pallette code is modified from seminar03, courtesy of Dean Attali
	colour_scheme <- "BuPu"
	palette <- colorRampPalette(rev(brewer.pal(n = 9, colour_scheme)))
	paletteSize <- 256
	cols <- palette(paletteSize)
	#heatmap
	annotation <- metadata[legend] #get the legend
	pheatmap(x, color = cols,
					 cluster_rows = row_dendrogram, 
					 cluster_cols = col_dendrogram, # turn on/off dendrogram
					 annotation = annotation,
					 fontsize_row = size,
					 fontsize_col = size,
					 show_rownames = names,
					 show_colnames = names,
					 main = title)
}


########################################################
##### plotScatter
########################################################

plotScatter <- function(df, x, y, gp, title){
	p <- ggplot(df, aes_string(x, y, color = gp)) + 
			 	geom_point(position = position_jitter(width = 0.05)) +
			 	stat_summary(aes_string(group = gp), fun.y = mean, geom ="line") +
			 	ggtitle(title)
	return(p)
}
