
# Folder content 

This folder contains analysis that aggregate all top tables and find important probes of FDR that pass a given threshold and find overlap of those probes for the next step in functional enrichment analysis. 

This folder contains 
- significant_genes.md
	- [top tables](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/tree/master/data/topTables) generated from previous analyses are loaded here and probes with FDR less than a given threshold are stored inside the R data 
	- Venn diagram showing overlapping probes is generated using this report. 
	- Feel free to modify parameters in the report, such as threshold, to explore the data further. 
- R data 
	- cool_genes_e5.RData
		- probes from all four top tables that have passed certain threshold. To see how they are generated, refer to the markdown
		- it's a list of length four containing probes from each of the four comparisons 
	- intersect_genes
		- probes that overlap between each two-group comparison. 
		- a list of length 4 

- To see how the R data are generated and how to read from them, please refer to the markdown 

