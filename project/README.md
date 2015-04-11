# project
Rashedul  
Friday, April 10, 2015  

#### Individual report - Rashedul (Team-1: Colorectal cancer)

In this study we performed differential expression analysis of microarray data collected from different stages of colorectal cancer patients e.g, adenoma, cancer as well as from normal patient samples. The group members have put much more efforts to the project that I initially anticipated. In our project, we followed the work procedure described on Seminar08 and therefore, our individual role was divided according to each step: 

**Data QC:**  Beryl Zhuang, Ka Ming Nip

**Exploratory analysis:** Santina Lin, Eva Yap

**Differential Methylation analysis:** Eva Yap, Beryl Zhuang, Ka Ming Nip, Santina Lin

My role in this project was the downstream analysis of differential methylation results. Since our intention for the project was to attempt to identify what methylation changes might be present in different stages of colorectal cancer. I looked mostly at gene function and their pathways using GREAT bioinformatics tools. For functional enrichment analysis, I did the pre-processing of data using [Rscripts1](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/blob/master/rscripts/DMR_positions_at_chr_pre-process.R) and [Rscript2](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/blob/master/rscripts/DMR_positions_at_chr_plot.R). Along with it, I looked at the distribution of significant (FDR 1e-4) CGIs at chromosome level by plotting them at each chromosome except ‘X’.

**My specific role:**


We did differential methylation analysis in three pairs of comparison e.g., ;;;;. We found that only 34 CGIs are common in all three groups. In my analysis with GREAT, I found that these 34 CGIs are associated with 54 genes. On the other hand, Bioconductor package found 21 genes are associated with those 34 CGIs. Interesting 17 genes out of 21 overlapped with those 54 genes of GREAT. However, GREAT did not find any hit for functional enrichment terms. The possible reasons could be the number of CGIs (34) is very few for GREAT to do functional analysis. To overcome this, I took more CGIs to fed in GREAT from each pair-wise differential methylation analysis (stored in topTable). Since our focus is on promoter regions, I looked at proximal regions of the genes (1500bp upstream and 300bp downstream of the genes). Then I found enriched terms and pathway information that are closely related to cancer and sometimes specifically with the colorectal cancer. The details of this results and statistics can be found [here](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/blob/master/analysis_reports/06positions_at_DMR/Positions_of_DMRs.md).

I used the instruction of Seminar08 to plot the distribution of CGIs along chromosomes. In our study we discarded X chromosome of its hypermethylation in female. But we kept Y chromosome in our study. It is interesting to look at the chromosome level distribution of CGIs that in our study Y chromosome has very lower amount (I would say negligible) of differential methylation. This result showed that the rational of not excluding Y chromosome in our study and Y chromosome is less likely to bias our results between sexes in colorectal cancer patients. For this analysis Rscript is [here](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/blob/master/rscripts/GREAT_enrichment_analysis.r) and the detailed resutls are [here](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/tree/master/analysis_reports/05functional_enrichment).

My minor contribution is to help in developing project proposal. 

**Conclusion:**  
I think there are many future directions for this project e.g., to look at the caveats of the differential analysis results and to give more focus on biology to find out the key regulators of the colorectal cancer. Also we can integrate and compare our results with other publicly available datasets to make a robust prediction of the results and functional validity. 
