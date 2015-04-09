
# K mean clustering analysis on beta values of CGI 

This folder contains some __exploratory analysis__ and quality check on the methylation data: 
- beta normalized values for CGI (CpG islands)


The rationale is to see how well the samples can cluster in their own group given just the methylation values, and compare the results to that in the previous [k mean clustering analysis](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/tree/master/analysis_reports/02kmeans)

# Analysis and results 
- With k = 4 in k-mean clustering analysis, adenoma and cancer tend to cluster in groups, sometimes with normal samples. Normal samples (normal-H and normal-C) tend to cluster together but their groups can contain cancer or adenoma samples. 
- sum of square measurement was used to see the noise in each cluster as k, number of clusters, increase. This is a common way to determine how many cluster, or the value of K, to used in a k-mean clustering analysis
	- the sharpest drop was observed at k = 3, which could be the most optimal k value. This would make sense if normal-H and normal-C are very similar in methylation levels. 
	- clustering analysis with k = 3 is performed on the samples. While two groups have exclusively adenoma and cancer samples (equal mix), the third group has samples from all four groups. 
