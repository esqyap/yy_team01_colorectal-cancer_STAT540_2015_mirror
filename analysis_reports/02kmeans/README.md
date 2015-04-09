
# K mean clustering analysis 

This folder contains some __exploratory analysis__ and quality check on the methylation data: 
- raw filtered data
- normalized beta value 
- normalized M values

The rationale is to see how well the samples can cluster in their own group given just the methylation values. 

# Analysis and results 
	- With k = 4 in k-mean clustering analysis, adenoma and cancer tend to cluster in groups and normal-H and normal-C cluster in two other groups. 
	- sum of square measurement was used to see the noise in each cluster as k, number of clusters, increase. This is a common way to determine how many cluster, or the value of K, to used in a k-mean clustering analysis
		- the sharpest drop was observed at k = 3, which could be the most optimal k value. This would make sense if normal-H and normal-C are very similar in methylation levels. 
