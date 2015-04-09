
# Density plot, CPA, and heatmap 

This folder contains some exploratory analysis and quality check on the methylation data: 
	- raw filtered data
	- normalized beta value 
	- normalized M values
	- normalized beta value at CGI 
	- normalized M values at CGI 

# Analysis and results

- Density plot: shows that normalizing the data allows us to work with data with the same distribution of beta values (CGI) across all groups. 
	 - Higher resolution of the density plot before and after normalization is also in the root folder "figures" and was used for the poster. 
- Heatmap: display correlation of raw beta values and normalized beta values across sample groups and in different colon regions. We do not observe obvious correlation for samples grouped by colon regions. 
- Principal component analysis: PCA shows that most information in our datasets can be captured with the first two components. Samples from different groups tend to cluster together, though the adenoma and cancer samples are more spread out. 
