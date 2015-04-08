
# Content 

Here's the folder for hierarchical clustering. Exploratory data analysis to look at whether samples within each group have methylation levels that are similar enough that hierarchical clustering using ward method can cluster them into the same group. 

- hierarchicalClustering: 
	- clustering analysis on filtered raw data 

- hierarchicalClustering_normalized.cgi 
	- clustering analysis on normalized data (beta value) in which CGI probes are averaged. 


The clustering analysis for each data set includes 
- War's minimum variance method  'ward'
- Furtherest neighbor / compact  'complete' 

Each data set also has three different version: 
- original
- with incomplete rows removed (NA removed)
- with NA values set to zero

Therefore, there are six analyses for each of the two dataset. 


# Result 

Several data processing (removing NA or setting NA to zeros) and two different unsupervised hierarchical clustering analysis (ward and complete) were used. Ward method on original data (without removing NA or setting NA to zeros) were used for better clustering results and visualization on the poster. 

See the results in the folder called "figures" at the root directory. Figures that are used in the poster were generated from the markdown, though not shown in the markdown because images were directly saved to the "figures" folder instead of being displayed in the markdown. 

On filtered raw data (i.e. X chromosome probes removed) : adenoma and cancer are mostly clustered into two different groups. Some cancer samples are mixed with normal-H and normal-C in two other groups. Some adenoma and cancer are in two groups that are more homogenous, and the rest of the samples are more mixed together. 

On data with CGI Beta values: the analysis shows that normal samples are clustered into a group by itself whereas adenoma and cancer are clustered into three other groups. In short, good separation of normal samples and non-normal samples. 

Interpretation: 
This exploratory analysis helps to show that: 
1) methylation levels do differ among groups. 
2) normalization and using averaged beta value (for only CpG islands) do improve the signals that will help differentiate groups 
