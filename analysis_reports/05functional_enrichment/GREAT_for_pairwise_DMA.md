# GREAT ANALYSIS
Rashedul  
Thursday, April 09, 2015  

**Scripting and interpretation:** Rashedul Islam

[GREAT](http://bejerano.stanford.edu/great/public/html/) is a bioinformatics tool that assings funcitons of non-coding regulatory regions. Here we did enrichment analysis for three differential methylation analyses (DMA) at FDR cutoff 1e-4. In this enrichment analysis we did not include between normal_Cancer vs normal_healthy tissues because they have only 8 CGIs at 1e-4. Other three groups are:

1. normal_vs_adenoma
2. normal_vs_cancer
3. adenoma_vs_cancer

**Methods and dataset:** We collected top hits from the topTable at threshold 1e-4 for each three pairs. There genomic cordinates were taken in all cases using this [Rscript](). Then these genomic coordinates were converted into bed format using text editor, as described [here](). Then the bed file was fed on GREAT and the distance of gene association was selected 1.5Kbp upstream and 300bp downstream because we are looking for promoter regions. GREAT gives functional enrichment terms from differnt functional databases. Here our interest is on pathways and MSigDB [2] perturbation database to see the funcitons asssociated with colorectal cancer. 

**Results:** In normal_vs_adenoma, the pathways predicted by GREAT is 'p38 MAPK' that means the genes involved with these CGIs plays role in 'p38 MAPK' pathway. This recent study [3] showed that 'p38 MAPK' pathway has an anti-apoptopic role and leads to cancer. Inhibition of the pathway has shown a reduced number of cell viability in human colon cancer. 

In normal_vs_cancer, we did not get any pathway and in adenoma_vs_cancer, we found a bery broad term for the genes involved in 'Generic Transcription pathway'.

MSigDB [2] perturbation database shows the possible chemical and genetic perturbations of the genes. Their results are almost similar in three pairs of comparison. Alteration in promoter marks and their involvement in cancers like acute promyelocytic leukemia and breast cancer.

Details of statistics and results are given in the table;



**Reference:**  1.

1. Cory Y McLean, Dave Bristor, Michael Hiller, Shoa L Clarke, Bruce T Schaar, Craig B Lowe, Aaron M Wenger, and Gill Bejerano. "GREAT improves functional interpretation of cis-regulatory regions". Nat. Biotechnol. 28(5):495-501, 2010. [PMID 20436461](http://www.ncbi.nlm.nih.gov/pubmed/20436461).
2. Tamayo, et al. [2005, PNAS 102, 15545-15550](http://www.pnas.org/cgi/content/abstract/102/43/15545) and Mootha, Lindgren, et al. [2003, Nat Genet 34, 267-273](http://www.nature.com/ng/journal/v34/n3/abs/ng1180.html).
3. Yang, Shi Yu, et al. "Inhibition of the p38 MAPK pathway sensitises human colon cancer cells to 5-fluorouracil treatment." International journal of oncology 38.6 (2011): 1695-1702. [PMID: 21424124](http://www.ncbi.nlm.nih.gov/pubmed/21424124).
