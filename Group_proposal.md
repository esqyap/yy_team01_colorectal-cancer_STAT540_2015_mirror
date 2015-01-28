
UBC STAT 540 Project Proposal - January 2015
================================================
Identifying differentially methylated regions underlying colorectal cancer progression
----------------

Project group members:
* Rashedul Islam
* Santina Lin
* Ka Ming Nip
* Eva Yap
* Beryl Zhuang

## Introduction: 
Colorectal cancer (CRC) initiation and step-wise progression are associated with the accumulation of genetic and epigenetic alterations. Epigenetic mechanisms are processes that can change gene expression without modifying the primary sequence of DNA. The best-characterized epigenetic process is DNA methylation, which regulates gene transcription by covalent addition of a methyl group at the 5-carbon position of cytosine within cytosine-guanine dinucleotides (CpG). Aberrant DNA methylation can lead to malignancy by hypermethylation of CpG islands, resulting in transcriptional silencing of tumour suppressor genes. CpG islands are sequences with high CpG ratios located within gene promoters. Methylation at these sites promotes association of methyl-binding proteins (MBPs) and subsequent recruitment of transcriptional repressors.

Changes in DNA methylation is one of the early molecular events involved in CRC initiation and many studies __(cite here?)__ have identified these epigenetic abnormalities in precursor lesions such as aberrant crypt foci and adenomas. __(may be emphasize how early identification of methylation sites can help diagnosis? blah blah)__

* What's CRC (colorectal cancer) 
* DNA methylation
* how they're related (cancer and methylation, and CRC in specific)

## Part 1: about the data: which data we picked and why, and description of the data
For our project, we selected *blah blah blah***

* Why we pick this dataset 
- 2014 - recent 
* other justification 

__(maybe Rashedul can edit this part?)__

* Description of the dataset 
    - the three groups (normal ,adenoma, CRC) :41 normal colon tissue, 42 colon adenomas,
and 64 cancers
    - the types of data (probes, samples,  how it's generated generated from the array) 
    - details about data files 

## Part 2: Research Question(s)

Can we identify differential DNA methylation patterns that underlie the progression of CRC? If so, what are the genes that are differentially methylated in the three groups?

## Part 3: Proposed Work 
As revealed by the corresponding GEO Platform GPL13534, 485577 HumanMethylation450 probes were used in this study. Given our limited computational resources and time, we opt for downsizing the number of probes into a manageable size (~10000). We will survey the literature to compile a list of genes known to be correlated with CRC and select probes that hybridize to these genes. To identify novel genes that could be involved in the tumour progression, we will also include probes that are found on or near CpG islands, regions that are some associated with promoters. Therefore, we could infer what genes are downregulated or upregulated at the level of transcription of we identify any significant changes in methylation. To explore the data more throughly, we could also include some randomly selected probes in our analysis or consider looking through all the probes if we successfully built a automated analysis pipieline. 

At the start of the project, we want to perform some data cleaning and reformatting to ensure that we work with high-quality data. We would first to a exploratory analysis to see whether there are missing values or missing probes in some samples to see how consistent our data are and decide how we want to account for the discrepancy among each sample. We would also do clustering analysis on data from three different patient groups to see if data from the same group would cluster together. Not only these initial check could help us in designing our analysis steps, but also give us a sense of what the raw data look like. 

In order to do group comparison among the data of the three different tumour progression stages, we will perform several statistical tests, such as ANOVA and t-test as covered in this class. Depending on the remaining time and resource, we could also build prediction models based on our data to see if we could do prediction on the tumour progression based on the methylation data in early stages of colorectal cancer. 


## References
* blah (need to fill in)
