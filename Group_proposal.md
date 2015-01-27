<<<<<<< Updated upstream
UBC STAT 540 Project Proposal - January 2015
================================================
Differential Methylation Analysis of Colorectal Cancer
=======
---
output:
  html_document:
    keep_md: yes
---
UBC STAT 540 Project Proposal - January 2015
================================================
Identifying differentially methylated regions underlying colorectal cancer progression
>>>>>>> Stashed changes
----------------

Project group members:
* Rashedul Islam
* Santina Lin
* Ka Ming Nip
* Eva Yap
* Beryl Zhuang

## Introduction: 
**Colorectal cancer (CRC) initiation and step-wise progression are associated with the accumulation of genetic and epigenetic alterations. Epigenetic mechanisms are processes that can change gene expression without modifying the primary sequence of DNA. The best-characterized epigenetic process is DNA methylation, which regulates gene transcription by covalent addition of a methyl group at the 5-carbon position of cytosine within cytosine-guanine dinucleotides (CpG). Aberrant DNA methylation can lead to malignancy by hypermethylation of CpG islands resulting in transcriptional silencing of tumour suppressor genes. CpG islands are sequences with high CpG ratios located within gene promoters and methylation at these sites promotes association of methyl-binding proteins (MBPs) and subsequent recruitment of transcriptional repressors.**

**Changes in DNA methylation is one of the early molecular events involved in CRC initiation and many studies have identified these epigenetic abnormalities in precursor lesions such as aberrant crypt foci and adenomas. For our project, we selected *blah blah blah***

* What's CRC (colorectal cancer) 
* DNA methylation
* how they're related (cancer and methylation, and CRC in specific)

## Part 1: Why's 
* Why we pick this dataset 
- 2014 - recent 
* other justification 

## Part 2: About the data 
__(maybe Rashedul can edit this part?)__

* Description of the dataset 
    - the three groups (normal ,adenoma, CRC) 
    - the types of data (probes, samples,  how it's generated generated from the array) 
    - details about data files 

## Part 3: Research Question(s)
* Can we identify differential DNA methylation patterns that underlie the progression of CRC?

## Part 4: Proposed Work 
As revealed by the corresponding GEO Platform GPL13534, 485577 HumanMethylation450 probes were used in this study. Given our limited computational resources and time, we opt for downsizing the number of probes into a manageable size (~10000). We will survey the literature to compile a list of genes known to be correlated with CRC. Then, we will select probes that hybridize to these genes with a flanking window of 1000 bp **[we definitely need citation for this window size]**. The remaining probes are selected arbitrarily to meet our targeted number of probes.

* other things we might need to do
   - reformatting and cleaning data
   - build a pipeline

## References
* blah
