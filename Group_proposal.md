UBC STAT 540 Project Proposal - January 2015
================================================
Differential Methylation Analysis of Colorectal Cancer
----------------

Project group members:
* Rashedul Islam
* Santina Lin
* Ka Ming Nip
* Eva Yap
* Beryl Zhuang

## Introduction: 
__(Eva can put the intro in here)__

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
