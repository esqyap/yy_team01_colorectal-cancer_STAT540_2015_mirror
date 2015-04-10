## Project title
Identifying aberrant methylation patterns underlying colorectal cancer progression

## Project Summary 
Aberrant DNA methylation can lead to malignancy by hyper-methylation of CpG islands (CGIs) resulting in transcriptional silencing of tumor suppressor genes. CGIs are sequences with high CpG fractions (>50%) located within gene promoters and methylation at these sites promotes association of methyl- binding proteins and subsequent recruitment of transcriptional repressors **(1,2)**. Colorectal cancer (CRC) accounts for the second highest cancer-related mortality among men and third among women in Canada, and it progresses from precursor lesions such as adenomas **(1,3)**. For our project, we compared methylation patterns between normal mucosa, adenoma, and colorectal tumor. By identifying differentially methylated (DM) CGIs between these three groups, we hope to determine aberrant methylation underlying CRC progression.

The workflow of our analysis is summarized by this [flow chart](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/blob/master/figures/workflow.png). Briefly, we performed the necessary data cleaning, processing, quality checking, and normalization on DNA methylation array data from [Series GSE48684](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48684) to generate M-values for each CGIs for all samples. This served as the input for our differential methylation (DM) analysis in which we performed pairwise comparisons between sample groups (i.e. normal vs. adenoma, adenoma vs. cancer, and normal vs. adenoma) using the R package `limma`. We identified 34 common differentially methylated CGIs between the three pairwise comparisons and these were selected for functional enrichment analysis. Literature review showed that the enriched functions are implicated in progression of various cancer types including CRC **(4,5,6)**. Given that these DM CGIs are also identified in adenomas compared to normal mucosa, these aberrant methylation patterns could be markers for early detection, risk assessment, and disease monitoring of CRC.

**References**
1. Li, J. et al. "Epigenetic Biomarkers: Potential Applications in Gastrointestinal Cancers." ISRN gastroenterology. Mar 6; 2014:464015.
2. Luo, Y. et al. "Differences in DNA methylation signatures reveal multiple pathways of progression from adenoma to colorectal cancer." Gastroenterology. 2014 Aug; 147(2):418-29.e8.
3. Canadian Cancer Society’s Advisory Committee on Cancer Statistics. Canadian Cancer Statistics 2014. Toronto, ON: Canadian Cancer Society; 2014.
4. Ostman, A. et al. "Protein-tyrosine phosphatases and cancer." Nat Rev Cancer. 2006 Apr; 6(4):307-20.
5. Derynck, R. et al. “ TGF-β signaling in tumor suppression and cancer progression.” Nat Genet. 2001 Oct; 29(2): 117-29.
6. Garcia-Lora, A. et al. “MHC Class I Antigens, Immune Surveillance, and Tumour Immune Escape.” J Cell Physiol. 2003 Jun; 195(3): 346-55.

## Group Members

Member	| Graduate Program |	Lab Group | Interest/Expertise | Gemstone |
------------- | -------------|------------- |------------- |------------- |
Rashedul Islam @Rashedul	|Bioinformatics Training Program| Dr. Martin Hirst (CHiBi) | *De novo* Sequence Assembly, Genetics and Molecular Biology, RNA-seq and ChIP-seq Data  Analysis | sapphire |
Santina Lin @Santina  |Bioinformatics Training Program| Dr. Steven Jones (BCGSC) | Data Visualization, Data Analysis, Data Mining  |	quartz |
Ka Ming Nip @kmnip	|Bioinformatics| Dr. Inanc Birol (BCGSC) | Data Visualization, Sequence Assembly, Structural Variation Analyses, Pathogenomics| amethyst |
Eva Yap	@evayap|Experimental Medicine|	Dr. Aly Karsan (BCGSC) | Biochemistry | emerald |
Beryl Zhuang @BerylZhuang	|Bioinformatics Training Program| Dr. William Hsiao (BCCDC) | Genetics, Application of Bioinformatics tools |	beryl |
- CHiBi = Centre for High-Throughput Biology
- BCGSC = Canada's Michael Smith Genome Sciences Centre		
- BCCDC = BC Centre for Disease Control

## Table of contents
- [initial proposal](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/blob/master/initial_project_summary.md)
- [final proposal](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/blob/master/Group_proposal.md)
- [work flow]
- [analysis scripts and results]
- [Wiki Page](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/wiki) for references and resources.

