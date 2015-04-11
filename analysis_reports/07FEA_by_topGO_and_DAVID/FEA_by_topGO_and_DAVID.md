# Functional enrichment analysis by topGO and DAVID
Beryl  
2015-04-10  
Goal: identify relevant biological processes or molecular functions that are enriched With the 34 candidate CGIs. 




```r
if(file.exists("../../data/FEA/1e-04/enrichment_table.tsv")){
	FEA_tb <- read.table("../../data/FEA/1e-04/enrichment_table.tsv", header = T)
} else {
	source("../../rscripts/FEA_topGO_analysis.R")
	FEA_tb <- read.table("../../data/FEA/1e-04/enrichment_table.tsv", header = T)
	}
```

The top 20 enriched terms:


```r
library(knitr)
kable(FEA_tb, format = "markdown", align ="l")
```



|GO.ID      |Term                                        |Annotated |Significant |Expected |Rank.in.classicFisher |classicFisher |classicKS |elimKS |
|:----------|:-------------------------------------------|:---------|:-----------|:--------|:---------------------|:-------------|:---------|:------|
|GO:0003674 |molecular_function                          |15366     |11          |9.33     |106                   |0.330         |< 1e-30   |<1e-30 |
|GO:0005488 |binding                                     |13091     |10          |7.95     |100                   |0.272         |< 1e-30   |<1e-30 |
|GO:0005515 |protein binding                             |9229      |8           |5.60     |94                    |0.199         |< 1e-30   |<1e-30 |
|GO:0046872 |metal ion binding                           |4485      |3           |2.72     |113                   |0.514         |2.3e-19   |<1e-30 |
|GO:0003824 |catalytic activity                          |5310      |4           |3.22     |109                   |0.404         |3.0e-21   |<1e-30 |
|GO:0097159 |organic cyclic compound binding             |6239      |7           |3.79     |69                    |0.086         |2.0e-28   |<1e-30 |
|GO:0042605 |peptide antigen binding                     |18        |0           |0.01     |122                   |1.000         |7.7e-05   |<1e-30 |
|GO:1901363 |heterocyclic compound binding               |6182      |7           |3.75     |66                    |0.083         |7.7e-28   |<1e-30 |
|GO:0005001 |transmembrane receptor protein tyrosine ... |98        |0           |0.06     |123                   |1.000         |1.1e-11   |<1e-30 |
|GO:0003677 |DNA binding                                 |2833      |1           |1.72     |121                   |0.823         |2.4e-14   |<1e-30 |
|GO:0046332 |SMAD binding                                |128       |0           |0.08     |124                   |1.000         |1.8e-10   |<1e-30 |
|GO:0003676 |nucleic acid binding                        |4144      |2           |2.52     |120                   |0.720         |4.9e-19   |<1e-30 |
|GO:0043565 |sequence-specific DNA binding               |1162      |0           |0.71     |125                   |1.000         |8.9e-08   |<1e-30 |
|GO:0042288 |MHC class I protein binding                 |10        |0           |0.01     |126                   |1.000         |0.0094    |<1e-30 |
|GO:0046978 |TAP1 binding                                |4         |0           |0.00     |127                   |1.000         |0.0239    |<1e-30 |
|GO:0005005 |transmembrane-ephrin receptor activity      |15        |0           |0.01     |128                   |1.000         |4.6e-07   |<1e-30 |
|GO:0032395 |MHC class II receptor activity              |4         |0           |0.00     |129                   |1.000         |0.0249    |<1e-30 |
|GO:0001106 |RNA polymerase II transcription corepres... |45        |0           |0.03     |130                   |1.000         |1.1e-07   |<1e-30 |
|GO:0032041 |NAD-dependent histone deacetylase activi... |24        |0           |0.01     |131                   |1.000         |2.5e-07   |<1e-30 |
|GO:0046970 |NAD-dependent histone deacetylase activi... |24        |0           |0.01     |132                   |1.000         |2.5e-07   |<1e-30 |

There are a few very broad terms, such as molecular function, binding, protein binding, that are not informative.

Among the other top terms, such as transmembrane receptor protein tyrosine phosphatase activity, SMAD binding, and MHC class I protein binding are related to tumour suppressor genes, which have been found to be involved in CRC progression and other cancer types too[1, 2, 3]. The enrichment results show that the candidate genes are relevant to CRC. However, there are no new discoveries among the top enriched terms.

####The top 5 nodes from Fisher test
![img](https://raw.githubusercontent.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/master/data/FEA/1e-04/Fisher_top5nodes.png?token=AIR6N8h1hFe2UbsotEXVJDbGtMlvAx9mks5VMezTwA%3D%3D)


#### Top 15 nodes from K-S test
![img](https://raw.githubusercontent.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/master/data/FEA/1e-04/KS_top15nodes.png?token=AIR6NwpTrEQAQUAoLVWqdhPYX7XZG_bnks5VMe0swA%3D%3D)


### Enrichment analysis by DAVID
See detail and results [here](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/tree/master/data/DAVID)

####DAVID results: Functional annotation clustering
![img](https://raw.githubusercontent.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/master/data/DAVID/DAVID_anno_clustering_normal_adenoma_cancer.png?token=AIR6N1P2a5HlvA1RKG56rTar2i6eToeWks5VMfE-wA%3D%3D)


####DAVID results: Functional annotation chart
![img](https://raw.githubusercontent.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/master/data/DAVID/DAVID_anno_chart_normal_adenoma_cancer.png?token=AIR6N9vYuN6Bb-zUngobNfveriZMbKhuks5VMfGswA%3D%3D)


References:

1. [Ostman, A. et al. "Protein-tyrosine phosphatases and cancer." Nat Rev Cancer. 2006 Apr; 6(4):307-20.](http://www.ncbi.nlm.nih.gov/pubmed/16557282)
2. [Derynck, R. et al. “ TGF-ß signaling in tumor suppression and cancer progression.” Nat Genet. 2001 Oct; 29(2): 117-29.](http://www.ncbi.nlm.nih.gov/pubmed/11586292)
3. [Garcia-Lora, A. et al. “MHC Class I Antigens, Immune Surveillance, and Tumour Immune Escape.” J Cell Physiol. 2003 Jun; 195(3): 346-55.](http://www.ncbi.nlm.nih.gov/pubmed/12704644)
