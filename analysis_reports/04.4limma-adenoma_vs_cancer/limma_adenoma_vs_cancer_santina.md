# Two Group: Adenoma and Cancer
Santina  
Sunday, April 05, 2015  

# Goal 
Beryl and Eva have been trying to do limma on normal samples (normal-H + normal-C) and cancer samples. They are getting different results so I'll try to see if I can reproduce one of their results. 



# Subset Data 

I'll use normalized cgi M value. 

Load data: 

```r
load("../../data/GSE48684_raw_filtered.m.norm.cgi.Rdata") # M.norm.CGI
load('../../data/metadata.Rdata') # metadata
```

Inspect 

```r
head(M.norm.CGI[, 1:5])
```

```
##                          GSM1183439  GSM1183440 GSM1183441 GSM1183442
## chr1:10003165-10003585   -2.6552923 -2.59856415 -2.5445386 -2.5589572
## chr1:1002663-1005318     -0.4691256 -0.09143954 -0.2533638 -0.2153825
## chr1:100315420-100316009 -3.2051837 -3.28066606 -3.1908783 -3.3171275
## chr1:100435297-100436070 -2.6509063 -2.53842308 -2.3740902 -2.6930074
## chr1:100503482-100504404 -2.5912106 -2.61793251 -2.5333289 -2.6985588
## chr1:10057121-10058108   -1.8793069 -1.60397303 -1.9452079 -1.4107345
##                           GSM1183443
## chr1:10003165-10003585   -2.55737653
## chr1:1002663-1005318     -0.09556418
## chr1:100315420-100316009 -3.44138495
## chr1:100435297-100436070 -2.58548225
## chr1:100503482-100504404 -2.77076081
## chr1:10057121-10058108   -1.89517507
```

```r
head(metadata)
```

```
##               group                                title geo_accession
## GSM1183439 normal-H Genomic DNA from normal individual 1    GSM1183439
## GSM1183440 normal-H Genomic DNA from normal individual 2    GSM1183440
## GSM1183441 normal-H Genomic DNA from normal individual 3    GSM1183441
## GSM1183442 normal-H Genomic DNA from normal individual 4    GSM1183442
## GSM1183443 normal-H Genomic DNA from normal individual 5    GSM1183443
## GSM1183444 normal-H Genomic DNA from normal individual 6    GSM1183444
##                       tissue colon_region gender stage
## GSM1183439 colorectal mucosa        colon   male  <NA>
## GSM1183440 colorectal mucosa        colon   male  <NA>
## GSM1183441 colorectal mucosa        colon female  <NA>
## GSM1183442 colorectal mucosa        colon   male  <NA>
## GSM1183443 colorectal mucosa        colon   male  <NA>
## GSM1183444 colorectal mucosa        colon   male  <NA>
```

```r
#review what we have in our samples 
table(metadata$group)
```

```
## 
##  adenoma   cancer normal-C normal-H 
##       42       64       24       17
```


Subset 

We'll use metadata to know which columns in the `M.norm.CGI` to keep.

```r
# sample_names <- rownames(metadata[metadata$group %in% c("normal-H", "normal-C", "cancer"), ]) 
# data_sub <- M.norm.CGI[, sample_names]
# I won't do this so I can keep things in order 

# get column names for "adenoma"
adeoma_samples <- rownames(metadata[metadata$group == "adenoma", ])

# get column names for "cancer"
cancer_samples <- rownames(metadata[metadata$group == "cancer", ])

adenoma_data <- M.norm.CGI[, adeoma_samples]
ncol(adenoma_data) 
```

```
## [1] 42
```

```r
cancer_data <- M.norm.CGI[, cancer_samples]
ncol(cancer_data)  
```

```
## [1] 64
```

```r
data_sub <- cbind(adenoma_data, cancer_data)
ncol(data_sub) # instead of 147
```

```
## [1] 106
```

Data cleaning: removing rows with NA values 

```r
(a <- nrow(data_sub))
```

```
## [1] 26403
```

```r
data_sub <- data_sub[complete.cases(data_sub), ]
(a1 <- nrow(data_sub))
```

```
## [1] 26366
```

So we have removed 37 rows.

# Limma 

Time to do the real stuff ! 

Design matrix: 

```r
des <- data.frame(replicate = colnames(data_sub), 
                                    condition = factor(c(rep("adenoma", ncol(adenoma_data)), rep("cancer", ncol(cancer_data))))) # we can do this because we know the colnames are in order 


des <- model.matrix(~condition, des)
head(des)
```

```
##   (Intercept) conditioncancer
## 1           1               0
## 2           1               0
## 3           1               0
## 4           1               0
## 5           1               0
## 6           1               0
```



```r
library(limma)

limmaFit <- lmFit(data_sub, des)
EbFit <- eBayes(limmaFit)
EbFit_table <- topTable(EbFit, number = nrow(EbFit))
```

```
## Removing intercept from test coefficients
```

```r
head(EbFit_table)
```

```
##                               logFC    AveExpr         t      P.Value
## chr14:21924242-21924568   0.4306780 -2.5934715  8.817376 2.499918e-14
## chr19:49255778-49256495   0.8004475  0.2226608  8.647366 5.999543e-14
## chr19:7267384-7267688     0.3670212  3.3895078  8.503081 1.258678e-13
## chr18:77515844-77516144   0.5166213  3.6460875  8.446653 1.680831e-13
## chr5:140222122-140223108 -0.9016578  1.8883670 -8.437083 1.765281e-13
## chr3:170625831-170626841 -0.3099769 -0.4642330 -8.320590 3.203437e-13
##                             adj.P.Val        B
## chr14:21924242-21924568  6.591284e-10 22.05752
## chr19:49255778-49256495  7.909197e-10 21.22023
## chr19:7267384-7267688    9.308681e-10 20.51143
## chr18:77515844-77516144  9.308681e-10 20.23474
## chr5:140222122-140223108 9.308681e-10 20.18784
## chr3:170625831-170626841 1.407697e-09 19.61772
```


Just arrange the table a bit

```r
# keep  needed columns 
EbFit_table <- EbFit_table[, c("logFC", "t", "P.Value", "adj.P.Val")]
# rename column names 
colnames(EbFit_table)  <- c("log.fc", "test.stat", "p.value", "q.value")
head(EbFit_table)
```

```
##                              log.fc test.stat      p.value      q.value
## chr14:21924242-21924568   0.4306780  8.817376 2.499918e-14 6.591284e-10
## chr19:49255778-49256495   0.8004475  8.647366 5.999543e-14 7.909197e-10
## chr19:7267384-7267688     0.3670212  8.503081 1.258678e-13 9.308681e-10
## chr18:77515844-77516144   0.5166213  8.446653 1.680831e-13 9.308681e-10
## chr5:140222122-140223108 -0.9016578 -8.437083 1.765281e-13 9.308681e-10
## chr3:170625831-170626841 -0.3099769 -8.320590 3.203437e-13 1.407697e-09
```

Looks good, now I'll save it : 


```r
write.table(EbFit_table, file="normal_vs_cancer_santina.tsv", sep = "\t")
```


