# Two Group: Normal and Cancer
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

# get column names for "normal"
normal_samples <- rownames(metadata[metadata$group %in% c("normal-H", "normal-C"), ])

# get column names for "cancer"
cancer_samples <- rownames(metadata[metadata$group == "cancer", ])

normal_data <- M.norm.CGI[, normal_samples]
ncol(normal_data)  # 24 + 17 = 41
```

```
## [1] 41
```

```r
cancer_data <- M.norm.CGI[, cancer_samples]
ncol(cancer_data)  
```

```
## [1] 64
```

```r
data_sub <- cbind(normal_data, cancer_data)
ncol(data_sub) # instead of 147
```

```
## [1] 105
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
## [1] 26372
```

So we have removed 31 rows.

# Limma 

Time to do the real stuff ! 

Design matrix: 

```r
des <- data.frame(replicate = colnames(data_sub), 
                                    condition = factor(c(rep("normal", 41), rep("cancer", 64)))) # we can do this because we know the colnames are in order 


des <- model.matrix(~condition, des)
head(des)
```

```
##   (Intercept) conditionnormal
## 1           1               1
## 2           1               1
## 3           1               1
## 4           1               1
## 5           1               1
## 6           1               1
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
## chr8:145103285-145108027 -1.8737262 -0.9755083 -17.59271 3.731297e-33
## chr17:27939298-27940770  -1.2733401  2.3219560 -15.08254 4.324114e-28
## chr1:1476093-1476669     -1.4322686  1.1098149 -14.96461 7.622567e-28
## chr8:67344497-67344989   -2.2061031 -0.2268103 -14.61692 4.093653e-27
## chr2:210010093-210010367  0.8659652  1.1148044  14.41510 1.093036e-26
## chr21:10884732-10885087   0.6808281  1.1380020  14.11443 4.761984e-26
##                             adj.P.Val        B
## chr8:145103285-145108027 9.840177e-29 64.59087
## chr17:27939298-27940770  5.701776e-24 53.23153
## chr1:1476093-1476669     6.700744e-24 52.67768
## chr8:67344497-67344989   2.698946e-23 51.03475
## chr2:210010093-210010367 5.765107e-23 50.07433
## chr21:10884732-10885087  2.093051e-22 48.63449
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
## chr8:145103285-145108027 -1.8737262 -17.59271 3.731297e-33 9.840177e-29
## chr17:27939298-27940770  -1.2733401 -15.08254 4.324114e-28 5.701776e-24
## chr1:1476093-1476669     -1.4322686 -14.96461 7.622567e-28 6.700744e-24
## chr8:67344497-67344989   -2.2061031 -14.61692 4.093653e-27 2.698946e-23
## chr2:210010093-210010367  0.8659652  14.41510 1.093036e-26 5.765107e-23
## chr21:10884732-10885087   0.6808281  14.11443 4.761984e-26 2.093051e-22
```

Looks good, now I'll save it : 


```r
write.table(EbFit_table, file="normal_vs_cancer_santina.tsv", sep = "\t")
```


Normal and adenoma

normal-H and C 

adenoma and cancer 
