
# Top tables from limma 

tsv files storing tables generated using limma on filtered normalized M values of CGI. 

Top table results from limma analysis using the script in the r markdown in the folder `analysis_reports`.

- `adenoma_vs_cancer.tsv` is generated using `analysis_reports/04.1limma-normal-C_vs_normal-H/04.1limma-adenoma_vs_cancer_santina.md`. 
- `normalC_vs_normalH.tsv` is generated using `analysis_reports/04.1limma-normal-C_vs_normal-H/04.1limma-normalC_vs_normalH_santina.md`.
- `normal_vs_adenoma.tsv` is generated using `analysis_reports/04.1limma-normal-C_vs_normal-H/04.1limma-normal_vs_adenoma_santina.md`.
- `normal_vs_cancer.tsv` is generated using `analysis_reports/04.1limma-normal-C_vs_normal-H/04.1limma-normal_vs_cancer_santina.md`.

The input to the markdown, btw, is `GSE48684_raw_filtered.m.norm.cgi.Rdata` This R data isn't in our repo because it's too large, but it can be generated with the script `rscripts/02norm_and_aggregate.R`
