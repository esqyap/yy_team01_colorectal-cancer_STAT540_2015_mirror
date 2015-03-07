####Add brief description of the scripts here:

[get\_data.R](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/blob/master/rscripts/get_data.R)- Beryl
  - download the raw methylation data, save the methylation data and metadata as Rdata

[filter\_probes.R](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/blob/master/rscripts/filter_probes.R) - Ka Ming and Beryl
  - filter the raw data, select probes of CG islands and non ChrX.
  - the output __../data/raw\_data\_filter__ is for normalization and raw data clustering
  
[process\_metadata.R](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/blob/master/rscripts/process_metadata.R) - Beryl
  - select the relevant columns from the raw metadata, new columns include `("group", "title", "geo_accession", "tissue", "colon_region", "gender", "stage")`
  - output is __../data/metadata.Rdata__