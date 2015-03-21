# R scripts

### Add brief description of the scripts here:

[01-get\_data\_and\_filter.R](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/blob/master/rscripts/get_data_and_filter.R)- Beryl
  - download the raw methylation data, save the raw methylation data and metadata as Rdata
  - filter the raw data, select probes of CG islands and non ChrX.
	- select the relevant columns from the raw metadata, new columns include `("group", "title", "geo_accession", "tissue", "colon_region", "gender", "stage")`
	- beta value is not normalized
  
[02-norm\_and\_aggregate.R](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/blob/master/rscripts/norm_and_aggregate.R)- Ka Ming
  - normalize filtered beta-values
  - convert beta-values to M-values
  - aggregate both beta-values and M-values to CpG islands
  
[Misc](https://github.com/STAT540-UBC/yy_team01_colorectal-cancer_STAT540_2015/tree/master/rscripts/Misc)
  - miscellaneous/junk scripts
  - scripts that are not crucial to the work flow
  