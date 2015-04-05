#####################################################
# Ka Ming Nip


library(dplyr)
library(limma)

cutoff <- 1e-5

for (f in list.files("../data", pattern="topTable.*.Rdata")){
  load(paste("../data/", f, sep=""))
  
  print(f)
  print(nrow(subset(myTopTable, adj.P.Val < cutoff)))
}

#####################################################
# End of script
#####################################################