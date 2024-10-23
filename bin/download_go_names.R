#!/usr/bin/Rscript
library(GO.db)
go <- keys(GO.db, keytype="GOID")
df <- select(GO.db, columns=c("GOID","TERM"), keys=go, keytype="GOID")
write.table(df, file="GO_to_name", sep="\t", quote=F)
