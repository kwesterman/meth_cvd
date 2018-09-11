library(data.table)
system.time(a <- fread("../int/unprunedGenos.csv", sep=","))
ids <- a$subjID
a <- as.matrix(a[,-1])
rownames(a) <- ids
saveRDS(a, "../int/unprunedGenos.rds", compress=F)
