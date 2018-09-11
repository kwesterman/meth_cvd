library(data.table)
system.time(a <- fread("../int/prunedGenos.csv", sep=","))
ids <- a$subjID
a <- as.matrix(a[,-1])
saveRDS(a, "../int/genos.rds", compress=F)
