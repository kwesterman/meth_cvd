library(tidyverse)
library(igraph)

load("../int/Mvals.RData")
Mval_variances <- apply(Mvals, 1, var)
Mvals <- Mvals[Mval_variances>1,]

system.time(corMat <- cor(t(Mvals[1:5000,])))
corMat <- corMat[abs(corMat)>0.7]
