load("../int/rgSet.RData")
library(minfi)
estCellCounts <- estimateCellCounts(rgSet)
save("estCellCounts", file="../int/estCellCounts_minfi.RData")