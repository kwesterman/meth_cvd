load("../int/rgSet.RData")
library(shinyMethyl)
summary <- shinySummarize(rgSet)
save("summary", file="../int/rgSet_shinySummary.RData")