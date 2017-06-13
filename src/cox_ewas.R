suppressMessages(silent <- lapply(c("tidyverse","doParallel","itertools","survival","gtools"), 
                                  library, character.only=T))

# Load methylation and sample data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")

# Load all covariate data and prepare survival object
load("../int/phenoData.RData")
load("../int/sampleData.RData")
load("../int/eventData.RData")
load("../int/estCellCounts.RData")
load("../int/CP_PCs.RData")
load("../int/SVs.RData")

CP_PCs <- CP_PCs[rownames(estCellCounts),]  # B/c used rgSet as input, must trim samples that were filtered
colnames(SVs) <- paste("SV", 1:ncol(SVs), sep="_")
technicalCovars <- data.frame(cbind(estCellCounts, CP_PCs, SVs))
technicalCovars$sampleKey <- rownames(technicalCovars)

survData <- sampleData %>%
  left_join(phenoData, by="shareid") %>%
  left_join(eventData, by="shareid") %>%
  left_join(technicalCovars, by="sampleKey") %>%
  dplyr::mutate(event=!is.na(timeToEvent),  # event == TRUE when a specific time of event exists, otherwise FALSE
                time=na.replace(timeToEvent, 1000)) %>%  ## NOTE: this 800 needs to be switched for an *actual* number
  dplyr::slice(match(colnames(Mvals), sampleKey))  # Re-order so rows match methylation columns
survData$survObj <- Surv(time=survData$time, event=survData$event, type="right")  # all subjects without an event time get max. # days as time to censor

# Sanity checks
print("Dimensions of survData:")
dim(survData)
print("Dimensions of M-value matrix:")
dim(Mvals)
stopifnot(ncol(Mvals)==nrow(survData),  # Ensure same # of samples in methylation and covariate data
          all(colnames(Mvals)==survData$sampleKey))  # Ensure identical order of samples for methylation and covariate data

# Formulas for models of interest
model_list <- list(
  unadjusted="survObj~meth",
  basic="survObj~meth+sex+age+smoking_now+bmi",
  basic_wbc="survObj~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran",
  basic_wbc_10PC=paste0("survObj~meth+sex+age+smoking_now+bmi+", paste0("PC",1:10,"_cp",collapse="+")),
  basic_wbc_10SVs=paste0("survObj~meth+sex+age+smoking_now+bmi+", paste0(colnames(SVs)[1:10], collapse="+")),
  basic_wbc_50SVs=paste0("survObj~meth+sex+age+smoking_now+bmi+", paste0(colnames(SVs)[1:50], collapse="+")),
  all=paste0("survObj~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+",
                   paste0("PC",1:10,"_cp",collapse="+"), "+", paste0(colnames(SVs)[1:10], collapse="+"))
)

run_cox <- function(probeData, model_spec) {
  # Given an index (corresponding to a methylation site), bind that methylation data
  # to the covariate data and run Cox proportional hazards regression
  # Return vector of NAs if the regression fails
  CpG <- rownames(probeData)
  modelData <- cbind(survData, meth=as.numeric(probeData))
  returnVal <- tryCatch({
    cox.fit <- coxph(as.formula(model_spec), data=modelData)
    c(CpG=CpG, summary(cox.fit)$coef['meth',c('coef','z','Pr(>|z|)')])
  }, error=function(e) {
    print(e)
    return(c(CpG, NA, NA, NA))
  }, warning=function(w) {
    print(w)
    return(c(CpG, NA, NA, NA))
  })
  returnVal
} 

Mvals <- Mvals[1:100000,]
print(paste(detectCores(), "cores detected. Running regressions..."))
cl <- makePSOCKcluster(detectCores(), outfile="progress.txt")
registerDoParallel(cl)
res_qqCompare <- foreach(mod=model_list, .final=function(l) setNames(l, names(model_list))) %:%
  foreach(meth=enumerate(iter(Mvals, by="row")), .combine=rbind, .packages="survival") %dopar% {
    if(meth$index %% 10000 == 0) print(paste0("Model ", mod, ", CpG ", meth$index, " of ", nrow(Mvals), "."))
    run_cox(meth$value, mod)
  }
stopImplicitCluster()

print("Saving results...")
save("res_qqCompare", file=paste0("../output/res_qqCompare.RData"))