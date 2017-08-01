suppressMessages(silent <- lapply(c("tidyverse","doParallel","itertools","survival","gtools","coxme"), 
                                  library, character.only=T))

# Load methylation and sample data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")

# Load non-methylation (event + covariate) data and sample similarity matrix
load("../int/nonMethData.RData")
nonMethData$survObj <- Surv(time=nonMethData$timeToEvent, event=nonMethData$event, type="right")

# Sanity checks
final_id_set <- intersect(nonMethData$sampleKey, colnames(Mvals))
nonMethData <- nonMethData[match(final_id_set, nonMethData$sampleKey),]
Mvals <- Mvals[,final_id_set]
print(paste("Dimensions of event/covariate matrix:", dim(nonMethData)))
print(paste("Dimensions of M-value matrix:", dim(Mvals)))
stopifnot(ncol(Mvals)==nrow(nonMethData),  # Ensure same # of samples in methylation and covariate data
          all(colnames(Mvals)==nonMethData$sampleKey))  # Ensure identical order of samples for methylation and covariate data

## Run regressions
# Mvals <- Mvals[sample(1:nrow(Mvals), 50000),]
cl <- makePSOCKcluster(detectCores(), outfile="progress.txt")
registerDoParallel(cl)

myTry <- function(expr, CpG) {
  # Captures model failures and returns successful result or vector of NAs
  tryCatch(expr,
           error = function(e) {print(e); return(c(CpG, rep(NA, 2)))},
           warning = function(w) {print(w); return(c(CpG, rep(NA, 2)))})
}

# Formulas for fixed-effects only models of interest
run_cox <- function(probeData, model_spec) {
  # Given a row of the M-value matrix (corresponding to a CpG site), bind that methylation
  # data to the covariate data and run Cox proportional hazards regression
  CpG <- rownames(probeData)
  modelData <- cbind(nonMethData, meth=as.numeric(probeData))
  myTry({
    cox.fit <- coxph(as.formula(model_spec), data=modelData)
    c(CpG=CpG, summary(cox.fit)$coef['meth',c('coef','Pr(>|z|)')])
  }, CpG)
} 
model_list <- list(
  basic_wbc_20PC=paste0("survObj~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                        paste0("PC",1:20,collapse="+"))
)
res <- foreach(mod=model_list, .final=function(l) setNames(l, names(model_list))) %:%
  foreach(meth=enumerate(iter(Mvals, by="row")), .combine=rbind, .packages=c("survival")) %dopar% {
    if(meth$index %% 10000 == 0) print(paste0("Model ", mod, ", CpG "))
    run_cox(meth$value, mod)
  }
stopImplicitCluster()

print("Saving results...")
save("res", file=paste0("../int/ewasRes.RData"))