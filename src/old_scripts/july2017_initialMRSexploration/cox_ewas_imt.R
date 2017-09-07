suppressMessages(silent <- lapply(c("tidyverse","doParallel","itertools","survival","gtools","coxme"), 
                                  library, character.only=T))

# Load methylation and sample data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")

# Load non-methylation (event + covariate) data and sample similarity matrix
load("../int/nonMethData2.RData")

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
           error = function(e) {print(e); return(c(CpG, rep(NA, 3)))},
           warning = function(w) {print(w); return(c(CpG, rep(NA, 3)))})
}
model_list <- list(
  # glu_WBC=paste0("GLUCOSE~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
  #                         paste0("PC",1:20,collapse="+")),
                   # glu_noGran=paste0("GLUCOSE~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+", 
                   #        paste0("PC",1:20,collapse="+")),
                   imt=paste0("IMT_mean~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                          paste0("PC",1:20,collapse="+")))
run_lm <- function(probeData, mod) {
  # Given a row of the M-value matrix (corresponding to a CpG site), bind that methylation
  # data to the covariate data and run Cox proportional hazards regression
  CpG <- rownames(probeData)
  modelData <- cbind(nonMethData, meth=as.numeric(probeData))
  myTry({
    lm.fit <- lm(as.formula(mod), 
                     data=modelData)
    c(CpG=CpG, summary(lm.fit)$coefficients['meth',c('Estimate','Pr(>|t|)')])
  }, CpG)
} 
res <- foreach(mod=model_list, .final=function(l) setNames(l, names(model_list))) %:%
  foreach(meth=enumerate(iter(Mvals, by="row")), .combine=rbind) %dopar% {
    if(meth$index %% 10000 == 0) print(paste0("Model ", mod, ", CpG "))
    run_lm(meth$value, mod)
  }
stopImplicitCluster()

print("Saving results...")
save("res", file=paste0("../int/ewasRes_IMTonly.RData"))