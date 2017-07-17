suppressMessages(silent <- lapply(c("tidyverse","doParallel","minfi","itertools","survival","gtools","coxme"), 
                                  library, character.only=T))

# Load methylation and sample data
load("../int/chaoMvals.RData")  # Loads Mvals object
print("Methylation data loaded.")

# Load non-methylation (event + covariate) data and sample similarity matrix
load("../int/nonMethData.RData")
nonMethData$survObj <- Surv(time=nonMethData$timeToEvent, event=nonMethData$event, type="right")
load("../int/Kmat.RData")

# Sanity checks
final_shareid_set <- intersect(nonMethData$shareid, colnames(Mvals))
nonMethData <- nonMethData[match(final_shareid_set, nonMethData$shareid),]
Mvals <- Mvals[,final_shareid_set]
# Kmat <- Kmat[final_shareid_set, final_shareid_set]
print(paste("Dimensions of event/covariate matrix:", dim(nonMethData)))
print(paste("Dimensions of M-value matrix:", dim(Mvals)))
stopifnot(ncol(Mvals)==nrow(nonMethData),  # Ensure same # of samples in methylation and covariate data
          all(colnames(Mvals)==nonMethData$shareid))  # Ensure identical order of samples for methylation and covariate data

## Run regressions
Mvals <- ilogit2(Mvals)  # convert back to betas for this test
Mvals <- Mvals[1:10000,]
cl <- makePSOCKcluster(detectCores(), outfile="progress.txt")
registerDoParallel(cl)

myTry <- function(expr, CpG) {
  # Captures model failures and returns successful result or vector of NAs
  tryCatch(expr,
           error = function(e) {print(e); return(c(CpG, rep(NA, 3)))},
           warning = function(w) {print(w); return(c(CpG, rep(NA, 3)))})
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
  # unadjusted="survObj~meth",
  # basic="survObj~meth+sex+age+smoking_now+bmi",
  # basic_wbc="survObj~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran",
  # basic_wbc_20PC=paste0("survObj~meth+sex+age+smoking_now+bmi+", paste0("PC",1:20,collapse="+"))
  basic_wbc_20CPA=paste0("survObj~meth+sex+age+smoking_now+bmi+", paste0("PC",1:20,"_cp",collapse="+"))
  # basic_wbc_20PC_20CPA=paste0("survObj~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+",
  #                             paste0("PC",1:20,collapse="+"),"+",paste0("PC",1:20,"_cp",collapse="+"))
  # basic_wbc_10SVs=paste0("survObj~meth+sex+age+smoking_now+bmi+", paste0("SV",1:10,collapse="+")),
  # basic_wbc_lipids="survObj~meth+sex+age+smoking_now+bmi+SBP+TOT_CHOL+HDL_CHOL+TRIG+GLUCOSE"
  # all=paste0("survObj~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+",
  #                  paste0("PC",1:10,"_cp",collapse="+"), "+", paste0(colnames(SVs)[1:10], collapse="+"))
)
res_fixed <- foreach(mod=model_list, .final=function(l) setNames(l, names(model_list))) %:%
  foreach(meth=enumerate(iter(Mvals, by="row")), .combine=rbind, .packages=c("survival")) %dopar% {
    if(meth$index %% 10000 == 0) print(paste0("Model ", mod, ", CpG "))
    run_cox(meth$value, mod)
  }

# run_coxme <- function(probeData, model_spec, Kmat=NULL) {
#   # Given a row of the M-value matrix (corresponding to a CpG site), bind that methylation
#   # data to the covariate data and run Cox proportional hazards regression
#   CpG <- rownames(probeData)
#   modelData <- cbind(nonMethData, meth=as.numeric(probeData))
#   varlist <- if(is.null(Kmat)) NULL else coxmeMlist(Kmat)
#   myTry({
#     coxme.fit <- coxme(model_spec, data=modelData, varlist=varlist)
#     coef <- unname(coef(coxme.fit)['meth'])
#     se <- sqrt(diag(vcov(coxme.fit))[1])
#     p <- 1 - pchisq((coef/se)^2, 1)
#     c(CpG=CpG, coef=coef, p=p)
#   }, CpG)
# }
# me_model_list <- list(
#   # kmat=survObj~meth+sex+age+smoking_now+bmi+(1|sampleKey),  # Sim. matrix is labeled by sample key, not shareid
#   slideRowCol_RE=survObj~meth+sex+age+smoking_now+bmi+(1|Sentrix_ID)+(1|Sentrix_Row)+(1|Sentrix_Col)
# )
# res_mixed <- foreach(mod=me_model_list, .final=function(l) setNames(l, names(me_model_list))) %:%
#   foreach(meth=enumerate(iter(Mvals, by="row")), .combine=rbind, .packages="coxme") %dopar% {
#     if(meth$index %% 10000 == 0) print(paste0("Model ", mod, ", CpG ", meth$index))
#     run_coxme(meth$value, mod)
#   }
# 
# stopImplicitCluster()
res <- c(res_fixed)

print("Saving results...")
save("res", file=paste0("../int/chaoewasRes.RData"))