suppressMessages(silent <- lapply(c("tidyverse","doParallel","itertools","survival","gtools","coxme"), 
                                  library, character.only=T))

### CREATE MRS

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
res <- foreach(meth=enumerate(iter(Mvals, by="row")), .combine=rbind, .packages=c("survival")) %dopar% {
  if(meth$index%%20000==0) print(paste("Reached probe index", meth$index))
  run_cox(meth$value, paste0("survObj~meth+sex+age+smoking_now+bmi+", paste0("PC",1:20,collapse="+")))
}
stopImplicitCluster()
resDF <- data.frame(res, stringsAsFactors=F)
names(resDF) <- c("CpG","beta","p")
resDF$p <- as.numeric(resDF$p)
resDF$fdr <- p.adjust(resDF$p, method="BH")
keep_CpGs <- na.omit(resDF$CpG[resDF$fdr < 0.05])
print(head(keep_CpGs))
print(paste(length(keep_CpGs), "component CpGs."))
Mvals_keep <- Mvals[keep_CpGs,]
allModelData <- cbind(nonMethData, t(Mvals_keep))
model_spec <- paste0("survObj~sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                     paste0("PC",1:20,"_cp",collapse="+"),"+",paste0(keep_CpGs,collapse="+"))
mrs.fit <- coxph(as.formula(model_spec), data=allModelData)
coefs <- summary(mrs.fit)$coef[keep_CpGs,'coef']  # Extract regression coefficients

## Calculate MRS
source("helpers.R")
nonMethData$frs <- calc_FRS(nonMethData[,-which(colnames(nonMethData)=="survObj")])
nonMethData$mrs <- as.vector(as.matrix(allModelData[,keep_CpGs]) %*% coefs)
mrsDat <- nonMethData
save("mrsDat", file="../int/mrsDat.RData")

## Calculate Zhang MRS
zhang_MRS <- calc_zhang_mrs(ilogit2(Mvals))
mrsDat <- left_join(dplyr::select(mrsDat, -survObj), zhang_MRS, by="sampleKey")  # Adds Zhang MRS

## Calculate cumulative exposures

library(readxl)
fram_lipid_varNames <- data.frame(exam=1:7,
                                  chol=c("A9","B352","C429","D448","E667","F726","G704"),
                                  hdl=c("A10","B355","C431","D449","E668","F725","G703"),
                                  ldl=c("A12","B357","C442",NA,NA,NA,NA),
                                  tg=c("A13","B358","C433","D451","E670","F727","G706"),
                                  glu=c("A31","B737","C434","D452","E671","F724","G705"),
                                  sysBP=c("A55","B24","C184","D192","E485","F476","G271"))

exams <- apply(fram_lipid_varNames, 1, function(row) {
  fileName <- grep(paste0("c1\\.ex1_", row["exam"]), list.files("../data/phenotypes/"), value=T)
  phenDF <- read_excel(paste0("../data/phenotypes/", fileName), sheet=2)
  cols <- c(shareid="shareid", na.omit(row)[-1])
  setNames(phenDF[cols], names(cols))
})
names(exams) <- paste0("exam", 1:7)
allExams <- bind_rows(exams, .id="exam") %>%
  mutate(nonHDL=chol-hdl)
cumVals <- allExams %>%
  group_by(shareid) %>%
  summarise_at(vars(-exam), mean)
names(cumVals)[-1] <- paste0(names(cumVals)[-1], "_cum")


mrsDat <- inner_join(mrsDat, cumVals, by="shareid")
save("mrsDat", file="../int/mrsDat2.RData")



