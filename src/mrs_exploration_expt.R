suppressMessages(silent <- lapply(c("tidyverse","doParallel","itertools","survival","gtools","coxme","glmnet"), 
                                  library, character.only=T))
source("helpers.R")

### CREATE MRS

# Load methylation and sample data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")

# Load non-methylation (event + covariate) data and sample similarity matrix
load("../int/nonMethData2.RData")
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
# myTry <- function(expr, CpG) {
#   # Captures model failures and returns successful result or vector of NAs
#   tryCatch(expr,
#            error = function(e) {print(e); return(c(CpG, rep(NA, 3)))},
#            warning = function(w) {print(w); return(c(CpG, rep(NA, 3)))})
# }

# Formulas for fixed-effects only models of interest
# run_cox <- function(allProbeData, model_spec) {
#   # Given a row of the M-value matrix (corresponding to a CpG site), bind that methylation
#   # data to the covariate data and run Cox proportional hazards regression
#   CpG <- rownames(probeData)
#   modelData <- cbind(nonMethData, meth=as.numeric(probeData))
#     design_mat_form <- paste0("~sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
#                               paste0("PC",1:20,collapse="+"))
#     x <- model.frame(as.formula(design_mat_form), modelData, na.action=na.pass)
#     x <- apply(x, 2, scale)
#     y <- cbind(time=trainSet$timeToEvent, status=trainSet$event)
#     mrs.fit <- glmnet(x, y, family="cox", alpha=0.5)
#     # mrs.fit <- glmnet(x, y, family="binomial", alpha=0.5)
#     coefs <- mrs.fit$beta[,ncol(mrs.fit$beta)]
#     coefs[coefs!=0]
# } 
# cl <- makePSOCKcluster(detectCores(), outfile="progress.txt")
# registerDoParallel(cl)
# res <- foreach(meth=enumerate(iter(Mvals, by="row")), .combine=rbind, .packages=c("survival","glmnet")) %dopar% {
#   if(meth$index%%20000==0) print(paste("Reached probe index", meth$index))
#   run_cox(meth$value, paste0("survObj~meth+sex+age+smoking_now+bmi+", paste0("PC",1:20,collapse="+")))
# }
# stopImplicitCluster()
# resDF <- data.frame(res, stringsAsFactors=F)
# names(resDF) <- c("CpG","beta","p")
# resDF$p <- as.numeric(resDF$p)
# resDF$fdr <- p.adjust(resDF$p, method="BH")
# keep_CpGs <- na.omit(resDF$CpG[resDF$fdr < 0.05])
# print(head(keep_CpGs))
# print(paste(length(keep_CpGs), "component CpGs."))
# Mvals_keep <- Mvals[keep_CpGs,]
# allModelData <- cbind(nonMethData, t(Mvals_keep))
# model_spec <- paste0("survObj~sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
#                      paste0("PC",1:20,"_cp",collapse="+"),"+",paste0(keep_CpGs,collapse="+"))
# mrs.fit <- coxph(as.formula(model_spec), data=allModelData)
# coefs <- summary(mrs.fit)$coef[keep_CpGs,'coef']  # Extract regression coefficients

# vars <- apply(Mvals, 1, var)
# Mvals <- Mvals[vars>0.5,]
# Mvals <- Mvals[1:10000,]
nonMethData$frs <- calc_FRS(nonMethData[,-which(colnames(nonMethData)=="survObj")])

# modelData <- cbind(nonMethData, t(Mvals))
design_mat_form <- paste0("~sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                          paste0("PC",1:20,collapse="+"))
# design_mat_form <- paste0("~sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
#                           paste0("PC",1:20,collapse="+"),"+",paste(rownames(Mvals),collapse="+"))
covars <- model.frame(as.formula(design_mat_form), nonMethData, na.action=na.pass)
# x <- model.frame(as.formula(design_mat_form), modelData, na.action=na.pass)
x <- cbind(covars, t(Mvals))[complete.cases(covars),]
x <- apply(x, 2, scale)
cvd_surv <- cbind(time=nonMethData$timeToEvent, status=nonMethData$event)[complete.cases(covars),]
frs <- nonMethData$frs[complete.cases(covars)]
imt_mean <- nonMethData$IMT_mean

print("CVD:")
mrs.fit <- glmnet(x, cvd_surv, family="cox", alpha=0.5)
coefs_cvd <- mrs.fit$beta[,ncol(mrs.fit$beta)]
coefs_cvd <- coefs_cvd[coefs_cvd!=0]
keep_CpGs_cvd <- names(coefs_cvd)
print(paste(length(coefs_cvd), "nonzero coefficients"))
print(head(coefs_cvd, 25))

print("FRS:")
mrs.fit <- glmnet(x, frs, family="gaussian", alpha=0.5)
coefs_frs <- mrs.fit$beta[,ncol(mrs.fit$beta)]
coefs_frs <- coefs_frs[coefs_frs!=0]
keep_CpGs_frs <- names(coefs_frs)
print(paste(length(coefs_frs), "nonzero coefficients"))
print(head(coefs_frs, 25))

print("IMT:")
mrs.fit <- glmnet(x, imt_mean, family="gaussian", alpha=0.5)
coefs_imt <- mrs.fit$beta[,ncol(mrs.fit$beta)]
coefs_imt <- coefs_imt[coefs_imt!=0]
keep_CpGs_imt <- names(coefs_imt)
print(paste(length(coefs_imt), "nonzero coefficients"))
print(head(coefs_imt, 25))

## Calculate MRS
nonMethData$mrs_cvd <- as.vector(t(Mvals[keep_CpGs_cvd,]) %*% coefs_cvd)
nonMethData$mrs_frs <- as.vector(t(Mvals[keep_CpGs_frs,]) %*% coefs_frs)
nonMethData$mrs_imt <- as.vector(t(Mvals[keep_CpGs_imt,]) %*% coefs_imt)
mrsDat <- nonMethData
# save("mrsDat", file="../int/mrsDat.RData")

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

foodLib <- c(satfat="NUT_SATFAT", monfat="NUT_MONFAT", polyfat="NUT_POLY", beans="FFD60",
             alcohol="NUT_ALCO", sucrose="NUT_SUCR", sodium="NUT_SODIUM", n3="NUT_OMEGA",
             folate="NUT_FOLEQ", sugars="NUT_SUGTOT", transfat="NUT_TRN02", ssb="FFD145",
             proantho="NUT_PROMON", ecg="NUT_UECG", kale="FFD66", prcsdmeat="FFD78")

library(readxl)
dietDat_all <- read_excel("../data/diet/phs000007.v28.pht002350.v4.p10.c1.vr_ffreq_ex08_1_0615s.HMB-IRB-MDS_ex8_diet.xlsx", sheet=2)
dietDat <- dplyr::select(dietDat_all, shareid, foodLib)
names(dietDat) <- c("shareid", names(foodLib))

# load("../int/mrsDat2.RData")
dietDat <- left_join(dplyr::select(mrsDat, shareid, mrs), dietDat, by="shareid")
save("mrsDat", file="../int/mrsDat_penalized.RData")


# dietMRScors <- cor(dietDat[-1], use="pairwise.complete.obs")
# library(gplots)
# jpeg("../output/dietMRS.jpg")
# heatmap.2(dietMRScors, Rowv=F, Colv=F)




