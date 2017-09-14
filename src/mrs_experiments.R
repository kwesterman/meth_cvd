suppressMessages(silent <- lapply(c("tidyverse","survival","glmnet","minfi","caret","doParallel","itertools"), 
                                  library, character.only=T))

# Load methylation data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")
Mvals_full <- Mvals

# Load non-methylation (event + covariate) data
load("../int/nonMethData.RData")
nonMethData <- replace_na(nonMethData,  # Replace missing values (very few) with median
                          list(bmi=median(nonMethData$bmi, na.rm=T),
                               smoking_now=median(nonMethData$smoking_now, na.rm=T)))
# Make sure samples and their ordering are identical in methylation and metadata
Mvals <- Mvals[,match(nonMethData$sampleKey, colnames(Mvals))]
print(paste("Dimensions of event/covariate matrix:", paste(dim(nonMethData), collapse=" x ")))
print(paste("Dimensions of M-value matrix:", paste(dim(Mvals), collapse=" x ")))
stopifnot(all(colnames(Mvals)==nonMethData$sampleKey))  # Ensure identical number and order of samples for methylation and covariate data

# For continued use in testing variants on MRS model
trainTestMRS <- function(mod, meta, methMat, trainset, testset, alpha=0.5) {
  covars_frame <- model.frame(as.formula(mod), meta, na.action=na.pass)  # To allow NAs to pass through
  covars_mat <- model.matrix(as.formula(mod), covars_frame)[,-1]  # Converts factors to numeric
  design_mat <- cbind(covars_mat, t(methMat))
  cvd_surv <- Surv(time=meta$timeToEvent, event=meta$event)  # Outcome
  
  mrs.fit <- glmnet(design_mat[trainset,], cvd_surv[trainset], family="cox", alpha=alpha,  # Train MRS model
                    penalty.factor=ifelse(grepl("^cg", colnames(design_mat)), 1, 0))
  coefs <- mrs.fit$beta[,ncol(mrs.fit$beta)]  # Extract coefficients
  coefs_meth <- coefs[coefs!=0 & grepl("cg", names(coefs))]  # Methylation coefficients only
  
  mrs_calc <- as.vector(design_mat[,names(coefs_meth)] %*% coefs_meth)  # Calculate the MRS for all samples
  mrs.test <- coxph(cvd_surv~mrs, data=cbind(cvd_surv=cvd_surv, mrs=mrs_calc,  # Test the predictivity of the MRS using a basic Cox model
                                             meta), subset=testset)
  
  list(components=coefs_meth, res=summary(mrs.test)$coef)
}

# Define train/test split
set.seed(10)
trainset <- createDataPartition(factor(nonMethData$event), p=0.6)[[1]]  # Set up train/test partition
testset <- setdiff(1:nrow(nonMethData), trainset)
train_shareids <- nonMethData$shareid[trainset]

# Basic set of covariates for the MRS model
basic_model <- "~sex+age+smoking_now+CD4T+NK+Bcell+Mono+Gran+PC1_cp"

# Precompute M-value variances
Mval_variances <- apply(Mvals[,trainset], 1, var)


# # Compare variance thresholds
# print("Comparison of variance thresholds...")
# for (thresh in seq(0.2, 0.8, by=0.05)) {
#   print(paste("Using M-value variance threshold of", thresh))
#   res <- trainTestMRS(basic_model, nonMethData, Mvals[Mval_variances>thresh,],
#                       trainset, testset)
#   print(res$res)
# }


# # Compare logistic regression approach
# trainTestMRS_logistic <- function(mod, meta, methMat, trainset, testset) {
#   covars_frame <- model.frame(as.formula(mod), meta, na.action=na.pass)  # To allow NAs to pass through
#   covars_mat <- model.matrix(as.formula(mod), covars_frame)[,-1]  # Converts factors to numeric
#   design_mat <- cbind(covars_mat, t(methMat))
#   cvd_event <- meta$event  # Outcome
# 
#   mrs.fit <- glmnet(design_mat[trainset,], cvd_event[trainset], family="binomial", alpha=0.5,  # Train MRS model
#                     penalty.factor=ifelse(grepl("^cg", colnames(design_mat)), 1, 0))
#   coefs <- mrs.fit$beta[,ncol(mrs.fit$beta)]  # Extract coefficients
#   coefs_meth <- coefs[coefs!=0 & grepl("cg", names(coefs))]  # Methylation coefficients only
# 
#   mrs_calc <- as.vector(design_mat[,names(coefs_meth)] %*% coefs_meth)  # Calculate the MRS for all samples
#   mrs.test <- glm(cvd_event~mrs, data=cbind(cvd_event=cvd_event, mrs=mrs_calc,  # Test the predictivity of the MRS using a basic Cox model
#                                              meta), subset=testset)
# 
#   list(components=coefs_meth, res=summary(mrs.test)$coefficients)
# }
# 
# print("Logistic regression approach using M-value variance threshold of 0.4:")
# res <- trainTestMRS_logistic(basic_model, nonMethData, Mvals[Mval_variances>0.4,],
#                              trainset, testset)
# print(res$res)


# # Compare EWAS-based approach
# print("EWAS approach...")
# cl <- makePSOCKcluster(detectCores())
# registerDoParallel(cl)
# 
# myTry <- function(expr, CpG) {
#   # Captures model failures and returns successful result or vector of NAs
#   tryCatch(expr,
#            error = function(e) {print(e); return(c(CpG, rep(NA, 3)))},
#            warning = function(w) {print(w); return(c(CpG, rep(NA, 3)))})
# }
# 
# # Formulas for fixed-effects only models of interest
# run_cox <- function(probeData, covarData, model_spec, subset) {
#   # Given a row of the M-value matrix (corresponding to a CpG site), bind that methylation
#   # data to the covariate data and run Cox proportional hazards regression
#   CpG <- rownames(probeData)
#   modelData <- cbind(covarData, meth=as.numeric(probeData),
#                      survObj=Surv(time=covarData$timeToEvent, event=covarData$event))
#   myTry({
#     cox.fit <- coxph(as.formula(model_spec), data=modelData, subset=subset)
#     c(CpG=CpG, summary(cox.fit)$coef['meth',c('coef','z','Pr(>|z|)')])
#   }, CpG)
# }
# 
# res <- foreach(meth=iter(Mvals[Mval_variances>0.4,], by="row"),
#                .combine=rbind, .packages="survival") %dopar%
#   run_cox(meth, nonMethData, paste0("survObj",basic_model,"+meth"), trainset)
# 
# stopCluster(cl)
# 
# resDF <- res %>%
#   data.frame(stringsAsFactors=F) %>%
#   dplyr::rename(p=Pr...z..) %>%
#   mutate_at(c("coef","z","p"), as.numeric) %>%
#   mutate(fdr=p.adjust(p, method="BH")) %>%
#   arrange(p)
# 
# print("Using top 100 EWAS CpGs:")
# mrs_coefs_100 <- setNames(resDF[1:100,"coef"], resDF[1:100,"CpG"])
# mrs_100 <- as.vector(t(Mvals)[,names(mrs_coefs_100)] %*% mrs_coefs_100)
# mrs.test <- coxph(Surv(time=nonMethData$timeToEvent, event=nonMethData$event)~mrs_100,
#                   data=cbind(nonMethData, mrs_100=mrs_100), subset=testset)
# print(summary(mrs.test)$coef)
# print("...then in training set to verify:")
# mrs.test <- coxph(Surv(time=nonMethData$timeToEvent, event=nonMethData$event)~mrs_100,
#                   data=cbind(nonMethData, mrs_100=mrs_100), subset=trainset)
# print(summary(mrs.test)$coef)
# 
# print("Using top 500 EWAS CpGs:")
# mrs_coefs_500 <- setNames(resDF[1:500,"coef"], resDF[1:500,"CpG"])
# mrs_500 <- as.vector(t(Mvals)[,names(mrs_coefs_500)] %*% mrs_coefs_500)
# mrs.test <- coxph(Surv(time=nonMethData$timeToEvent, event=nonMethData$event)~mrs_500,
#                   data=cbind(nonMethData, mrs_500=mrs_500), subset=testset)
# print(summary(mrs.test)$coef)
# 
# print("Using top 1000 EWAS CpGs:")
# mrs_coefs_1000 <- setNames(resDF[1:1000,"coef"], resDF[1:1000,"CpG"])
# mrs_1000 <- as.vector(t(Mvals)[,names(mrs_coefs_1000)] %*% mrs_coefs_1000)
# mrs.test <- coxph(Surv(time=nonMethData$timeToEvent, event=nonMethData$event)~mrs_1000,
#                   data=cbind(nonMethData, mrs_1000=mrs_1000), subset=testset)
# print(summary(mrs.test)$coef)
# 
# print("Using top 2000 EWAS CpGs:")
# mrs_coefs_2000 <- setNames(resDF[1:2000,"coef"], resDF[1:2000,"CpG"])
# mrs_2000 <- as.vector(t(Mvals)[,names(mrs_coefs_2000)] %*% mrs_coefs_2000)
# mrs.test <- coxph(Surv(time=nonMethData$timeToEvent, event=nonMethData$event)~mrs_2000,
#                   data=cbind(nonMethData, mrs_2000=mrs_2000), subset=testset)
# print(summary(mrs.test)$coef)
# 
# print("Top 1000 EWAS CpGs in elastic net regression:")
# print(trainTestMRS(basic_model, nonMethData, Mvals[names(mrs_coefs_1000),], trainset, testset)$res)


# # Taking past events into account
# print("Analysis taking past events into account (variance threshold of 0.4)...")
# nonMethData <- replace_na(nonMethData, list(pastEvent=F))
# print("Adjust for past events in training only:")
# res <- trainTestMRS(paste0(basic_model,"+pastEvent"), nonMethData, Mvals[Mval_variances>0.4,],
#                     trainset, testset)
# print(res$res)
# print("Remove everyone with past events completely:")
# npeTrainSet <- trainset[nonMethData$pastEvent[trainset]==F]
# npeTestSet <- testset[nonMethData$pastEvent[testset]==F]
# print(paste("...", sum(nonMethData$event[npeTrainSet]), "events in train set."))
# print(paste("...", sum(nonMethData$event[npeTestSet]), "events in test set."))
# res <- trainTestMRS(basic_model, nonMethData, Mvals[Mval_variances>0.4,], npeTrainSet, npeTestSet)
# print(res$res)
# print("Remove everyone with past events from only the train set:")
# print(paste("...", sum(nonMethData$event[npeTrainSet]), "events in train set."))
# print(paste("...", sum(nonMethData$event[testset]), "events in test set."))
# res <- trainTestMRS(basic_model, nonMethData, Mvals[Mval_variances>0.4,], npeTrainSet, testset)
# print(res$res)
# print("Remove everyone with past events from only the test set:")
# print(paste("...", sum(nonMethData$event[trainset]), "events in train set."))
# print(paste("...", sum(nonMethData$event[npeTestSet]), "events in test set."))
# res <- trainTestMRS(basic_model, nonMethData, Mvals[Mval_variances>0.4,], trainset, npeTestSet)
# print(res$res)


# # Interpersonal stability of the MRS
# print("ICC of the MRS (trained on entire dataset):")
# res <- trainTestMRS(basic_model, nonMethData, Mvals[Mval_variances>0.4,],
#                     trainset=seq(1,nrow(nonMethData)), testset=seq(1,nrow(nonMethData)))
# mrsValues <- t(Mvals_full[names(res$components),]) %*% res$components
# library(ICC)
# load("../int/metaData.RData")
# iccData <- data.frame(shareid=factor(metaData$shareid[match(colnames(Mvals_full), metaData$sampleKey)]),
#                       mrs=mrsValues)
# iccRes <- ICCest(x=shareid, y=mrs, data=iccData)
# print(iccRes$ICC)


# # See how pre-selection of specific CpG sets affects model performance
# print("Test pre-selected CpG sets with biological relevance...")
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# mon_to_mac <- read_csv("../data/literature/mon_to_mac_DMRs_Wallner2016.csv") %>%
#   filter(!is.na(Location)) %>%
#   mutate(chr=gsub(":.*", "", Location),
#          start=as.integer(gsub("-.*", "", gsub(".*:", "", Location))),
#          end=as.integer(gsub(".*-", "", Location))) %>%
#   select(chr, start, end)
# 
# mon_var_schroder <- read_csv("../data/literature/variable_monocyte_DMRs_Schroder2017.csv") %>%
#   mutate(chr=paste0("chr", gsub(":.*", "", location)),
#          start=as.integer(gsub("-.*", "", gsub(".*:", "", location))),
#          end=as.integer(gsub(".*-", "", location))) %>%
#   select(chr, start, end)
# 
# cpgAnnot <- Locations %>%
#   data.frame(stringsAsFactors=F) %>%
#   rownames_to_column(var="cpg")
# 
# mon_to_mac_cpgs <- inner_join(cpgAnnot, mon_to_mac, by="chr") %>%
#   filter(pos>=start, pos<=end)
# 
# mon_var_schroder_cpgs <- inner_join(cpgAnnot, mon_var_schroder, by="chr") %>%
#   filter(pos>=start, pos<=end) 
# 
# # Read in Ecker variable CpGs (more than one sheet, each separately)
# library(readxl)
# variable_monSpecific <- read_excel("../data/literature/variable_cpgs_Ecker2017.xlsx", sheet="Monocytes")
# variable_neutSpecific <- read_excel("../data/literature/variable_cpgs_Ecker2017.xlsx", sheet="Neutrophils")
# variable_tcellSpecific <- read_excel("../data/literature/variable_cpgs_Ecker2017.xlsx", sheet="T cells")
# variable_monNeutSpecific <- read_excel("../data/literature/variable_cpgs_Ecker2017.xlsx", 
#                                        sheet="Monocytes + neutrophils")
# 
# print("Monocyte-to-macrophage CpGs:")
# print(trainTestMRS(basic_model, nonMethData, Mvals[rownames(Mvals) %in% mon_to_mac_cpgs$cpg,], 
#                    trainset, testset)$res)
# print("Schroder CpGs:")
# print(trainTestMRS(basic_model, nonMethData, Mvals[rownames(Mvals) %in% mon_var_schroder_cpgs$cpg,], 
#                    trainset, testset)$res)
# print("Monocyte-specific variable CpGs:")
# print(trainTestMRS(basic_model, nonMethData, Mvals[rownames(Mvals) %in% variable_monSpecific$`Probe ID`,], 
#                    trainset, testset)$res)
# print("Neutrophil-specific variable CpGs:")
# print(trainTestMRS(basic_model, nonMethData, Mvals[rownames(Mvals) %in% variable_neutSpecific$`Probe ID`,], 
#                    trainset, testset)$res)
# print("T cell-specific variable CpGs:")
# print(trainTestMRS(basic_model, nonMethData, Mvals[rownames(Mvals) %in% variable_tcellSpecific$`Probe ID`,], 
#                    trainset, testset)$res)
# print("MonocyteAndNeutrophil-specific variable CpGs:")
# print(trainTestMRS(basic_model, nonMethData, Mvals[rownames(Mvals) %in% variable_monNeutSpecific$`Probe ID`,], 
#                    trainset, testset)$res)


# # Sensitivity of results to elastic net alpha
# print("Sensitivity of results to elastic net alpha parameter (using variance threshold of 0.4)...")
# for (alpha in c(0,0.25,0.5,0.75,1)) {
#   print(paste0("Alpha = ", alpha, ":"))
#   print(trainTestMRS(basic_model, nonMethData, Mvals[Mval_variances>0.4,], 
#                      trainset, testset, alpha=alpha)$res)
# }



