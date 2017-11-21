suppressMessages(silent <- lapply(c("tidyverse","survival","glmnet","minfi",
                                    "caret","doParallel","itertools","pROC","goseq"), library, character.only=T))

train_mrs <- function(mod, meta, meth, trainset, alpha=0.5) {
  covars_frame <- model.frame(as.formula(mod), meta, na.action=na.pass)  # To allow NAs to pass through
  covars_mat <- model.matrix(as.formula(mod), covars_frame)[,-1]  # Converts factors to numeric
  design_mat <- cbind(covars_mat, t(meth))
  cvd_surv <- Surv(time=meta$time, event=meta$event)  # Outcome
  mrs.fit <- glmnet(design_mat[trainset,], cvd_surv[trainset], family="cox", alpha=alpha,  # Train MRS model
                    penalty.factor=ifelse(grepl("^cg", colnames(design_mat)), 1, 0))
  coefs <- mrs.fit$beta[,ncol(mrs.fit$beta)]  # Extract coefficients
  list(fit=mrs.fit,  # Output is a list containing the model fit and the CpG model coefficients
       coefs=coefs[grepl("cg", names(coefs))])
}

calc_mrs <- function(coefs, meth) as.vector(t(meth[names(coefs),]) %*% coefs)

test_mrs <- function(meta, mrs, testset) {
  meta$mrs <- scale(mrs)
  meta$cvd_surv <- Surv(time=meta$time, event=meta$event)  # Outcome
  mrs.test <- coxph(cvd_surv~mrs, data=meta, subset=testset) # Test the MRS using a Cox model
  mrs.coefs <- summary(mrs.test)$coef
  summary(mrs.test)$coef[,c("exp(coef)","Pr(>|z|)")]
}

trainTestMRS <- function(mod, meta, meth, trainset, testset, alpha=0.5) {
  mrs_coefs <- train_mrs(mod, meta, meth, trainset, alpha=alpha)$coefs
  mrs_values <- calc_mrs(mrs_coefs, meth)
  mrs_test <- test_mrs(meta, mrs_values, testset)
  list(components=mrs_coefs, values=mrs_values, res=mrs_test)
}

betas <- readRDS("../int/betas.qc.norm.filt.rds")  # Loads beta values object

# Load non-methylation (event + covariate) data
nonMethData <- readRDS("../int/nonMethData.rds")
nonMethData <- replace_na(nonMethData,  # Replace missing values (very few) with median
                          list(bmi=median(nonMethData$bmi, na.rm=T),
                               smk_now=median(nonMethData$smk_now, na.rm=T),
                               ht_med=median(nonMethData$ht_med, na.rm=T),
                               lipid_med=median(nonMethData$lipid_med, na.rm=T),
                               dm_med=median(nonMethData$dm_med, na.rm=T)))

# # Load ComBat-specific PCA
# load("../int/PCA_combat.RData")
# combat_PCs <- setNames(PCs, c("sampleKey", paste0("PC",1:20,"_combat")))
# nonMethData <- inner_join(nonMethData, combat_PCs, by="sampleKey")

# Make sure samples and their ordering are identical in methylation and metadata
betas <- betas[,match(nonMethData$sampleKey, colnames(betas))]
print(paste("Dimensions of event/covariate matrix:", paste(dim(nonMethData), collapse=" x ")))
print(paste("Dimensions of Beta-value matrix:", paste(dim(betas), collapse=" x ")))
stopifnot(all(colnames(betas)==nonMethData$sampleKey))  # Ensure identical number and order of samples for methylation and covariate data

metaData <- readRDS("../int/metaData.rds")
fhsCPACOR <- readRDS("../int/cpacor.fit.fhs.rds")
fhsCPPCs <- setNames(cbind(metaData$sampleKey[metaData$study=="fhs"], data.frame(fhsCPACOR$x)),
                     c("sampleKey", paste0("cp", colnames(fhsCPACOR$x),"_indiv")))
whiCPACOR <- readRDS("../int/cpacor.fit.whi.rds")
whiCPPCs <- setNames(cbind(metaData$sampleKey[metaData$study=="whi"], data.frame(whiCPACOR$x)),
                     c("sampleKey", paste0("cp", colnames(whiCPACOR$x),"_indiv")))
allCPPCs <- rbind(fhsCPPCs, whiCPPCs)
nonMethData <- nonMethData %>% inner_join(allCPPCs, by="sampleKey")


cpacor_adj_model <- paste0("~",paste0("cpPC",1:10,"_indiv",collapse="+"))

whiSet <- which(nonMethData$study=="whi")
fhsSet <- which(nonMethData$study=="fhs")

train_mrs_withLambda <- function(mod, meta, meth, trainset, lambda, alpha=0.5) {
  covars_frame <- model.frame(as.formula(mod), meta, na.action=na.pass)  # To allow NAs to pass through
  covars_mat <- model.matrix(as.formula(mod), covars_frame)[,-1]  # Converts factors to numeric
  design_mat <- cbind(covars_mat, t(meth))
  cvd_surv <- Surv(time=meta$time, event=meta$event)  # Outcome
  mrs.fit <- glmnet(design_mat[trainset,], cvd_surv[trainset], family="cox", alpha=alpha,  # Train MRS model
                    penalty.factor=ifelse(grepl("^cg", colnames(design_mat)), 1, 0), lambda=lambda)
  coefs <- mrs.fit$beta[,ncol(mrs.fit$beta)]  # Extract coefficients
  list(fit=mrs.fit,  # Output is a list containing the model fit and the CpG model coefficients
       coefs=coefs[grepl("cg", names(coefs))])
}

trainTestMRS_withLambda <- function(mod, meta, meth, trainset, testset, lambda, alpha=0.5) {
  mrs_coefs <- train_mrs_withLambda(mod, meta, meth, trainset, lambda, alpha=alpha)$coefs
  mrs_values <- calc_mrs(mrs_coefs, meth)
  mrs_test <- test_mrs(meta, mrs_values, testset)
  list(components=mrs_coefs, values=mrs_values, res=mrs_test)
}


load("../cache/mrs_experiments/ewas-cpacorOnly_ca68f367b16930007fea044e45722b9c.RData")
top250 <- ewasResDF_cpacorOnly[1:250,"CpG"]



library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno450k <- data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19), stringsAsFactors=F)
cvd_gwas_associations <- read_tsv("../data/literature/gwas-association-downloaded_2017-09-18-cardiovascular disease.tsv")
cvd_gwas_cpgs_byPos <- cvd_gwas_associations %>%
  mutate(chr=paste0("chr", CHR_ID),
         snp_pos=CHR_POS) %>%
  select(chr, snp_pos) %>%
  inner_join(select(anno450k, Name, chr, pos), by="chr") %>%
  mutate(snp_pos=as.numeric(snp_pos)) %>%
  filter(pos>=snp_pos-1000, pos<=snp_pos+1000) %>%
  dplyr::rename(cpg=Name) %>%
  distinct(cpg)
cvd_cpgs <- cvd_gwas_cpgs_byPos$cpg
cvdRes <- trainTestMRS(cpacor_adj_model, nonMethData, betas[rownames(betas) %in% cvd_cpgs,], whiSet, fhsSet)$res


train_mrs <- function(mod, meta, meth, trainset, 
                      alpha=0.5, nfolds=5, lmr=0.1, lambda.type="lambda.min", parallel=T) {
  covars_frame <- model.frame(as.formula(mod), meta, na.action=na.pass)  # To allow NAs to pass through
  covars_mat <- model.matrix(as.formula(mod), covars_frame)[,-1]  # Converts factors to numeric
  design_mat <- cbind(covars_mat, t(meth))
  cvd_surv <- Surv(time=meta$time, event=meta$event)  # Outcome
  mrs.fit.cv <- cv.glmnet(design_mat[trainset,], cvd_surv[trainset], family="cox", alpha=alpha,  # Train MRS model
                          penalty.factor=ifelse(grepl("^cg", colnames(design_mat)), 1, 0), 
                          nfolds=nfolds, lambda.min.ratio=lmr, parallel=parallel)
  coefs <- coef(mrs.fit.cv, s=lambda.type)
  cpgCoefs <- coefs[grepl("^cg", rownames(coefs)),]
  list(fit=mrs.fit.cv,  # Output is a list containing the model fit and the CpG model coefficients
       allCoefs=coefs,
       cpgCoefs=cpgCoefs)
}

calc_mrs <- function(coefs, meth) as.vector(t(meth[names(coefs),]) %*% coefs)

test_mrs <- function(meta, mrs, testset) {
  meta$mrs <- scale(mrs)
  meta$cvd_surv <- Surv(time=meta$time, event=meta$event)  # Outcome
  mrs.test <- coxph(cvd_surv~mrs, data=meta, subset=testset) # Test the MRS using a Cox model
  summary(mrs.test)$coef[,c("exp(coef)","Pr(>|z|)")]
}

trainTestMRS <- function(mod, meta, meth, trainset, testset, alpha=0.5, nfolds=5, lmr=0.1, lt="lambda.min") {
  mrs_fit <- train_mrs(mod, meta, meth, trainset, alpha=alpha, nfolds=nfolds, lmr=lmr, lambda.type="lambda.min")
  mrs_values <- calc_mrs(mrs_fit$cpgCoefs, meth)
  mrs_test <- test_mrs(meta, mrs_values, testset)
  list(fit=mrs_fit, components=mrs_fit$cpgCoefs, values=mrs_values, res=mrs_test)
}

my_test <- function(mod, meta, meth, trainset, testset, alpha=0.5, nfolds=5, lmr=0.1, lt="lambda.min", parallel=T) {
  set.seed(1)
  mrs_fit <- train_mrs(mod, meta, meth, trainset, alpha=alpha, nfolds=nfolds, lmr=lmr, lambda.type=lt, parallel=parallel)
  print(paste("Lambda.min:", mrs_fit$fit$lambda.min))
  print(paste("Lambda.1se:", mrs_fit$fit$lambda.1se))
  mrs_values <- calc_mrs(mrs_fit$cpgCoefs, meth)
  mrs_test <- test_mrs(meta, mrs_values, testset)
  print(mrs_test)
  list(fit=mrs_fit, components=mrs_fit$cpgCoefs, values=mrs_values, res=mrs_test)
}

train_mrs_binary <- function(mod, meta, meth, trainset, 
                             alpha=0.5, nfolds=5, lmr=0.1, lambda.type="lambda.min", parallel=F) {
  covars_frame <- model.frame(as.formula(mod), meta, na.action=na.pass)  # To allow NAs to pass through
  covars_mat <- model.matrix(as.formula(mod), covars_frame)[,-1]  # Converts factors to numeric
  design_mat <- cbind(covars_mat, t(meth))
  eventVec <- factor(meta$event)
  mrs.fit.cv <- cv.glmnet(design_mat[trainset,], eventVec[trainset], family="binomial", alpha=alpha,
                          penalty.factor=ifelse(grepl("^cg", colnames(design_mat)), 1, 0), 
                          nfolds=nfolds, lambda.min.ratio=lmr, parallel=parallel)
  coefs <- coef(mrs.fit.cv, s=lambda.type)
  cpgCoefs <- coefs[grepl("^cg", rownames(coefs)),]
  list(fit=mrs.fit.cv,  # Output is a list containing the model fit and the CpG model coefficients
       allCoefs=coefs,
       cpgCoefs=cpgCoefs)
}

test_mrs_binary <- function(meta, mrs, testset) {
  meta$mrs <- scale(mrs)
  mrs.test <- glm(event~mrs, data=meta, family="binomial", subset=testset) # Test the MRS using a Cox model
  res <- summary(mrs.test)$coef["mrs",c("Estimate","Pr(>|z|)")]  # Outputs a vector
  tibble(OR_per_SD=exp(res["Estimate"]), p=res["Pr(>|z|)"])
}

trainTestMRS_binary <- function(mod, meta, meth, trainset, testset, alpha=0.5, nfolds=5, lmr=0.1, lt="lambda.min", parallel=F) {
  mrs_fit <- train_mrs_binary(mod, meta, meth, trainset, 
                              alpha=alpha, nfolds=nfolds, lmr=lmr, lambda.type=lt, parallel=parallel)
  mrs_values <- calc_mrs(mrs_fit$cpgCoefs, meth)
  mrs_test <- tryCatch(test_mrs_binary(meta, mrs_values, testset), error=function(e) matrix(NA, 1, 2))
  list(fit=mrs_fit, components=mrs_fit$cpgCoefs, values=mrs_values, res=mrs_test)
}

