suppressMessages(silent <- lapply(c("tidyverse","survival","glmnet","minfi","caret","doParallel","itertools",
                                    "randomForest","pROC","kableExtra"), library, character.only=T))
lazyLoad("../cache/incident_ewas/clean-data_fdb7874b9175a1162bb9b84a199905e6")
betas_whi <- readRDS("../int/betas.qc.norm.filt_whi.rds")
betas_fhs <- readRDS("../int/betas.qc.norm.filt_fhs.rds")
betas_whi <- betas_whi[,match(nonMethData$sampleKey[nonMethData$study=="whi"], colnames(betas_whi))]
betas_fhs <- betas_fhs[,match(nonMethData$sampleKey[nonMethData$study=="fhs"], colnames(betas_fhs))]

nmd_fhs <- filter(nonMethData, study=="fhs")
nmd_whi <- filter(nonMethData, study=="whi")
betaVariances_fhs <- apply(betas_fhs, 1, var)

meth_rf_proxy <- function(pheno, betas, nmd, betaVariances, cpgSubset=NULL) {
  rmSamples <- is.na(nmd[[pheno]])
  phenVec <- nmd[[pheno]][!rmSamples]
  cpgSet <- {
    if (!is.null(cpgSubset)) rownames(betas) %in% cpgSubset 
    else betaVariances_fhs>quantile(betaVariances_fhs,0.5)
  }
  relevantBetas <- betas[cpgSet,!rmSamples]
  glm.fit <- cv.glmnet(t(relevantBetas), phenVec, nfolds=3, parallel=T)
  # phen_pred <- predict(glm.fit, newx=t(relevantBetas), s="lambda.min")
  # cvd_surv <- Surv(time=nmd$time, event=nmd$event)[!rmSamples]
  # print(paste0("Predict CVD events from ", pheno, ":"))
  # print(coxph(cvd_surv ~ phenVec))
  # print(paste0("Predict CVD events from methylation-based approximation:"))
  # print(coxph(cvd_surv ~ phen_pred))
  coefs <- coef(glm.fit, s="lambda.min")
  setNames(as.vector(coefs), rownames(coefs))
}

estimate_from_coefs <- function(cpgCoefs, betas) {
  sharedCpGs <- intersect(names(cpgCoefs), rownames(betas))
  t(betas[sharedCpGs,]) %*% cpgCoefs[sharedCpGs]
}

p_from_z <- function(coxFit, predictor) 2*pnorm(-abs(summary(coxFit)$coef[predictor,"z"]))

phenos <- c("tg","ldl","hdl","smk_py","bmi","sbp","age")

# res <- lapply(phenos, function(pheno) {
#   print(pheno)
#   
#   phenoCoefs <- meth_rf_proxy(pheno, betas_fhs, filter(nonMethData, study=="fhs"), betaVariances_fhs)
#   
#   estimate_fhs <- estimate_from_coefs(phenoCoefs, betas_fhs)
#   estimate_whi <- estimate_from_coefs(phenoCoefs, betas_whi)
#   
#   rsq_test <- cor.test(estimate_whi, nmd_whi[[pheno]])$estimate^2
#   
#   cvd_surv_fhs <- Surv(time=nmd_fhs$time, event=nmd_fhs$event)
#   cvd_surv_whi <- Surv(time=nmd_whi$time, event=nmd_whi$event)
#   
#   rf_train <- coxph(as.formula(paste0("cvd_surv_fhs~",pheno)), data=nmd_fhs)
#   rf_test <- coxph(as.formula(paste0("cvd_surv_whi~",pheno)), data=nmd_whi)
#   
#   predRF_train <- coxph(cvd_surv_fhs~estimate_fhs)
#   predRF_test <- coxph(cvd_surv_whi~estimate_whi)
#   
#   c(rsq_test=rsq_test,
#     rf_train=p_from_z(rf_train), 
#     rfPred_train=p_from_z(predRF_train),
#     rf_test=p_from_z(rf_test),
#     rfPred_test=p_from_z(predRF_test),
#     nTest_withRF=sum(!is.na(nmd_whi[[pheno]])))
# })
# 
# fullNames <- c(rsq_test.cor="Rsq in test\n(var. expl.)", rf_train="RF -> CVD in train",
#                rfPred_train="Pred. RF -> CVD in train", rf_test="RF -> CVD in test",
#                rfPred_test="Pred. RF -> CVD in test", nTest_withRF="Sample size w/\nRF in test")  
# 
# resDF <- cbind(pheno=phenos, data.frame(do.call(rbind, res))) %>%
#   gather(field, value, rsq_test.cor:nTest_withRF) %>%
#   mutate(text=case_when(field=="rsq_test.cor" ~ as.character(round(value,2)),
#                         field=="nTest_withRF" ~ as.character(round(value)),
#                         TRUE ~ ""),
#          negLogP=ifelse(field %in% c("rsq_test.cor","nTest_withRF"), 0, -log10(as.numeric(value))),
#          negLogP=pmin(20,negLogP),
#          field=factor(fullNames[field], levels=unique(fullNames[field])))
#   
# 
# print(resDF)
# 
# jpeg("initialPlot.jpg")
# ggplot(resDF, aes(x=field, y=forcats::fct_rev(pheno))) +
#   scale_x_discrete(position="top") +
#   geom_tile(aes(fill=negLogP), color="#6475A4") +
#   scale_fill_gradient2(low='#0000FF',mid='#FFFFFF',high='#FF0000') +
#   geom_text(aes(label=text)) +
#   theme(axis.title=element_blank(), axis.text.x=element_text(angle=30))
# dev.off()

whi_estimates <- lapply(setNames(phenos, phenos), function(pheno) {
  print(pheno)
  phenoCoefs <- meth_rf_proxy(pheno, betas_fhs, filter(nonMethData, study=="fhs"), betaVariances_fhs)
  estimate_whi <- estimate_from_coefs(phenoCoefs, betas_whi)
  estimate_whi
})

saveRDS(whi_estimates, "../int/rf_estimates_whi.rds")


# ---------------------------
# 
# smk_cpgs <- scan(file="few_smoking_cpgs.txt", what=character())
# 
# betas_fhs <- readRDS("../int/betas.qc.norm.filt_fhs.rds")
# lazyLoad("../cache/incident_ewas/clean-data_fdb7874b9175a1162bb9b84a199905e6")
# nmd_fhs <- filter(nonMethData, study=="fhs")
# betas_fhs <- betas_fhs[,match(nonMethData$sampleKey[nonMethData$study=="fhs"], colnames(betas_fhs))]
# 
# betas_smk_fhs <- betas_fhs[rownames(betas_fhs) %in% smk_cpgs,]
# betas_smk_fhs_small <- betas_smk_fhs[1:10,]
# cvd_surv <- Surv(time=nmd_fhs$time, event=nmd_fhs$event)
# glm.fit <- cv.glmnet(t(betas_smk_fhs_small), nmd_fhs$tg)
# py_pred <- predict(glm.fit, newx=t(betas_smk_fhs_small), s="lambda.min")
# summary(coxph(cvd_surv~nmd_fhs$smk_py))
# summary(coxph(cvd_surv~py_pred))


