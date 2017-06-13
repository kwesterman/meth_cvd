library(tidyverse)
library(survival)
library(gtools)

resMat_files <- grep("resMats_[0-9]+.*", list.files("../int/"), value=T)
resMat_list <- lapply(resMat_files, function(f) {
  load(paste0("../int/", f))
  resMats[[length(resMats)]]
})
resMat <- do.call(rbind, resMat_list)
clean_resMat <- function(mat) {
  data.frame(mat, stringsAsFactors=F) %>%
    dplyr::rename(beta=coef, p=Pr...z..) %>%
    mutate_at(vars(one_of("beta","z","p")), as.numeric) %>%
    na.omit()
}
resMat <- clean_resMat(resMat)
resMat$fdr <- p.adjust(resMat$p, "BH")
keep_CpGs <- resMat$CpG[resMat$fdr < 0.01]
print(head(keep_CpGs))

load("../int/Mvals.RData")
print("M-values loaded.")
Mvals_keep <- data.frame(t(Mvals[keep_CpGs,]))
print(dim(Mvals_keep))

load("../int/phenoData.RData")
load("../int/sampleData.RData")
load("../int/eventData.RData")

survData <- sampleData %>%
  left_join(phenoData, by="shareid") %>%
  left_join(eventData, by="shareid") %>%
  dplyr::mutate(event=!is.na(timeToEvent),  # event == TRUE when a specific time of event exists, otherwise FALSE
                time=na.replace(timeToEvent, 1000)) %>%  ## NOTE: this 800 needs to be switched for an *actual* number
  dplyr::slice(match(colnames(Mvals), sampleKey))  # Re-order so rows match methylation columns
survData$survObj <- Surv(time=survData$time, event=survData$event, type="right")  # all subjects without an event time get max. # days as time to censor
print("Survival data prepared.")

smaller_model <- paste0("survObj~sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                paste0("PC",1:10,"_cp",collapse="+"),"+",paste0(keep_CpGs,collapse="+"))

modelData <- cbind(survData, Mvals_keep)
small.fit <- coxph(as.formula(smaller_model), data=modelData)
coef_mat <- summary(small.fit)$coef
CpG_coefs <- coef_mat[keep_CpGs,'coef']
methRiskScore <- apply(Mvals_keep, 1, function(r) sum(r*CpG_coefs))

source("FRS.R")
modelData <- cbind(survData, 
                   frs=calc_FRS(survData[,-which(colnames(survData)=="survObj")]),
                   mrs=methRiskScore)
FRSonly_model <- paste0("survObj~frs+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                        paste0("PC",1:10,"_cp",collapse="+"))
FRSonly.fit <- coxph(as.formula(FRSonly_model), data=modelData)

MRSonly_model <- paste0("survObj~mrs+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                        paste0("PC",1:10,"_cp",collapse="+"))
MRSonly.fit <- coxph(as.formula(MRSonly_model), data=modelData)

combined_model <- MRSonly_model <- paste0("survObj~mrs+frs+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                                          paste0("PC",1:10,"_cp",collapse="+"))
combined.fit <- coxph(as.formula(combined_model), data=modelData)

save("FRSonly.fit", "MRSonly.fit", "combined.fit", file="../int/tmp_coxfitresults.RData")

