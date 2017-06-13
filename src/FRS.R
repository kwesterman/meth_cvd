# Calculation of Framingham risk score based on D'Agostino 2008 
# doi: https://doi.org/10.1161/CIRCULATIONAHA.107.699579


calc_FRS <- function(pData) {
FRS_coefs <- list(male=c(logAge=2.32888, logTC=1.20904, logHDL=-0.70833, logSBPnontreat=2.76157, 
                         logSBPtreat=2.82263, smoking=0.52973, diabetes=0.69154),
                  female=c(logAge=3.06117, logTC=1.12370, logHDL=-0.93263, logSBPnontreat=1.93303,
                           logSBPtreat=1.99881, smoking=0.65451, diabetes=-0.57367))


FRS_data <- pData %>%
  mutate(logAge=log(age),
         logTC=log(TOT_CHOL),
         logHDL=log(HDL_CHOL),
         logSBP=log(SBP),
         smoking=smoking_now,
         diabetes=T2D_med) %>%
  select(shareid, sex, logAge, logTC, logHDL, HT_med, logSBP, smoking, diabetes)

## negative for diabetes??
attach(FRS_data)
frs <- ifelse(sex==1,
              2.32888*logAge + 1.20904*logTC - 0.70833*logHDL + 2.76157*logSBP*(1-HT_med) + 2.82263*logSBP*(HT_med) + 0.52973*smoking + 0.69154*diabetes,
              3.06117*logAge + 1.12370*logTC - 0.93263*logHDL + 1.93303*logSBP*(1-HT_med) + 1.99881*logSBP*(HT_med) + 0.65451*smoking - 0.57367*diabetes)
frs
}

# load("../int/phenoData.RData")
# calc_FRS(phenoData)
# 
# load("../int/eventData.RData")
# load("../int/sampleData.RData")
# survData <- sampleData %>%
#   left_join(phenoData, by="shareid") %>%
#   left_join(eventData, by="shareid") %>%
#   dplyr::mutate(event=!is.na(timeToEvent),  # event == TRUE when a specific time of event exists, otherwise FALSE
#                 time=na.replace(timeToEvent, 1000),  ## NOTE: this 1000 needs to be switched for an *actual* number
#                 frs=calc_FRS(.))  
# survData$survObj <- Surv(time=survData$time, event=survData$event, type="right") 
# coxph(survObj~frs, data=survData)


