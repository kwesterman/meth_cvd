lapply(c("tidyverse","minfi","wateRmelon"), library, character.only=T)
dataDir <- "~/Box Sync/Analysis/data/"  # At the moment, this only contains a few samples for development
targets <- read.metharray.sheet(dataDir, pattern="^sample_sheet4.csv$") %>%
  filter(Basename!="character(0)")
rgSet <- read.metharray.exp(targets=targets, extended=T)
boxplot(log(getGreen(rgSet)))
boxplot(log(getRed(rgSet)))
mSet.filt <- pfilter(mn=rgSet, pnthresh=0.05, perc=5, pthresh=5)  # Samples and probes must have <5% of representatives with detection P>0.05
  # Also outputs a MethylSet using preprocessRaw
mSet.filt.bmiq_matrix <- BMIQ(mSet.filt)
mSet.filt.bmiq <- RatioSet(mSet.filt.bmiq_matrix)

# Cell counts
data("FlowSorted.Blood.450k.compTable")
cellRef <- as.matrix(FlowSorted.Blood.450k.compTable)[,c('CD8T','CD4T','NK','Bcell','Mono','Gran')]
commonCpGs <- intersect(rownames(melon.filt.bmiq), rownames(cellRef))
estCellCounts <- projectMix(Y=mSet.filt.bmiq[commonCpGs,], Xmat=cellRef[commonCpGs,])

# All covariates
targets_withKey <- mutate(targets, Sample_pKey=paste(Slide, Sample_Well, sep="_"))
estCellCounts_df <- rownames_to_column(as.data.frame(estCellCounts), var="Sample_pKey")
covData <- inner_join(targets_withKey, estCellCounts_df, by="Sample_pKey")

# CVD event data
soe2015 <- read.delim("data/pheno_files/phs000007.v29.pht000309.v12.p10.c1.vr_soe_2015_a_0991s.HMB-IRB-MDS.txt", skip=10)
soe2015.clean <- soe2015 %>%
  mutate(numDays=DATE-37*365) %>%  # NOTE: specific values here will change
  filter(numDays>0) %>%  # Don't keep events that occurred prior to methylation collection
  group_by(shareid) %>%  # To trim to earliest relevant event per person (what about people who have already had events?)
  filter(EVENT %in% 1:26) %>%  # Only events of interest (CVD, stroke, etc. -- see data dictionary)
  slice(which.min(numDays)) %>%  # Only keep first relevant event per person
  ungroup()

# Survival analysis
survData <- left_join(targets, soe2015.clean, by=c("shareid"="Sample_Name")) %>%
  mutate(event=!is.na(numDays),  # event == TRUE when a specific time of event exists, otherwise FALSE
         time=na.replace(numDays, 800))  ## NOTE: this 800 needs to be switched for an *actual* number
survData$survObj <- Surv(time=survData$time, event=survData$event, type="right")  # all subjects without an event time get max. # days as time to censor
stopifnot(ncol(mSet.filt.bmiq)==nrow(survData))

cox.fit <- coxph(survObj~Pool_ID, data=survData)

run_cox <- function(i) {
  survData$meth <- mSet.filt.bmiq[i,]
  cox.fit <- coxph(survObj~meth+CD8T+CD4T+NK+Bcell+Mono+Gran)
  summary(cox.fit)$logtest['pvalue']
}

stopifnot(all(colnames(mSet.filt.bmiq)==survData$Sample_pKey))  # Ensure that order of samples for covariate data and methylation data are identical
mclapply(1:nrow(mSet.filt.bmiq), run_cox, mc.cores=4)


##############################
# Survival analysis
# Make up survival data for the moment:
survData <- data.frame(id=sample(targets$Sample_Name, 10, replace=F),  # Random set of samples "get CVD"
                       time=rnorm(10, 400, 50),
                       abc=rnorm(10, 20, 10)) %>%  # Arbitrary times to event
  right_join(targets, by=c("id"="Sample_Name")) %>%  # Join synthetic event data with sample metadata
  mutate(event=!is.na(time),  
         time=na.replace(time, 800))
survData$survObj <- Surv(time=survData$time, event=survData$event, type="right")  # all subjects without an event time get max. # days as time to censor
cox.fit <- coxph(survObj~Pool_ID+abc, data=survData)



