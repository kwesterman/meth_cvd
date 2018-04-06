### Prep

suppressMessages(silent <- lapply(c("tidyverse","survival","glmnet","minfi","caret","doParallel","itertools",
                                    "randomForest","pROC","kableExtra"), library, character.only=T))

betas_whi <- readRDS("../int/betas.qc.norm.filt_whi.rds")
betas_fhs <- readRDS("../int/betas.qc.norm.filt_fhs.rds")

# Load metadata
metaData <- readRDS("../int/metaData.rds")

# Load estimated cell counts
estCellCounts_whi <- readRDS("../int/estCellCounts_whi.rds")
estCellCounts_fhs <- readRDS("../int/estCellCounts_fhs.rds")

# Load CPACOR principal component adjustment factors
load("../int/CPACOR_whi.RData")
CP_PCs_whi <- CP_PCs
load("../int/CPACOR_fhs.RData")
CP_PCs_fhs <- CP_PCs

# Load PCA results
load("../int/PCA.fit_whi.RData")
PCs_whi <- PCs
load("../int/PCA.fit_fhs.RData")
PCs_fhs <- PCs

# Load Illumina 450k annotation
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno450k <- data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19), stringsAsFactors=F)

estCellCounts <- bind_rows(estCellCounts_whi, estCellCounts_fhs)
CP_PCs <- bind_rows(CP_PCs_whi, CP_PCs_fhs)
PCs <- bind_rows(PCs_whi, PCs_fhs)

### Create master covariate/event data frame
nonMethData <- Reduce(function(x,y) inner_join(x, y, by="sampleKey"), 
                      list(metaData, estCellCounts, CP_PCs, PCs)) %>%
  distinct(subjID, .keep_all=T)  # Hacky...better way of deciding which biological replicates to keep?
nonMethData <- replace_na(nonMethData, 
                          list(bmi=median(nonMethData$bmi, na.rm=T),
                               smk_now=0, smk_py=0, ht_med=0, lipid_med=0, dm_med=0))
nmd_whi <- filter(nonMethData, study=="whi")
nmd_fhs <- filter(nonMethData, study=="fhs")

betas_whi <- betas_whi[,match(nonMethData$sampleKey[nonMethData$study=="whi"], colnames(betas_whi))]
betas_fhs <- betas_fhs[,match(nonMethData$sampleKey[nonMethData$study=="fhs"], colnames(betas_fhs))]
stopifnot(all(colnames(betas_whi)==nonMethData$sampleKey[nonMethData$study=="whi"]),
          all(colnames(betas_fhs)==nonMethData$sampleKey[nonMethData$study=="fhs"]))

myTry <- function(expr, CpG) {
  # Captures model failures and returns successful result or vector of NAs
  tryCatch(expr,
           error = function(e) {print(e); return(c(CpG, rep(NA, 3)))},
           warning = function(w) {print(w); return(c(CpG, rep(NA, 3)))})
}

run_coxModel <- function(probeData, covarData, model_spec, subset) {
  # Given a row of the M-value matrix (corresponding to a CpG site), bind that methylation
  # data to the covariate data and run Cox proportional hazards regression
  CpG <- rownames(probeData)
  probeData <- as.vector(probeData)
  outlierTF <- probeData<quantile(probeData,0.25)-5*IQR(probeData) | 
    probeData>quantile(probeData,0.75)+5*IQR(probeData)
  modelData <- cbind(covarData, meth=as.numeric(probeData),
                     survObj=Surv(time=covarData$time, event=covarData$event))
  myTry({
    cox.fit <- coxph(as.formula(model_spec), data=modelData, subset=!outlierTF)
    c(CpG=CpG, summary(cox.fit)$coef['meth',c('coef','z','Pr(>|z|)')])
  }, CpG)
}


### Actual approach

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno450k <- data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19), stringsAsFactors=F)
geneGroupListDF <- anno450k %>% 
  dplyr::rename(cpg=Name, gene=UCSC_RefGene_Name, group=UCSC_RefGene_Group) %>%
  dplyr::select(cpg, gene, group) %>%
  separate_rows(gene, group, sep=";") %>% 
  distinct(cpg, gene, group) %>%
  filter(cpg %in% rownames(betas_whi),  # NOTE THAT I DID THIS -- MAYBE FIX LATER
         cpg %in% rownames(betas_fhs)) %>%
  group_by(gene, group) %>% 
  summarise(n=n(), cpgs=list(cpg)) %>%
  filter(gene!="") %>%
  filter(n>1) %>%
  mutate(label=paste(gene, group, sep="_"))

geneGroup_to_eigenCpG <- function(cpgSet) {
  whiIntersect <- intersect(cpgSet, rownames(betas_whi)) 
  # if (length(whiIntersect)==0) return (NA)
  # else if (length(whiIntersect)==1) return (scale(betas_whi[whiIntersect,]))
  prcomp(t(betas_whi[whiIntersect,]), scale.=T)$rotation[,"PC1"]
}

calc_eigencpg <- function(betas, weights) {
  includeCpgs <- intersect(rownames(betas), names(weights))
  as.vector(scale(t(betas[includeCpgs,])) %*% weights[includeCpgs])
}

pc1weights <- purrr::map(geneGroupListDF$cpgs, geneGroup_to_eigenCpG)
names(pc1weights) <- geneGroupListDF$label
saveRDS(pc1weights, file="../int/regionApproach_pc1weights.rds")

betas_whi_split <- lapply(1:nrow(geneGroupListDF), 
                          function(idx) list(locus=geneGroupListDF$label[idx],
                                             betas=betas_whi[geneGroupListDF$cpgs[[idx]],],
                                             weights=pc1weights[[idx]])) 
betas_fhs_split <- lapply(1:nrow(geneGroupListDF), 
                          function(idx) list(locus=geneGroupListDF$label[idx],
                                             betas=betas_fhs[geneGroupListDF$cpgs[[idx]],],
                                             weights=pc1weights[[idx]])) 

cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)

ewasModel_whi <- paste0("survObj~meth+age+race+bmi+smk_now+smk_py+CD4T+CD8T+Bcell+NK+Mono+Gran+", 
                        paste0("cpPC",1:5,collapse="+"))
eigenProbeResList_whi <- foreach(loc=iter(betas_whi_split), .packages="survival") %dopar% {
  eigencpg <- calc_eigencpg(loc$betas, loc$weights)
  run_coxModel(eigencpg, nmd_whi, ewasModel_whi)
}
names(eigenProbeResList_whi) <- names(pc1weights)
eigenProbeRes_whi <- do.call(rbind, eigenProbeResList_whi) %>%
  data.frame(stringsAsFactors=F, check.names=F) %>%
  rownames_to_column(var="cpgSet") %>%
  dplyr::rename(p=`Pr(>|z|)`) %>%
  mutate_at(c("coef","z","p"), as.numeric) %>%
  mutate(fdr=p.adjust(p, method="BH")) %>%
  arrange(p)
saveRDS(eigenProbeRes_whi, file="../int/regionApproach_eigencpgRes_whi.rds")

ewasModel_fhs <- paste0("survObj~meth+sex+age+bmi+smk_now+smk_py+CD4T+CD8T+Bcell+NK+Mono+Gran+center+",
                        paste0("cpPC",1:7,collapse="+"))
eigenProbeResList_fhs <- foreach(loc=iter(betas_fhs_split), .packages="survival") %dopar% {
  eigencpg <- calc_eigencpg(loc$betas, loc$weights)
  run_coxModel(eigencpg, nmd_fhs, ewasModel_fhs)
}
names(eigenProbeResList_fhs) <- names(pc1weights)
eigenProbeRes_fhs <- do.call(rbind, eigenProbeResList_fhs) %>%
  data.frame(stringsAsFactors=F, check.names=F) %>%
  rownames_to_column(var="cpgSet") %>%
  dplyr::rename(p=`Pr(>|z|)`) %>%
  mutate_at(c("coef","z","p"), as.numeric) %>%
  mutate(fdr=p.adjust(p, method="BH")) %>%
  arrange(p)
saveRDS(eigenProbeResList_fhs, file="../int/regionApproach_eigencpgRes_fhs.rds")

# stopCluster(cl)

make_qqplot <- function(pVec, plotTitle="Title") {
  pVec <- pVec[!is.na(pVec)]
  qqplot(-log10(1:length(pVec)/length(pVec)), -log10(pVec), pch=".", main=plotTitle, xlab="Expected (-logP)", ylab="Observed (-logP)")
  abline(0,1,col="red")
}
gControl <- function(pVals) {
  # See van Iterson 2017 methods and/or Lehne 2015 code for details on genomic control for EWAS
  # Below is modeled after Lehne 2015
  lambda <- median(qchisq(pVals, df=1, lower.tail=F), na.rm=T)/qchisq(0.5, df=1)
  round(lambda, 2)
}

make_qqplot(eigenProbeRes_whi$p, plotTitle=paste0("lambda = ", gControl(eigenProbeRes_whi$p)))


make_qqplot(eigenProbeRes_fhs$p, plotTitle=paste0("lambda = ", gControl(eigenProbeRes_fhs$p)))


bothNominalSites <- na.omit(intersect(eigenProbeRes_whi$cpgSet[eigenProbeRes_whi$p<0.01], eigenProbeRes_fhs$cpgSet[eigenProbeRes_fhs$p<0.01]))
bothSuggestiveSites <- na.omit(intersect(eigenProbeRes_whi$cpgSet[eigenProbeRes_whi$p<0.001], eigenProbeRes_fhs$cpgSet[eigenProbeRes_fhs$p<0.001]))

whiCoefs <- setNames(eigenProbeRes_whi$coef, eigenProbeRes_whi$cpgSet)
fhsCoefs <- setNames(eigenProbeRes_fhs$coef, eigenProbeRes_fhs$cpgSet)


commonRes <- inner_join(data.frame(CpG=bothSuggestiveSites, 
                                   coef_direction=sign(whiCoefs[bothSuggestiveSites])),
                        select(anno450k, Name, UCSC_RefGene_Name, UCSC_RefGene_Group), 
                        by=c("CpG"="Name"))
