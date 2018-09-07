library(tidyverse)
library(WGCNA)
library(survival)

betas_whi <- readRDS("../int/betas.qc.norm.filt_whi.rds")
beta_variances <- apply(betas_whi, 1, var)

enableWGCNAThreads(8)  # use enableWGCNAThreads() instead of allow...()


# Subset of individuals for threshold choice and pre-clustering
set.seed(123)
sample_subset <- sample(1:ncol(betas_whi), 100)  
meth_dat_subset <- t(betas_whi[beta_variances > quantile(beta_variances, 0.5),
                               sample_subset])

# Analysis for power choice
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
system.time(sft <- pickSoftThreshold(meth_dat_subset, 
                                     powerVector=powers, verbose=5))
## soft threshold choice on 100 samples for top 50% variable genes took 3.5 hrs
saveRDS(sft, file="../int/wgcna/whiSoftThreshold.rds")

# Code from tutorial to assess proper soft threshold choice
par(mfrow = c(1, 2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab="Soft Threshold (power)", 
     ylab="Scale Free Topology Model Fit,signed R^2", 
     type="n", main=paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels=powers, cex=cex1, col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90, col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main=paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     labels=powers, cex=cex1, col="red")
## --> based on visual analysis, preliminary choice of beta=6

# Pre-clustering using projective K-means
system.time(pkm <- projectiveKMeans(t(betas_whi[, sample_subset])))
## pkm on 100 samples for top 50% variable genes took a couple of hours
saveRDS(pkm, file="../int/wgcna/whiPreclustering.rds")


# Run WGCNA on full matrix of WHI beta values
system.time(bwnet <- blockwiseModules(t(betas_whi),
                                      blocks=pkm$clusters,
                                      power=6,
                                      networkType="unsigned",
                                      saveTOMs=T,
                                      saveTOMFileBase="../int/wgcna/whiTOM",
                                      verbose=2))
## blockwise network construction took 4.5 hrs when given pkm clustering blocks
saveRDS(bwnet, file="../int/wgcna/whiNet.rds")


# Calculate module eigengenes
system.time(MEs <- moduleEigengenes(t(betas_whi), bwnet$colors, excludeGrey=T))
save(MEs, file="../int/wgcna/whiMEs.rds")
## eigengene calculation took 0.5 hrs


# Module preservation in FHS
multi_expr <- list(whi=list(data=t(betas_whi)),
                  fhs=list(data=t(betas_fhs[(
                    rownames(betas_fhs) %in% wgcnaInputCpgs), ])))
multi_color <- list(whi=whiNet$colors)
system.time(mp <- modulePreservation(multi_expr, 
                                     multi_color,
                                     networkType="unsigned",
                                     referenceNetworks=1,
                                     nPermutations=5,
                                     randomSeed=1,
                                     quickCor=1,
                                     verbose=3))
## modulePreservation takes ~40 mins. for full set of markers
saveRDS(mp, file="../int/wgcna/fhsMP.rds")


# eg_model_res_list <- lapply(MEs$eigengenes, function(eg) {
#   cox.fit <- coxph(Surv(time, event) ~ eg + dnaPull, data=nmd_whi)
#   summary(cox.fit)$coef["eg", c("coef", "z", "Pr(>|z|)")]
# })
# eg_model_res_df <- do.call(rbind, eg_model_res_list) %>%
#   as.data.frame() %>%
#   rownames_to_column(var="module") %>%
#   dplyr::rename(p=`Pr(>|z|)`) %>%
#   mutate(module=gsub("^ME", "", module)) %>%
#   arrange(p)
# sigModules <- eg_model_res_df$module[eg_model_res_df$p<0.05]

library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
enrichments <- lapply(setNames(sigModules,sigModules), function(m) {
  cpgSet <- rownames(betas_whi)[bwnet$colors==m]
  gsaTbl <- missMethyl::gometh(cpgSet, collection="GO")
  gsaTopResTbl <- arrange(filter(gsaTbl, FDR<0.25), FDR)
})

# egModelRes.adj.list <- lapply(MEs$eigengenes, function(eg) {
#   cox.fit <- coxph(Surv(time,event)~eg+dnaPull+age+bmi+race, data=nmd_whi)
#   summary(cox.fit)$coef["eg",c("coef","z","Pr(>|z|)")]
# })
# egModelRes.adj.df <- do.call(rbind, egModelRes.adj.list) %>%
#   as.data.frame() %>%
#   rownames_to_column(var="module") %>%
#   arrange(`Pr(>|z|)`)

egModelRes.adj.list <- lapply(MEs$eigengenes, function(eg) {
  cox.fit <- coxph(Surv(time,event)~eg+dnaPull+age+bmi+race+CD4T+CD8T+Bcell+Mono+NK+Gran, data=nmd_whi)
  summary(cox.fit)$coef["eg",c("coef","z","Pr(>|z|)")]
})
egModelRes.adj.df <- do.call(rbind, egModelRes.adj.list) %>%
  as.data.frame() %>%
  rownames_to_column(var="module") %>%
  dplyr::rename(p=`Pr(>|z|)`) %>%
  mutate(module=gsub("^ME", "", module)) %>%
  arrange(p)
sigModules.adj <- egModelRes.adj.df$module[egModelRes.adj.df$p<0.05]
enrichments.adj <- lapply(setNames(sigModules.adj,sigModules.adj), function(m) {
  cpgSet <- rownames(betas_whi)[bwnet$colors==m]
  gsaTbl <- missMethyl::gometh(cpgSet, collection="GO")
  gsaTopResTbl <- arrange(filter(gsaTbl, FDR<0.25), FDR)
})

# fhsEigengenes <- lapply(setNames(sigModules,sigModules), function(col) {
allSigModules <- unique(c(sigModules, sigModules.adj))
wgcnaInputCpgs <- rownames(betas_whi)
whiNet <- bwnet
fhsEigengenes <- lapply(setNames(allSigModules, allSigModules), function(col) {
  print(col)
  pc1rots <- prcomp(t(betas_whi[wgcnaInputCpgs[whiNet$colors==col],]), scale.=T)$rot[,"PC1"]
  commonCpgs <- intersect(rownames(betas_fhs), names(pc1rots))
  as.vector(t(betas_fhs[commonCpgs,]) %*% pc1rots[commonCpgs])
})
save("fhsEigengenes", file="../int/wgcna/fhsEGs.RData")
# coxph(Surv(time,event)~fhsEigengenes$darkturquoise, data=nmd_fhs)

fhseg_model_res_list <- lapply(fhsEigengenes, function(eg) {
  cox.fit <- coxph(Surv(time,event)~eg+center, data=nmd_fhs)
  summary(cox.fit)$coef["eg",c("coef","z","Pr(>|z|)")]
})
fhseg_model_res_df <- do.call(rbind, fhseg_model_res_list) %>%
  as.data.frame() %>%
  rownames_to_column(var="module") %>%
  dplyr::rename(p=`Pr(>|z|)`) %>%
  mutate(module=gsub("^ME", "", module)) %>%
  arrange(p)

fhsEgModelRes.adj.list <- lapply(fhsEigengenes, function(eg) {
  cox.fit <- coxph(Surv(time,event)~eg+center+cpPC1+cpPC2+cpPC3+cpPC4+cpPC5+cpPC6+cpPC7+age+bmi+sex+smk_now+CD4T+CD8T+Bcell+Mono+NK+Gran, data=nmd_fhs)
  summary(cox.fit)$coef["eg",c("coef","z","Pr(>|z|)")]
})
fhsEgModelRes.adj.df <- do.call(rbind, fhsEgModelRes.adj.list) %>%
  as.data.frame() %>%
  rownames_to_column(var="module") %>%
  dplyr::rename(p=`Pr(>|z|)`) %>%
  mutate(module=gsub("^ME", "", module)) %>%
  arrange(p)


wbcpca.fit <- prcomp(select(nmd_fhs, CD8T, CD4T, NK, Bcell, Mono, Gran), scale.=T)
nmd_fhs$WBC_PC1 <- wbcpca.fit$x[,"PC1"]
nmd_fhs$WBC_PC2 <- wbcpca.fit$x[,"PC2"]

corPrepDF <- cbind(select(nmd_fhs, age, bmi, ldl, hdl, tg, sbp, WBC_PC1, WBC_PC2), 
                   data.frame(fhsEigengenes))
corDF <- cor(corPrepDF, use="pairwise.complete.obs") %>%
  data.frame() %>%
  rownames_to_column(var="module") %>%
  filter(module %in% names(fhsEigengenes)) %>%
  gather(key=riskFactor, value=correlation, -module) %>%
  filter(!(riskFactor %in% names(fhsEigengenes)))
rfCorPlt <- ggplot(corDF, aes(x=riskFactor, y=module, fill=correlation)) +
  geom_tile() +
  scale_fill_gradient2(low="blue",mid="white",high="red")

pDF_whi <- select(corDF, module) %>%
  inner_join(select(eg_model_res_df, module, p), by="module") %>%
  inner_join(select(egModelRes.adj.df, module, p), by="module") %>%
  dplyr::rename(p.raw=p.x, p.adj=p.y) %>%
  gather(key=model, value=p, p.raw:p.adj) %>%
  mutate(negLogP=-log10(p))
whiPPlt <- ggplot(pDF_whi, aes(x=model, y=module, fill=negLogP)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="black") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank())



fhsPres.df <- rownames_to_column(mp$preservation$Z$ref.whi$inColumnsAlsoPresentIn.fhs, var="module") %>%
  dplyr::rename(Zpres=Zsummary.pres)

pDF_fhs <- select(corDF, module) %>%
  # inner_join(select(fhsPres.df, module, moduleSize, Zpres), by="module") %>%
  inner_join(select(fhseg_model_res_df, module, p), by="module") %>%
  inner_join(select(fhsEgModelRes.adj.df, module, p), by="module") %>%
  dplyr::rename(p.raw=p.x, p.adj=p.y) %>%
  gather(key=model, value=p, -module, -moduleSize) %>%
  mutate(negLogP=-log10(p))
fhsPPlt <- ggplot(pDF_fhs, aes(x=model, y=module, fill=negLogP)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="black") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank())

plot_grid(rfCorPlt, whiPPlt, fhsPPlt, nrow=1, rel_widths=c(3,1,1))

#                   fhs=t(betas_fhs[rownames(betas_fhs) %in% colnames(methDatVar),]))
# multi_color <- list(whi=whiNet$colors)
# system.time(mp <- modulePreservation(multi_expr, multi_color,
#                                      networkType="unsigned",
#                                      referenceNetworks=1,
#                                      nPermutations=50,
#                                      randomSeed=1,
#                                      quickCor=1,
#                                      verbose=3))
# save("mp", file="../whiToFhsMP.RData")






################################################################
## ALTER BELOW -- CALCULATE CELL COUNT-ADJUSTED CORRELATIONS WITH INCIDENT CVD ##
# whiEGWAS.ccAdj.list <- lapply(MEs$eigengenes, function(eg) {
#   cox.fit <- coxph(Surv(time,event)~eg+dnaPull+CD4T+CD8T+Bcell+Mono+NK+Gran, data=nmd_whi)
#   summary(cox.fit)$coef["eg",c("coef","z","Pr(>|z|)")]
# })
# whiEGWAS.ccAdj.df <- do.call(rbind, whiEGWAS.ccAdj.list) %>%
#   as.data.frame() %>%
#   rownames_to_column(var="module") %>%
#   dplyr::rename(p=`Pr(>|z|)`) %>%
#   mutate(module=gsub("^ME", "", module)) %>%
#   arrange(p)
# sigModules.ccAdj <- whiEGWAS.ccAdj.df$module[whiEGWAS.ccAdj.df$p<0.05]
# enrichments.ccAdj <- lapply(setNames(sigModules.ccAdj,sigModules.ccAdj), function(m) {
#   cpgSet <- rownames(betas_whi)[bwnet$colors==m]
#   gsaTbl <- missMethyl::gometh(cpgSet, collection="GO")
#   gsaTopResTbl <- arrange(filter(gsaTbl, FDR<0.25), FDR)
# })


# fhsEigengenes.ccAdj <- lapply(setNames(sigModules.ccAdj, sigModules.ccAdj), function(m) {
#   pca.fit <- prcomp(t(betas_whi[wgcnaInputCpgs[whiNet$colors==m],]), scale.=T)
#   pc1rots <- pca.fit$rot[,"PC1"]
#   if (sign(cor(pca.fit$x[,"PC1"], whiMEs$eigengenes[[m]]))==1) pc1rots <- -pc1rots
#   commonCpgs <- intersect(rownames(betas_fhs), names(pc1rots))
#   as.vector(t(betas_fhs[commonCpgs,]) %*% pc1rots[commonCpgs])
# })
# save("fhsEigengenes", file="../int/wgcna/fhsEGs.RData")
# # coxph(Surv(time,event)~fhsEigengenes$darkturquoise, data=nmd_fhs)
# 
# fhsEgModelRes.ccAdj.list <- lapply(fhsEigengenes.ccAdj, function(eg) {
#   nmd_fhs$eg <- eg
#   cox.fit <- coxph(Surv(time,event)~eg+age+center+cpPC1+cpPC2+cpPC3+cpPC4+cpPC5+cpPC6+cpPC7+CD4T+CD8T+Bcell+Mono+NK+Gran, data=nmd_fhs)
#   # cox.fit <- coxph(Surv(time,event)~eg+center+CD4T+CD8T+Bcell+Mono+NK+Gran, data=nmd_fhs)
#   summary(cox.fit)$coef["eg",c("coef","z","Pr(>|z|)")]
# })
# fhsEgModelRes.ccAdj.df <- do.call(rbind, fhsEgModelRes.ccAdj.list) %>%
#   as.data.frame() %>%
#   rownames_to_column(var="module") %>%
#   dplyr::rename(p=`Pr(>|z|)`) %>%
#   mutate(module=gsub("^ME", "", module)) %>%
#   arrange(p)
# 
# fhsEgModelRes.ccAdj.list <- lapply(fhsEigengenes.ccAdj, function(eg) {
#   cox.fit <- coxph(Surv(time,event)~eg+center+cpPC1+cpPC2+cpPC3+cpPC4+cpPC5+cpPC6+cpPC7+age+bmi+sex+smk_now+CD4T+CD8T+Bcell+Mono+NK+Gran, data=nmd_fhs)
#   summary(cox.fit)$coef["eg",c("coef","z","Pr(>|z|)")]
# })
# fhsEgModelRes.adj.df <- do.call(rbind, fhsEgModelRes.adj.list) %>%
#   as.data.frame() %>%
#   rownames_to_column(var="module") %>%
#   dplyr::rename(p=`Pr(>|z|)`) %>%
#   mutate(module=gsub("^ME", "", module)) %>%
#   arrange(p)

### FORMALLY EXPLORE THE SEX DIFFERENCES IN REPLICATION OF THESE MODULES ###



sigByRF <- expand.grid(module=sigModules.ccAdj, rf=c("age","bmi","ldl","hdl","tg","sbp","diabetes","smk_now"))
sigByRF$cor <- apply(sigByRF, 1, function(row) {
  cor.test(MEs$eigengenes[[paste0("ME",row["module"])]], as.numeric(nmd_whi[[row["rf"]]]))$estimate
})
sigByRF$moduleLabel=factor(sigByRF$module, levels=sigModules.ccAdj)

rfPlt <- ggplot(sigByRF, aes(x=rf, y=module, fill=cor)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low='#0000FF',mid='#FFFFFF',high='#FF0000') +
  # ggtitle("Correlations of DMPs with traditional risk factors") +
  theme(legend.title=element_blank(), axis.title=element_blank())

library(qgraph)
corInput.df <- cbind(data.frame(MEs$eigengenes[paste0("ME",sigModules.ccAdj)]), 
                     select(nmd_whi, unique(as.character(sigByRF$rf))))
corMat <- cor(corInput.df, use="pairwise.complete.obs")
qgraph(corMat, graph="pcor", threshold="bonferroni", sampleSize=2000, layout="spring")
qgraph(corMat, graph="cor", threshold="bonferroni", sampleSize=2000, layout="spring")

##############################################

## Basic EWAS of ccAdj-significant modules
a <- whiRes[whiRes$CpG %in% rownames(betas_whi)[bwnet$colors %in% sigModules.ccAdj],]
b <- fhsRes[fhsRes$CpG %in% a$CpG,]


#################3

## PhenoAge stuff
phenoAgeWeights <- read_csv("../data/literature/phenoAge.xls")
paCpgsInWhi <- intersect(phenoAgeWeights$CpG, rownames(betas_whi))
whiPA <- as.vector(t(betas_whi[paCpgsInWhi,]) %*% phenoAgeWeights$Weight[match(paCpgsInWhi,phenoAgeWeights$CpG)])


### Does membership in a module mean better replication across to FHS?
cpgToModule <- data.frame(CpG=rownames(betas_whi), module=bwnet$colors)
compareDF <- whiResFull %>%
  inner_join(cpgToModule, by="CpG") %>%
  mutate(inSig=module %in% sigModules.ccAdj) %>%
  inner_join(select(fhsResFull, CpG, coef, p), by="CpG", suffix=c(".whi",".fhs"))
repDF <- compareDF %>%
  filter(p.whi<0.05) %>%
  mutate(replicated=(sign(coef.fhs)==sign(coef.whi) & p.fhs<0.05)) %>%
  group_by(inSig, replicated) %>%
  summarise(n=n()) %>%
  arrange(desc(inSig), desc(replicated))
repDF$n[1]/repDF$n[2]
repDF$n[3]/repDF$n[4]
chisq.test(matrix(repDF$n, 2, 2, byrow=T))




######################################################

anno450kRanges <- GRanges(anno450k$chr, IRanges(start=anno450k$pos, end=anno450k$pos))

read_epiAnno <- function(broadPeaksFile) {
  # Given a path to a broadPeaks annotation file, read in the set of regions and output a GRanges object
  epiAnno.df <- suppressMessages(read_tsv(broadPeaksFile, col_names=F, progress=F)) %>%
    select(1, 2, 3, ncol(.)) %>%
    setNames(c("chr","start","end","negLogQ")) %>%
    filter(negLogQ>2)
  GRanges(seqnames=epiAnno.df$chr, IRanges(start=epiAnno.df$start, end=epiAnno.df$end))
}

# test_epiAnno <- function(epiAnno.gr, ewasRes) {
#   # Given a set of peak calls/regions of interest for an epigenomic annotation of interest (as GRanges 
#   # object), return the -log(p) testing whether a set of EWAS results are associated with the annotation
#   cpgsInAnno <- findOverlaps(anno450kRanges, epiAnno.gr)
#   annoWithMembership <- mutate(anno450k, member=1:nrow(anno450k) %in% queryHits(cpgsInAnno))
#   annoWithMembershipAndPvals <- inner_join(ewasRes, annoWithMembership, by=c("CpG"="Name"))
#   tryCatch(unname(t.test(negLogP~member, data=annoWithMembershipAndPvals)$statistic),
#            error=function(e) NA)
# }

test_epiAnno <- function(epiAnno.gr, ewasRes) {
  # Given a set of peak calls/regions of interest for an epigenomic annotation of interest (as GRanges
  # object), return the -log(p) testing whether a set of EWAS results are associated with the annotation
  cpgsInAnno <- findOverlaps(anno450kRanges, epiAnno.gr)
  annoWithMembership <- mutate(anno450k, member=1:nrow(anno450k) %in% queryHits(cpgsInAnno))
  annoWithMembershipAndPvals <- inner_join(ewasRes, annoWithMembership, by=c("CpG"="Name"))
  returnVal <- tryCatch({
    p <- phyper(sum(annoWithMembershipAndPvals$member*(annoWithMembershipAndPvals$negLogP>-log10(0.05))),
                sum(annoWithMembershipAndPvals$member),  # Number of sites in peaks
                sum(!annoWithMembershipAndPvals$member),  # Number of sites not in peaks
                sum(annoWithMembershipAndPvals$negLogP>-log10(0.05)))  # Number of sites with p<0.05
    -log10(min(p, 1-p))
  },
  error=function(e) NA)
  if (is.na(returnVal) || is.infinite(returnVal)) NA else returnVal
}

assemble_epiAnno <- function(annoType, cellType, ewasRes) {
  # Given an epigenomic annotation type (e.g. H3K4me1), output a named list of GRanges objects
  annoDir <- paste0("../data/literature/roadmap/", annoType)
  annoFile <- grep(cellType, list.files(annoDir, full.names=T), value=T)
  if (length(annoFile)==0) return(NA) else epiAnno.gr <- read_epiAnno(annoFile)
  tStat <- test_epiAnno(epiAnno.gr, ewasRes)
  data.frame(EpiFeature=annoType, CellTypeCode=cellType, tStat=as.numeric(tStat))
}

epiFeatures <- list.files("../data/literature/roadmap/", pattern="^[DH].*")
cellTypeCodes <- unique(gsub("^.*/|-.*", "", list.files("../data/literature/roadmap/", recursive=T)))
roadmapMeta <- read_csv("../data/literature/roadmap/jul2013.roadmapData.qc\ -\ Consolidated_EpigenomeIDs_QC.csv", col_types=cols())
# relevant_tissues <- grep("fetal|blood|liver|endoth|cardi|heart", epiAnnoSummary$CellType, ignore.case=T, value=T)
featureCellGrid <- expand.grid(EpiFeature=epiFeatures, CellTypeCode=cellTypeCodes) %>%
  mutate(CellType=roadmapMeta$`Standardised epigenome name`[match(CellTypeCode,roadmapMeta$EID)]) %>%
  filter(grepl("fetal|blood|liver|endoth|cardi|heart", CellType, ignore.case=T))
mergedEwasResByModule <- inner_join(whiResFull, fhsResFull, by="CpG", suffix=c(".whi",".fhs")) %>%
  mutate(negLogP=-(log10(p.whi)+log10(p.fhs))/2) %>%
  select(CpG, negLogP) %>%
  inner_join(cpgToModule, by="CpG")

mergedBlue <- filter(mergedEwasResByModule, module=="blue")
mergedBrown <- filter(mergedEwasResByModule, module=="brown")

mergedResModules <- lapply(setNames(sigModules.ccAdj,sigModules.ccAdj), 
                           function(m) filter(mergedEwasResByModule, module==m))
  
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)
tStats <- lapply(mergedResModules, function(mrm) {
  foreach(anno=iter(featureCellGrid, by="row"), .combine=rbind,
          .packages=c("tidyverse","GenomicRanges")) %do%
    assemble_epiAnno(anno$EpiFeature, anno$CellTypeCode, mrm) 
})
# tStats_blue <- foreach(anno=iter(filter(featureCellGrid,EpiFeature=="H3K27ac"), by="row"), .combine=rbind,
#                   .packages=c("tidyverse","GenomicRanges")) %do%
#   assemble_epiAnno(anno$EpiFeature, anno$CellTypeCode, mergedBlue) 
# tStats_brown <- foreach(anno=iter(filter(featureCellGrid,EpiFeature=="H3K27ac"), by="row"), .combine=rbind,
#                        .packages=c("tidyverse","GenomicRanges")) %do%
#   assemble_epiAnno(anno$EpiFeature, anno$CellTypeCode, mergedBrown) 
stopCluster(cl)

epiAnnoSummaries <- lapply(tStats, function(m) {
  m %>%
    na.omit() %>%
    mutate(CellType=roadmapMeta$`Standardised epigenome name`[match(CellTypeCode,roadmapMeta$EID)]) %>%
    group_by(CellType) %>%
    # filter(all(epiFeatures %in% EpiFeature)) %>%
    summarise(CellTypeAbsMean=sum(abs(tStat))) %>%
    arrange(desc(CellTypeAbsMean)) %>%
    mutate(CellType=factor(CellType, levels=unique(CellType)))
})

epiAnnoSummary <- tStats$firebrick4 %>%
  na.omit() %>%
  mutate(CellType=roadmapMeta$`Standardised epigenome name`[match(CellTypeCode,roadmapMeta$EID)]) %>%
  group_by(CellType) %>%
  # filter(all(epiFeatures %in% EpiFeature)) %>%
  summarise(CellTypeAbsMean=sum(abs(tStat))) %>%
  arrange(desc(CellTypeAbsMean)) %>%
  mutate(CellType=factor(CellType, levels=unique(CellType)))



