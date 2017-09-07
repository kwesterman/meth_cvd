cardio_gxe <- read_delim("../data/literature/CM_GxEdata.txt", delim="\t") %>%
  filter(`Interaction significance`=="significant interaction")

library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19, quietly=T)
annot <- data.frame(cbind(Locations, Other, Manifest, SNPs.Illumina, Islands.UCSC))
allAnnotGenes <- unique(unlist(strsplit(annot$UCSC_RefGene_Name, split=";")))
intersectGenes <- intersect(toupper(allAnnotGenes), cardio_gxe$Gene)
annot_gxeGenes <- annot %>%
  separate_rows(UCSC_RefGene_Name, sep=";") %>%
  filter(UCSC_RefGene_Name %in% intersectGenes)
# alter above to include strsplit

load("../int/Mvals.RData")
gxeMvals <- Mvals[match(annot_gxeGenes$Name, rownames(Mvals)),]
save("annot_gxeGenes", "gxeMvals", file="../int/gxeMvals.RData")

################


# Diet
foodLib <- c(SFA="NUT_SATFAT", MUFA="NUT_MONFAT", PUFA="NUT_POLY", omega3="NUT_OMEGA",
             alcohol="NUT_ALCO", sucrose="NUT_SUCR", caffeine="NUT_CAFF")
library(readxl)
dietData_all <- read_excel("../data/diet/phs000007.v28.pht002350.v4.p10.c1.vr_ffreq_ex08_1_0615s.HMB-IRB-MDS_ex8_diet.xlsx", sheet=2)
dietData <- dplyr::select(dietData_all, shareid, foodLib)
names(dietData) <- c("shareid", names(foodLib))

# Phenotypes
load("../int/nonMethData.RData")
phenos <- inner_join(nonMethData, dietData, by="shareid")
survObj <- Surv(phenos$timeToEvent, even=phenos$event)

# Methylation
load("../int/gxeMvals.RData")
Mvals <- gxeMvals[,match(phenos$sampleKey, colnames(gxeMvals))]

cardio_gxe <- read_delim("../data/literature/CM_GxEdata.txt", delim="\t") %>%
  filter(`Interaction significance`=="significant interaction")
regInfo <- inner_join(cardio_gxe, annot_gxeGenes, by=c("Gene"="UCSC_RefGene_Name"))


trigsSFA <- filter(regInfo, 
                   # Phenotype=="HDL-C", 
                   # `Environmental factor`=="PUFA intake",
                   Name %in% rownames(Mvals))
trigsSFA <- sample_n()

funca <- function(dietVar) {
  # trigsSFA$interactionPval <- sapply(trigsSFA$Name, function(cpg) {
  a <- sapply(unique(trigsSFA$Name)[1:2000], function(cpg) {
    form <- as.formula("survObj~meth+age+sex+bmi+smoking_now+pastEvent+Gran+Mono+Bcell+NK+CD8T+PC1_cp")
    cox1 <- coxph(form, data=cbind(phenos, meth=Mvals[cpg,]))
    # summary(cox1)$coef[paste0('meth:',dietVar),'Pr(>|z|)']
    summary(cox1)$coef["meth","Pr(>|z|)"]
  })
  hist(a)
}
funca("PUFA")

load("../int/ewasRes7_pastEventsAdjustment.RData")
ewasRes <- res$wbcPastEvent1CPA
a <- data.frame(ewasRes, stringsAsFactors=F) %>%
  dplyr::mutate(p=as.numeric(Pr...z..)) %>%
  inner_join(annot, by=c("CpG"="Name")) %>%
  separate_rows(UCSC_RefGene_Name, sep=";") %>%
  mutate(cgeGene=UCSC_RefGene_Name %in% cardio_gxe$Gene) %>%
  distinct(CpG, .keep_all=T) %>%
  mutate(negLogP=-log(p))
ggplot(data=a, aes(x=DHS, y=negLogP)) + geom_boxplot()
hist(a$p)




load("../int/mrsData_1CPA_noPastEvents.RData")
mrsData$survObj <- Surv(time=mrsData$timeToEvent, event=mrsData$event)
coxph(survObj~mrs_half, data=mrsData, subset=!mrsData$inTrainSet)





