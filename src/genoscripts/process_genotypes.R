library(tidyverse)

metaData <- readRDS("../int/metaData.rds")

recode_vcf <- function(vcfMat) {
  # Take in a VCF matrix (genotypes only, no metadata in initial columns) and recode as 0/1/2
  apply(vcfMat, 2, function(callSet) {  # Need apply to avoid "long vector" issue when subsetting large matrices
    callSet[callSet=="0/0"] <- 0
    callSet[callSet=="0/1"] <- 1
    callSet[callSet=="1/1"] <- 2
    callSet[callSet=="./."] <- NA
    callSet
  })
}

## FHS

Affy500kMarkerInfo <- read_csv("../data/fhs/gen/Affy500kMarkerInfo.csv", skip=23, col_types=cols(`# Chr`="c"))
fhsCallsDF_c1 <- read_tsv("../data/fhs/gen/FHS_SHARe_Affy500K_subjects_c1.vcf", skip=30)
fhsCallsDF_c2 <- read_tsv("../data/fhs/gen/FHS_SHARe_Affy500K_subjects_c2.vcf", skip=30)

fhsCallsMat <- as.matrix(cbind(as.matrix(fhsCallsDF_c1[,-(1:9)]), as.matrix(fhsCallsDF_c2[,-(1:9)])))
rownames(fhsCallsMat) <- Affy500kMarkerInfo$`Rs#`[match(fhsCallsDF_c1$ID, Affy500kMarkerInfo$`Ss#`)]
colnames(fhsCallsMat) <- gsub(".*_", "", colnames(fhsCallsMat))
fhsCallsFinal <- recode_vcf(fhsCallsMat[,intersect(colnames(fhsCallsMat),metaData$subjID)])
rm(fhsCallsDF_c1, fhsCallsDF_c2, fhsCallsMat)

## WHI

Affy6MarkerInfo <- read_csv("../data/whi/gen/Affy6MarkerInfo.txt", skip=22,
                            col_types=cols_only(`Probe Set ID`="c", `dbSNP RS ID`="c"))
whiShareCallsDF_c1 <- read_tsv("../data/whi/gen/NHLBI_SHARE_WHI_subject_level_c1.vcf", skip=31)
whiShareCallsDF_c2 <- read_tsv("../data/whi/gen/NHLBI_SHARE_WHI_subject_level_c2.vcf", skip=31)

whiShareCallsMat <- as.matrix(cbind(whiShareCallsDF_c1[,-(1:9)], whiShareCallsDF_c2[,-(1:9)]))
rownames(whiShareCallsMat) <- Affy6MarkerInfo$`dbSNP RS ID`[match(whiShareCallsDF_c1$ID, 
                                                                  Affy6MarkerInfo$`Probe Set ID`)]
colnames(whiShareCallsMat) <- gsub(".*_", "", colnames(whiShareCallsMat))
whiShareCallsFinal <- recode_vcf(whiShareCallsMat[,intersect(colnames(whiShareCallsMat),metaData$subjID)])
rm(whiShareCallsDF_c1, whiShareCallsDF_c2, whiShareCallsMat)

# Omni1QuadMarkerInfo <- read_csv("../int/Omni1QuadMarkerInfo.csv", skip=7)
whiGarnetCallsDF_c1 <- read_tsv("../data/whi/gen/GARNET_WHI_TOP_subject_level_c1.vcf", skip=31)
whiGarnetCallsDF_c2 <- read_tsv("../data/whi/gen/GARNET_WHI_TOP_subject_level_c2.vcf", skip=31)

whiGarnetCallsMat <- as.matrix(cbind(whiGarnetCallsDF_c1[,-(1:9)], whiGarnetCallsDF_c2[,-(1:9)]))
rownames(whiGarnetCallsMat) <- whiGarnetCallsDF_c1$ID
whiGarnetSampleInfo <- read_tsv("../data/whi/gen/garnetSampleInfo.txt", skip=15)
colnames(whiGarnetCallsMat) <- whiGarnetSampleInfo$ShareID[match(gsub(".*_","",colnames(whiGarnetCallsMat)),
                                                                 whiGarnetSampleInfo$SubjectID)]
whiGarnetCallsFinal <- recode_vcf(whiGarnetCallsMat[,intersect(colnames(whiGarnetCallsMat),metaData$subjID)])
rm(whiGarnetCallsDF_c1, whiGarnetCallsDF_c2, whiGarnetCallsMat)

commonSnps <- Reduce(intersect, list(rownames(fhsCallsFinal), 
                                     rownames(whiShareCallsFinal), rownames(whiGarnetCallsFinal)))
allGenotypes <- cbind(fhsCallsFinal[commonSnps,],
                      whiShareCallsFinal[commonSnps,],
                      whiGarnetCallsFinal[commonSnps,])
saveRDS(allGenotypes, file="../int/genotypes.rds", compress=F)

