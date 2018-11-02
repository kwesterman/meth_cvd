# Miscellaneous helper functions conaining functionality from 
# relevant papers/methods


run_CPA <- function(RGset) {
  # Performs control probe adjustment (CPA) by running PCA on positive control probes from 450k array
  # Code adapted from Lehne et al. 2015
  # Args: 
  #   RGset: unnormalized minfi RGChannelSet
  
  # Type II probes
  TypeII.Name <- getProbeInfo(RGset, type = "II")$Name
  TypeII.Green <- getGreen(RGset)[getProbeInfo(RGset, type = "II")$AddressA,]
  TypeII.Red <- getRed(RGset)[getProbeInfo(RGset, type = "II")$AddressA,]
  rownames(TypeII.Red) <- TypeII.Name
  colnames(TypeII.Red) <- sampleNames(RGset)
  rownames(TypeII.Green) <- TypeII.Name
  colnames(TypeII.Green) <- sampleNames(RGset)
  
  # Type I probes, split into green and red channels
  TypeI.Green.Name <- getProbeInfo(RGset, type = "I-Green")$Name
  TypeI.Green.M <- getGreen(RGset)[getProbeInfo(RGset, type = "I-Green")$AddressB,]
  rownames(TypeI.Green.M) <- TypeI.Green.Name
  colnames(TypeI.Green.M) <- sampleNames(RGset)
  TypeI.Green.U <- getGreen(RGset)[getProbeInfo(RGset, type = "I-Green")$AddressA,]
  rownames(TypeI.Green.U) <- TypeI.Green.Name
  colnames(TypeI.Green.U) <- sampleNames(RGset)
  
  TypeI.Red.Name <- getProbeInfo(RGset, type = "I-Red")$Name
  TypeI.Red.M <- getRed(RGset)[getProbeInfo(RGset, type = "I-Red")$AddressB,]
  rownames(TypeI.Red.M) <- TypeI.Red.Name
  colnames(TypeI.Red.M) <- sampleNames(RGset)
  TypeI.Red.U <- getRed(RGset)[getProbeInfo(RGset, type = "I-Red")$AddressA,]
  rownames(TypeI.Red.U) <- TypeI.Red.Name
  colnames(TypeI.Red.U) <- sampleNames(RGset)
  
  #BSC1 control probes
  BSCI.Green.Name =getProbeInfo(RGset, type = "Control")[16:18,]$ExtendedType
  BSCI.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(BSCI.Green.Name), dimnames = list(BSCI.Green.Name, sampleNames(RGset)))
  BSCI.Green[BSCI.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[16:18,]$Address,]
  BSCI.Red.Name =getProbeInfo(RGset, type = "Control")[22:24,]$ExtendedType
  BSCI.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(BSCI.Red.Name), dimnames = list(BSCI.Red.Name, sampleNames(RGset)))
  BSCI.Red[BSCI.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[22:24,]$Address,]
  
  #BSC2 control probes
  BSCII.Red.Name =getProbeInfo(RGset, type = "Control")[28:31,]$ExtendedType
  BSCII.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(BSCII.Red.Name), dimnames = list(BSCII.Red.Name, sampleNames(RGset)))
  BSCII.Red[BSCII.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[28:31,]$Address,]
  
  #STAINING
  stain.Red.Name =getProbeInfo(RGset, type = "Control")[2,]$ExtendedType
  stain.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(stain.Red.Name), dimnames = list(stain.Red.Name, sampleNames(RGset)))
  stain.Red[stain.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[2,]$Address,]
  stain.Green.Name =getProbeInfo(RGset, type = "Control")[4,]$ExtendedType
  stain.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(stain.Green.Name), dimnames = list(stain.Green.Name, sampleNames(RGset)))
  stain.Green[stain.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[4,]$Address,]
  
  #EXTENSION
  extensionA.Red.Name =getProbeInfo(RGset, type = "Control")[7,]$ExtendedType
  extensionA.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(extensionA.Red.Name), dimnames = list(extensionA.Red.Name, sampleNames(RGset)))
  extensionA.Red[extensionA.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[7,]$Address,]
  extensionT.Red.Name =getProbeInfo(RGset, type = "Control")[8,]$ExtendedType
  extensionT.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(extensionT.Red.Name), dimnames = list(extensionT.Red.Name, sampleNames(RGset)))
  extensionT.Red[extensionT.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[8,]$Address,]
  extensionC.Green.Name =getProbeInfo(RGset, type = "Control")[9,]$ExtendedType
  extensionC.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(extensionC.Green.Name), dimnames = list(extensionC.Green.Name, sampleNames(RGset)))
  extensionC.Green[extensionC.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[9,]$Address,]
  extensionG.Green.Name =getProbeInfo(RGset, type = "Control")[10,]$ExtendedType
  extensionG.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(extensionG.Green.Name), dimnames = list(extensionG.Green.Name, sampleNames(RGset)))
  extensionG.Green[extensionG.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[10,]$Address,]
  
  #HYBRIDISATION
  hybridH.Green.Name =getProbeInfo(RGset, type = "Control")[11,]$ExtendedType
  hybridH.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(hybridH.Green.Name), dimnames = list(hybridH.Green.Name, sampleNames(RGset)))
  hybridH.Green[hybridH.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[11,]$Address,]
  hybridM.Green.Name =getProbeInfo(RGset, type = "Control")[12,]$ExtendedType
  hybridM.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(hybridM.Green.Name), dimnames = list(hybridM.Green.Name, sampleNames(RGset)))
  hybridM.Green[hybridM.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[12,]$Address,]
  hybridL.Green.Name =getProbeInfo(RGset, type = "Control")[13,]$ExtendedType
  hybridL.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(hybridL.Green.Name), dimnames = list(hybridL.Green.Name, sampleNames(RGset)))
  hybridL.Green[hybridL.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[13,]$Address,]
  
  #TARGET REMOVAL
  target.Green.Name =getProbeInfo(RGset, type = "Control")[14:15,]$ExtendedType
  target.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(target.Green.Name), dimnames = list(target.Green.Name, sampleNames(RGset)))
  target.Green[target.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[14:15,]$Address,]
  
  #Specificity I
  specI.Green.Name =getProbeInfo(RGset, type = "Control")[32:34,]$ExtendedType
  specI.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(specI.Green.Name), dimnames = list(specI.Green.Name, sampleNames(RGset)))
  specI.Green[specI.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[32:34,]$Address,]
  specI.Red.Name =getProbeInfo(RGset, type = "Control")[38:40,]$ExtendedType
  specI.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(specI.Red.Name), dimnames = list(specI.Red.Name, sampleNames(RGset)))
  specI.Red[specI.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[38:40,]$Address,]
  
  #Specificity II
  specII.Red.Name =getProbeInfo(RGset, type = "Control")[44:46,]$ExtendedType
  specII.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(specII.Red.Name), dimnames = list(specII.Red.Name, sampleNames(RGset)))
  specII.Red[specII.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[44:46,]$Address,]
  
  #NON POLYMORPHIC
  np.Red.Name =getProbeInfo(RGset, type = "Control")[47:48,]$ExtendedType
  np.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(np.Red.Name), dimnames = list(np.Red.Name, sampleNames(RGset)))
  np.Red[np.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[47:48,]$Address,]
  np.Green.Name =getProbeInfo(RGset, type = "Control")[49:50,]$ExtendedType
  np.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(np.Green.Name), dimnames = list(np.Green.Name, sampleNames(RGset)))
  np.Green[np.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[49:50,]$Address,]
  
  #Normalisation
  control=getProbeInfo(RGset, type = "Control")
  normC.Green.Name=control[control[,2]=='NORM_C',4]
  normC.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(normC.Green.Name), dimnames = list(normC.Green.Name, sampleNames(RGset)))
  normC.Green[normC.Green.Name,] <- getGreen(RGset)[control[control[,2]=='NORM_C',1],]
  normG.Green.Name=control[control[,2]=='NORM_G',4]
  normG.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(normG.Green.Name), dimnames = list(normG.Green.Name, sampleNames(RGset)))
  normG.Green[normG.Green.Name,] <- getGreen(RGset)[control[control[,2]=='NORM_G',1],]
  normA.Red.Name=control[control[,2]=='NORM_A',4]
  normA.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(normA.Red.Name), dimnames = list(normA.Red.Name, sampleNames(RGset)))
  normA.Red[normA.Red.Name,] <- getRed(RGset)[control[control[,2]=='NORM_A',1],]
  normT.Red.Name=control[control[,2]=='NORM_T',4]
  normT.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(normT.Red.Name), dimnames = list(normT.Red.Name, sampleNames(RGset)))
  normT.Red[normT.Red.Name,] <- getRed(RGset)[control[control[,2]=='NORM_T',1],]
  
  #combine ctrl probe intensities
  ctrl <- rbind(as.matrix(BSCI.Green), as.matrix(BSCI.Red), as.matrix(BSCII.Red), (stain.Red), (stain.Green), 
                (extensionA.Red), (extensionT.Red), (extensionC.Green), (extensionG.Green), (hybridH.Green), (hybridM.Green), 
                (hybridL.Green),as.matrix(target.Green),as.matrix(specI.Green),as.matrix(specI.Red), 
                as.matrix(specII.Red),(np.Red[1,]),(np.Red[2,]),(np.Green[1,]),(np.Green[2,]),
                as.matrix(normC.Green),as.matrix(normG.Green), as.matrix(normA.Red),as.matrix(normT.Red))
  pca.fit <- prcomp(na.omit(t(ctrl)))
  pca.fit
}


refactor <- function(
  betavals, k, covarfile=NULL, t=500, numcomp=NULL, 
  stdth=0.02, out="refactor") {
  ## Performs sparse PCA to return k components theoretically representing cell counts or other substructure
  ## Adapted from R code based on Rahmani et al. 2016
  
  ranked_filename = paste(out, ".out.rankedlist.txt", sep="")
  components_filename = paste(out, ".out.components.txt", sep="")
  
  print('Starting ReFACTor v1.0...');
  
  print('Reading input files...');
  
  # O = as.matrix(read.table(data_file))
  # O = Mvals
  # sample_id <- O[1, -1] # extract samples ID
  # O <- O[-1,] # remove sample ID from matrix
  # cpgnames <- O[, 1] ## set rownames
  # O <- O[, -1] 
  sample_id <- colnames(betavals)
  cpgnames <- rownames(betavals)
  O <- betavals
  O = matrix(as.numeric(O),nrow=nrow(O),ncol=ncol(O))
  
  print(paste("Excluding sites with low variance (std < ", stdth, ")..."), sep="")
  sds = apply(t(O), 2, sd)
  m_before = length(sds)
  include = which(sds >= stdth)
  O = O[include,]
  cpgnames = cpgnames[include]
  print(paste((m_before - length(which(sds >= stdth))), " sites were excluded due to low variance...", sep=""))
  
  if (is.null(numcomp) || is.na(numcomp)) 
  {
    numcomp = k
  }
  
  # Adjust the data for the covariates
  if (!is.null(covarfile))
  {
    covs = as.matrix(read.table(covarfile))
    sample_id2 <- covs[, 1]
    if (!all(sample_id == sample_id2)){
      print("ERROR: The order of the samples in the covariates file must be the same as the order in the data file")
      quit()
    }
    covs <- covs[,-1]
    if (length(covs) > dim(O)[2])
    {
      covs = matrix(as.numeric(covs),nrow=nrow(covs),ncol=ncol(covs))
    }else{
      covs = as.numeric(covs)
    }    
    O_adj = O
    for (site in 1:nrow(O))
    {
      model <- lm(O[site,] ~  covs)
      O_adj[site,] = residuals(model)
    }
    O = O_adj
  }
  
  print('Running a standard PCA...')
  pcs = prcomp(scale(t(O)));
  
  coeff = pcs$rotation
  score = pcs$x
  
  print('Compute a low rank approximation of input data and rank sites...')
  x = score[,1:k]%*%t(coeff[,1:k]);
  An = scale(t(O),center=T,scale=F)
  Bn = scale(x,center=T,scale=F)
  An = t(t(An)*(1/sqrt(apply(An^2,2,sum))))
  Bn = t(t(Bn)*(1/sqrt(apply(Bn^2,2,sum))))
  
  
  # Find the distance of each site from its low rank approximation.
  distances = apply((An-Bn)^2,2,sum)^0.5 ;
  dsort = sort(distances,index.return=T);
  ranked_list = dsort$ix
  
  print('Compute ReFACTor components...')
  sites = ranked_list[1:t];
  pcs = prcomp(scale(t(O[sites,])));
  first_score <- score[,1:k];
  score = pcs$x
  
  # print('Saving a ranked list of the data features...');
  # write(t(cpgnames[ranked_list]),file=ranked_filename,ncol=1)
  # #write(t(cbind(ranked_list,cpgnames[ranked_list])),file=ranked_filename,ncol=2)
  # 
  # print('Saving the ReFACTor components...');
  # write(t(score[,1:numcomp]), file=components_filename, ncol=numcomp)
  
  print('ReFACTor is Done');
  result <- list(refactor_components=score[,1:numcomp], ranked_list=ranked_list, standard_pca=first_score) 
  return(result$refactor_components)
}


calc_FRS <- function(pData) {
  # Calculation of Framingham risk score based on D'Agostino 2008 
  # doi: https://doi.org/10.1161/CIRCULATIONAHA.107.699579
  
  FRS_data <- pData %>%
    mutate(age=pmin(pmax(age, 30), 74),  # Constrain continuous values to within specific ranges
           chol=pmin(pmax(chol, 100), 405),
           hdl=pmin(pmax(hdl, 10), 100),
           sbp=pmin(pmax(sbp, 90), 200)) %>%
    mutate(logAge=log(age),  # Take logs of continuous values
           logTC=log(chol),
           logHDL=log(hdl),
           logSBP=log(sbp),
           smoking=smk_now) %>%
    select(sex, logAge, logTC, logHDL, ht_med, logSBP, smoking, diabetes)

  with(FRS_data, {
       weightedSum <- ifelse(
         sex=="M",
         (3.06117 * logAge + 1.12370 * logTC - 0.93263 * logHDL + 
            1.93303 * logSBP * (1 - ht_med) + 1.99881 * logSBP * (ht_med) + 
            0.65451 * smoking + 0.57367 * diabetes),
         (2.32888 * logAge + 1.20904 * logTC - 0.70833 * logHDL + 
            2.76157 * logSBP * (1 - ht_med) + 2.82263 * logSBP * (ht_med) + 
            0.52973 * smoking + 0.69154 * diabetes))
       frs <- ifelse(
         sex=="M",
         1 - 0.88936 ^ (exp(weightedSum - 23.9802)),
         1 - 0.95012 ^ (exp(weightedSum - 26.1931)))
       frs
  })
}


calc_zhang_mrs <- function(betas) {
  # Calculate mortality risk score from Zhang et al. 2017
  
  zhang_coefs <- data.frame(cpg=c("cg01612140","cg05575921","cg06126421","cg08362785","cg10321156",
                                  "cg14975410","cgcg19572487","cg23665802","cg24704287","cg25983901"),
                            coef=c(-0.38253,-0.92224,-1.70129,2.71749,-0.02073,-0.04156,
                                   -0.28069,-0.89440,-2.98637,-1.80325))
  zhang_subset <- betas[zhang_coefs$cpg,]
  mrs_vals <- t(zhang_subset) %*% zhang_coefs$coef
  data.frame(sampleKey=colnames(betas), zhang_mrs=mrs_vals, stringsAsFactors=F)
}


calc_khera2016grs <- function(grs_vcf) {
  
  riskAlleleString <- c("ATAAGTTCGTTTGCTGACCCCTCCCATGGCGACCTGTTACTCTTTGCTCA")
  khera2016grs <- tibble(snp=c("rs599839","rs17114036","rs11206510","rs4845625","rs17465637","rs1561198",
                               "rs6544713","rs515135","rs2252641","rs6725887","rs9818870",
                               "rs1878406","rs7692387","rs273909","rs10947789","rs17609940",
                               "rs12526453","rs12190287","rs2048327","rs3798220","rs10455872",
                               "rs4252120","rs2023938","rs10953541","rs11556924","rs2954029",
                               "rs3217992","rs4977574","rs579459","rs2505083","rs2047009","rs501120",
                               "rs2246833","rs12413409","rs974819","rs964184","rs2259816","rs3184504",
                               "rs9319428","rs4773144","rs9515203","rs2895811","rs3825807","rs7173743",
                               "rs17514846","rs12936587","rs216172","rs46522","rs1122608","rs9982601"),
                         riskAllele=strsplit(riskAlleleString, "")[[1]],
                         riskEstimate=c(1.11,1.11,1.08,1.04,1.14,1.05,1.06,1.08,1.04,1.12,1.07,1.06,
                                        1.06,1.09,1.06,1.07,1.1,1.07,1.06,1.51,1.45,1.06,1.07,
                                        1.08,1.09,1.04,1.16,1.29,1.07,1.06,1.05,1.07,1.06,1.12,1.07,
                                        1.13,1.08,1.07,1.05,1.07,1.08,1.06,1.08,1.07,1.05,1.06,1.07,
                                        1.06,1.1,1.13))
  
  manualRiskAnnotations <- c(rs9818870="T", rs1878406="T", rs12190287="C", rs17514846="A", rs12413409="G", 
                             rs501120="T", rs4773144="G", rs3825807="A", rs1122608="G", rs4845625="T", 
                             rs17465637="C", rs11206510="T", rs17114036="A", rs9982601="T", 
                             rs2252641="C", rs273909="G", rs12526453="C", rs10455872="G", rs579459="T")
  
  vcf_anno_portion <- dplyr::select(grs_vcf, ID, REF, ALT)
  grs_calculator.df <- inner_join(vcf_anno_portion, khera2016grs, by=c("ID"="snp")) %>%
    mutate(riskAlleleCorrected=case_when(riskAllele==REF | riskAllele==ALT ~ riskAllele,
                                         ID %in% names(manualRiskAnnotations) ~ 
                                           manualRiskAnnotations[ID],
                                         TRUE ~ as.character(NA)),
           weight=ifelse(riskAlleleCorrected==ALT, log(riskEstimate), -log(riskEstimate)))
  
  
  geno_mat <- t(as.matrix(grs_vcf[,-(1:9)]))
  colnames(geno_mat) <- grs_vcf$ID
  geno_mat <- ifelse(geno_mat=="0/0", 0, ifelse(geno_mat=="0/1", 1, ifelse(geno_mat=="1/1", 1, NA)))
  grs <- as.vector(geno_mat[,grs_calculator.df$ID] %*% grs_calculator.df$weight)
  grs
  # fix_risk_allele <- function(rsid, riskAllele, annotatedAlleles) {
  #   case_when(length(strsplit(annotatedAlleles, "/")[[1]])>2 ~ manualRiskAnnotations[rsid],
  #             grepl(riskAllele, annotatedAlleles) ~ riskAllele,
  #             TRUE ~ c(A="T", C="G", G="C", `T`="A")[riskAllele])
  # }
  
  
  
  #   mutate(riskAlleleCorrected=pmap_chr(list(snp, riskAllele, allele), fix_risk_allele),
  #          loc=paste0(chr_name, ":", chrom_start),
  #          lnOR=log(riskEstimate)) %>%
  #   dplyr::rename(alleleOptions=allele) %>%
  #   dplyr::select(snp, loc, alleleOptions, riskAlleleCorrected, lnOR)
  # 
  # 
  # # library(biomaRt)
  # # ensembl <- useMart("ENSEMBL_MART_SNP", host="grch37.ensembl.org", dataset="hsapiens_snp")
  # # snpLocs <- getBM(attributes=c("refsnp_id","chr_name","chrom_start","allele"),
  # #                  filters="snp_filter",
  # #                  values=khera2016grs$snp,
  # #                  mart=ensembl)
  # 
  # grs_calculator.df <- inner_join(khera2016grs, snpLocs, by=c("snp"="refsnp_id")) %>%
  #   mutate(riskAlleleCorrected=pmap_chr(list(snp, riskAllele, allele), fix_risk_allele),
  #          loc=paste0(chr_name, ":", chrom_start),
  #          lnOR=log(riskEstimate)) %>%
  #   dplyr::rename(alleleOptions=allele) %>%
  #   dplyr::select(snp, loc, alleleOptions, riskAlleleCorrected, lnOR)
  # 
  # genotypeAnno <- tibble(locAllele=colnames(genotypes)) %>%
  #   separate(locAllele, into=c("loc","allele"), sep="_") %>%
  #   inner_join(grs_calculator.df, by="loc") %>%
  #   # mutate(genoInOptions=map2_lgl(allele, alleleOptions, function(x,y) x %in% strsplit(y, "/")[[1]]))
  #   mutate(needsFlip=ifelse(allele==riskAlleleCorrected, F, T),
  #          grsCoef=ifelse(needsFlip, -lnOR, lnOR),
  #          col_name=paste(loc, allele, sep="_"))
  # 
  # grs_weights <- setNames(genotypeAnno$grsCoef, genotypeAnno$col_name)[colnames(genotypes)]
  # 
  # grs <- as.vector(genotypes %*% grs_weights)
  # grs
}