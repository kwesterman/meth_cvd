#!/bin/bash

FHS_SHARE_C1_PREFIX=~/kw/ncbi/dbGaP-60682_FHS_gen_IRBMDS/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v17.p10.c1.HMB-IRB-MDS/GenotypeFiles/phg000006.v9.FHS_SHARe_Affy500K.genotype-calls-matrixfmt.c1/subject_level_PLINK_sets/FHS_SHARe_Affy500K_subjects_c1
FHS_SHARE_C2_PREFIX=~/kw/ncbi/dbGaP-60683_FHS_gen_IRBNPUMDS/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v17.p10.c2.HMB-IRB-NPU-MDS/GenotypeFiles/phg000006.v9.FHS_SHARe_Affy500K.genotype-calls-matrixfmt.c2/subject_level_PLINK_sets/FHS_SHARe_Affy500K_subjects_c2
WHI_SHARE_C1_PREFIX=~/kw/ncbi/dbGaP-60684_WHI_gen_IRB/PhenoGenotypeFiles/ChildStudyConsentSet_phs000386.WHI.v7.p3.c1.HMB-IRB/GenotypeFiles/phg000061.v2.p1.WHI_NHLBI.genotype-calls-matrixfmt.c1/subject_level_filtered_PLINK_sets/NHLBI_SHARE_WHI_subject_level_c1
WHI_SHARE_C2_PREFIX=~/kw/ncbi/dbGaP-60685_WHI_gen_IRBNPU/PhenoGenotypeFiles/ChildStudyConsentSet_phs000386.WHI.v7.p3.c2.HMB-IRB-NPU/GenotypeFiles/phg000061.v2.p1.WHI_NHLBI.genotype-calls-matrixfmt.c2/subject_level_filtered_PLINK_sets/NHLBI_SHARE_WHI_subject_level_c2
WHI_GARNET_C1_PREFIX=~/kw/ncbi/dbGaP-60684_WHI_gen_IRB/PhenoGenotypeFiles/ChildStudyConsentSet_phs000315.WHI.v7.p3.c1.HMB-IRB/GenotypeFiles/phg000139.v1.GARNET_WHI.genotype-calls-matrixfmt.c1/subject_level_filtered_PLINK_sets/GARNET_WHI_TOP_subject_level_c1
WHI_GARNET_C2_PREFIX=~/kw/ncbi/dbGaP-60685_WHI_gen_IRBNPU/PhenoGenotypeFiles/ChildStudyConsentSet_phs000315.WHI.v7.p3.c2.HMB-IRB-NPU/GenotypeFiles/phg000139.v1.GARNET_WHI.genotype-calls-matrixfmt.c2/subject_level_filtered_PLINK_sets/GARNET_WHI_TOP_subject_level_c2

# Tables for conversion of variant IDs to RS IDs
cat /cluster/kappa/90-days-archive/ordovaslab/ncbi/dbGaP-60682_FHS_gen_IRBMDS/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v17.p10.c1.HMB-IRB-MDS/GenotypeFiles/phg000006.v9.FHS_SHARe_Affy500K.marker-info.MULTI/Original_marker_info/phg000006.FHS.genotype-calls.Affy500K.v3.p3.marker-info | tail -n +25 | cut -d ',' -f 4,5 | sed 's/,/\t/g' > ../int/affy500k_ssToRs.txt
cat /cluster/kappa/90-days-archive/ordovaslab/ncbi/dbGaP-60684_WHI_gen_IRB/PhenoGenotypeFiles/ChildStudyConsentSet_phs000386.WHI.v7.p3.c1.HMB-IRB/GenotypeFiles/phg000061.v2.p1.WHI_NHLBI.marker-info.MULTI/GenomeWideSNP_6.na30.annot.csv | tail -n +24 | cut -d ',' -f 1,2 | awk '{gsub(/\"/,"")};1' | sed 's/,/\t/g' > ../int/affy6_subjToRs.txt

# Table for conversion of GARNET sample IDs to shareids
garnet_subjIDtoShareid=$(cat /cluster/kappa/90-days-archive/ordovaslab/ncbi/dbGaP-60684_WHI_gen_IRB/PhenoGenotypeFiles/ChildStudyConsentSet_phs000315.WHI.v7.p3.c1.HMB-IRB/GenotypeFiles/phg000139.v1.GARNET_WHI.sample-info.MULTI/phg000139.v1_genotype_released_files_manifest.txt | tail -n +17 | cut -f 2,16)
GARNET_ALLFAM=$(cat ${WHI_GARNET_C1_PREFIX}.fam ${WHI_GARNET_C2_PREFIX}.fam)
join -1 2 -2 1 -o 1.1,1.2,1.1,2.2 <(echo "$GARNET_ALLFAM" | sort -k 2b,2) <(echo "$garnet_subjIDtoShareid" | sort -k 1b,1 | uniq) > ../int/garnet_famUpdate.txt

#join -1 2 -2 1 -o 1.1,1.2,1.1,2.2 <(sort ${WHI_GARNET_C1_PREFIX}.fam -k 2b,2) <(echo "$garnet_subjIDtoShareid" | sort -k 1b,1) > garnet_famUpdate_c1.txt
#join -1 2 -2 1 -o 1.1,1.2,1.1,2.2 <(sort ${WHI_GARNET_C2_PREFIX}.fam -k 2b,2) <(echo "$garnet_subjIDtoShareid" | sort -k 1b,1) > garnet_famUpdate_c2.txt

# Create new bed/bim/fam files with rsIDs and shareids

#FHS_SHARE_PREFIXES=($FHS_SHARE_C1_PREFIX, $FHS_SHARE_C2_PREFIX)
#for p in $FHS_SHARE_PREFIXES; do
#	plink --bfile $p --update-name ../int/affy500k_ssToRs.txt --make-bed --out ../int/fhsShareUpdate
#done
plink --bfile $FHS_SHARE_C1_PREFIX --bmerge $FHS_SHARE_C2_PREFIX --update-name ../int/affy500k_ssToRs.txt --make-bed --out ../int/fhsShareUpdate

#WHI_SHARE_PREFIXES=($WHI_SHARE_C1_PREFIX, $WHI_SHARE_C2_PREFIX)
#for p in $WHI_SHARE_PREFIXES; do
#	plink --bfile $p --update-name ../int/affy6_subjToRs.txt --make-bed --out ../int/whiShareUpdate
#done
plink --bfile $WHI_SHARE_C1_PREFIX --bmerge $WHI_SHARE_C2_PREFIX --update-name ../int/affy6_subjToRs.txt --make-bed --out ../int/whiShareUpdate

plink --bfile $WHI_GARNET_C1_PREFIX --bmerge $WHI_GARNET_C2_PREFIX --update-ids ../int/garnet_famUpdate.txt --make-bed --out ../int/whiGarnetUpdate
#plink --bfile $WHI_GARNET_C2_PREFIX --update-ids ../int/garnet_famUpdate_c2.txt --make-bed --out ../int/whiGarnetUpdate

ALL_SHAREIDS=$(cat ../int/metaData.csv | tail -n +2 | cut -d ',' -f 2)
join -1 2 -2 1 -o 1.1,1.2 <(cat ../int/fhsShareUpdate.fam ../int/whiShareUpdate.fam ../int/whiGarnetUpdate.fam | sort -k 2b,2) <(echo "$ALL_SHAREIDS" | sort) | cut -f 1,2 | uniq > ../int/allShareIDs.txt

printf "../int/fhsShareUpdate\n../int/whiShareUpdate\n../int/whiGarnetUpdate" > ../int/mergeList.txt

plink --merge-list ../int/mergeList.txt --make-bed --out ../int/mergeTest

plink --bfile ../int/fhsShareUpdate --exclude ../int/mergeTest-merge.missnp --make-bed --out ../int/fhsShareUpdate_tmp
plink --bfile ../int/whiShareUpdate --exclude ../int/mergeTest-merge.missnp --make-bed --out ../int/whiShareUpdate_tmp
plink --bfile ../int/whiGarnetUpdate --exclude ../int/mergeTest-merge.missnp --make-bed --out ../int/whiGarnetUpdate_tmp

printf "../int/fhsShareUpdate_tmp\n../int/whiShareUpdate_tmp\n../int/whiGarnetUpdate_tmp" > ../int/mergeList_tmp.txt

plink --merge-list ../int/mergeList_tmp.txt --keep ../int/allShareIDs.txt  --recode A --out test
