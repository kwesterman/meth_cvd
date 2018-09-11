#!/bin/bash

### Convert dbGaP genotype imputation files from Mach into usable additive genotype matrices (0/1/2)


module load gcta
module load plink

# FHS
fhs_c1=/cluster/home/kweste01/kw/ncbi/dbGaP-60162_FHS_gen_IRBMDS/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v17.p10.c1.HMB-IRB-MDS/GenotypeFiles/\
phg000679.v1.FHS_SHARe_Imputed_1000G.genotype-imputed-data.c1
fhs_c2=/cluster/home/kweste01/kw/ncbi/dbGaP-60163_FHS_gen_IRBNPUMDS/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v17.p10.c2.HMB-IRB-NPU-MDS/GenotypeFiles/\
phg000679.v1.FHS_SHARe_Imputed_1000G.genotype-imputed-data.c2

for ((i=1; i<=22; i++)); do
	gunzip < $fhs_c1/imputed-metrics/machout.chr$i.info.gz > tmp_fhs_chr$i.info
	gunzip < $fhs_c1/machout.chr$i.dose_GRU.gz > tmp_fhs_c1_chr$i.dose
	gunzip < $fhs_c2/machout.chr$i.dose_NPU.gz > tmp_fhs_c2_chr$i.dose
	cat tmp_fhs_c1_chr$i.dose tmp_fhs_c2_chr$i.dose > tmp_fhs_chr$i.dose
	awk '{print $1"->"$1}' tmp_fhs_chr$i.dose > tmp_fhs_idcol.txt
	paste tmp_fhs_idcol.txt <(cut -f 2- tmp_fhs_chr$i.dose) > tmp_fhs_chr$i.dosefix
	~/kw/opt/dose2plink -dose tmp_fhs_chr$i.dosefix \
		-info tmp_fhs_chr$i.info \
		-out ../int/plinksets/fhs_chr$i
done

# List of SHARe IDs with methylation to keep
awk -F, '{print $1,$1}' ../int/metaData_fhs.csv > ../int/meth_ids.txt 

#for ((i=1; i<=22; i++)); do 
#	# Use GCTA to read in zipped dose and info files from imputation, filter for imputation quality and relevant subjects, and output a plink set
#	gcta64 --dosage-mach-gz ${fhs_c1}/machout.chr${i}.dose_GRU.gz ${fhs_c1}/imputed-metrics/machout.chr${i}.info.gz --keep ../int/meth_ids.txt --imput-rsq 0.3 \
#		--make-bed --out ../int/plink_imputed/fhs_c1_chr${i}
#	gcta64 --dosage-mach-gz ${fhs_c2}/machout.chr${i}.dose_NPU.gz ${fhs_c2}/imputed-metrics/machout.chr${i}.info.gz --keep ../int/meth_ids.txt --imput-rsq 0.3 \
#		--make-bed --out ../int/plink_imputed/fhs_c2_chr${i}
#done

SNPANNO=../int/snp_annotations/snp_annot_hg19_nodups.txt  # Contains 1000 Genomes SNP locations, rsIDs, and alleles from Ensembl

#ls ../int/plink_imputed/fhs_c*.bed | sed 's/.bed//g' > ../int/mergeListFhs.txt  # All chromosome plinksets for FHS
#
#plink --merge-list ../int/mergeListFhs.txt \
#--snps-only just-acgt \
#--update-name $SNPANNO 3 6 \
#--a2-allele $SNPANNO 4 3 \
#--geno 0.01 \
#--maf 0.01 \
#--make-bed \
#--out ../int/plink_imputed/fhs_genos
#
#plink --bfile ../int/plink_imputed/fhs_genos \
#--update-chr $SNPANNO 1 3 \
#--update-map $SNPANNO  2 3 \
#--make-bed --out ../int/plink_imputed/fhs_genos
#rm ../int/plink_imputed/*~

# Loop through bfiles, removing SNPs with identical alleles and trimming to only SNPs w/ ACGT notation
#ls ../int/plink_imputed/*.bed | sed 's/.bed//g' > ../int/mergeList.txt
#mkdir ../int/plink_imputed_noIdentical
#for pref in $(cat ../int/mergeList.txt); do
#	prefOnly=$(echo $pref | cut -d '/' -f 4)
#	awk '$5==$6 {print $2}' ${pref}.bim > ../int/identicalAlleles.txt
#	plink --bfile $pref --exclude ../int/identicalAlleles.txt --alleleACGT --make-bed --out ../int/plink_imputed_noIdentical/tmp
#	plink --bfile ../int/plink_imputed_noIdentical/tmp --snps-only just-acgt --make-bed --out ../int/plink_imputed_noIdentical/${prefOnly}
#	rm ../int/plink_imputed_noIdentical/tmp.*
#done

# ls ../int/plink_imputed_noIdentical/*.bed | sed 's/.bed//g' > ../int/mergeList_noIdentical.txt
# #join -1 2 -2 1 -o 1.1,1.1, <(cat ${whi_c1_root}/phg000592.v1.WHI_Imputation.sample-info.MULTI/phg000592.v1_release_manifest.txt | tail -n +17 | cut -f 1,2 | sort -k 2,2) <(cat ../int/meth_ids.txt | sort -k 1,1) | awk '{print $1,$1}' > ../int/whi_sampIDs.txt
# plink --merge-list ../int/mergeList_noIdentical.txt --make-bed --out ../int/plink_imputed_noIdentical/noIdenticalMerge  # Merges all datasets into single plink set with only relevant subjects and performs variant QC
# mkdir ../int/plink_imputed_noMultiAllele
# for pref in $(cat ../int/mergeList_noIdentical.txt); do
#         newPref=../int/plink_imputed_noMultiAllele/$(echo $pref | cut -d '/' -f 4)
#         plink --bfile $pref --exclude ../int/plink_imputed_noIdentical/noIdenticalMerge-merge.missnp --make-bed --out $newPref
# 	cat ${newPref}.bim > tmpBim
# 	cut -f 2 ${newPref}.bim | awk -F ':' '{print $1,$2}' > tmpLoc
# 	paste <(awk '{print $1}' tmpLoc) <(cut -f 2 tmpBim) <(awk '{print 0,$2}' OFS='\t' tmpLoc) <(awk '{print $5,$6}' OFS='\t' tmpBim) > newBim && mv newBim ${newPref}.bim
# done
#ls ../int/plink_imputed_noMultiAllele/*.bed | sed 's/.bed//g' > ../int/mergeList_noMultiAllele.txt
#mkdir -p ../int/finalGenos
#cat ${whi_c1_root}/phg000592.v1.WHI_Imputation.sample-info.MULTI/phg000592.v1_release_manifest.txt | tail -n +17 | cut -f 1,2 | awk '{print $1,$1,$2,$2}' | sort -u -k 3 > ../int/idUpdate.txt
#plink --merge-list ../int/mergeList_noMultiAllele.txt --geno 0.01 --maf 0.01 --update-ids ../int/idUpdate.txt --indep 50 5 1.5 --make-bed --out ../int/finalGenos/unpruned  # Generates filtered plink set and list of SNPs to prune
#plink --bfile ../int/finalGenos/unpruned --extract ../int/finalGenos/unpruned.prune.in --recode A --out ../int/finalGenos/pruned  # Generates a .raw file with 0/1/2 genotype calls
