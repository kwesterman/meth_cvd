#!/bin/bash

module load gcta
module load plink2

GENODIR=../int/plinksets

# Use dose2plink to convert minimac output to plink 1 dosage format
fhs_c1=/cluster/home/kweste01/kw/ncbi/dbGaP-60162_FHS_gen_IRBMDS/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v17.p10.c1.HMB-IRB-MDS/GenotypeFiles/\
phg000679.v1.FHS_SHARe_Imputed_1000G.genotype-imputed-data.c1
fhs_c2=/cluster/home/kweste01/kw/ncbi/dbGaP-60163_FHS_gen_IRBNPUMDS/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v17.p10.c2.HMB-IRB-NPU-MDS/GenotypeFiles/\
phg000679.v1.FHS_SHARe_Imputed_1000G.genotype-imputed-data.c2

#for ((i=1; i<=22; i++)); do
#	gunzip < $fhs_c1/imputed-metrics/machout.chr$i.info.gz > tmp_fhs_chr$i.info
#	awk '{if ($7 < 0.3) print $1}' tmp_fhs_chr$i.info > $GENODIR/fhs_chr${i}_low_qual_snps.txt
#	gunzip < $fhs_c1/machout.chr$i.dose_GRU.gz > tmp_fhs_chr${i}_c1.dose
#	gunzip < $fhs_c2/machout.chr$i.dose_NPU.gz > tmp_fhs_chr${i}_c2.dose
#	cat tmp_fhs_chr${i}_c1.dose tmp_fhs_chr${i}_c2.dose > tmp_fhs_chr$i.dose
#	awk '{print $1"->"$0}' tmp_fhs_chr$i.dose > tmp_fhs_chr$i.dosefix
#	~/kw/opt/dose2plink -dose tmp_fhs_chr$i.dosefix \
#		-info tmp_fhs_chr$i.info \
#		-out $GENODIR/fhs_chr$i
#	rm tmp_fhs*
#done

# "Custom" merge of chromosomes
#gunzip < $GENODIR/fhs_chr1.pdat.gz > $GENODIR/fhs.pdat
#for ((i=2; i<=22; i++)); do
#	gunzip < $GENODIR/fhs_chr$i.pdat.gz | tail -n +2 >> $GENODIR/fhs.pdat
#done

# Housekeeping prep for plink2 filters
cat $GENODIR/fhs_chr*_low_qual_snps.txt > $GENODIR/fhs_low_qual_snps.txt
awk -F, '{print $2,$2}' ../int/metaData.csv > ../int/meth_ids.txt
SNPANNO=../int/snp_annotations/snp_annot_hg19_nodups.txt  # Contains 1000 Genomes SNP locations, rsIDs, and alleles from Ensembl

# Preprocessing and filtering of dosages in plink2
plink2 --import-dosage $GENODIR/fhs.pdat \
	--psam $GENODIR/fhs_chr1.pfam \
	--exclude $GENODIR/fhs_low_qual_snps.txt \
	--keep ../int/meth_ids.txt \
	--make-pgen \
	--out $GENODIR/tmp_fhs

plink2 --pfile $GENODIR/tmp_fhs \
	--update-name $SNPANNO 3 6 \
	--ref-allele $SNPANNO 4 3 \
	--make-pgen \
	--out $GENODIR/fhs

rm $GENODIR/tmp_fhs*
