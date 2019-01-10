#!/bin/bash


module load plink2

METHDIR=../data/fhs/meth
GENODIR=../int/plinksets

# Gather share IDs of FHS Offspring participants with methylation data available
cat $METHDIR/sample_attributes_c1.txt | tail -n +12 | awk '{print $2,$2}' | sed 's/_724//g' > meth_ids.txt
cat $METHDIR/sample_attributes_c2.txt | tail -n +12 | awk '{print $2,$2}' | sed 's/_724//g' >> meth_ids.txt

# Collect genotypes and keep track of low-quality imputed SNPs by r2
fhs_c1=/cluster/home/kweste01/kw/ncbi/dbGaP-60162_FHS_gen_IRBMDS/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v17.p10.c1.HMB-IRB-MDS/GenotypeFiles/\
phg000679.v1.FHS_SHARe_Imputed_1000G.genotype-imputed-data.c1
fhs_c2=/cluster/home/kweste01/kw/ncbi/dbGaP-60163_FHS_gen_IRBNPUMDS/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v17.p10.c2.HMB-IRB-NPU-MDS/GenotypeFiles/\
phg000679.v1.FHS_SHARe_Imputed_1000G.genotype-imputed-data.c2

for ((i=1; i<=22; i++)); do
	# For each chromosome, keep track of low-quality imputations and use dose2plink
	# to convert dbGaP dosage files into plink-readable format
	gunzip < $fhs_c1/imputed-metrics/machout.chr$i.info.gz > tmp_fhs_chr$i.info
	awk '{if ($7 < 0.9) print $1}' tmp_fhs_chr$i.info > fhs_chr${i}_low_qual_snps.txt
	#gunzip < $fhs_c1/machout.chr$i.dose_GRU.gz > tmp_fhs_chr${i}_c1.dose
	#gunzip < $fhs_c2/machout.chr$i.dose_NPU.gz > tmp_fhs_chr${i}_c2.dose
	#cat tmp_fhs_chr${i}_c1.dose tmp_fhs_chr${i}_c2.dose > tmp_fhs_chr$i.dose
	#awk '{print $1"->"$0}' tmp_fhs_chr$i.dose > tmp_fhs_chr$i.dosefix
	#~/kw/opt/dose2plink -dose tmp_fhs_chr$i.dosefix \
	#	-info tmp_fhs_chr$i.info \
	#	-out $GENODIR/fhs_chr$i
	rm tmp_fhs*
done

# "Custom" merge of chromosomes
#gunzip < $GENODIR/fhs_chr1.pdat.gz > $GENODIR/fhs.pdat
#for ((i=2; i<=22; i++)); do
#	gunzip < $GENODIR/fhs_chr$i.pdat.gz | tail -n +2 >> $GENODIR/fhs.pdat
#done

# Housekeeping prep for plink2 filters
cat fhs_chr*_low_qual_snps.txt > fhs_low_qual_snps.txt
SNPANNO=../int/snp_annotations/snp_annot_hg19_nodups.txt  # Contains 1000 Genomes SNP locations, rsIDs, and alleles from Ensembl

# Processing of dosages using plink2: 
# exclude SNPs w/ r2 < 0.9 gathered above, 
# update SNP annotations from dbSNP file, 
# keep only subjects with methylation data available,
# convert to hard-calls (dosage distance to integer of <0.2, otherwise set missing)
# output a plink1 fileset (bed/bim/fam)
plink2 --import-dosage $GENODIR/fhs.pdat \
	--psam $GENODIR/fhs_chr1.pfam \
	--exclude fhs_low_qual_snps.txt \
	--keep meth_ids.txt \
	--hard-call-threshold 0.2 \
	--make-bed \
	--out tmp_fhs_with_meth

plink2 --bfile tmp_fhs_with_meth \
	--update-name $SNPANNO 3 6 \
	--ref-allele $SNPANNO 4 3 \
	--make-bed \
	--out fhs_with_meth
# Updating of SNP annotation must come in a separate command to allow --exclude to work properly
