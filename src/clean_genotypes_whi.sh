#!/bin/bash

### Convert dbGaP genotype imputation files from Mach into usable additive genotype matrices (0/1/2)

module load gcta
module load plink
module load plink2

GENODIR=../int/plinksets

# WHI
whi_c1_root=/cluster/home/kweste01/kw/ncbi/dbGaP-60168_WHI_gen_IRB/PhenoGenotypeFiles/ChildStudyConsentSet_phs000746.WHI.v2.p3.c1.HMB-IRB/GenotypeFiles
whi_c2_root=/cluster/home/kweste01/kw/ncbi/dbGaP-60170_WHI_gen_IRBNPU/PhenoGenotypeFiles/ChildStudyConsentSet_phs000746.WHI.v2.p3.c2.HMB-IRB-NPU/GenotypeFiles

#use_dose2plink_whi () {
#	local substudy=$1	
#	local dir_c1=$2
#	local dir_c2=$3
#	for ((i=1; i<=22; i++)); do
#		#gunzip < $dir_c1/machout.chr$i.info.gz > \
#		#	tmp_whi_${substudy}_chr$i.info
#		gunzip < $dir_c1/machout.chr$i.dose_GRU.gz > tmp_whi_${substudy}_chr${i}_c1.dose
#		gunzip < $dir_c2/machout.chr$i.dose_NPU.gz > tmp_whi_${substudy}_chr${i}_c2.dose
#		cat tmp_whi_${substudy}_chr${i}_c1.dose \
#			tmp_whi_${substudy}_chr${i}_c2.dose > tmp_whi_${substudy}_chr$i.dose
#		#awk '{print $1"->"$0}' tmp_whi_${substudy}_chr$i.dose > \
#		#	tmp_whi_${substudy}_chr$i.dosefix
#		~/kw/opt/dose2plink -dose tmp_whi_${substudy}_chr$i.dosefix \
#			-info tmp_whi_${substudy}_chr$i.info \
#			-out $GENODIR/whi_chr$i
#		rm tmp_whi_${substudy}*
#	done
#}	
merge_chromosomes () {
	local group=$1
	#echo "Merging $group"

	#gunzip < $GENODIR/whi_${group}_chr1.pdat.gz > $GENODIR/whi_$group.pdat
	#for ((i=2; i<=22; i++)); do
	#	gunzip < $GENODIR/whi_${group}_chr$i.pdat.gz | tail -n +2 >> $GENODIR/whi_$group.pdat
	#done
	#cat $GENODIR/whi_${group}_chr1.pfam > $GENODIR/whi_${group}.pfam

	#cat $GENODIR/whi_${group}_chr*_low_qual_snps.txt > $GENODIR/whi_${group}_low_qual_snps.txt

	# Read-in, and pfile output using plink2
	plink2 --import-dosage $GENODIR/whi_$group.pdat \
		--psam $GENODIR/whi_$group.pfam \
		--make-pgen \
		--out $GENODIR/whi_${group}_tmp

	## ...while plink2 --alleleACGT is not yet operational
	mv $GENODIR/whi_${group}_tmp.pvar $GENODIR/whi_${group}_tmp.pvartmp
	awk '{gsub(1,"A",$4); gsub(2,"C",$4); gsub(3,"G",$4); gsub(4,"T",$4); \
		gsub(1,"A",$5); gsub(2,"C",$5); gsub(3,"G",$5); gsub(4,"T",$5); print}' \
		< $GENODIR/whi_${group}_tmp.pvartmp \
		> $GENODIR/whi_${group}_tmp.pvar

	# Update rsIDs and reference alleles (from dbSNP) and apply R2 filter
	plink2 --pfile $GENODIR/whi_${group}_tmp \
		--exclude $GENODIR/whi_${group}_low_qual_snps.txt \
		--update-name $SNPANNO 3 6 \
		--ref-allele $SNPANNO 4 3 \
		--make-pgen \
		--out $GENODIR/whi_${group}

	rm -f $GENODIR/whi_${group}_chr*  # delete chromosome-specific files after manual merge
	rm -f $GENODIR/whi_${group}_tmp.*
}


awk -F, '{print $2,$2}' ../int/metaData.csv > ../int/meth_ids.txt
whi_sample_info=${whi_c1_root}/phg000592.v1.WHI_Imputation.sample-info.MULTI/phg000592.v1_release_manifest.txt 
cat $whi_sample_info | tail -n +17 | cut -f 1,2 | awk '{print $1,$1,$2,$2}' | sort -u -k 3 \
	> ../int/whi_sample_to_subject.txt
join -1 3 -2 1 -o 1.1 1.1 \
	../int/whi_sample_to_subject.txt \
	<(cat ../int/meth_ids.txt | sort -k 1,1) \
	> ../int/whi_sample_ids.txt
SNPANNO=../int/snp_annotations/snp_annot_hg19_nodups.txt

# Convert .dose files from mach into plink-usable format
#whi_c1_as264=${whi_c1_root}/phg000592.v1.WHI_AS264.genotype-imputed-data.c1
#whi_c2_as264=${whi_c2_root}/phg000592.v1.WHI_AS264.genotype-imputed-data.c2
#for ((i=1; i<=22; i++)); do
#	gunzip < $whi_c1_as264/AS264chr${i}c.info.gz > tmp_whi_as264_chr$i.info
#	awk '{if ($7 < 0.3) print $1}' tmp_whi_as264_chr$i.info > $GENODIR/whi_as264_chr${i}_low_qual_snps.txt
#	gunzip < $whi_c1_as264/AS264chr${i}c.dose.c1.gz > tmp_whi_as264_chr${i}_c1.dose
#	gunzip < $whi_c2_as264/AS264chr${i}c.dose.c2.gz > tmp_whi_as264_chr${i}_c2.dose
#	cat tmp_whi_as264_chr${i}_c1.dose tmp_whi_as264_chr${i}_c2.dose > tmp_whi_as264_chr$i.dose
#	~/kw/opt/dose2plink -dose tmp_whi_as264_chr$i.dose \
#		-info tmp_whi_as264_chr$i.info \
#		-out $GENODIR/whi_as264_chr$i
#	rm tmp_whi_as264*
#done
#merge_chromosomes as264
#
#whi_c1_garnet=${whi_c1_root}/phg000592.v1.WHI_GARNET.genotype-imputed-data.c1
#whi_c2_garnet=${whi_c2_root}/phg000592.v1.WHI_GARNET.genotype-imputed-data.c2
#for ((i=1; i<=22; i++)); do
#	gunzip < $whi_c1_garnet/GARNETchr${i}c.info.gz > tmp_whi_garnet_chr$i.info
#	awk '{if ($7 < 0.3) print $1}' tmp_whi_garnet_chr$i.info > $GENODIR/whi_garnet_chr${i}_low_qual_snps.txt
#	#gunzip < $whi_c1_garnet/GARNETchr${i}c.dose.c1.gz > tmp_whi_garnet_chr${i}_c1.dose
#	#gunzip < $whi_c2_garnet/GARNETchr${i}c.dose.c2.gz > tmp_whi_garnet_chr${i}_c2.dose
#	#cat tmp_whi_garnet_chr${i}_c1.dose tmp_whi_garnet_chr${i}_c2.dose > tmp_whi_garnet_chr$i.dose
#	#~/kw/opt/dose2plink -dose tmp_whi_garnet_chr$i.dose \
#	#	-info $whi_c1_garnet/GARNETchr${i}c.info.gz \
#	#	-out $GENODIR/whi_garnet_chr$i
#	rm tmp_whi_garnet*
#done
#merge_chromosomes garnet
#
#whi_c1_gecco_cyto=${whi_c1_root}/phg000592.v1.WHI_GECCO_cyto.genotype-imputed-data.c1
#whi_c2_gecco_cyto=${whi_c2_root}/phg000592.v1.WHI_GECCO_cyto.genotype-imputed-data.c2
#for ((i=1; i<=22; i++)); do
#	gunzip < $whi_c1_gecco_cyto/GECCOchr${i}cyto.info.gz > tmp_whi_gecco_cyto_chr$i.info
#	awk '{if ($7 < 0.3) print $1}' tmp_whi_gecco_cyto_chr$i.info > $GENODIR/whi_gecco_cyto_chr${i}_low_qual_snps.txt
#	#gunzip < $whi_c1_gecco_cyto/GECCOchr${i}cyto.dose.c1.gz > tmp_whi_gecco_cyto_chr${i}_c1.dose
#	#gunzip < $whi_c2_gecco_cyto/GECCOchr${i}cyto.dose.c2.gz > tmp_whi_gecco_cyto_chr${i}_c2.dose
#	#cat tmp_whi_gecco_cyto_chr${i}_c1.dose tmp_whi_gecco_cyto_chr${i}_c2.dose > tmp_whi_gecco_cyto_chr$i.dose
#	#~/kw/opt/dose2plink -dose tmp_whi_gecco_cyto_chr$i.dose \
#	#	-info $whi_c1_gecco_cyto/GECCOchr${i}cyto.info.gz \
#	#	-out $GENODIR/whi_gecco_cyto_chr$i
#	rm tmp_whi_gecco_cyto*
#done
#merge_chromosomes gecco_cyto
#
#whi_c1_gecco_init=${whi_c1_root}/phg000592.v1.WHI_GECCO_init.genotype-imputed-data.c1
#whi_c2_gecco_init=${whi_c2_root}/phg000592.v1.WHI_GECCO_init.genotype-imputed-data.c2
#for ((i=1; i<=22; i++)); do
#	gunzip < $whi_c1_gecco_init/GECCOchr${i}init.info.gz > tmp_whi_gecco_init_chr$i.info
#	awk '{if ($7 < 0.3) print $1}' tmp_whi_gecco_init_chr$i.info > $GENODIR/whi_gecco_init_chr${i}_low_qual_snps.txt
#	#gunzip < $whi_c1_gecco_init/GECCOchr${i}init.dose.c1.gz > tmp_whi_gecco_init_chr${i}_c1.dose
#	#gunzip < $whi_c2_gecco_init/GECCOchr${i}init.dose.c2.gz > tmp_whi_gecco_init_chr${i}_c2.dose
#	#cat tmp_whi_gecco_init_chr${i}_c1.dose tmp_whi_gecco_init_chr${i}_c2.dose > tmp_whi_gecco_init_chr$i.dose
#	#~/kw/opt/dose2plink -dose tmp_whi_gecco_init_chr$i.dose \
#	#	-info $whi_c1_gecco_init/GECCOchr${i}init.info.gz \
#	#	-out $GENODIR/whi_gecco_init_chr$i
#	rm tmp_whi_gecco_init*
#done
#merge_chromosomes gecco_init
#
#whi_c1_hipfx=${whi_c1_root}/phg000592.v1.WHI_HIPFX.genotype-imputed-data.c1
#whi_c2_hipfx=${whi_c2_root}/phg000592.v1.WHI_HIPFX.genotype-imputed-data.c2
#for ((i=1; i<=22; i++)); do
#	gunzip < $whi_c1_hipfx/HIPFXchr${i}c.info.gz > tmp_whi_hipfx_chr$i.info
#	awk '{if ($7 < 0.3) print $1}' tmp_whi_hipfx_chr$i.info > $GENODIR/whi_hipfx_chr${i}_low_qual_snps.txt
#	#gunzip < $whi_c1_hipfx/HIPFXchr${i}c.dose.c1.gz > tmp_whi_hipfx_chr${i}_c1.dose
#	#gunzip < $whi_c2_hipfx/HIPFXchr${i}c.dose.c2.gz > tmp_whi_hipfx_chr${i}_c2.dose
#	#cat tmp_whi_hipfx_chr${i}_c1.dose tmp_whi_hipfx_chr${i}_c2.dose > tmp_whi_hipfx_chr$i.dose
#	#~/kw/opt/dose2plink -dose tmp_whi_hipfx_chr$i.dose \
#	#	-info $whi_c1_hipfx/HIPFXchr${i}c.info.gz \
#	#	-out $GENODIR/whi_hipfx_chr$i
#	rm tmp_whi_hipfx*
#done
#merge_chromosomes hipfx
#
#whi_c1_whims=${whi_c1_root}/phg000592.v1.WHIMS.genotype-imputed-data.c1
#whi_c2_whims=${whi_c2_root}/phg000592.v1.WHIMS.genotype-imputed-data.c2
#for ((i=1; i<=22; i++)); do
#	gunzip < $whi_c1_whims/WHIMSchr${i}c.info.gz > tmp_whi_whims_chr$i.info
#	awk '{if ($7 < 0.3) print $1}' tmp_whi_whims_chr$i.info > $GENODIR/whi_whims_chr${i}_low_qual_snps.txt
#	#gunzip < $whi_c1_whims/WHIMSchr${i}c.dose.c1.gz > tmp_whi_whims_chr${i}_c1.dose
#	#gunzip < $whi_c2_whims/WHIMSchr${i}c.dose.c2.gz > tmp_whi_whims_chr${i}_c2.dose
#	#cat tmp_whi_whims_chr${i}_c1.dose tmp_whi_whims_chr${i}_c2.dose > tmp_whi_whims_chr$i.dose
#	#~/kw/opt/dose2plink -dose tmp_whi_whims_chr$i.dose \
#	#	-info $whi_c1_whims/WHIMSchr${i}c.info.gz \
#	#	-out $GENODIR/whi_whims_chr$i
#	rm tmp_whi_whims*
#done
#merge_chromosomes whims
#
#whi_c1_share_aa=${whi_c1_root}/phg000592.v1.WHI_SHARE_aa.genotype-imputed-data.c1
#whi_c2_share_aa=${whi_c2_root}/phg000592.v1.WHI_SHARE_aa.genotype-imputed-data.c2
#for ((i=1; i<=22; i++)); do
#	gunzip < $whi_c1_share_aa/SHAREchr${i}aa.info.gz > tmp_whi_share_aa_chr$i.info
#	awk '{if ($7 < 0.3) print $1}' tmp_whi_share_aa_chr$i.info > $GENODIR/whi_share_aa_chr${i}_low_qual_snps.txt
#	#gunzip < $whi_c1_share_aa/SHAREchr${i}aa.dose.c1.gz > tmp_whi_share_aa_chr${i}_c1.dose
#	#gunzip < $whi_c2_share_aa/SHAREchr${i}aa.dose.c2.gz > tmp_whi_share_aa_chr${i}_c2.dose
#	#cat tmp_whi_share_aa_chr${i}_c1.dose tmp_whi_share_aa_chr${i}_c2.dose > tmp_whi_share_aa_chr$i.dose
#	#~/kw/opt/dose2plink -dose tmp_whi_share_aa_chr$i.dose \
#	#	-info $whi_c1_share_aa/SHAREchr${i}aa.info.gz \
#	#	-out $GENODIR/whi_share_aa_chr$i
#	rm tmp_whi_share_aa*
#done
#merge_chromosomes share_aa
#
#whi_c1_share_ha=${whi_c1_root}/phg000592.v1.WHI_SHARE_ha.genotype-imputed-data.c1
#whi_c2_share_ha=${whi_c2_root}/phg000592.v1.WHI_SHARE_ha.genotype-imputed-data.c2
#for ((i=1; i<=22; i++)); do
#	gunzip < $whi_c1_share_ha/SHAREchr${i}ha.info.gz > tmp_whi_share_ha_chr$i.info
#	awk '{if ($7 < 0.3) print $1}' tmp_whi_share_ha_chr$i.info > $GENODIR/whi_share_ha_chr${i}_low_qual_snps.txt
#	#gunzip < $whi_c1_share_ha/SHAREchr${i}ha.dose.c1.gz > tmp_whi_share_ha_chr${i}_c1.dose
#	#gunzip < $whi_c2_share_ha/SHAREchr${i}ha.dose.c2.gz > tmp_whi_share_ha_chr${i}_c2.dose
#	#cat tmp_whi_share_ha_chr${i}_c1.dose tmp_whi_share_ha_chr${i}_c2.dose > tmp_whi_share_ha_chr$i.dose
#	#~/kw/opt/dose2plink -dose tmp_whi_share_ha_chr$i.dose \
#	#	-info $whi_c1_share_ha/SHAREchr${i}ha.info.gz \
#	#	-out $GENODIR/whi_share_ha_chr$i
#	rm tmp_whi_share_ha*
#done
#merge_chromosomes share_ha

declare -a groups=("as264" "garnet" "gecco_cyto" "gecco_init" "hipfx" "whims" "share_aa" "share_ha")

# Merge chromosomes and create plink2 sets
for group in ${groups[@]}; do
	merge_chromosomes $group
done

# Merge chromosomes and create plink2 sets 
#for group in ${groups[@]}; do
#	plink2 --pfile $GENODIR/whi_$group --make-bed --out $GENODIR/whi_filtered_$group
#done
#ls $GENODIR/whi_filtered_*.bed | sed 's/.bed//g' > $GENODIR/whi_merge_list.txt

#plink --merge-list $GENODIR/whi_merge_list.txt \
#	--make-bed --out $GENODIR/whi
#
#for group in ${groups[@]}; do
#	plink --bfile $GENODIR/whi_filtered_$group \
#		--exclude $GENODIR/whi-merge.missnp \
#		--make-bed --out $GENODIR/whi_filtered_noMA_$group
#done
#
#ls $GENODIR/whi_filtered_noMA_*.bed | sed 's/.bed//g' > $GENODIR/whi_merge_list_noMA.txt
#
#plink --merge-list $GENODIR/whi_merge_list_noMA.txt \
#	--update-name $SNPANNO 3 6 \
#	--alleleACGT \
#	--update-ids ../int/whi_sample_to_subject.txt \
#       	--make-bed --out $GENODIR/whi
#plink2 --bfile $GENODIR/whi --make-pgen --out $GENODIR/whi
