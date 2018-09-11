# Use plink and vcftools to clean and properly format genotype data

GEN_PREFIX_C1=$1
GEN_PREFIX_C2=$2
GEN_SETNAME=$3

~/kw/multiomics/src/add_fam_phenotypes.sh $GEN_PREFIX_C1
~/kw/multiomics/src/add_fam_phenotypes.sh $GEN_PREFIX_C2

module load plink
plink --bfile $GEN_PREFIX_C1 --bmerge $GEN_PREFIX_C2 --geno 0.05 --maf 0.01 --hwe 0.0001 --make-bed --out ~/kw/multiomics/int/intFile
plink --bfile ~/kw/multiomics/int/intFile --mind 0.05 --recode A --out ~/kw/multiomics/int/${GEN_SETNAME}
# plink --bfile $GEN_PREFIX_C1 --bmerge $GEN_PREFIX_C2 --snps-only --geno 0.05 --mind 0.05 --maf 0.01 --hwe 0.0001 --recode A --out ~/kw/multiomics/int/${GEN_SETNAME}

#PHENOFILE=~/kw/multiomics/int/metaData.csv
#PHENOSORT=$(awk -F ',' '{if ($25=="TRUE") {print $2,1} else {print $2,0}}' $PHENOFILE | sort -k 1)
#
#GEN_PREFIX_C1=$1
#FAMFILE_C1=${GEN_PREFIX}.fam
#FAMSORT_C1=$(sort -k 2 $FAMFILE_C1)
#JOINFILE_C1=$(join -1 2 -2 1 -a 1 -e -9 -o 1.1,1.2,1.3,1.4,1.5,2.2 <(echo "$FAMSORT_C1") <(echo "$PHENOSORT"))
#cp $FAMFILE_C1 ${GEN_PREFIX_C1}_orig.fam
#echo "$JOINFILE" > ${GEN_PREFIX}.fam
#
#GEN_PREFIX_C2=$1
#FAMFILE_C2=${GEN_PREFIX}.fam
#FAMSORT_C2=$(sort -k 2 $FAMFILE_C2)
#JOINFILE_C2=$(join -1 2 -2 1 -a 1 -e -9 -o 1.1,1.2,1.3,1.4,1.5,2.2 <(echo "$FAMSORT_C2") <(echo "$PHENOSORT"))
#cp $FAMFILE_C2 ${GEN_PREFIX_C2}_orig.fam
#echo "$JOINFILE" > ${GEN_PREFIX}.fam
