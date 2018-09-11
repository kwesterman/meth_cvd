GEN_PREFIX=$1  # bed/bim/fam file prefix


if [ ! -f ${GEN_PREFIX}_orig.fam ]; then  # store original .fam file as backup
        cp ${GEN_PREFIX}.fam $FAMFILE ${GEN_PREFIX}_orig.fam
fi

FAMFILE=${GEN_PREFIX}_orig.fam
FAMSORT=$(sort -k 2b,2 $FAMFILE)  # sort .fam file in preparation for join

PHENOFILE=~/kw/multiomics/int/metaData.csv
PHENOSORT=$(awk -F ',' '{if ($25=="TRUE") {print $2,2} else {print $2,1}}' $PHENOFILE | sort -k 1b,1)  # extract shareids and case/control status (coded as 1/2) from metadata file

if grep -q GARNET <<< $GEN_PREFIX; then
	ANNOFILE=$(dirname ${GEN_PREFIX}.fam)/../../phg000139.v1.GARNET_WHI.sample-info.MULTI/phg000139.v1_genotype_released_files_manifest.txt
	ANNOSORT=$(cat $ANNOFILE | tail -n +17 | cut -f2,16 | sort -k 1b,1 | uniq)
	PHENOSORT=$(join -1 2 -2 1 -o 1.1,2.2 <(echo "$ANNOSORT" | sort -k 2b,2) <(echo "$PHENOSORT" | sort -k 1b,1) | sort -k 1b,1)
	# FAMSORT=$(join -1 2 -2 1 -a 1 -e -1 -o 1.1,2.2,1.3,1.4,1.5,1.6 <(echo "$FAMSORT") <(echo "$ANNOSORT") | sort -k 2b,2)
fi

JOINFILE=$(join -1 2 -2 1 -a 1 -e -9 -o 1.1,1.2,1.3,1.4,1.5,2.2 <(echo "$FAMSORT") <(echo "$PHENOSORT"))  # join on shareids (within-subject ID from .fam file) to create updated .fam w/ phenotypes
echo "$JOINFILE" > ${GEN_PREFIX}.fam
