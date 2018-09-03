CPV_DIR=~/Downloads/combined-pvalues-master/cpv
STEP_SIZE=50
P_COL=4
DATASET=$1
RAW_BED=${DATASET}Res.bed

cd ../int

python2.7 $CPV_DIR/acf.py -d 1:1000:$STEP_SIZE -c $P_COL $RAW_BED > acf_$DATASET.txt
python2.7 $CPV_DIR/slk.py --acf acf_${DATASET}.txt -c $P_COL ${DATASET}Res.bed > acf_${DATASET}.bed
python2.7 $CPV_DIR/peaks.py --dist 500 --seed 0.1 acf_${DATASET}.bed > regions_${DATASET}.bed
python2.7 $CPV_DIR/region_p.py -p $RAW_BED -r regions_${DATASET}.bed -s $STEP_SIZE -c $P_COL > regions.sig_${DATASET}.bed

cd ../src
