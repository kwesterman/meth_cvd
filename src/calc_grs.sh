#!/bin/bash

module load python/3.6.0
module load plink2

GENODIR=../int/plinksets

python << EOF
import pandas as pd
grs_weights = (pd.read_csv("../data/literature/CoronaryArteryDisease_PRS_LDpred_rho0.001_v3.txt", 
			   delim_whitespace=True, skiprows=15)
	       .loc[:, ["chr", "position_hg19", "effect_allele", "effect_weight"]])
snp_annot = (pd.read_csv("../int/snp_annotations/snp_annot_hg19_nodups.txt",
			 delim_whitespace=True, header=None, names=["chr", "position_hg19", "id"],
			 usecols=[0, 1, 2])
	     .query('chr not in ["X", "Y"]'))
snp_annot['chr'] = snp_annot['chr'].astype('int64')
grs_weights = grs_weights.merge(snp_annot, on=["chr", "position_hg19"])	
grs_weights = grs_weights[["id", "effect_allele", "effect_weight"]]
grs_weights.to_csv("../int/grs_weights.txt", sep="\t", header=False, index=False)
EOF

# WHI
plink2 --pfile $GENODIR/whi \
       --score ../int/grs_weights.txt 1 2 3 \
       --out ../int/whi_grs

# FHS
plink2 --pfile $GENODIR/fhs \
       --score ../int/grs_weights.txt 1 2 3 \
       --out ../int/fhs_grs
