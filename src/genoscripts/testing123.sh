module load glibc


FHS_INFO=/cluster/kappa/90-days-archive/ordovaslab/ncbi/dbGaP-60163_FHS_gen_IRBNPUMDS/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v17.p10.c2.HMB-IRB-NPU-MDS/GenotypeFiles/phg000679.v1.FHS_SHARe_Imputed_1000G.genotype-imputed-data.c2/imputed-metrics/machout.chr18.info.gz
FHS_DOSE=/cluster/kappa/90-days-archive/ordovaslab/ncbi/dbGaP-60163_FHS_gen_IRBNPUMDS/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v17.p10.c2.HMB-IRB-NPU-MDS/GenotypeFiles/phg000679.v1.FHS_SHARe_Imputed_1000G.genotype-imputed-data.c2/machout.chr18.dose_NPU.gz

#awk 'BEGIN {FS="\t"; OFS="\t"; print $1,$2,$3,$4,$5,$6,$7}' $FHS_INFO > temp.txt

#~/kw/opt/fcgene-1.0.7/fcgene --mach-mlprob $FHS_DOSE --mach-mlinfo temp.txt --rsq 0.3 --oformat plink-dosage

module load gcta

head -100 /cluster/kappa/90-days-archive/ordovaslab/ncbi/dbGaP-60163_FHS_gen_IRBNPUMDS/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v17.p10.c2.HMB-IRB-NPU-MDS/GenotypeFiles/phg000679.v1.FHS_SHARe_Imputed_1000G.genotype-imputed-data.c2/test/machout.chr18.dose_NPU | awk 'BEGIN {FS="\t"} {print $1}' > myids.txt

gcta64 --dosage-mach-gz $FHS_DOSE $FHS_INFO --imput-rsq 0.3 --keep myids.txt --thread-num 4 --make-bed --out mytest
