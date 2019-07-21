#!/bin/bash

cd ../int/snp_annotations

rsync -avP rsync://ftp.ensembl.org/ensembl/pub/grch37/release-93/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz .
cat 1000GENOMES-phase_3.vcf.gz | gunzip | tail -n +43 | cut -f 1-5 > snp_annot_hg19.txt
paste snp_annot_hg19.txt <(awk 'BEGIN {OFS=":"} {print $1,$2}' snp_annot_hg19.txt) > snp_annot_hg19_withLoci.txt
cat snp_annot_hg19_withLoci.txt | sort -u -k3,3 | sort -u -k6,6 -V > snp_annot_hg19_nodups.txt  # Unique rsID and chr:pos locus
rm snp_annot_hg19.txt snp_annot_hg19_withLoci.txt
