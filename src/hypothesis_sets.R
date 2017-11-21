library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

anno450k <- data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19), stringsAsFactors=F) %>%
  mutate(cpg=Name)

# Wallner 2016 monoyte-to-macrophage transition DMRs
mon_to_mac <- suppressWarnings(read_csv("../data/literature/mon_to_mac_DMRs_Wallner2016.csv")) %>%
  filter(!is.na(Location)) %>%
  mutate(chr=gsub(":.*", "", Location),
         start=as.integer(gsub("-.*", "", gsub(".*:", "", Location))),
         end=as.integer(gsub(".*-", "", Location))) %>%
  select(chr, start, end) %>%
  inner_join(anno450k, by="chr") %>%
  filter(pos>=start, pos<=end)

# Schroder 2017 variable regions in monocytes from WGBS
mon_var_schroder <- read_csv("../data/literature/variable_monocyte_DMRs_Schroder2017.csv") %>%
  mutate(chr=paste0("chr", gsub(":.*", "", location)),
         start=as.integer(gsub("-.*", "", gsub(".*:", "", location))),
         end=as.integer(gsub(".*-", "", location))) %>%
  select(chr, start, end) %>%
  inner_join(anno450k, by="chr") %>%
  filter(pos>=start, pos<=end)

# Ecker 2017 variable CpGs (more than one sheet, each separately)
library(readxl)
variable_monSpecific <- read_excel("../data/literature/variable_cpgs_Ecker2017.xlsx", sheet="Monocytes")
variable_neutSpecific <- read_excel("../data/literature/variable_cpgs_Ecker2017.xlsx", sheet="Neutrophils")
variable_tcellSpecific <- read_excel("../data/literature/variable_cpgs_Ecker2017.xlsx", sheet="T cells")
variable_monNeutSpecific <- read_excel("../data/literature/variable_cpgs_Ecker2017.xlsx",
                                       sheet="Monocytes + neutrophils")

# CVD GWAS SNP-based
cvd_gwas_associations <- read_tsv("../data/literature/gwas-association-downloaded_2017-09-18-cardiovascular disease.tsv")
cvd_gwas_genes <- unique(unlist(strsplit(cvd_gwas_associations$`REPORTED GENE(S)`, split=", ")))

cvd_gwas_cpgs_byPos <- cvd_gwas_associations %>%
  mutate(chr=paste0("chr", CHR_ID),
         snp_pos=CHR_POS) %>%
  select(chr, snp_pos) %>%
  inner_join(select(anno450k, Name, chr, pos), by="chr") %>%
  mutate(snp_pos=as.numeric(snp_pos)) %>%
  filter(pos>=snp_pos-1000, pos<=snp_pos+1000) %>%
  dplyr::rename(cpg=Name) %>%
  distinct(cpg)

genes_to_cpgs_df <- anno450k %>%
  dplyr::rename(gene=UCSC_RefGene_Name) %>%
  select(cpg, gene) %>%
  separate_rows(gene, sep=";") %>%
  filter(gene!="")
genes_to_cpgs <- split(genes_to_cpgs_df$cpg, genes_to_cpgs_df$gene)

cvd_gwas_cpgs_byGene <- unique(unlist(genes_to_cpgs[cvd_gwas_genes]))

# Chia-PET-based
javierre_cpData <- read_tsv("../data/literature/javierre2016_PCHiC_peak_matrix_cutoff5.tsv") 
system.time(cvdPromoterAssoc_contacts <- javierre_cpData %>%
  filter(Mon>5) %>%
  select(baitName, oeChr, oeStart, oeEnd) %>%
  separate_rows(baitName, sep=";") %>%
  filter(baitName %in% cvd_gwas_genes, !is.na(baitName)) %>%
  mutate(chr=paste0("chr", oeChr)) %>%
  mutate(chunk=rep(1:50, each=ceiling(nrow(.)/50), length.out=nrow(.))) %>%
  nest(-chunk) %>%  # In chunks so size doesn't blow up
  mutate(data=map(data, function(ch) {
    inner_join(ch, select(anno450k, Name, chr, pos), by="chr") %>%
      filter(pos>=oeStart, pos<=oeEnd)
  })) %>% 
  unnest() %>%
  distinct(cpg))

## Final list of CpGs
hypothesis_sets <- list(mon_to_mac=mon_to_mac$cpg,
                        mon_var_schroder=mon_var_schroder$cpg,
                        mon_var_ecker=variable_monSpecific$`Probe ID`,
                        neut_var_ecker=variable_neutSpecific$`Probe ID`,
                        tcell_var_ecker=variable_tcellSpecific$`Probe ID`,
                        monNeut_var_ecker=variable_monNeutSpecific$`Probe ID`,
                        gwas_byPos=cvd_gwas_cpgs_byPos$cpg,
                        gwas_byGene=cvd_gwas_cpgs_byGene,
                        mon_cvdGene_contacts=cvdPromoterAssoc_contacts$cpg)
