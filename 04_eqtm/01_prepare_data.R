library(dplyr)
library(data.table)

# Set up parameters

output.eqtm.pre <- "~/bio/datasets/eQTM/"

src.pheno.data.pre <- "~/bio/datasets/pheno/"
src.snps.data.pre  <-"~/bio/datasets/snps/"
# src.pheno.data.pre <- "/binder/mgp/datasets/2020_DexStim_Array_Human/pheno/"
pheno.fn           <- paste0(src.pheno.data.pre, "pheno_full_for_kimono.csv")

# Load data

db_gene_cg <- fread("~/bio/datasets/kimono/mapping/mapping_cpg_gene_ensg_full.csv")
db_ilmn_gene <- fread("~/bio/datasets/kimono/mapping/mapping_ilmn_ensg_gene.csv")
db_gene_snp <- fread("~/bio/datasets/kimono/mapping/mapping_snp_gene_distance.csv")

snp.kimono.mtrx    <- fread("~/bio/datasets/kimono/input/snp.csv")
snp.bim            <- fread("~/bio/datasets/snps/Dex_genoData_SNPs.bim")
pheno              <- fread(pheno.fn, na.strings = c('#N/A', NA), dec = ",") %>% setDT()
methyl.mtrx        <- readRDS("~/bio/datasets/methylation/10_final_qc_data/dex_methyl_beta_combat_mtrx.rds")
gex.mtrx.veh       <- fread("~/bio/datasets/kimono/input/gex_veh.csv")
gex.mtrx.dex       <- fread("~/bio/datasets/kimono/input/gex_dex.csv")


# Prepare and save the SNP matrix

snp <- melt(snp.kimono.mtrx, id="RNA_ID")
head(snp)

snp <- left_join(snp, pheno[, c("DNA_ID", "RNA_ID")]) %>% mutate(SNP = variable) %>% select(DNA_ID, SNP, value)
head(snp)

fwrite(snp, 
       paste0(src.snps.data.pre, "snp_for_eQTL.csv"),
       quote = F, row.names = F, sep = ";")


if (dim(table((colnames(snp) == c("DNA_ID", "SNP", "value")))) > 1)
    colnames(snp) <- c("DNA_ID", "SNP", "value")
    
snp.mtrx <- dcast(snp, formula = SNP ~ DNA_ID)

fwrite(snp.mtrx, 
       paste0(src.snps.data.pre, "snp_mtrx.csv"),
       quote = F, row.names = F, sep = ";")


snp.mtrx <- fread( paste0(src.snps.data.pre, "snp_mtrx.csv"))

# Prepare methylation data

# Extract only baseline samples:
veh.ids <- pheno[Dex == 0 & !is.na(DNAm_ID), .(DNA_ID, DNAm_ID)]
dex.ids <- pheno[Dex == 1 & !is.na(DNAm_ID), .(DNA_ID, DNAm_ID)]

methyl.mtrx.veh <- methyl.mtrx[, colnames(methyl.mtrx) %in% veh.ids$DNAm_ID]
methyl.mtrx.dex <- methyl.mtrx[, colnames(methyl.mtrx) %in% dex.ids$DNAm_ID]

# Match colnames from methylation beta mtrx to values in "veh.ids" table
methyl.mtrx.veh           <- methyl.mtrx.veh[,match(veh.ids$DNAm_ID, colnames(methyl.mtrx.veh))]
colnames(methyl.mtrx.veh) <- veh.ids$DNA_ID
methyl.mtrx.veh           <- methyl.mtrx.veh %>% data.frame()
methyl.mtrx.veh["CpG_ID"] <- rownames(methyl.mtrx.veh)
methyl.mtrx.veh           <- methyl.mtrx.veh %>% dplyr::select(CpG_ID, everything())

methyl.mtrx.dex           <- methyl.mtrx.dex[,match(dex.ids$DNAm_ID, colnames(methyl.mtrx.dex))]
colnames(methyl.mtrx.dex) <- dex.ids$DNA_ID
methyl.mtrx.dex           <- methyl.mtrx.dex %>% data.frame()
methyl.mtrx.dex["CpG_ID"] <- rownames(methyl.mtrx.dex)
methyl.mtrx.dex           <- methyl.mtrx.dex %>% dplyr::select(CpG_ID, everything())

fwrite(methyl.mtrx.veh, 
       paste0(output.eqtm.pre, "methyl_beta_mtrx_veh_mad_25.csv"),
       quote = F, row.names = F, sep = ";")

fwrite(methyl.mtrx.dex, 
       paste0(output.eqtm.pre, "methyl_beta_mtrx_dex_mad_25.csv"),
       quote = F, row.names = F, sep = ";")

# Get significant CpGs:
dmps.fn      <- "~/bio/datasets/methylation/20_DMA/02_dmp/dmps_significant_with_beta_stat_bcc_pcs_beta_0_p_90.txt"
dmps.df      <- fread(dmps.fn, sep = "\t")
dmp.sign.ids <- dmps.df$Probe_Id

methyl.mtrx.dex <- fread(paste0(output.eqtm.pre, "methyl_beta_mtrx_dex.csv"))
methyl.mtrx.veh <- fread(paste0(output.eqtm.pre, "methyl_beta_mtrx_veh.csv"))

methyl.mtrx.dex.sign <- methyl.mtrx.dex[V1 %in% dmp.sign.ids]
methyl.mtrx.veh.sign <- methyl.mtrx.veh[CpG_ID %in% dmp.sign.ids]

fwrite(methyl.mtrx.veh.sign, 
       paste0(output.eqtm.pre, "methyl_beta_mtrx_veh_beta_0_p_90.csv"),
       quote = F, row.names = F, sep = ";")

fwrite(methyl.mtrx.dex.sign, 
       paste0(output.eqtm.pre, "methyl_beta_mtrx_dex_beta_0_p_90.csv"),
       quote = F, row.names = F, sep = ";")
    
# Preapare GEX data
veh.ids <- pheno[Dex == 0 & !is.na(DNAm_ID), .(DNA_ID, RNA_ID)]
dex.ids <- pheno[Dex == 1 & !is.na(DNAm_ID), .(DNA_ID, RNA_ID)]

colnames(gex.mtrx.veh)[1] <- "RNA_ID"
colnames(gex.mtrx.dex)[1] <- "RNA_ID"

gex.mtrx.veh <- left_join(gex.mtrx.veh, veh.ids) %>% dplyr::select(-RNA_ID)
gex.mtrx.dex <- left_join(gex.mtrx.dex, dex.ids) %>% dplyr::select(-RNA_ID)

gex.mtrx.veh.t <- t(gex.mtrx.veh) %>% data.frame()
colnames(gex.mtrx.veh.t) <- gex.mtrx.veh$DNA_ID
gex.mtrx.veh.t <- gex.mtrx.veh.t[, match(veh.ids$DNA_ID, colnames(gex.mtrx.veh.t))]
gex.mtrx.veh.t["ENSG_ID"] <- rownames(gex.mtrx.veh.t)
gex.mtrx.veh.t <- gex.mtrx.veh.t %>% dplyr::select(ENSG_ID, everything())

gex.mtrx.dex.t <- t(gex.mtrx.dex) %>% data.frame()
colnames(gex.mtrx.dex.t) <- gex.mtrx.dex$DNA_ID
gex.mtrx.dex.t <- gex.mtrx.dex.t[,match(dex.ids$DNA_ID, colnames(gex.mtrx.dex.t))]
gex.mtrx.dex.t["ENSG_ID"] <- rownames(gex.mtrx.dex.t)
gex.mtrx.dex.t <- gex.mtrx.dex.t %>% dplyr::select(ENSG_ID, everything())

fwrite(gex.mtrx.veh.t, 
       paste0(output.eqtm.pre, "gex_mtrx_veh.csv"),
       quote = F, row.names = F, sep = ";")

fwrite(gex.mtrx.dex.t, 
       paste0(output.eqtm.pre, "gex_mtrx_dex.csv"),
       quote = F, row.names = F, sep = ";")

# Prepare file with SNP coordinates
snp.loc <- snp.bim %>% dplyr::select(V2, V1, V4) %>% unique()
colnames(snp.loc) <- c("SNP", "chr", "pos")

fwrite(snp.loc, 
       paste0(output.eqtm.pre, "snp_locations.csv"),
       quote = F, row.names = F, sep = ";")

# Prepare file with SNP coordinates
snp.loc <- snp.bim %>% dplyr::select(V2, V1, V4) %>% unique()
colnames(snp.loc) <- c("SNP", "chr", "pos")

fwrite(snp.loc, 
       paste0(output.eqtm.pre, "snp_locations.csv"),
       quote = F, row.names = F, sep = ";")

# Prepare file with CpG coordinates

cpg.anno <- fread("~/bio/datasets/methylation/20_DMA/02_dmp/dmp_bcc_pcs_anno.csv")
colnames(cpg.anno)

cpg.loc           <- cpg.anno %>% dplyr::select(Name, chr, pos) %>% unique()
colnames(cpg.loc) <- c("CpG_ID", "chr", "pos")
cpg.loc$chr       <- gsub("[^0-9.]", "",  cpg.loc$chr)

fwrite(cpg.loc, 
       paste0(output.eqtm.pre, "cpg_locations.csv"),
       quote = F, row.names = F, sep = ";")

# Prepare file with ENSG coordinates

ilmn.anno <- fread("~/bio/datasets/gene_expression/v3_v4_sharedContent_QC.txt")
ilmn.loc  <- ilmn.anno %>% dplyr::select(`X.PROBE_ID`, Chr, P_start, P_end) %>% unique()
colnames(ilmn.loc) <- c("Illumina_ID", "chr", "P_start", "P_end")

map.ilmn.ensg <- fread("~/bio/datasets/kimono/mapping/mapping_ilmn_ens.csv")
ensg.loc <- inner_join(ilmn.loc, map.ilmn.ensg) %>% dplyr::select(ENSG_ID = Ensemble_ID, chr, P_start, P_end) %>% unique() 

fwrite(ensg.loc, 
       paste0(output.eqtm.pre, "ensg_locations.csv"),
       quote = F, row.names = F, sep = ";")

# Prepare mapping file with snp ids, cpgs, ilmn ids, ennsg ids

map.cpg.ensg.gene <- fread("~/bio/datasets/kimono/mapping/mapping_cpg_gene_ensg_full.csv")
map.cpg.gene.ensg.ilmn <- inner_join(map.cpg.ensg.gene, map.ilmn.ensg)

fwrite(map.cpg.gene.ensg.ilmn, 
       paste0(output.eqtm.pre, "mapping_cpg_gene_ensg_ilmn_full.csv"),
       quote = F, row.names = F, sep = ";")

fwrite(map.cpg.gene.ensg.ilmn, 
       "~/bio/datasets/kimono/mapping/mapping_cpg_gene_ensg_ilmn_full.csv",
       quote = F, row.names = F, sep = ";")

# Prepare bio data

pheno <- na.omit(pheno)
covariates <- colnames(pheno)

cov.list <- c("DNA_ID",
              "Dex", "Sex", "Status", "Age", "BMI_D1",
              "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

# VEH
bio.mtrx <- pheno[Dex == 0 & !is.na(DNAm_ID), ] %>% dplyr::select(cov.list)

bio.mtrx.t <- dcast(melt(bio.mtrx, id.vars = "DNA_ID"), variable ~ DNA_ID)
colnames(bio.mtrx.t)[1] <-  "Feature"

fwrite(bio.mtrx.t, 
       paste0(output.eqtm.pre, "bio_mtrx_veh.csv"),
       quote = F, row.names = F, sep = ";")

# DEX
bio.mtrx <- pheno[Dex == 1 & !is.na(DNAm_ID), ] %>% dplyr::select(cov.list)

bio.mtrx.t <- dcast(melt(bio.mtrx, id.vars = "DNA_ID"), variable ~ DNA_ID)
colnames(bio.mtrx.t)[1] <-  "Feature"

fwrite(bio.mtrx.t, 
       paste0(output.eqtm.pre, "bio_mtrx_dex.csv"),
       quote = F, row.names = F, sep = ";")

