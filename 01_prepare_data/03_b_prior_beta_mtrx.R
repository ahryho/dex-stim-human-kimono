library(optparse)

# CLI parsing
option_list = list(
  make_option(c("-r", "--directory"),
              action = "store",
              type = "character",
              default = "~/bio/datasets/methylation/20_DMA/",
              help = "a directory where the methylation data are saved",
              metavar = "character"),
  make_option(c("-d", "--data"),
              type = "character",
              default = "02_dmp/dmps_significant_with_beta_stat_bcc_pcs_beta_10_p_5.txt",
              help = "an input data filename, significant CpGs, result of filtered lmer",
              metavar = "character"),
  make_option(c("-b", "--betamtrx"),
              type = "character",
              default = "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data/dex_methyl_beta_combat_mtrx.rds",
              help = "QC methylation beta matrix",
              metavar = "character"),
  make_option(c("-m", "--mapping"),
              type = "character",
              default = "/mapping/mapping_cpg_gene_ensg_full.csv",
              help = "QC methylation beta matrix",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = "~/bio/datasets/kimono/",
              help = "an output directory",
              metavar = "character")
)

opt <- parse_args(OptionParser(option_list=option_list))

library(biomaRt)
library(dplyr)
library(data.table)

# Set up parameters

# src.data.pre       <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/20_DMA/" 
# output.data.pre <-"/binder/mgp/datasets/2020_DexStim_Array_Human/kimono/"

# Load dmps and dmrs

# dmps.fn <- paste0(src.data.pre, "03_dmr/dmps_unique_bcc_pcs_beta_10_p_5.txt")
# dmps.df <- read.csv(dmps.fn, sep = "\t")
# 
# dmrs.fn <- paste0(src.data.pre, "03_dmr/comb-p/dex.anno.hg19.bed")
# dmrs.df <- read.csv(dmrs.fn, sep = "\t") # dmr > 0.05

# Load just dmps
system(paste0("ls -lh ", opt$directory, "02_dmp"))
dmps.fn <- paste0(opt$directory, "02_dmp/dmps_significant_with_beta_stat.txt")
dmps.fn <- paste0(opt$directory, opt$data)
dmps.df <- fread(dmps.fn, sep = "\t")

# Load mapping table
mapping.tbl <- read.csv2(paste0(opt$out, opt$mapping))

# Export only "CpG_ID to ENSG_ID" mapping table for kimono
map.cpg.ensg.tbl <- mapping.tbl[, c("CpG_ID", "Ensemble_ID")]
map.cpg.ensg.tbl <- map.cpg.ensg.tbl[!duplicated(map.cpg.ensg.tbl),]
map.cpg.ensg.tbl <- na.omit(map.cpg.ensg.tbl)

system(paste0("ls -lh ", opt$out))
write.csv2(map.cpg.ensg.tbl, 
           paste0(opt$out, "/input/prior_cpg_ensg.csv"), 
           quote = F, row.names = F)


# Preapare DMPs beta mtrx for kimono

opt$betamtrx <- "~/bio/datasets/methylation/10_final_qc_data/dex_methyl_beta_combat_mtrx.rds"
beta.mtrx    <- readRDS(opt$betamtrx)

cpg.sign.ids <- dmps.df$Probe_Id

beta.kimono.mtrx <- beta.mtrx[rownames(beta.mtrx) %in% cpg.sign.ids, ]
# beta.kimono.mtrx <- beta.mtrx
beta.kimono.mtrx <- as.data.frame(t(beta.kimono.mtrx))

fwrite(beta.kimono.mtrx, 
       paste0(opt$out, "/input/methylation_beta_mtrx.csv"), 
       quote = F, row.names = T)  

