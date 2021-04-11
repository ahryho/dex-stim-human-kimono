library(dplyr)
library(data.table)

# Set up parameters

snp.src.dir <- "~/bio/datasets/snps/"
# snp.src.dir <- "/home/ahryhorzhevska/mpip/datasets/2020_DexStim_Array_Human/snps/"

src.pheno.data.pre <- "~/bio/datasets/pheno/"
# src.pheno.data.pre <- "/binder/mgp/datasets/2020_DexStim_Array_Human/pheno/"

kimono.data.pre <- "~/bio/datasets/kimono/"
# kimono.data.pre <- "/binder/mgp/datasets/2020_DexStim_Array_Human/kimono/"

pheno.fn             <- paste0(src.pheno.data.pre, "pheno_full_for_kimono.csv")
pheno.dex.kimono.fn  <- paste0(kimono.data.pre, "input/pheno_dex.csv")
gex.kimono.fn        <- paste0(kimono.data.pre, "input/gex.csv")
map.snp.gene.fn      <- paste0(kimono.data.pre, "/mapping/mapping_snp_gene_distance.csv")
snp.mtrx.fn          <- paste0(kimono.data.pre, "mapping/snp_mtrx.csv")

# Load data:  GEX for kimono matrix (gene IDs in the columns)
pheno           <- fread(pheno.fn, na.strings = c('#N/A', ''))
pheno.dex       <- fread(pheno.dex.kimono.fn, na.strings = c('#N/A', ''))
gex             <- fread(gex.kimono.fn)
map.snp.gene.df <- fread(map.snp.gene.fn)
snp.mtrx        <- fread(snp.mtrx.fn)

# Prepare snp mtrx for kimono

# Ensemble IDs in GEX matrix
gex.ensg.ids <- colnames(gex)[-1]

# Distance threshold
dist.trsh       <- 1e4

# Extract SNP matrix for kimono
snps.ids.subset <- map.snp.gene.df[Ensemble_ID %in% gex.ensg.ids,][Distance < dist.trsh][, SNP] %>% unique(.)
snp.submtrx     <- snp.mtrx[sample %in% snps.ids.subset, ]

snp.kimono.tbl  <- as.data.frame(t(snp.submtrx[, -c(1:2)]))
colnames(snp.kimono.tbl) <- snp.submtrx$sample
snp.kimono.tbl$DNA_ID    <- rownames(snp.kimono.tbl)

# Order the samples in the snp tbl based on the order in the pheno data
snp.kimono.tbl <- left_join(snp.kimono.tbl, pheno[Dex == 1, c("DNA_ID", "RNA_ID")], by = "DNA_ID") %>% select(-DNA_ID) %>% setDT
snp.kimono.tbl <- snp.kimono.tbl[order(match(RNA_ID, pheno.dex$RNA_ID))]
snp.kimono.tbl <- snp.kimono.tbl %>% select(RNA_ID, everything())

# Save snp kimono table

fwrite(snp.kimono.tbl, 
       paste0(kimono.data.pre, "/input/snp.csv"), 
       quote = F, row.names = F)  

# Create gex - snp prior

prior.snp.ensg.tbl <- map.snp.gene.df[Ensemble_ID %in% gex.ensg.ids,][Distance < dist.trsh][, .(SNP, Ensemble_ID)] %>% unique()

fwrite(prior.snp.ensg.tbl, 
       paste0(kimono.data.pre, "/input/prior_snp_ensg.csv"), 
       quote = F, row.names = F) 
