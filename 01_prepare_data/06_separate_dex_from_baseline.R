# Separate dex from baseline

# Load libraries

library(data.table)
library(dplyr)

# Set parameters

output.data.pre <- "~/bio/datasets/kimono/"

beta.kimono.fn   <- paste0(output.data.pre, "input/methylation_beta_mtrx.csv")
gex.kimono.fn    <- paste0(output.data.pre, "input/gex.csv")
pheno.fn         <- paste0(output.data.pre, "input/pheno.csv")
pheno.full.fn    <- paste0("~/bio/datasets/pheno/pheno_full_for_kimono.csv")  

methyl.veh.kimono.fn <- paste0(output.data.pre, "input/methylation_beta_mtrx_veh.csv")
methyl.dex.kimono.fn <- paste0(output.data.pre, "input/methylation_beta_mtrx_dex.csv")

gex.veh.kimono.fn    <- paste0(output.data.pre, "input/gex_veh.csv")
gex.dex.kimono.fn    <- paste0(output.data.pre, "input/gex_dex.csv")

pheno.veh.kimono.fn  <- paste0(output.data.pre, "input/pheno_veh.csv")
pheno.dex.kimono.fn  <- paste0(output.data.pre, "input/pheno_dex.csv")

# load layers
beta.kimono.mtrx <- fread(beta.kimono.fn)
gex.kimono.tbl   <- fread(gex.kimono.fn)
pheno            <- fread(pheno.fn, na.strings = c('#N/A', ''))

# Load full pheno tbl
pheno.full <- fread(pheno.full.fn)

# Check the number of samples in each table are the same

# methyl.sample.ids <- as.data.frame(rownames(beta.kimono.mtrx))
# colnames(methyl.sample.ids) <- "V1"
# beta.kimono.mtrx["V1"] <- rownames(beta.kimono.mtrx)
# beta.kimono.mtrx <- beta.kimono.mtrx %>% dplyr::select(V1, everything())

methyl.sample.ids <- beta.kimono.mtrx[, 1]
gex.sample.ids    <- gex.kimono.tbl[, 1]

pheno.full <- pheno.full[pheno.full$DNAm_ID %in% methyl.sample.ids$V1, ]
pheno.full <- pheno.full[pheno.full$RNA_ID %in% gex.sample.ids$V1, ]

beta.kimono.mtrx <- beta.kimono.mtrx[methyl.sample.ids$V1 %in% pheno.full$DNAm_ID]
gex.kimono.tbl   <- gex.kimono.tbl[gex.sample.ids$V1 %in% pheno.full$RNA_ID, ]
pheno            <- pheno[pheno$DNAm_ID %in% methyl.sample.ids$V1, ] 

# Save
fwrite(beta.kimono.mtrx, 
       beta.kimono.fn, 
       quote = F, row.names = F) 

fwrite(gex.kimono.tbl, 
       gex.kimono.fn, 
       quote = F, row.names = F)  

fwrite(pheno, 
       pheno.fn, 
       quote = F, row.names = F, sep = ";") 

# Separate dex from baseline

# Methylation beta mtrx
veh.methyl.ids <- pheno.full %>% filter(Dex == 0) %>% dplyr::select(DNAm_ID)
dex.methyl.ids <- pheno.full %>% filter(Dex == 1) %>% dplyr::select(DNAm_ID)

methyl.veh.kimono.tbl <- beta.kimono.mtrx[methyl.sample.ids$V1 %in% veh.methyl.ids$DNAm_ID,]
methyl.dex.kimono.tbl <- beta.kimono.mtrx[methyl.sample.ids$V1 %in% dex.methyl.ids$DNAm_ID,]

fwrite(methyl.veh.kimono.tbl, 
       methyl.veh.kimono.fn, 
       quote = F, row.names = F, sep = ";")  

fwrite(methyl.dex.kimono.tbl, 
       methyl.dex.kimono.fn, 
       quote = F, row.names = F, sep = ";")  

# Expression data
veh.gex.ids <- pheno.full %>% filter(Dex == 0) %>% dplyr::select(RNA_ID)
dex.gex.ids <- pheno.full %>% filter(Dex == 1) %>% dplyr::select(RNA_ID)

gex.veh.kimono.tbl <- gex.kimono.tbl[gex.kimono.tbl$V1 %in% veh.gex.ids$RNA_ID,]
gex.dex.kimono.tbl <- gex.kimono.tbl[gex.kimono.tbl$V1 %in% dex.gex.ids$RNA_ID,]

fwrite(gex.veh.kimono.tbl, 
       gex.veh.kimono.fn, 
       quote = F, row.names = F, sep = ";")  

fwrite(gex.dex.kimono.tbl, 
       gex.dex.kimono.fn, 
       quote = F, row.names = F, sep = ";")  

# Pheno data
pheno.veh.kimono.tbl <- pheno[pheno$Dex == 0, ]
pheno.dex.kimono.tbl <- pheno[pheno$Dex == 1, ]

fwrite(pheno.veh.kimono.tbl, 
       pheno.veh.kimono.fn, 
       quote = F, row.names = F, sep = ";")  

fwrite(pheno.dex.kimono.tbl, 
       pheno.dex.kimono.fn, 
       quote = F, row.names = F, sep = ";") 
