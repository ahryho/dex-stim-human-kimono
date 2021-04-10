# Map Gene IDs to phenotype variable

library(data.table)
library(tidyr)

output.data.pre    <- "~/bio/datasets/kimono/"
src.pheno.data.pre <- "~/bio/datasets/pheno/" 

pheno.fn        <- paste0(src.pheno.data.pre, "pheno_full_for_kimono.csv")
gex.kimono.fn   <- paste0(output.data.pre, "input/gex.csv")
output.map.fn   <- paste0(output.data.pre, "input/prior_gene_pheno.csv")
output.pheno.fn <- paste0(output.data.pre, "input/pheno.csv")

pheno.veh.kimono.fn  <- paste0(output.data.pre, "input/pheno_veh.csv")
pheno.dex.kimono.fn  <- paste0(output.data.pre, "input/pheno_dex.csv")

# Load pheno data and  GEX for kimono matrix (gene IDs in the columns)
pheno <- fread(pheno.fn, na.strings = c('#N/A', ''))
gex   <- fread(gex.kimono.fn)

# Extract all gene IDs from GEX matrix and required features from pheno table
gene.ids   <- colnames(gex)[-1]
covariates <- colnames(pheno)

cov.kimono <- c("Dex", "Sex", "Status", "Age", 
                # "Batch",
                "V1", "V2", "V3",
                "PC1", "PC2")

# Crerate pheno table as an input layer for kimono
pheno.kimono <- data.frame(pheno)[, c("DNAm_ID", "RNA_ID", cov.kimono) ]

pheno.kimono <- pheno.kimono[!is.na(pheno.kimono$DNAm_ID),]

# Save pheno tbl
fwrite(pheno.kimono, 
       output.pheno.fn, 
       quote = F, row.names = F, sep = ";") 

# Create all pairs of pheno and gene
mapping.tbl <- crossing("Gene" = gene.ids, "Cov" = cov.kimono)

# Save final mapping
fwrite(mapping.tbl, 
       output.map.fn, 
       quote = F, row.names = F, sep = ";") 

# Separate dex from baseline

pheno.veh.kimono.tbl <- pheno.kimono %>% filter(Dex == 0)
pheno.dex.kimono.tbl <- pheno.kimono %>% filter(Dex == 1)

fwrite(pheno.veh.kimono.tbl, 
       pheno.veh.kimono.fn, 
       quote = F, row.names = F, sep = ";")  

fwrite(pheno.dex.kimono.tbl, 
       pheno.dex.kimono.fn, 
       quote = F, row.names = F, sep = ";") 



