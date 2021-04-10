
library(data.table)
library(dplyr)
library(minfi)
# Set up parameters

src.data.pre    <- "~/bio/datasets/gene_expression/" 
# src.data.pre   <- "/binder/mgp/datasets/2020_DexStim_Array_Human/gene_expression
output.data.pre <- "~/bio/datasets/kimono/"
mapping.tbl.fn  <- paste0(output.data.pre, "/mapping/mapping_ilmn_ens.csv")

# Load gene-expression data

gex.fn <- paste0(src.data.pre, "ESet_dex_sva.rda")
x      <- load(gex.fn)
gex    <- get(x)
rm(x)

# Extract gex table for kimono

gex.tbl  <- as.data.frame(exprs(gex))
# dim(gex.tbl)
probeids <- rownames(gex.tbl)
gex.tbl$Probe_Id <- probeids

# Create kimono format

# check if the mapping file exists in the folder, if no, create a mapping file for ILMN -> ENS IDs
if (file.exists(mapping.tbl.fn))
    mapping.tbl <- read.csv2(paste0(output.data.pre, "/mapping/mapping_ilmn_ens.csv")) else{
  
  library(biomaRt)      
  ensembl     <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
  mapping.tbl <- getBM(attributes = c('illumina_humanht_12_v3', 'ensembl_gene_id', 'external_gene_name'),
                       filters = 'illumina_humanht_12_v3', 
                       values = probeids, mart = ensembl)
  colnames(mapping.tbl) <- c("Illumina_ID", "Ensemble_ID", "Gene_ID")
  rownames(mapping.tbl) <- mapping.tbl$Illumina_ID
  
  write.csv2(mapping.tbl, 
             paste0(output.data.pre, "/mapping/mapping_ilmn_ens.csv"), 
             quote = F, row.names = F)
}

# Transform GEX data
system(paste0("ls -lh ", src.data.pre, "20_DEA/02_sign_gex/"))

gex.sign.df.fn <- paste0(src.data.pre, "20_DEA/02_sign_gex/gex_significant_with_beta_20_p_1.txt")
gex.sign.df    <- read.csv(gex.sign.df.fn, sep = "\t")
gex.sign.df    <- gex.sign.df[order(gex.sign.df$pFDR),]

gex.sign.ids   <- gex.sign.df$Probe_Id
gex.sign.eset.tbl <- gex.tbl[gex.sign.ids,]

# Mapping table for all ILMN ids that are present in the input gex table
map.sign.tbl <- mapping.tbl[mapping.tbl$Illumina_ID %in% gex.sign.ids, ]

# Get ENS ids
map.tbl.tmp <- right_join(gex.sign.df, map.sign.tbl, by = c("Probe_Id" = "Illumina_ID"))[, c("Probe_Id", "Ensemble_ID", "pFDR")] 
map.tbl.tmp <- map.tbl.tmp[!duplicated(map.tbl.tmp$Ensemble_ID), c("Probe_Id", "Ensemble_ID")]
mapping.tbl <- map.tbl.tmp
rm(map.tbl.tmp, map.sign.tbl)

# Prepare input format for kimono
gex.kimono.tbl <- inner_join(gex.tbl, mapping.tbl, by = "Probe_Id") %>%
  dplyr::select(-c("Probe_Id")) # %>% select(Ensemble_ID, everything())

rownames(gex.kimono.tbl) <- gex.kimono.tbl$Ensemble_ID
gex.kimono.tbl <- gex.kimono.tbl %>% dplyr::select(-c("Ensemble_ID"))
gex.kimono.tbl <- as.data.frame(t(gex.kimono.tbl))

system(paste0("ls -lh ", output.data.pre, "input"))
system(paste0("cd ", output.data.pre, "; cp -r input input_01; ls -lh"))

fwrite(gex.kimono.tbl, 
       paste0(output.data.pre, "/input/gex.csv"), 
       quote = F, row.names = T)  

