# SNP

library(dplyr)
library(data.table)
library(biomaRt)

# Set up parameters

snp.src.dir     <- "~/bio/datasets/snps/"
output.data.pre <- "~/bio/datasets/kimono/"
snp.annovar.fn  <- paste0(snp.src.dir, "Dex_genoData_SNPs_forAnovar.txt.variant_function")

# Load data

snp.db <- fread(snp.annovar.fn, header=T) %>% setDT
str(snp.db)
head(snp.db)

# Extract nearbyGenes

nearby.gene        <- strsplit(as.character(snp.db$nearByGen), ",")
nearby.gene        <- rapply(nearby.gene, function(x) gsub("dist=", "", x), how = "replace") 
names(nearby.gene) <- snp.db$SNP

# Create mapping table SNP to Gene(distance)
map.snp.tbl <- as.data.table(data.frame("SNP"= rep(as.character(snp.db$SNP), times = lapply(nearby.gene, length)),
                                    "Gene_with_dist" = unlist(nearby.gene, use.names = FALSE)))

map.snp.tbl[, c("Gene_ID", "Distance") := tstrsplit(Gene_with_dist, "(", fixed = TRUE)] %>% 
  .[, Distance := gsub(")", "", Distance)] %>% .[, Gene_with_dist := NULL]
colnames(map.snp.tbl) <- c("SNP", "Gene_ID", "Distance")


# Set distances and replace zero distances with NA
map.snp.tbl$Distance <- as.integer(map.snp.tbl$Distance)
map.snp.tbl[Gene_ID =="NONE", Gene_ID := NA]

# Merge map.snp.tbl with original SNP db
map.snp.tbl <- inner_join(map.snp.tbl, snp.db[, .(SNP, Location, CHR, BP1)], by = "SNP")
map.snp.tbl <- na.omit(map.snp.tbl)

# These are positions where distance is quite confusing - will omit them for this
snp.db[SNP %in% map.snp.tbl[which(is.na(map.snp.tbl$distance)),]$SNP,]

# Check structure of data base. Distance, CHR and BP1 should be integer, rest variables - character 

str(map.snp.tbl)

# Set colnames
colnames(map.snp.tbl) <- c("SNP", "Gene_ID", "Distance", "Location", "Chromosome", "BP1")
map.snp.tbl           <- unique(map.snp.tbl)

# Create SNP - Gene_ID - ENSG_ID mapping table
map.snp.tbl <- fread(paste0(output.data.pre, "/mapping/mapping_snp_gene_distance.csv"))
# gene.ids <- unique(unlist(strsplit(map.snp.tbl$UCSC_RefGene_Name, ";")))
gene.ids <- unique(map.snp.tbl$Gene_ID)

ensembl     <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
mapping.tbl <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id'),
                     filters = 'external_gene_name', 
                     values = gene.ids, mart = ensembl)

colnames(mapping.tbl) <- c("Gene_ID", "Ensemble_ID")

# Join "CpG_ID to Gene_ID" and "Gene_ID to ENSG_ID"

map.snp.gene.ensg.tbl <- full_join(map.snp.tbl, mapping.tbl)
map.snp.gene.ensg.tbl <- na.omit(map.snp.gene.ensg.tbl)

write.csv2(map.snp.gene.ensg.tbl, 
           paste0(output.data.pre, "/mapping/mapping_snp_gene_distance.csv"), 
           quote = F, row.names = F)
