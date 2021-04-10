library(optparse)
library(plyr)
library(biomaRt)

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
              default = "02_dmp/dmp_bcc_pcs_anno.csv",
              help = "an input data filename, result of lmer + annotation",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = "~/bio/datasets/kimono/",
              help = "an ourput directory",
              metavar = "character")
)

opt <- parse_args(OptionParser(option_list=option_list))

# DMPS annotated

dmps.anno.fn <- paste0(opt$directory, opt$data)
dmps.anno.df <- fread(dmps.anno.fn, sep = "\t")

# Approach 1

# Extract Gene symbol and map to ENSG

# map.tbl <- dmps.anno.df[, c("Name", "GencodeCompV12_Accession", "Pval_adj")]
map.tbl <- dmps.anno.df[, c("Name", "UCSC_RefGene_Name", "Pval_adj")]
nr.na <- nrow(map.tbl[map.tbl$GencodeCompV12_Accession == "", ])

# Create CpG_ID - Gene_ID mapping table

gene.ids.list <- strsplit(map.tbl$UCSC_RefGene_Name, ";")
names(gene.ids.list) <- map.tbl$Name
map.cpg.gene.tbl <- ldply(gene.ids.list, data.frame)
map.cpg.gene.tbl <- map.cpg.gene.tbl[!duplicated(map.cpg.gene.tbl),] 
colnames(map.cpg.gene.tbl) <- c("CpG_ID", "Gene_ID")

# Create Gene_ID - ENSG_ID mapping table

# gene.ids <- unique(unlist(strsplit(map.tbl$UCSC_RefGene_Name, ";")))
gene.ids <- unique(map.cpg.gene.tbl$Gene_ID)

ensembl     <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
mapping.tbl <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id'),
                     filters = 'external_gene_name', 
                     values = gene.ids, mart = ensembl)

colnames(mapping.tbl) <- c("Gene_ID", "Ensemble_ID")

# Join "CpG_ID to Gene_ID" and "Gene_ID to ENSG_ID"

map.cpg.gene.ensg.tbl <- full_join(map.cpg.gene.tbl, mapping.tbl)

write.csv2(map.cpg.gene.ensg.tbl, 
           paste0(opt$out, "/mapping/mapping_cpg_gene_ensg.csv"), 
           quote = F, row.names = F)