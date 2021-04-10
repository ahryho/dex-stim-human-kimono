# Prepare BioGrid

library(biomaRt)
library(dplyr)
library(data.table)

output.data.pre <- "~/bio/datasets/kimono/"
src.data.pre    <- "~/bio/datasets/biogrid/"
species.fn      <- "BIOGRID-ORGANISM-4.3.195.mitab/BIOGRID-ORGANISM-Homo_sapiens-4.3.195.mitab.txt"

biogrid <- read.csv(paste0(src.data.pre, species.fn), sep = '\t')

# Extract etrez IDs for A and B

entrez.gene.a.ids <- sub('.*:\\s*', '', biogrid$X.ID.Interactor.A) # extract only ids
length(unique(entrez.gene.a.ids)) # 17893

entrez.gene.b.ids <- sub('.*:\\s*', '', biogrid$ID.Interactor.B)
length(unique(entrez.gene.b.ids)) # 23906

entrez.genes.ab.tbl           <- as.data.frame(cbind(entrez.gene.a.ids, entrez.gene.b.ids))
colnames(entrez.genes.ab.tbl) <- c("Entrez_A", "Entrez_B")
entrez.genes.ab.tbl$Entrez_A  <- as.numeric(entrez.genes.ab.tbl$Entrez_A)
entrez.genes.ab.tbl$Entrez_B  <- as.numeric(entrez.genes.ab.tbl$Entrez_B)
entrez.genes.ab.tbl           <- entrez.genes.ab.tbl[!duplicated(entrez.genes.ab.tbl), ]

entrez.ids <- unique(c(entrez.gene.a.ids, entrez.gene.b.ids))
# x <- c("nsp11", "ORF9b")

# Map unique Entrez IDs to Ensemble IDs using biomaRt

ensembl     <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 

biomart.tbl <- getBM(attributes = c('entrezgene_id', 'ensembl_gene_id'),
                           filters = 'entrezgene_id', 
                           values = unique(entrez.ids), mart = ensembl) # 20812

colnames(biomart.tbl) <- c("ENTREZ_ID", "ENSG_ID")

# Map Entrez AB to Enseml AB and export for kimono as a prior

map.tmp.tbl <- left_join(entrez.genes.ab.tbl, biomart.tbl, by = c("Entrez_A" = "ENTREZ_ID")) %>%
  mutate(ENSG_A = ENTREZ_ID) %>%
  dplyr::select(-ENTREZ_ID)

map.tmp.tbl <- left_join(map.tmp.tbl, biomart.tbl, by = c("Entrez_B" = "ENTREZ_ID")) %>%
  mutate(ENSG_B = ENTREZ_ID) %>%
  dplyr::select(-ENTREZ_ID)

mapping.kimono.tbl <- na.omit(map.tmp.tbl) %>% 
                        dplyr::select(ENSG_A, ENSG_B)
mapping.kimono.tbl <- mapping.kimono.tbl[!duplicated(mapping.kimono.tbl),]

# Save results

# Full list of Entrez ID to Ensemble ID
fwrite(biomart.tbl, 
       paste0(output.data.pre, "mapping/mapping_entrez_ensg.csv"), 
       quote = F, row.names = F, sep = ";")  

# Full prior for kimono: ENSG_A to ENSG_B
fwrite(mapping.kimono.tbl, 
       paste0(output.data.pre, "mapping/prior_biogrid_ensg_full.csv"), 
       quote = F, row.names = F, sep = ";") 
