# Check what id from BioGrid to use for mapping to Ensembl ID

library(biomaRt)

src.data.pre <- "~/bio/datasets/biogrid/BIOGRID-ORGANISM-4.3.195.mitab/"
species.fn   <- "BIOGRID-ORGANISM-Homo_sapiens-4.3.195.mitab.txt"

biogrid <- read.csv(paste0(src.data.pre, species.fn), sep = '\t')
colnames(biogrid)
biogrid[1, 1]

# Extract etrez IDs for A and B

entrez.gene.a.ids <- sub('.*:\\s*', '', biogrid$X.ID.Interactor.A) # extract only ids
length(unique(entrez.gene.a.ids)) # 17893

biogrid.a.ids <- gsub('.*biogrid:(.+)\\|entrez.*', '\\1', biogrid$Alt.IDs.Interactor.A)
length(unique(biogrid.a.ids)) # 17936

uniprot.a.ids <- sub('.*swiss-prot:(.+?)\\|.*', '\\1', biogrid$Alt.IDs.Interactor.A)# extract only ids
length(unique(uniprot.a.ids)) # 17829

# Look for pattern

en2sembl     <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
searchAttributes(mart = ensembl, pattern = "entrez") # entrezgene_id
searchAttributes(mart = ensembl, pattern = "biogrid") # biogrid
searchAttributes(mart = ensembl, pattern = "uniprot") # uniprot_gn_id

# Map Entrez ID to Ensemble IDs using biomaRt

mapping.entrez.a.tbl <- getBM(attributes = c('entrezgene_id', 'ensembl_gene_id'),
                              filters = 'entrezgene_id', 
                              values = unique(entrez.gene.a.ids), mart = ensembl) # 16909

# Map bioGrid ID to Ensemble IDs using biomaRt

mapping.biogrid.a.tbl <- getBM(attributes = c('biogrid', 'ensembl_gene_id'),
                               filters = 'biogrid', 
                               values = unique(biogrid.a.ids), mart = ensembl) # 10407

# Map Uniprot ID to Ensemble IDs using biomaRt

mapping.uiprot.a.tbl <- getBM(attributes = c('uniprot_gn_id', 'ensembl_gene_id'),
                              filters = 'uniprot_gn_id', 
                              values = unique(uniprot.a.ids), mart = ensembl) # 15836
