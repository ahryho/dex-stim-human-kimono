
# Extract only those Ensemble IDs which are in my GEX matrix

library(data.table)

output.data.pre  <- "~/bio/datasets/kimono/"
gex.kimono.fn    <- paste0(output.data.pre, "input/gex.csv")
prior.biogrid.fn <- paste0(output.data.pre, "input/prior_biogrid_ensg.csv")

# Load GEX matrix

gex.kimono.tbl <- fread(gex.kimono.fn)

# Load full prior for kimono: ENSG_A to ENSG_B

mapping.kimono.tbl <- fread(paste0(output.data.pre, "mapping/prior_biogrid_ensg_full.csv")) 

# Extract required ENSG IDs

ensg.ids <- colnames(gex.kimono.tbl)[-1]
mapping.final.kimono.tbl <- mapping.kimono.tbl[mapping.kimono.tbl$ENSG_A %in% ensg.ids, ]

# Save prior for kimono: ENSG_A to ENSG_B
fwrite(mapping.final.kimono.tbl, 
       prior.biogrid.fn,
       quote = F, row.names = F, sep = ";") 
