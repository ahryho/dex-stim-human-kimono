
src.pheno.data.pre <- "~/bio/datasets/pheno/" # /binder/mgp/datasets/2020_DexStim_Array_Human/pheno
src.data.pre       <- "~/bio/datasets/gene_expression/"

# Mapping table with sample ids

samplesheet.fn <- paste0(src.pheno.data.pre, "mars_dna_gex_meth_ids.csv")
samplesheet    <- read.csv(samplesheet.fn, sep = ";")

# Gene expression data

gex.fn <- paste0(src.data.pre, "ESet_dex_sva.rda")
x   <- load(gex.fn)
gex <- get(x)
rm(x)

# Create common pheno table

gex.pheno  <- pData(gex)
meth.pheno <- read.csv(paste0("~/bio/datasets/methylation/00_sample_sheets/pheno_with_pcs.csv"), sep = ";", header = T)

pheno <- left_join(samplesheet, gex.pheno, by = "RNA_ID", suffix = c("", ".gex"))
pheno <- full_join(pheno, meth.pheno, by = c("DNAm_ID" = "Sample_Name"), suffix = c("", ".meth"))

pheno <- pheno %>% select(- c(Sample_ID.gex, Dex.gex, DNA_ID.gex, Individual, Status.meth, Sex.meth, Real_Age, BMI))

write.csv2(pheno, 
           paste0(src.pheno.data.pre, "pheno_full_for_kimono.csv"), 
           quote = F, row.names = F)






