library(dplyr)
library(data.table)

# Set up parameters

snp.src.dir <- "~/bio/datasets/snps/"
snp.src.dir <- "/home/ahryhorzhevska/mpip/datasets/2020_DexStim_Array_Human/snps/"

src.pheno.data.pre <- "~/bio/datasets/pheno/"
src.pheno.data.pre <- "/binder/mgp/datasets/2020_DexStim_Array_Human/pheno/"

kimono.data.pre <- "~/bio/datasets/kimono/"
kimono.data.pre <- "/binder/mgp/datasets/2020_DexStim_Array_Human/kimono/"

pheno.fn        <- paste0(src.pheno.data.pre, "pheno_full_for_kimono.csv")
gex.kimono.fn   <- paste0(kimono.data.pre, "input/gex.csv")
map.snp.gene.fn <- paste0(kimono.data.pre, "/mapping/mapping_snp_gene_distance.csv")

# Load data:  GEX for kimono matrix (gene IDs in the columns)
pheno <- fread(pheno.fn, na.strings = c('#N/A', ''))
gex   <- fread(gex.kimono.fn)
map.snp.gene.df <- fread(map.snp.gene.fn)

# Ensemble IDs in GEX matrix
gex.ensg.ids <- colnames(gex)[-1]


exprilmn <- colnames(expressiondt)[-c(1:3)] # which ILMN do I have?
subilmn <- db_ilmn_gene[ILMN %in% exprilmn,] # which ILMN can I map

# these are the SNPs that I will be working with
# only in mapping and in distance below 0.5 Mega bp
snpsafterfilter <- db_gene_snp[ gene %in% subilmn$Gene,][distance< 1e4, ][, SNP] %>%  unique(.)
length(snpsafterfilter)

####################
idorder <- fread(paste0(workspacedir,"idorder.txt"))
colnames(idorder) <- gsub("biological.", "",colnames(idorder))
new_x_2 <- fread(file=paste0(workspacedir, "new_x_2.txt"))

########### filter for SNPs
sum(new_x_2$sample %in% snpsafterfilter)
fil_new_x_2 <- new_x_2[sample %in% snpsafterfilter,] # take on SNPs that passed filter
dim(fil_new_x_2)

mynames <- fil_new_x_2$sample

# transpose all but the first column (name)
t_new_x_2 <- as.data.frame(t(fil_new_x_2[,-1])) # transpose table
colnames(t_new_x_2) <- mynames
t_new_x_2$DNA_ID <- rownames(t_new_x_2)

# snpsdt <- merge(idorder,t_new_x_2 , by="DNA_ID") # merge DNA and RNA_ID
# snpsdt <- snpsdt[match(idorder$RNA_ID, snpsdt$RNA_ID),]
snpsdt <- t_new_x_2[match(unique(idorder$DNA_ID), t_new_x_2$DNA_ID),]
snpsdt <- setDT(snpsdt, keep.rownames = TRUE)[] %>%  rename(DNA_ID=rn)
ind <-(which(colnames(snpsdt)=="DNA_ID")[2])
snpsdt <- snpsdt[,!ind, with=F]



fwrite(snpsdt, file=paste0(workspacedir, "snpsdt.txt"))


dim(snpsdt);snpsdt[1:5,1:5]

db_gene_snp_fil <- db_gene_snp[(SNP %in% snpsafterfilter) & (distance < 1e4),]
db_gene_snp_fil2 <-merge(db_gene_snp_fil, db_ilmn_gene, by.x="gene", by.y="Gene", all=F,allow.cartesian=TRUE)
prior_mrna_SNP <- db_gene_snp_fil2[,.(ILMN, SNP)]
fwrite(prior_mrna_SNP, file=paste0(workspacedir, "inputMONI/prior_mrna_SNP.csv"))
