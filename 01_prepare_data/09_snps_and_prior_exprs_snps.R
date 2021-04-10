
library(dplyr)
library(data.table)

# Set up parameters

snp.src.dir <- "~/bio/datasets/snps/"
snp.src.dir <- "/home/ahryhorzhevska/mpip/datasets/2020_DexStim_Array_Human/snps/"

src.pheno.data.pre <- "~/bio/datasets/pheno/"
src.pheno.data.pre <- "/binder/mgp/datasets/2020_DexStim_Array_Human/pheno/"

pheno.fn        <- paste0(src.pheno.data.pre, "pheno_full_for_kimono.csv")

kimono.data.pre <- "~/bio/datasets/kimono/"
kimono.data.pre <- "/binder/mgp/datasets/2020_DexStim_Array_Human/kimono/"
# Load data

x <- fread(paste0(snp.src.dir, "plinkgen_X.txt"))
y <- fread(paste0(snp.src.dir, "Dex_genoData_SNPs.sample"), header=F)
z <- fread(paste0(snp.src.dir, "SNPposition.txt"), header = F)

pheno <- fread(pheno.fn, na.strings = c('#N/A', ''))

# Substract one sub-matrix from another: first number of each patient - second number 

#  1 indicates homozygote for reference
# -1 indicates heterozygote
#  0 inidicates homozygote for alternative

numsam <- (nrow(y) - 2) * 3
new_x  <- x[, seq(1, numsam, by = 3), with = F] - x[, seq(1, numsam, by = 3) + 1, with = F]
new_x[1:5]

# Code the numbers so that 0: homo-reference; 1: hetero; 2: homo-alternative
new_x[new_x == 1]  <-"ref"
new_x[new_x == -1] <-"hetero"
new_x[new_x == 0]  <-"alt"

new_x[new_x == "ref"]    <- 0
new_x[new_x == "hetero"] <- 1
new_x[new_x == "alt"]    <- 2

z$ID <- rownames(z)
colnames(z)[1] <- "SNP"
z$SNP           <- as.character(z$SNP)

new_x_1 <- cbind(z$ID, new_x)

# Add DNA_ID to the mat
colnames(new_x_1) <- c("sample", y$V1[3:length(y$V1)])


# Keep only columns (samples) that are in the rest of the data

dna.ids   <- unique(pheno[!is.na(pheno$DNAm_ID), "DNA_ID"])
whichcols <- c("sample", colnames(new_x_1)[colnames(new_x_1) %in% dna.ids$DNA_ID])

# extract them
new_x_2 <- new_x_1[,..whichcols, with = F] 

# as integer
changeCols <- whichcols[-1]
new_x_2[,(changeCols):= lapply(.SD, as.integer), .SDcols = changeCols]
str(new_x_2)

new_x_3 <- cbind(z$SNP, new_x_2)

fwrite(new_x_3, 
       paste0(kimono.data.pre, "mapping/snp_mtrx.csv"),
       quote = F, row.names = F, sep = ";")

####################################################
####################################################
####################################################
# this section deals with SNPdata for my benchmarking as opposed to Christoph benchmarking

# SNPdata <- melt(new_x_2, id="sample")
# print("melted data")
# 
# 
# #save
# write.table(SNPdata, file="./SNPdata.txt", quote = F, row.names = F)
# save(SNPdata, file="./SNPdata.rda")
####################################################
####################################################
####################################################

workspacedir="/storage/groups/ccm01/workspace/Benchmark_MONI_MDD/Benchmarking/"

db_gene_snp <- fread(paste0(workspacedir,"Gene_SNP_Distance.txt"))
db_gene_cg <- fread(paste0(workspacedir,"Gene_CGsite_Distance.txt"))
db_ilmn_gene <- fread(paste0(workspacedir,"ILMN_to_GenesID.txt"))

expressiondt <- fread(file=paste0(workspacedir, "inputMONI/expressiondt_moni.csv"))

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
