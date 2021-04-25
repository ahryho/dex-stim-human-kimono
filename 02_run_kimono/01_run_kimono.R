
###############################################
### 1 libraries
###############################################
library(kimono)
library(data.table)

###############################################
data.dir.pre  <- "/binder/mgp/datasets/2020_DexStim_Array_Human/kimono/"
# data.dir.pre  <- "~/bio/datasets/kimono/"

args      <- commandArgs(T)
treatment <- as.character(args[1]) #"dex"

treatment <- "dex"

kimono.res.fn <- paste0(data.dir.pre, "output/kimono_res_", treatment, ".csv")
###############################################
### 2 Data
###############################################
# read in data

print("read data")

layer.gex    <- fread(paste0(data.dir.pre, "input/gex_", treatment, ".csv"))
layer.methyl <- fread(paste0(data.dir.pre, "input/methylation_beta_mtrx_", treatment, ".csv"))
layer.snp    <- fread(paste0(data.dir.pre, "input/snp_", treatment, ".csv"))
layer.pheno  <- fread(paste0(data.dir.pre, "input/pheno_", treatment, ".csv"), dec = ",")

# read in mapping (prior)
prior.gex       <- fread(paste0(data.dir.pre, "input/prior_biogrid_ensg.csv"))
prior.cpg.gex   <- fread(paste0(data.dir.pre, "input/prior_cpg_ensg.csv"))
prior.gex.pheno <- fread(paste0(data.dir.pre, "input/prior_gene_pheno.csv"))
prior.snp.gex   <- fread(paste0(data.dir.pre, "input/prior_snp_ensg.csv"))

print("data read successfully")

# make sure samples are in the same order
idorder      <- layer.pheno$RNA_ID
layer.gex    <- layer.gex[match(idorder, layer.gex$V1),]
layer.snp    <- layer.snp[match(idorder, layer.snp$RNA_ID),]

idorder      <- layer.pheno$DNAm_ID
layer.methyl <- layer.methyl[match(idorder, layer.methyl$V1),]

# remove sample column
layer.gex    <- layer.gex[, -1]
layer.methyl <- layer.methyl[, -1]
layer.snp    <- layer.snp[, -1]
layer.pheno  <- layer.pheno[, -c(1,2)]

###############################################
# 3 Assemble into lists
###############################################
print("assemble data")

# data list
input_list <- list(
  as.data.table(layer.gex),
  as.data.table(layer.methyl),
  as.data.table(layer.snp),
  as.data.table(layer.pheno)
)
names(input_list) <- c('expr',
                       'methyl',
                       'snp',
                       'pheno')
#########################
# mapping list
mapping_list <- list(
  as.data.table(prior.gex),
  as.data.table(prior.cpg.gex),
  as.data.table(prior.snp.gex),
  as.data.table(prior.gex.pheno)
)
#########################
# meta info
metainfo <- data.frame('ID' = c('prior_expr', 'methyl_expr', 'snp_expr', 'expr_pheno'),
                       'main_to'   =  c(1, 2, 3, 4)
)

print("data assembled")

###############################################
# 4 Run MONI
###############################################
print("run kimono")
start.time <- Sys.time(); 
start.time

kimono.res <- kimono(input_list, mapping_list, metainfo,  main_layer = 1, min_features = 2)#, stab_sel = F)

end.time <- Sys.time()
print("kimono done")
end.time - start.time

###############################################
# 5 Save result
###############################################
print("write results out")
fwrite(kimono.res, 
       kimono.res.fn, 
       quote = F, row.names = F, sep = ";") 

print("result written out")
