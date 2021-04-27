
###############################################
# data.dir.pre  <- "/binder/mgp/datasets/2020_DexStim_Array_Human/kimono/"
data.dir.pre  <- "~/bio/datasets/kimono/"
code.dir.pre <- "~/bio/code/mpip/dex-stim-human-kimono/02_run_kimono/kimono/R/"

args      <- commandArgs(T)
treatment <- as.character(args[1]) #"dex"
startnode <- as.numeric(args[2])

# treatment <- "dex"

kimono.res.fn <- paste0(data.dir.pre, "output/kimono_res_", treatment, "_chunk_", startnode, ".csv")

###############################################
### 1 libraries
###############################################

libraries <- c("data.table", "tidyr", "ggplot2", "reshape2", "oem", "ramify", "stringr",
               "magrittr", "foreach", "doParallel", "dplyr", "furrr", "purrr")
lapply(libraries, require, character.only = TRUE)

source(paste0(code.dir.pre,"infer_sgl_model.R"))
source(paste0(code.dir.pre, "utility_functions.R"))
source(paste0(code.dir.pre, "kimono.R"))

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

options(future.globals.maxSize = 100000 * 1024^2) # more memory for each thread = 100Gb
plan(multisession, workers = 12)

endnode <- min(startnode + 10000, length(colnames(input_list$methyl)))
node_list <- colnames(input_list$methyl)[startnode:endnode]
results <- future_map(node_list, run_kimono_para, stab_sel = TRUE, niterations = 20, .options = furrr_options(seed = TRUE)) #add one more layer of parallelization in kimono.R (seeds)

end.time <- Sys.time()
end.time
print("kimono done")
end.time - start.time

results <- do.call(rbind, results) # make one big data table



###############################################
# 5 Save result
###############################################
print("write results out")
fwrite(kimono.res, 
       kimono.res.fn, 
       quote = F, row.names = F, sep = ";") 

print("result written out")
