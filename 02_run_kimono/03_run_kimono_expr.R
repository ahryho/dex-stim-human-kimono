
###############################################
### 1 libraries
###############################################
library(kimono)
library(data.table)

###############################################
data.dir.pre  <- "/binder/mgp/datasets/2020_DexStim_Array_Human/kimono/"
# data.dir.pre  <- "~/bio/datasets/kimono/"

kimono.res.fn <- paste0(data.dir.pre, "output/kimono_res.csv")
###############################################
### 2 Data
###############################################
# read in data

print("read data")

layer.gex    <- fread(paste0(data.dir.pre, "input/gex_veh.csv"))
layer.pheno  <- fread(paste0(data.dir.pre, "input/pheno_veh.csv"))

# read in mapping (prior)
prior.gex       <- fread(paste0(data.dir.pre, "input/prior_biogrid_ensg.csv"))
prior.gex.pheno <- fread(paste0(data.dir.pre, "input/prior_gene_pheno.csv"))

print("data read successfully")

# make sure samples are in the same order
idorder      <- layer.pheno$RNA_ID
layer.gex    <- layer.gex[match(idorder, layer.gex$X),]

# remove sample column
layer.gex    <- layer.gex[, -1]
layer.pheno  <- layer.pheno %>% select(Sex, Age, Status)

###############################################
# 3 Assemble into lists
###############################################
print("assemble data")

# data list
input_list <- list(
  as.data.table(layer.gex),
  as.data.table(layer.pheno)
)
names(input_list) <- c('expr',
                       'pheno')
#########################
# mapping list
mapping_list <- list(
  as.data.table(prior.gex),
  as.data.table(prior.gex.pheno)
)
#########################
# meta info
metainfo <- data.frame('ID' = c('prior_expr', 'expr_pheno'),
                       'main_to'   =  c(1, 2)
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
