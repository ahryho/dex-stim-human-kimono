library(MatrixEQTL)
library(data.table)

args <- commandArgs(T)
treatment <- as.character(args[1]) #"veh"
eqtm.pre <- as.character(args[2]) 
# eqtm.pre <- "~/bio/datasets/eQTM/"

eqtm.in.pre  <- paste0(eqtm.pre, "input/")
eqtm.res.pre <- paste0(eqtm.pre, "result/")

cpg.loc.fn <- paste0(eqtm.in.pre, "cpg_locations.csv")
ensg.loc.fn <- paste0(eqtm.in.pre, "ensg_locations.csv")
snp.loc.fn <- paste0(eqtm.in.pre, "snp_locations.csv")

gex.layer.fn <- paste0(eqtm.in.pre, "gex_mtrx_", treatment, ".csv")
methyl.layer.fn <- paste0(eqtm.in.pre, "methyl_beta_mtrx_", treatment, ".csv")
bio.layer.fn  <- paste0(eqtm.in.pre, "bio_mtrx_", treatment, ".csv")

eqtm.cis.result.fn <- paste(eqtm.res.pre, "result", "eqtm_cis_result.csv")
eqtm.trans.result.fn <- paste(eqtm.res.pre, "result", "eqtm_trans_result.csv")

# Load data

cpg.loc  <- fread(cpg.loc.fn)
ensg.loc <- fread(ensg.loc.fn)
snp.loc <- fread(snp.loc.fn)

# gex.layer <- fread(gex.layer.fn) 
# methyl.layer <- fread(methyl.layer.fn)
# bio.layer <- fread(bio.layer.fn)

# Check the colnames of layers are in the same order

# gex.layer[, colnames(gex.layer) %in% colnames(methyl.layer)]
# bio.layer[, colnames(bio.layer) %in% colnames(methyl.layer)]

run_matrix_eqtl <- function(SNP_file, expression_file, cov_file, outfile_cis, outfile_trans, cis, trans){
  # 1. Set general parameters
  # Use linear model (other models: modelANOVA, modelLINEAR, or modelLINEAR_CROSS)
  useModel <- modelLINEAR
  # We need to set up the p-value cutoff and the error model (in this case assuming independent errors)
  # Only associations significant at this level will be saved
  pvOutputThreshold_cis = cis  # cis-eQTLS cutoff: 0.05 = 5e-2
  pvOutputThreshold_tra = trans  # trans-eQTLs cutoff: 0.01 = 1e-2 (0 means no trans-eQTLs)
  # Error covariance matrix: set to numeric() for identity
  errorCovariance <- numeric()
  # Distance for local (cis) gene-SNP pairs
  cisDist = 1e6  # cis window of 1Mb = 1e6, 100kb = 1e5
  
  # 2. Set Data Up
  # Now we need to set up the snp and gene expression data in the special format required by the MatrixEQTL package
  # 2.1 SNP Data
  snps = SlicedData$new()
  snps$fileDelimiter = ";"     # the TAB character
  snps$fileOmitCharacters = "NA" # denote missing values;
  snps$fileSkipRows = 1          # one row of column labels
  snps$fileSkipColumns = 1       # one column of row labels
  snps$fileSliceSize = 100000     # read file in pieces of 100,000 rows
  snps$LoadFile(SNP_file)
  
  # 2.2 Gene Expression Data
  gene = SlicedData$new()
  gene$fileDelimiter = ";"      # the TAB character
  gene$fileOmitCharacters = "NA" # denote missing values;
  gene$fileSkipRows = 1          # one row of column labels
  gene$fileSkipColumns = 1      # one column of row labels
  gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
  gene$LoadFile(expression_file)
  
  # 2.3 Covariates
  cvrt = SlicedData$new()
  cvrt$fileDelimiter = ";"      # the TAB character
  cvrt$fileOmitCharacters = "NA" # denote missing values;
  cvrt$fileSkipRows = 2          # one row of column labels + one row of treatment
  cvrt$fileSkipColumns = 1       # one column of row labels
  
  if(length(cov_file) > 0){
    cvrt$LoadFile(cov_file)  # read file if given
  }
  
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = outfile_trans,
    pvOutputThreshold = pvOutputThreshold_tra , #pvOutputThreshold_tra only cis
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = outfile_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = cpg.loc,
    genepos = ensg.loc,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)
  
  returm (me)
}

# Run matrixEQTL only for cis and save results into an RData file
me.all <- run_matrix_eqtl(SNP_file = methyl.layer.fn, 
                          expression_file = gex.layer.fn, 
                          cov_file = bio.layer.fn, 
                          outfile_cis = eqtm.cis.result.fn, 
                          outfile_trans = eqtm.trans.result.fn, 
                          cis = 5e-2, trans = 0)

saveRDS(me.all, file =  paste0(eqtm.res.pre, Sys.Date(), "eQTM_matrx.RDS"))

# base.dir = find.package("MatrixEQTL");
# useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS
# SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
# expression_file_name = paste(base.dir, "/data/GE.txt", sep="");

