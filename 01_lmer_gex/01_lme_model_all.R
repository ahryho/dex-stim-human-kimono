require(lme4)
require(foreign)

library(parallel)
library(foreach)
library(doParallel)

library(Biobase)

# 3. with age, sex, BMI, and MDD status, cell types or SV1-5 as covariates

# 1. Load data

# input.parameters.fn <- "input_parameters.csv"

src.data.pre       <- "~/bio/datasets/gene_expression/"
src.data.pre       <- "/binder/mgp/datasets/2020_DexStim_Array_Human/gene_expression/"
# lmer.res.out.fn <- paste0(src.data.pre, "lmer_all_plus_sv_1_5.csv")

args                <- commandArgs(T)
lmer.res.out.fn     <- as.character(args[1])

# Load gene-expression data

gex.fn <- paste0(src.data.pre, "ESet_dex_sva.rda")
x      <- load(gex.fn)
gex    <- get(x)
rm(x)

gex.tbl  <- exprs(gex)
pheno    <- pData(gex)

# Prepre the pheno data

pheno$Sample_ID <- as.factor(pheno$Sample_ID)
pheno$Sex       <- as.factor(pheno$Sex)
pheno$Age       <- as.numeric(pheno$Age)
pheno$BMI_D1    <- as.numeric(pheno$BMI_D1)
pheno$Status    <- as.factor(pheno$Status)
pheno$Dex       <- as.factor(pheno$Dex)

# Surrogative variables

sv.list <- colnames(pheno)[c(20:24)]#ncol(pheno))]
nr.sv   <- length(sv.list)
form <- foreach(sv.iter = 1:nr.sv, .combine = cbind) %dopar% {
  # tmp <-sv.list[sv.iter] 
  # data <- paste("pheno$", sv.list[sv.iter], sep = "")
  # assign(tmp, eval(parse(text = data)))
  sv.list[sv.iter]
}
form <- paste(form, collapse = "+")

# Cell types:
pheno$Cort_D1   <- as.numeric(pheno$Cort_D1)
pheno$ACTH_D1   <- as.numeric(pheno$ACTH_D1)
pheno$Leukos_D1 <- as.numeric(pheno$Leukos_D1)
pheno$Gran_D1   <- as.numeric(pheno$Gran_D1)
pheno$Mono_D1   <- as.numeric(pheno$Mono_D1)
pheno$Lymph_D1  <- as.numeric(pheno$Lymph_D1)

samples.veh.ids <- pheno$RNA_ID[pheno$Dex == 0]
samples.dex.ids <- pheno$RNA_ID[pheno$Dex == 1]

# 2. Making sure about samples in pheno and and betas matrix in the same order

table(colnames(gex.tbl) %in% pheno$RNA_ID)
all(pheno$RNA_ID == colnames(gex.tbl))


# 3. Build model

no.cores <- detectCores() - 1
cl <- makeCluster(no.cores)
registerDoParallel(cl)

res <- foreach(gex.iter = 2:nrow(gex.tbl), .combine = rbind, .packages = 'lme4') %dopar% {
  
  # for (gex.iter in 1:nrow(gex.tbl)){
  sample <- as.numeric(as.factor(pheno$Sample_ID))
  
  formula.null <- as.formula(paste0("gex.tbl[gex.iter, ] ~ 1 + Sex + Age + BMI_D1 + Status + ", 
                                    "Cort_D1 + ACTH_D1 + Leukos_D1 + Gran_D1 + Mono_D1 + Lymph_D1 + ",
                                    # form,
                                    "+ (1|sample)"))
  lmer.null  <- lmer(formula.null, REML = F, data = pheno)
  
  formula.model <- as.formula(paste0("gex.tbl[gex.iter, ] ~ 1 + Sex + Age + BMI_D1 + Status + ",
                                     "Cort_D1 + ACTH_D1 + Leukos_D1 + Gran_D1 + Mono_D1 + Lymph_D1",
                                     # form,
                                     " + Dex + (1|sample)"))
  lmer.model <- lmer(formula.model, REML = F, data = pheno)
  # lmer.model <- lmer(gex.tbl[gex.iter, ] ~ pheno$Sex + pheno$Age + pheno$BMI_D1 + pheno$Status + 
  #                     pheno$Cort_D1 + pheno$ACTH_D1 + pheno$Leukos_D1 + pheno$Gran_D1 + pheno$Mono_D1 + pheno$Lymph_D1 +
  #                     form + 
  #                     pheno$Dex + (1|sample), REML = F)
  # 
  res.anova  <- anova(lmer.null, lmer.model)
  
  
  lmer.pheno <- lmer(gex.tbl[gex.iter, ] ~ pheno$Dex + (1|sample), REML = F)
  coefs      <- data.frame(coef(summary(lmer.pheno)))
  coefs$p.z  <- 2 * (1 - pnorm(abs(coefs$t.value)))
  pval.group <- round(coefs$p.z[2], 3)

  mean.veh <- round(mean(gex.tbl[gex.iter, samples.veh.ids]), 3)
  mean.dex <- round(mean(gex.tbl[gex.iter, samples.dex.ids]), 3)
  sd.veh   <- round(sd(gex.tbl[gex.iter, samples.veh.ids]), 3)
  sd.dex   <- round(sd(gex.tbl[gex.iter, samples.dex.ids]), 3)
  var.veh  <- round(var(gex.tbl[gex.iter, samples.veh.ids]), 3)
  var.dex  <- round(var(gex.tbl[gex.iter, samples.dex.ids]), 3)
  
  fc       <- round(mean.dex - mean.veh, 3)
  var.all  <- round(var(gex.tbl[gex.iter, ]), 3)
  
  cbind(as.character(rownames(gex.tbl)[gex.iter]), 
        signif(res.anova$"Pr(>Chisq)"[2], 3),
        round(res.anova$Chisq[2], 3),
        mean.veh, sd.veh, signif(var.veh),
        mean.dex, sd.dex, signif(var.dex),
        fc, signif(var.all),
        signif(pval.group))
}

stopImplicitCluster()

res <- as.data.frame(res)
res.colnames  <- cbind("Probe_Id", "pVal", "Chisq", 
                      "Mean_Veh", "Sd_Veh", "Var_Veh", 
                      "Mean_Dex", "Sd_Dex", "Var_Dex",
                      "FC", "Var_Total_pheno", "pVal_pheno")
colnames(res) <- res.colnames

# 4. Adjust p values
mthd    <- "fdr" 
res$pFDR <- p.adjust(res$pVal, method = mthd) 

write.table(res,
            file = lmer.res.out.fn, row.names = F, quote = F, sep = "\t", col.names = T, append = F)