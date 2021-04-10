
# Set up parameters

src.data.pre   <- "~/bio/datasets/gene_expression/20_DEA/" 
# src.data.pre <- "/binder/mgp/datasets/2020_DexStim_Array_Human/gene_expression/20_DEA/"

lmer.mdl.fn <- paste0(src.data.pre, "01_lme_models/lmer_all_plus_sv_1_5.csv")
lmer.mdl.fn <- paste0(src.data.pre, "01_lme_models/lmer_all_plus_bcc.csv")
lmer.mdl    <- read.csv(lmer.mdl.fn, sep = "\t")

summary(lmer.mdl)
boxplot(lmer.mdl$FC)

delta.beta <- 0.2
p.thr      <- 0.01

gex.sign.df  <- lmer.mdl[abs(lmer.mdl$FC) >= delta.beta & lmer.mdl$pFDR <= p.thr, ]
gex.sign.df  <- gex.sign.df[order(gex.sign.df$pFDR),]

gex.sign.df.fn <- paste0(src.data.pre, "02_sign_gex/gex_significant_with")

write.table(gex.sign.df,
            file = paste0(gex.sign.df.fn, 
                          "_beta_", delta.beta * 100, 
                          "_p_", p.thr * 100, 
                          ".txt"), row.names = F, quote = F, sep = "\t", col.names = T)
