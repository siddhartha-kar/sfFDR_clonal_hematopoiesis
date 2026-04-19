### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(sffdr)

### Load data:
overall_bcx <- vroom("overall_bcx_wbc_rbc_platelet_harmonised_combined_sffdr.tsv.gz")
overall_clumped <- vroom("overall_clumped.valid.snp", delim = "\t", col_names = FALSE)

### Select p-value columns:
colnames(overall_bcx)
overall_bcx2 <- overall_bcx %>%
  select(-c(SNP, rsid, CHR, BP))
overall_bcx3 <- as.data.frame(overall_bcx2)
dim(overall_bcx3)
head(overall_bcx3)

### Set seed:
set.seed(123)

### Create variables:
p <- overall_bcx3$pval_overall
z <- as.matrix(overall_bcx3[,-1])

### Specify independent SNPs:
overall_clumped <- overall_clumped %>%
  rename(rsid = X1)
indep_snps_rsid <- overall_clumped$rsid
indep_snps <- overall_bcx$rsid %in% indep_snps_rsid

### Create model: 
mpi0 <- pi0_model(z = z,
                  indep_snps = indep_snps,
                  knots = c(0.01, 0.025, 0.05, 0.1))

### Estimate fpi0 using design matrix from mpi0:
fpi0 <- fpi0est(p = p,
                z = mpi0$zt,
                indep_snps = indep_snps,
                pi0_model = mpi0$fmod)

### Estimate FDR quantities and functional p-value:
sffdr_out <- sffdr(p,
                   fpi0 = fpi0$fpi0,
                   indep_snps = indep_snps)

### Plot significant results:
pdf("overall_ch_bcx_sffdr_significance_plot.pdf")
plot(sffdr_out, rng = c(0, 5e-8))
dev.off()

### Determine functional pval, qval and local fdr values:
fp <- sffdr_out$fpvalues
fq <- sffdr_out$fqvalues
flfdr <- sffdr_out$flfdr
fpi0 <- sffdr_out$fpi0

### Combine results and save:
rsid <- overall_bcx %>% select(rsid, CHR, BP, pval_overall)
sffdr_results <- data.frame(SNP = rsid$rsid,
                            CHR = rsid$CHR,
                            BP = rsid$BP,
                            Overall_original_pvalue = rsid$pval_overall,
                            Functional_pvalue = fp,
                            Functional_qvalue = fq,
                            Local_FDR = flfdr,
                            fpi0 = fpi0)
head(sffdr_results)
write.table(sffdr_results, file = "overall_ch_bcx_eur_sffdr_results.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)