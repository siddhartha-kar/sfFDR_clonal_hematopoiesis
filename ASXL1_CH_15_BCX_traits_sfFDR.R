### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(sffdr)

### Load data:
asxl1_bcx <- vroom("asxl1_bcx_wbc_rbc_platelet_harmonised_combined_sffdr.tsv.gz")
asxl1_clumped <- vroom("asxl1_clumped.valid.snp", delim = "\t", col_names = FALSE)

### Select p-value columns:
colnames(asxl1_bcx)
asxl1_bcx2 <- asxl1_bcx %>%
  select(-c(SNP, rsid, CHR, BP))
asxl1_bcx3 <- as.data.frame(asxl1_bcx2)
dim(asxl1_bcx3)
head(asxl1_bcx3)

### Set seed:
set.seed(123)

### Create variables:
p <- asxl1_bcx3$pval_asxl1
z <- as.matrix(asxl1_bcx3[,-1])

### Specify independent SNPs:
asxl1_clumped <- asxl1_clumped %>%
  rename(rsid = X1)
indep_snps_rsid <- asxl1_clumped$rsid
indep_snps <- asxl1_bcx$rsid %in% indep_snps_rsid

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
pdf("asxl1_ch_bcx_sffdr_significance_plot.pdf")
plot(sffdr_out, rng = c(0, 5e-8))
dev.off()

### Determine functional pval, qval and local fdr values:
fp <- sffdr_out$fpvalues 
fq <- sffdr_out$fqvalues
flfdr <- sffdr_out$flfdr
fpi0 <- sffdr_out$fpi0

### Combine results and save:
rsid <- asxl1_bcx %>% select(rsid, CHR, BP, pval_asxl1)
sffdr_results <- data.frame(SNP = rsid$rsid,
                            CHR = rsid$CHR,
                            BP = rsid$BP,
                            asxl1_original_pvalue = rsid$pval_asxl1,
                            Functional_pvalue = fp,
                            Functional_qvalue = fq,
                            Local_FDR = flfdr,
                            fpi0 = fpi0)
head(sffdr_results)
write.table(sffdr_results, file = "asxl1_ch_bcx_eur_sffdr_results.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)