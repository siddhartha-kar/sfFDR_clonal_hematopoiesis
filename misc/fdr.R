
dat <- read.table(
"fdr_input.txt",
header = TRUE,
sep = "\t",
stringsAsFactors = FALSE,
check.names = FALSE
)
dat$pval_META_onesided <- ifelse(
dat$Direction == "SAME",
dat$pval_META / 2,
ifelse(
dat$Direction == "DIFF",
1 - dat$pval_META / 2,
NA
)
)
dat$FDR_BH <- p.adjust(dat$pval_META_onesided, method = "BH")
dat$FDR_sig <- dat$FDR_BH < 0.05
dat[, c(
"CH_phenotype",
"rsID",
"pval_META",
"Direction",
"pval_META_onesided",
"FDR_BH",
"FDR_sig"
)]
write.table(
dat,
"fdr_results.txt",
sep = "\t",
quote = FALSE,
row.names = FALSE
)