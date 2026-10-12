args <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[1] else "meta_overall_CH.txt"
outfile <- if (length(args) >= 2) args[2] else sub("\\.[^.]+$", "_with_meta_het.txt", infile)

dat <- read.delim(infile, check.names = FALSE, stringsAsFactors = FALSE)

beta_cols <- grep("^beta_", names(dat), value = TRUE)
beta_cols <- setdiff(beta_cols, "beta_META")
datasets <- sub("^beta_", "", beta_cols)
se_cols <- paste0("se_", datasets)

keep <- se_cols %in% names(dat)
beta_cols <- beta_cols[keep]
datasets <- datasets[keep]
se_cols <- se_cols[keep]

fixed_effect_meta <- function(row_vals) {
  betas <- as.numeric(row_vals[beta_cols])
  ses   <- as.numeric(row_vals[se_cols])

  present <- !(is.na(betas) | is.na(ses))

  direction <- ifelse(!present, "?",
                      ifelse(betas > 0, "+",
                             ifelse(betas < 0, "-", "0")))
  direction <- paste(direction, collapse = "")

  if (!any(present)) {
    return(c(beta_META = NA, se_META = NA, or_META = NA,
             l95cl_META = NA, u95cl_META = NA, pval_META = NA,
             pval_het_META = NA, direction_META = direction))
  }

  betas_use <- betas[present]
  ses_use   <- ses[present]
  w <- 1 / (ses_use^2)

  beta_meta <- sum(w * betas_use) / sum(w)
  se_meta   <- sqrt(1 / sum(w))
  z_meta    <- beta_meta / se_meta
  pval_meta <- 2 * pnorm(abs(z_meta), lower.tail = FALSE)

  or_meta    <- exp(beta_meta)
  l95cl_meta <- exp(beta_meta - 1.96 * se_meta)
  u95cl_meta <- exp(beta_meta + 1.96 * se_meta)

  if (sum(present) >= 2) {
    Q <- sum(w * (betas_use - beta_meta)^2)
    df_het <- sum(present) - 1
    pval_het <- pchisq(Q, df = df_het, lower.tail = FALSE)
  } else {
    pval_het <- NA
  }

  c(beta_META = beta_meta,
    se_META = se_meta,
    or_META = or_meta,
    l95cl_META = l95cl_meta,
    u95cl_META = u95cl_meta,
    pval_META = pval_meta,
    pval_het_META = pval_het,
    direction_META = direction)
}

meta_res <- t(apply(dat, 1, fixed_effect_meta))
meta_res <- as.data.frame(meta_res, stringsAsFactors = FALSE)

num_cols <- c("beta_META", "se_META", "or_META", "l95cl_META", "u95cl_META", "pval_META", "pval_het_META")
for (nm in num_cols) meta_res[[nm]] <- as.numeric(meta_res[[nm]])

existing_meta_cols <- intersect(c("beta_META", "se_META", "or_META", "l95cl_META",
                                  "u95cl_META", "pval_META", "pval_het_META",
                                  "direction_META"), names(dat))
if (length(existing_meta_cols) > 0) {
  dat <- dat[ , !(names(dat) %in% existing_meta_cols), drop = FALSE]
}

out <- cbind(dat, meta_res)

write.table(out, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat("Wrote meta-analysis results to:", outfile, "\n")
cat("Direction string order:\n")
cat(paste(datasets, collapse = ", "), "\n")
