
if (!requireNamespace("MendelianRandomization", quietly = TRUE)) {
  stop("Install MendelianRandomization first: install.packages('MendelianRandomization')")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 first: install.packages('ggplot2')")
}
suppressPackageStartupMessages(library(MendelianRandomization))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Usage: Rscript mr_analysis.R [input_file] [output_directory]")

if (length(args) >= 1L) {
  input_file <- args[[1L]]
} else if (file.exists("ch_aml_mr.txt")) {
  input_file <- "ch_aml_mr.txt"
} else {
  input_file <- file.path("upload", "ch_aml_mr.txt")
}
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

output_dir <- if (length(args) == 2L) args[[2L]] else "mr_results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Could not create output directory: ", output_dir)

dat <- read.delim(input_file, sep = "\t", header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE, na.strings = c("", "NA"))
required <- c("rsid", "beta.x", "se.x", "beta.y", "se.y")
if (!all(required %in% names(dat))) {
  stop("Missing columns: ", paste(setdiff(required, names(dat)), collapse = ", "))
}
if (nrow(dat) < 3L) stop("At least three SNPs are needed for the sensitivity methods.")
dat$rsid <- trimws(as.character(dat$rsid))
if (anyNA(dat$rsid) || any(!nzchar(dat$rsid))) stop("Missing SNP identifier(s).")
if (anyDuplicated(dat$rsid)) stop("Duplicate SNP identifiers found.")

for (col in required[-1L]) {
  dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
  if (anyNA(dat[[col]]) || any(!is.finite(dat[[col]]))) {
    stop("Missing, nonnumeric, or nonfinite value in column ", col)
  }
}
if (any(dat$se.x <= 0) || any(dat$se.y <= 0)) {
  stop("Standard errors must be strictly positive.")
}
if (any(dat$beta.x == 0)) stop("beta.x contains zero; Wald ratios are undefined.")

exposure_label <- "Exposure"
outcome_label <- "Outcome"
to_mr_input <- function(d) {
  mr_input(bx = d$beta.x, bxse = d$se.x,
           by = d$beta.y, byse = d$se.y,
           snps = d$rsid, exposure = exposure_label, outcome = outcome_label)
}
mr_data <- to_mr_input(dat)

snp_checks <- data.frame(
  rsid = dat$rsid,
  F_approx = (dat$beta.x / dat$se.x)^2,
  Wald_ratio = dat$beta.y / dat$beta.x,
  Wald_SE_first_order = dat$se.y / abs(dat$beta.x)
)
write.csv(snp_checks, file.path(output_dir, "snp_checks.csv"), row.names = FALSE)

ivw_random <- mr_ivw(mr_data, model = "random")
ivw_fixed <- mr_ivw(mr_data, model = "fixed")
egger <- mr_egger(mr_data)
weighted_median <- mr_median(mr_data, weighting = "weighted",
                             iterations = 10000, seed = 314159265)
simple_median <- mr_median(mr_data, weighting = "simple",
                           iterations = 10000, seed = 314159265)

result_row <- function(label, fit, egger_slope = FALSE) {
  if (egger_slope) {
    data.frame(method = label, n_snps = as.numeric(fit@SNPs),
               estimate = as.numeric(fit@Estimate),
               standard_error = as.numeric(fit@StdError.Est),
               ci_95_lower = as.numeric(fit@CILower.Est),
               ci_95_upper = as.numeric(fit@CIUpper.Est),
               p_value = as.numeric(fit@Pvalue.Est))
  } else {
    data.frame(method = label, n_snps = as.numeric(fit@SNPs),
               estimate = as.numeric(fit@Estimate),
               standard_error = as.numeric(fit@StdError),
               ci_95_lower = as.numeric(fit@CILower),
               ci_95_upper = as.numeric(fit@CIUpper),
               p_value = as.numeric(fit@Pvalue))
  }
}
estimates <- rbind(
  result_row("IVW (multiplicative random effects)", ivw_random),
  result_row("IVW (fixed effect)", ivw_fixed),
  result_row("Weighted median", weighted_median),
  result_row("Simple median", simple_median),
  result_row("MR-Egger slope", egger, egger_slope = TRUE)
)
write.csv(estimates, file.path(output_dir, "mr_estimates.csv"), row.names = FALSE)

ivw_q <- as.numeric(unlist(ivw_random@Heter.Stat, use.names = FALSE))
egger_q <- as.numeric(unlist(egger@Heter.Stat, use.names = FALSE))
if (length(ivw_q) < 2L || length(egger_q) < 2L) {
  stop("The package did not return the expected Q statistics.")
}
diagnostics <- data.frame(
  statistic = c("IVW Cochran Q", "IVW Q df", "IVW Q p-value",
                "MR-Egger residual Q", "MR-Egger Q df", "MR-Egger Q p-value",
                "MR-Egger intercept", "MR-Egger intercept SE",
                "MR-Egger intercept p-value", "MR-Egger I2_GX"),
  value = c(ivw_q[[1L]], nrow(dat) - 1L, ivw_q[[2L]],
            egger_q[[1L]], nrow(dat) - 2L, egger_q[[2L]],
            as.numeric(egger@Intercept), as.numeric(egger@StdError.Int),
            as.numeric(egger@Pvalue.Int), as.numeric(egger@I.sq))
)
write.csv(diagnostics, file.path(output_dir, "mr_diagnostics.csv"), row.names = FALSE)

leave_one_out <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
  fit <- mr_ivw(to_mr_input(dat[-i, , drop = FALSE]), model = "random")
  data.frame(omitted_rsid = dat$rsid[[i]],
             estimate = as.numeric(fit@Estimate),
             standard_error = as.numeric(fit@StdError),
             ci_95_lower = as.numeric(fit@CILower),
             ci_95_upper = as.numeric(fit@CIUpper),
             p_value = as.numeric(fit@Pvalue))
}))
write.csv(leave_one_out, file.path(output_dir, "mr_leave_one_out.csv"),
          row.names = FALSE)

grDevices::pdf(file.path(output_dir, "mr_scatter.pdf"), width = 9, height = 7)
tryCatch(
  mr_plot(mr_data, line = "ivw", orientate = TRUE,
          interactive = FALSE, labels = TRUE),
  finally = grDevices::dev.off()
)
ggplot2::ggsave(file.path(output_dir, "mr_forest.pdf"),
                plot = mr_forest(mr_data, ordered = TRUE),
                width = 9, height = 9, units = "in")
ggplot2::ggsave(file.path(output_dir, "mr_funnel.pdf"),
                plot = mr_funnel(mr_data),
                width = 8, height = 7, units = "in")
ggplot2::ggsave(file.path(output_dir, "mr_leave_one_out.pdf"),
                plot = mr_loo(mr_data),
                width = 9, height = 9, units = "in")

capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
message("Analyzed ", nrow(dat), " SNPs. Results: ",
        normalizePath(output_dir, winslash = "/"))
message("Minimum per-SNP F approximation: ",
        signif(min(snp_checks$F_approx), 4),
        "; number below 10: ", sum(snp_checks$F_approx < 10))
print(estimates, row.names = FALSE)
