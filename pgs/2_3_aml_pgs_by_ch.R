analysis_df <- analysis_df_saved_for_later

required_pkgs <- c(required_pkgs, "survival")
install_if_missing(required_pkgs)

suppressPackageStartupMessages({
  library(survival)
})

# =========================
# User-defined settings
# =========================

ref_lower_pct <- 0
ref_upper_pct <- 2
comp_lower_pct <- 98
comp_upper_pct <- 100

# =========================
# Case-control counts
# =========================

counts <- table(analysis_df$aml)

# =========================
# Main Cox model
# =========================

cox_formula <- as.formula(
  Surv(time = age, time2 = dxage, event = aml) ~
    pgs_std + sex +
    PC1 + PC2 + PC3 + PC4 + PC5 +
    PC6 + PC7 + PC8 + PC9 + PC10
)

cox_model <- coxph(cox_formula, data = analysis_df)
cox_model_summary <- summary(cox_model)

hr_pgs <- tibble(
  term = "pgs_std",
  beta = unname(coef(cox_model)["pgs_std"]),
  HR = unname(exp(coef(cox_model)["pgs_std"])),
  lower_95 = unname(exp(confint(cox_model)["pgs_std", 1])),
  upper_95 = unname(exp(confint(cox_model)["pgs_std", 2])),
  p_value = unname(cox_model_summary$coefficients["pgs_std", "Pr(>|z|)"])
)

ph_test <- cox.zph(cox_model)

cat(
  "HR per SD increase in pgs_std: ",
  round(hr_pgs$HR, 3),
  " (95% CI ",
  round(hr_pgs$lower_95, 3),
  " to ",
  round(hr_pgs$upper_95, 3),
  "), p = ",
  signif(hr_pgs$p_value, 3),
  "\n",
  sep = ""
)

# =========================
# Decile analysis
# =========================

analysis_df <- analysis_df %>%
  mutate(
    pgs_decile = ntile(pgs_std, 10),
    pgs_decile = factor(pgs_decile, levels = 1:10, labels = paste0("D", 1:10))
  )

cox_model_deciles <- coxph(
  Surv(time = age, time2 = dxage, event = aml) ~
    pgs_decile + sex +
    PC1 + PC2 + PC3 + PC4 + PC5 +
    PC6 + PC7 + PC8 + PC9 + PC10,
  data = analysis_df
)

cox_model_deciles_summary <- summary(cox_model_deciles)

decile_results <- cox_model_deciles_summary$coefficients %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  filter(str_detect(term, "^pgs_decile")) %>%
  transmute(
    decile = str_remove(term, "^pgs_decile"),
    HR = exp(coef),
    lower_95 = exp(coef - 1.96 * `se(coef)`),
    upper_95 = exp(coef + 1.96 * `se(coef)`),
    p_value = `Pr(>|z|)`
  ) %>%
  bind_rows(
    tibble(
      decile = "D1",
      HR = 1,
      lower_95 = NA_real_,
      upper_95 = NA_real_,
      p_value = NA_real_
    )
  ) %>%
  mutate(decile = factor(decile, levels = paste0("D", 1:10))) %>%
  arrange(decile)

# =========================
# Null model + C-index + LRT
# =========================

cox_model_null <- coxph(
  Surv(time = age, time2 = dxage, event = aml) ~ sex,
  data = analysis_df
)

lp_full <- predict(cox_model, type = "lp")
lp_null <- predict(cox_model_null, type = "lp")

cidx_obj <- concordance(
  Surv(time = age, time2 = dxage, event = aml) ~ lp_full + lp_null,
  data = analysis_df,
  reverse = TRUE
)

c_values <- cidx_obj$concordance
c_var <- cidx_obj$var

c_full <- c_values[1]
c_null <- c_values[2]

se_full <- sqrt(c_var[1, 1])
se_null <- sqrt(c_var[2, 2])

cindex_results <- tibble(
  model = c("cox_model", "cox_model_null"),
  c_index = c(c_full, c_null),
  lower_95 = c(c_full - 1.96 * se_full, c_null - 1.96 * se_null),
  upper_95 = c(c_full + 1.96 * se_full, c_null + 1.96 * se_null)
)

delta_c <- c_full - c_null
se_delta <- sqrt(c_var[1, 1] + c_var[2, 2] - 2 * c_var[1, 2])
lower_delta <- delta_c - 1.96 * se_delta
upper_delta <- delta_c + 1.96 * se_delta
p_delta <- 2 * pnorm(abs(delta_c / se_delta), lower.tail = FALSE)

delta_cindex_results <- tibble(
  comparison = "cox_model - cox_model_null",
  delta_c_index = delta_c,
  lower_95 = lower_delta,
  upper_95 = upper_delta,
  p_value = p_delta
)

lrt_results <- anova(cox_model_null, cox_model, test = "LRT")

# =========================
# Extreme percentile analysis
# =========================

ref_label <- paste0("P", ref_lower_pct, "-<", ref_upper_pct)
comp_label <- paste0("P", comp_lower_pct, "-", comp_upper_pct)

analysis_df_pct <- analysis_df %>%
  mutate(
    pgs_percentile = percent_rank(pgs_std) * 100,
    pgs_group = case_when(
      pgs_percentile >= ref_lower_pct & pgs_percentile < ref_upper_pct ~ ref_label,
      pgs_percentile >= comp_lower_pct & pgs_percentile <= comp_upper_pct ~ comp_label
    )
  ) %>%
  filter(!is.na(pgs_group)) %>%
  mutate(pgs_group = factor(pgs_group, levels = c(ref_label, comp_label)))

cox_model_pct <- coxph(
  Surv(time = age, time2 = dxage, event = aml) ~
    pgs_group + sex +
    PC1 + PC2 + PC3 + PC4 + PC5 +
    PC6 + PC7 + PC8 + PC9 + PC10,
  data = analysis_df_pct
)

pct_term <- paste0("pgs_group", comp_label)

hr_pct <- tibble(
  comparison = paste0(comp_label, " vs ", ref_label),
  HR = exp(coef(cox_model_pct)[pct_term]),
  lower_95 = exp(confint(cox_model_pct)[pct_term, 1]),
  upper_95 = exp(confint(cox_model_pct)[pct_term, 2]),
  p_value = summary(cox_model_pct)$coefficients[pct_term, "Pr(>|z|)"]
)

# =========================
# Collect outputs
# =========================

analysis_results <- list(
  models = list(
    cox_model = cox_model,
    cox_model_deciles = cox_model_deciles,
    cox_model_null = cox_model_null,
    cox_model_pct = cox_model_pct
  ),
  results = list(
    hr_pgs = hr_pgs,
    ph_test = ph_test,
    decile_results = decile_results,
    cindex_results = cindex_results,
    delta_cindex_results = delta_cindex_results,
    lrt_results = lrt_results,
    hr_pct = hr_pct,
    counts = counts
  )
)

names(analysis_results)
analysis_results$results