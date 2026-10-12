# Commands to run in the Terminal:
# dx download sffdr_ch_gwas/pgs_165.txt
# dx download sffdr_ch_gwas/eid_pgs.txt
# dx download sffdr_ch_gwas/age_sex_10pcs_dor_participant.tsv
# dx download sffdr_ch_gwas/ch_binary_participant.tsv
# dx download sffdr_ch_gwas/ltl_pgs.txt
# dx download sffdr_ch_gwas/eid_ltl_participant.tsv
# dx download sffdr_ch_gwas/ltl_eid_pgs.txt

required_pkgs <- c("tidyverse")

install_if_missing <- function(pkgs) {
  missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs)
  }
}

install_if_missing(required_pkgs)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

if (!requireNamespace("ukbrapR", quietly = TRUE)) {
  remotes::install_github("lcpilling/ukbrapR@*release")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(ukbrapR)
})

eid_cov <- read_tsv("age_sex_10pcs_dor_participant.tsv", show_col_types = FALSE) %>%
  rename(
    age = p21022,
    sex = p31,
    PC1 = p22009_a1,
    PC2 = p22009_a2,
    PC3 = p22009_a3,
    PC4 = p22009_a4,
    PC5 = p22009_a5,
    PC6 = p22009_a6,
    PC7 = p22009_a7,
    PC8 = p22009_a8,
    PC9 = p22009_a9,
    PC10 = p22009_a10,
    recdate = p53_i0
  )

eid_cov_model <- eid_cov %>%
  select(-recdate)

pgs_variants <- read_tsv("pgs_165.txt", show_col_types = FALSE)

ltl_pgs_variants <- read_tsv("ltl_pgs.txt", show_col_types = FALSE)

pgs_weights <- pgs_variants %>%
  transmute(
    rsid = rsid,
    chr = chr,
    pos = pos,
    effect_allele = ea,
    other_allele = oa,
    beta = beta,
    af_effect = eaf
  )

ltl_pgs_weights <- ltl_pgs_variants %>%
  transmute(
    rsid = rsid,
    chr = chr,
    pos = pos,
    effect_allele = ea,
    other_allele = oa,
    beta = beta,
    af_effect = eaf
  )

# eid_pgs <- create_pgs(
#   in_file = pgs_weights,
#   out_file = "ch.pgs",
#   pgs_name = "pgs"
# )

# ltl_eid_pgs <- create_pgs(
#    in_file = ltl_pgs_weights,
#    out_file = "ltl.pgs",
#   pgs_name = "ltl_pgs"
# )

# write_tsv(ltl_eid_pgs, file = "ltl_eid_pgs.txt")

ltl_eid_pgs <- read_tsv("ltl_eid_pgs.txt", show_col_types = FALSE)

eid_ltl <- read_tsv("eid_ltl_participant.tsv", show_col_types = FALSE)

eid_ltl <- eid_ltl %>% rename(ltl = p22192_i0)

eid_pgs <- read_tsv("eid_pgs.txt", show_col_types = FALSE)

eid_ltl_pgs <- eid_ltl %>%
  inner_join(eid_pgs, by = "eid")

eid_ltl2_pgs <- eid_ltl_pgs %>%
  inner_join(ltl_eid_pgs, by = "eid")

control_stats <- eid_pgs %>%
  summarise(
    control_mean = mean(pgs, na.rm = TRUE),
    control_sd = sd(pgs, na.rm = TRUE)
  )

analysis_df <- eid_ltl2_pgs %>%
  mutate(
    pgs_std = (pgs - control_stats$control_mean) / control_stats$control_sd
  ) %>%
  inner_join(eid_cov_model, by = "eid") %>%
  select(
    eid, ltl, ltl_pgs, pgs_std, age, sex,
    PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10
  ) %>%
  filter(
    !is.na(ltl),
    !is.na(ltl_pgs),
    !is.na(pgs_std),
    !is.na(age),
    !is.na(sex),
    !if_any(starts_with("PC"), is.na)
  ) %>%
  mutate(
    sex = factor(sex)
  )

eid_ch <- read_tsv("ch_binary_participant.tsv")

eid_ch <- rename(eid_ch, ch = p30105)

eid_ch <- eid_ch %>%
    filter(!is.na(ch))

analysis_df <- analysis_df %>%
  inner_join(eid_ch, by = "eid")

analysis_df$ch_bin <- ifelse(is.na(analysis_df$ch), NA_integer_,
                     ifelse(analysis_df$ch >= 1, 1L, 0L))

analysis_df_saved_for_later <- analysis_df

model1 <- analysis_df %>%
    lm(ltl ~ pgs_std + age + I(age^2) + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = .)

model2 <- analysis_df %>%
    lm(ltl ~ pgs_std + ltl_pgs + age + I(age^2) + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = .)

model3 <- analysis_df %>%
    lm(ltl ~ pgs_std + ltl_pgs + age + I(age^2) + sex + ch_bin + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = .)

model4 <- analysis_df %>%
    filter(ch_bin == 0) %>%
    lm(ltl ~ pgs_std + ltl_pgs + age + I(age^2) + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = .)

model5 <- analysis_df %>%
    filter(ch_bin == 0) %>%
    lm(ltl ~ pgs_std + age + I(age^2) + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = .)

library(broom)

tidy(model1, conf.int = TRUE)
tidy(model2, conf.int = TRUE)
tidy(model3, conf.int = TRUE)
tidy(model4, conf.int = TRUE)
tidy(model5, conf.int = TRUE)