# Commands to run in the Terminal:
# dx download sffdr_ch_gwas/pgs_165.txt
# dx download sffdr_ch_gwas/eid_pgs.txt
# dx download sffdr_ch_gwas/age_sex_10pcs_dor_participant.tsv
# dx download sffdr_ch_gwas/diagnosis_df_mds.txt
# dx download sffdr_ch_gwas/diagnosis_df_all_mds.txt
# dx download sffdr_ch_gwas/ch_binary_participant.tsv

required_pkgs <- c("tidyverse", "lubridate")

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
  library(lubridate)
  library(ukbrapR)
})

mds_icd10 <- c("C931", "C946", "D46", "D460", "D461", "D462", "D464", "D465", "D466", "D467", "D469")
mds_icd9  <- c("2061","2384", "2385")

all_icd10 <- c(
  "C92", "C920", "C922", "C923", "C924", "C925", "C926", 
  "C927", "C928", "C929", "C93", "C930", "C931", "C932", 
  "C94", "C940", "C942", "C946", "D45", "D450", "D46", 
  "D460", "D461", "D462", "D464", "D465", "D466", "D467", 
  "D469", "D47", "D470", "D471", "D473", "D474"
)

all_icd9 <- c(
  "205", "2050", "2052", "2053", "2058", "2059", "206", "2060", 
  "2061", "2062", "207", "2070", "2072", "2384", "2385", "2387"
)

make_codes_df <- function(condition, icd10_codes, icd9_codes) {
  bind_rows(
    tibble(condition = condition, vocab_id = "ICD10", code = icd10_codes),
    tibble(condition = condition, vocab_id = "ICD9",  code = icd9_codes)
  )
}

codes_df_mds <- make_codes_df("mds", mds_icd10, mds_icd9)
codes_df_all <- make_codes_df("all", all_icd10, all_icd9)

ukbrapr_paths = data.frame(
  object=c("death","death_cause","hesin","hesin_diag","hesin_oper","gp_clinical","gp_scripts","selfrep_illness","cancer_registry","baseline_dates"),
  path=c(
    "ukbrapr_data/death.tsv",
    "ukbrapr_data/death_cause.tsv",
    "ukbrapr_data/hesin.tsv",
    "ukbrapr_data/hesin_diag.tsv",
    "ukbrapr_data/hesin_oper.tsv",
    "ukbrapr_data/gp_clinical.tsv",
    "ukbrapr_data/gp_scripts.tsv",
    "ukbrapr_data/selfrep_illness.tsv",
    "ukbrapr_data/cancer_registry.tsv",
    "ukbrapr_data/baseline_dates.tsv"
  )
)
ukbrapr_paths

get_diagnosis_df <- function(codes_df) {
  codes_df %>%
    get_diagnoses() %>%
    get_df()
}

diagnosis_df_mds <- get_diagnosis_df(codes_df_mds)
diagnosis_df_all <- get_diagnosis_df(codes_df_all)

# write_tsv(diagnosis_df_mds, file = "diagnosis_df_mds.txt")
# write_tsv(diagnosis_df_all, file = "diagnosis_df_all_mds.txt")

diagnosis_df_mds <- read_tsv("diagnosis_df_mds.txt", show_col_types = FALSE)
diagnosis_df_all <- read_tsv("diagnosis_df_all_mds.txt", show_col_types = FALSE)

diagnosis_df_mds_inci <- diagnosis_df_mds %>%
  select(eid, df, bin, bin_prev) %>%
  filter(bin == 1, bin_prev == 0) %>%
  rename(dxdate = df)

diagnosis_df_control <- diagnosis_df_all %>%
  select(eid, df, bin, bin_prev) %>%
  filter(bin == 0) %>%
  rename(dxdate = df)

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

eid_recdate <- eid_cov %>%
  select(eid, age, recdate)

eid_cov_model <- eid_cov %>%
  select(-recdate)

derive_dxage <- function(df, rec_df) {
  df %>%
    inner_join(rec_df, by = "eid") %>%
    mutate(
      years_diff = floor(time_length(interval(recdate, dxdate), "years")),
      dxage = age + years_diff
    ) %>%
    select(eid, bin, bin_prev, dxage)
}

diagnosis_df_mds_inci <- derive_dxage(diagnosis_df_mds_inci, eid_recdate)
diagnosis_df_control  <- derive_dxage(diagnosis_df_control, eid_recdate)

mds_surv_status <- bind_rows(diagnosis_df_mds_inci, diagnosis_df_control) %>%
  transmute(
    eid,
    mds = bin,
    dxage
  )

pgs_variants <- read_tsv("pgs_165.txt", show_col_types = FALSE)

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

# eid_pgs <- create_pgs(
#   in_file = pgs_weights,
#   out_file = "ch.pgs",
#   pgs_name = "pgs"
# )

# write_tsv(eid_pgs, file = "eid_pgs.txt")

eid_pgs <- read_tsv("eid_pgs.txt", show_col_types = FALSE)

analysis_df <- mds_surv_status %>%
  inner_join(eid_pgs, by = "eid")

control_stats <- eid_pgs %>%
  summarise(
    control_mean = mean(pgs, na.rm = TRUE),
    control_sd = sd(pgs, na.rm = TRUE)
  )

analysis_df <- analysis_df %>%
  mutate(
    pgs_std = (pgs - control_stats$control_mean) / control_stats$control_sd
  ) %>%
  inner_join(eid_cov_model, by = "eid") %>%
  select(
    eid, mds, pgs_std, age, sex,
    PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
    dxage
  ) %>%
  filter(
    !is.na(mds),
    !is.na(pgs_std),
    !is.na(age),
    !is.na(dxage),
    !is.na(sex),
    !if_any(starts_with("PC"), is.na)
  ) %>%
  mutate(
    mds = as.integer(mds),
    sex = factor(sex)
  ) %>%
  filter(dxage > age)

eid_ch <- read_tsv("ch_binary_participant.tsv")

eid_ch <- rename(eid_ch, ch = p30105)

eid_ch <- eid_ch %>%
    filter(!is.na(ch))

analysis_df <- analysis_df %>%
  inner_join(eid_ch, by = "eid")

analysis_df$ch_bin <- ifelse(is.na(analysis_df$ch), NA_integer_,
                     ifelse(analysis_df$ch >= 1, 1L, 0L))

analysis_df_saved_for_later <- analysis_df