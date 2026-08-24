

packages <- c("readxl","dplyr","tibble","modelsummary","fixest","sandwich","lmtest")
to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(readxl); library(dplyr); library(tibble)
library(modelsummary); library(fixest); library(sandwich); library(lmtest)

# 1. Load dataset

file_path <- "MASTER_PANEL_DATA_GLOBAL_RISK_corrected_finstress2_flagged.xlsx"

df_all <- read_excel(file_path, sheet = "Sheet1") %>%
  rename_with(tolower) %>%
  mutate(country = as.character(country), year = as.integer(year)) %>%
  filter(!is.na(usd_share_pct))

# 2. Variable names

dv         <- "usd_share_pct"
country_var <- "country"
year_var    <- "year"
trade      <- "log_us_trade_gdp"
size       <- "size_ratio_us_to_holder"
equity     <- "log_equity_turnover_gdp"
stress     <- "fin_stress_3part_raw_flexible"
stress_imp <- "fin_stress_3part_raw_flexible_any_imp"
rel_stress     <- "fin_stress_3part_rel_raw_flexible"
rel_stress_imp <- "fin_stress_3part_rel_raw_flexible_any_imp"
equity_imp <- "equity_turnover_to_gdp_pct_miss"
nato   <- "nato_member_dummy"
secexp <- "security_exposure_no_nato_flexible"
trump  <- "trump_2025"

required_cols <- c(dv, country_var, year_var, trade, size, equity,
                   stress, stress_imp, rel_stress, rel_stress_imp,
                   equity_imp, nato, secexp, trump)

missing_cols <- setdiff(required_cols, names(df_all))
if (length(missing_cols) > 0)
  stop(paste0("Missing columns: ", paste(missing_cols, collapse = ", ")))

# 3. Build samples

old_m3_vars   <- c(dv, trade, size, equity, stress, secexp, trump, country_var, year_var)
df_old_m3     <- df_all %>% filter(if_all(all_of(old_m3_vars), ~ !is.na(.x)))

core_vars_abs <- c(dv, trade, size, equity, stress, secexp, trump, nato, country_var, year_var)
df_core <- df_all %>%
  filter(.data[[stress_imp]] == 0, .data[[equity_imp]] == 0) %>%
  filter(if_all(all_of(core_vars_abs), ~ !is.na(.x))) %>%
  mutate(country = as.character(country))

core_vars_rel <- c(dv, trade, size, equity, rel_stress, secexp, trump, nato, country_var, year_var)
df_rel_core <- df_all %>%
  filter(.data[[rel_stress_imp]] == 0, .data[[equity_imp]] == 0) %>%
  filter(if_all(all_of(core_vars_rel), ~ !is.na(.x))) %>%
  mutate(country = as.character(country))

sample_check <- tibble(
  item = c("Dependent-variable sample",
           "Old FE complete-case sample before flag exclusions",
           "Old FE rows with imputed financial-stress components",
           "Old FE rows with imputed equity turnover",
           "Strict absolute-stress FE sample",
           "Strict relative-stress FE sample"),
  n = c(nrow(df_all), nrow(df_old_m3),
        sum(df_old_m3[[stress_imp]] == 1, na.rm = TRUE),
        sum(df_old_m3[[equity_imp]] == 1, na.rm = TRUE),
        nrow(df_core), nrow(df_rel_core))
)
print(sample_check)
write.csv(sample_check, "FE_OLS_strict_sample_check.csv", row.names = FALSE)
cat("\nStrict FE sample: N =", nrow(df_core),
    "| clusters =", n_distinct(df_core[[country_var]]),
    "| imputed stress =", sum(df_core[[stress_imp]] == 1, na.rm = TRUE),
    "| imputed equity =", sum(df_core[[equity_imp]] == 1, na.rm = TRUE), "\n\n")

# 4. Model helper
# feols() for estimation (correct within-estimator, works with modelsummary).
# lm() with factor(country) fitted on the SAME data to obtain HC1 vcov.
# The lm vcov is subset to substantive regressors only (dropping country dummies)
# and passed to modelsummary, which overrides feols default SEs.

run_fe <- function(formula_fe, data, cluster_var = "country") {
  
  needed_vars <- unique(c(all.vars(formula_fe), cluster_var))
  model_data  <- data %>%
    select(all_of(needed_vars[needed_vars %in% names(data)])) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  
  # feols fit for correct within-estimator
  fit_fe <- feols(formula_fe, data = model_data,
                  vcov = as.formula(paste0("~", cluster_var)))
  
  # Extract substantive RHS terms (before the | country part)
  fe_str     <- deparse(formula_fe)
  rhs_part   <- sub("\\|.*$", "", sub(".*~", "", fe_str))
  rhs_terms  <- trimws(unlist(strsplit(rhs_part, "\\+")))
  
  # Build equivalent lm formula with factor(country) dummies
  lm_formula <- as.formula(
    paste(dv, "~", paste(c(rhs_terms, "factor(country)"), collapse = " + "))
  )
  
  fit_lm <- lm(lm_formula, data = model_data)
  
  # HC1 vcov from lm
  cl_vcov_full <- sandwich::vcovCL(fit_lm,
                                   cluster = model_data[[cluster_var]],
                                   type = "HC1")
  

  fe_coef_names <- names(coef(fit_fe))
  keep <- intersect(fe_coef_names, rownames(cl_vcov_full))
  cl_vcov <- cl_vcov_full[keep, keep, drop = FALSE]
  
  list(
    fit        = fit_fe,
    vcov       = cl_vcov,
    n_obs      = nobs(fit_fe),
    n_clusters = n_distinct(model_data[[cluster_var]]),
    data       = model_data
  )
}

make_fe_formula <- function(rhs_terms) {
  as.formula(paste(dv, "~", paste(rhs_terms, collapse = " + "), "| country"))
}

# 5. Core FE specifications

f_FE1 <- make_fe_formula(c(trade, size, equity, stress, trump))
f_FE2 <- make_fe_formula(c(trade, size, secexp, trump))
f_FE3 <- make_fe_formula(c(trade, size, equity, stress, secexp, trump))
f_FE4 <- make_fe_formula(c(trade, size, equity, stress, secexp, trump,
                           paste0(trump, ":", secexp)))
f_FE5 <- make_fe_formula(c(trade, size, equity, stress, secexp, trump,
                           paste0(trump, ":", stress)))
f_FE6 <- make_fe_formula(c(trade, size, equity, stress, secexp, trump,
                           paste0(trump, ":", nato)))

FE1 <- run_fe(f_FE1, df_core)
FE2 <- run_fe(f_FE2, df_core)
FE3 <- run_fe(f_FE3, df_core)
FE4 <- run_fe(f_FE4, df_core)
FE5 <- run_fe(f_FE5, df_core)
FE6 <- run_fe(f_FE6, df_core)

# Verify SE sanity on FE3 trade coefficient
cat("FE3 trade SE check:\n")
cat("  feols default SE:", sqrt(FE3$fit$se["log_us_trade_gdp"]), "\n")
cat("  HC1 override SE: ", sqrt(FE3$vcov["log_us_trade_gdp","log_us_trade_gdp"]), "\n\n")

models_fe_core <- list(
  "FE1 Mercury"        = FE1$fit,
  "FE2 Mars"           = FE2$fit,
  "FE3 Full"           = FE3$fit,
  "FE4 Trump x SecExp" = FE4$fit,
  "FE5 Trump x Stress" = FE5$fit,
  "FE6 Trump x NATO"   = FE6$fit
)

vcovs_fe_core <- list(FE1$vcov, FE2$vcov, FE3$vcov, FE4$vcov, FE5$vcov, FE6$vcov)

coef_map_fe <- c(
  "log_us_trade_gdp" = "Log US trade/GDP",
  "size_ratio_us_to_holder" = "US-to-holder GDP size ratio",
  "log_equity_turnover_gdp" = "Log equity turnover/GDP",
  "fin_stress_3part_raw_flexible" = "Financial stress (3-part raw flexible)",
  "security_exposure_no_nato_flexible" = "Security exposure (ex-NATO)",
  "trump_2025" = "Trump 2025",
  
  "security_exposure_no_nato_flexible:trump_2025" = "Trump 2025 x Security exposure",
  "trump_2025:security_exposure_no_nato_flexible" = "Trump 2025 x Security exposure",
  
  "fin_stress_3part_raw_flexible:trump_2025" = "Trump 2025 x Financial stress",
  "trump_2025:fin_stress_3part_raw_flexible" = "Trump 2025 x Financial stress",
  
  "nato_member_dummy:trump_2025" = "Trump 2025 x NATO",
  "trump_2025:nato_member_dummy" = "Trump 2025 x NATO"
)
extra_rows_fe_core <- tibble(
  term = c("Country FE", "Standalone NATO dummy", "Strict observed sample",
           "SE method", "Financial stress imputed rows", "Equity imputed rows", "N clusters"),
  `FE1 Mercury`        = c("Yes","Excluded","Yes","HC1","0","0", FE1$n_clusters),
  `FE2 Mars`           = c("Yes","Excluded","Yes","HC1","0","0", FE2$n_clusters),
  `FE3 Full`           = c("Yes","Excluded","Yes","HC1","0","0", FE3$n_clusters),
  `FE4 Trump x SecExp` = c("Yes","Excluded","Yes","HC1","0","0", FE4$n_clusters),
  `FE5 Trump x Stress` = c("Yes","Excluded","Yes","HC1","0","0", FE5$n_clusters),
  `FE6 Trump x NATO`   = c("Yes","Excluded","Yes","HC1","0","0", FE6$n_clusters)
)

fe_note <- paste(
  "Country fixed effects included. Clustered standard errors (HC1) by country in parentheses,",
  "consistent with pooled OLS and wild cluster bootstrap. All FE core models use a common",
  "strict observed-data sample excluding rows where the 3-part flexible financial-stress",
  "index or equity-turnover variable relies on imputation. The standalone NATO dummy is",
  "excluded because it is mostly time-invariant and absorbed by country fixed effects."
)

modelsummary(models_fe_core, vcov = vcovs_fe_core, coef_map = coef_map_fe,
             statistic = "({std.error})",
             stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
             gof_omit = "AIC|BIC|Log.Lik.|RMSE|Std.Errors|FE",
             add_rows = extra_rows_fe_core, notes = fe_note,
             output = "Table9_FE_core_STRICT.html")

modelsummary(models_fe_core, vcov = vcovs_fe_core, coef_map = coef_map_fe,
             statistic = "({std.error})",
             stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
             gof_omit = "AIC|BIC|Log.Lik.|RMSE|Std.Errors|FE",
             add_rows = extra_rows_fe_core, notes = fe_note,
             output = "Table9_FE_core_STRICT.tex")

cat("\nTable 9: FE core models\n")
modelsummary(models_fe_core, vcov = vcovs_fe_core, coef_map = coef_map_fe,
             statistic = "({std.error})",
             stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
             gof_omit = "AIC|BIC|Log.Lik.|RMSE|Std.Errors|FE",
             add_rows = extra_rows_fe_core, output = "markdown")

# 6. FE robustness specifications

df_noHK      <- df_core %>% filter(country != "Hong Kong")
df_noPOL     <- df_core %>% filter(country != "Poland")
df_noSWE     <- df_core %>% filter(country != "Sweden")
df_noMEX     <- df_core %>% filter(country != "Mexico")
df_noHK_MEX  <- df_core %>% filter(!country %in% c("Hong Kong","Mexico"))

f_FE3_rel <- make_fe_formula(c(trade, size, equity, rel_stress, secexp, trump))

FE3_full     <- run_fe(f_FE3,     df_core)
FE3_noHK     <- run_fe(f_FE3,     df_noHK)
FE3_noPOL    <- run_fe(f_FE3,     df_noPOL)
FE3_noSWE    <- run_fe(f_FE3,     df_noSWE)
FE3_noMEX    <- run_fe(f_FE3,     df_noMEX)
FE3_noHK_MEX <- run_fe(f_FE3,     df_noHK_MEX)
FE3_rel      <- run_fe(f_FE3_rel, df_rel_core)

models_fe_robust <- list(
  "FE3 Full"       = FE3_full$fit,
  "FE3 ex-HK"      = FE3_noHK$fit,
  "FE3 ex-Poland"  = FE3_noPOL$fit,
  "FE3 ex-Sweden"  = FE3_noSWE$fit,
  "FE3 ex-MEX"     = FE3_noMEX$fit,
  "FE3 ex-HK/MEX"  = FE3_noHK_MEX$fit,
  "FE3 Rel Stress" = FE3_rel$fit
)

vcovs_fe_robust <- list(FE3_full$vcov, FE3_noHK$vcov, FE3_noPOL$vcov,
                        FE3_noSWE$vcov, FE3_noMEX$vcov,
                        FE3_noHK_MEX$vcov, FE3_rel$vcov)

coef_map_fe_robust <- c(
  "log_us_trade_gdp"                   = "Log US trade/GDP",
  "size_ratio_us_to_holder"            = "US-to-holder GDP size ratio",
  "log_equity_turnover_gdp"            = "Log equity turnover/GDP",
  "fin_stress_3part_raw_flexible"      = "Financial stress (3-part raw flexible)",
  "fin_stress_3part_rel_raw_flexible"  = "Rel. financial stress (3-part raw flexible)",
  "security_exposure_no_nato_flexible" = "Security exposure (ex-NATO)",
  "trump_2025"                         = "Trump 2025"
)

extra_rows_fe_robust <- tibble(
  term = c("Country FE","Standalone NATO dummy","Strict observed sample",
           "SE method","Excluded country/countries","N clusters"),
  `FE3 Full`       = c("Yes","Excluded","Yes","HC1","None",              FE3_full$n_clusters),
  `FE3 ex-HK`      = c("Yes","Excluded","Yes","HC1","Hong Kong",         FE3_noHK$n_clusters),
  `FE3 ex-Poland`  = c("Yes","Excluded","Yes","HC1","Poland",            FE3_noPOL$n_clusters),
  `FE3 ex-Sweden`  = c("Yes","Excluded","Yes","HC1","Sweden",            FE3_noSWE$n_clusters),
  `FE3 ex-MEX`     = c("Yes","Excluded","Yes","HC1","Mexico",            FE3_noMEX$n_clusters),
  `FE3 ex-HK/MEX`  = c("Yes","Excluded","Yes","HC1","Hong Kong; Mexico", FE3_noHK_MEX$n_clusters),
  `FE3 Rel Stress` = c("Yes","Excluded","Yes","HC1","None",              FE3_rel$n_clusters)
)

robust_note <- paste(
  "Country fixed effects included. Clustered standard errors (HC1) by country in parentheses.",
  "Robustness specifications use strict observed-data samples.",
  "The relative-stress model uses the relative financial-stress imputation flag.",
  "The standalone NATO dummy is excluded because it is mostly time-invariant",
  "and absorbed by country fixed effects."
)

modelsummary(models_fe_robust, vcov = vcovs_fe_robust, coef_map = coef_map_fe_robust,
             statistic = "({std.error})",
             stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
             gof_omit = "AIC|BIC|Log.Lik.|RMSE|Std.Errors|FE",
             add_rows = extra_rows_fe_robust, notes = robust_note,
             output = "Table10_FE_robustness_STRICT.html")

modelsummary(models_fe_robust, vcov = vcovs_fe_robust, coef_map = coef_map_fe_robust,
             statistic = "({std.error})",
             stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
             gof_omit = "AIC|BIC|Log.Lik.|RMSE|Std.Errors|FE",
             add_rows = extra_rows_fe_robust, notes = robust_note,
             output = "Table10_FE_robustness_STRICT.tex")

cat("\nTable 10: FE robustness models\n")
modelsummary(models_fe_robust, vcov = vcovs_fe_robust, coef_map = coef_map_fe_robust,
             statistic = "({std.error})",
             stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
             gof_omit = "AIC|BIC|Log.Lik.|RMSE|Std.Errors|FE",
             add_rows = extra_rows_fe_robust, output = "markdown")

# 7. Within-country variation summary

within_variation_vars <- c(trade, size, equity, stress, rel_stress, secexp, nato, trump)
within_variation <- df_core %>%
  group_by(country) %>%
  summarise(across(all_of(within_variation_vars),
                   list(sd = ~ sd(.x, na.rm=TRUE),
                        min = ~ min(.x, na.rm=TRUE),
                        max = ~ max(.x, na.rm=TRUE)),
                   .names = "{.col}_{.fn}"),
            n_obs = n(), .groups = "drop")

write.csv(within_variation, "FE_within_country_variation_STRICT.csv", row.names = FALSE)
print(within_variation)

# 8. Save model objects

saveRDS(
  list(sample_check = sample_check, df_core = df_core, df_rel_core = df_rel_core,
       models_fe_core = list(FE1=FE1, FE2=FE2, FE3=FE3, FE4=FE4, FE5=FE5, FE6=FE6),
       models_fe_robust = list(FE3_full=FE3_full, FE3_noHK=FE3_noHK,
                               FE3_noPOL=FE3_noPOL, FE3_noSWE=FE3_noSWE,
                               FE3_noMEX=FE3_noMEX, FE3_noHK_MEX=FE3_noHK_MEX,
                               FE3_rel=FE3_rel),
       within_variation = within_variation),
  "FE_OLS_STRICT_results.rds"
)

cat("\nDONE.\n")
cat("Estimation: feols() within-estimator\n")
cat("SE method: sandwich HC1 from lm() with factor(country) on same data\n")
cat("Key outputs:\n")
cat("- FE_OLS_strict_sample_check.csv\n")
cat("- Table9_FE_core_STRICT.html / .tex\n")
cat("- Table10_FE_robustness_STRICT.html / .tex\n")
cat("- FE_within_country_variation_STRICT.csv\n")
cat("- FE_OLS_STRICT_results.rds\n")