

library(readxl)
library(dplyr)
library(lmtest)
library(sandwich)
library(modelsummary)
library(AER)

# 0. Load data

df <- read_excel(
  "MASTER_PANEL_DATA_GLOBAL_RISK_corrected_finstress2_flagged.xlsx"
) %>%
  rename_with(tolower) %>%
  filter(!is.na(usd_share_pct))

# 1. Variable definitions

dv      <- "usd_share_pct"
trade   <- "log_us_trade_gdp"
size    <- "size_ratio_us_to_holder"
equity  <- "log_equity_turnover_gdp"
stress  <- "fin_stress_3part_raw_flexible"
nato    <- "nato_member_dummy"
secexp  <- "security_exposure_no_nato_flexible"
trump   <- "trump_2025"
gpr     <- "log_gpr_country"

# Derived interaction columns
df <- df %>%
  mutate(
    trump_x_secexp = trump_2025 * security_exposure_no_nato_flexible,
    trump_x_gpr    = trump_2025 * log_gpr_country
  )

# 2. IV sample — drop countries missing GPR
# Ireland, Morocco, Serbia, Uruguay have no GPR data
df_iv <- df %>%
  filter(!is.na(log_gpr_country))

cat(sprintf("IV sample: %d obs, %d countries\n",
            nrow(df_iv), n_distinct(df_iv$country)))
cat("Countries excluded (no GPR):",
    paste(setdiff(unique(df$country), unique(df_iv$country)), collapse = ", "), "\n\n")

# 3. Helper functions

# OLS with clustered SEs
run_ols <- function(formula, data, cluster_var = "country") {
  fit     <- lm(formula, data = data)
  cl_vcov <- vcovCL(fit, cluster = data[[cluster_var]], type = "HC1")
  list(fit = fit, cl_vcov = cl_vcov,
       n_clusters = length(unique(data[[cluster_var]])),
       n_obs = nobs(fit))
}

# 2SLS via ivreg with clustered SEs
run_iv <- function(formula, data, cluster_var = "country") {
  fit     <- ivreg(formula, data = data)
  cl_vcov <- vcovCL(fit, cluster = data[[cluster_var]], type = "HC1")
  list(fit = fit, cl_vcov = cl_vcov,
       n_clusters = length(unique(data[[cluster_var]])),
       n_obs = nobs(fit))
}

# First-stage F statistic (clustered Wald test)
first_stage_F <- function(endog, instrument, controls, data) {
  f_fs <- as.formula(paste(
    endog, "~", instrument, "+",
    paste(controls, collapse = " + ")
  ))
  fit_fs  <- lm(f_fs, data = data)
  cl_vcov <- vcovCL(fit_fs, cluster = data[["country"]], type = "HC1")
  # Wald test on instrument coefficient
  wt <- coeftest(fit_fs, vcov = cl_vcov)
  t_stat <- wt[instrument, "t value"]
  cat(sprintf("  First-stage: coef=%.4f, t=%.3f, F=%.2f\n",
              coef(fit_fs)[instrument],
              t_stat, t_stat^2))
  return(t_stat^2)
}

# 4. Mercury + Mars controls (excluding endogenous variable)
controls_base <- c(trade, size, equity, nato, trump)
controls_stress <- c(trade, size, equity, stress, nato, trump)

# 5. First-stage diagnostics

cat("=== FIRST-STAGE DIAGNOSTICS ===\n")
cat("\nIV1/IV2 — instrument: log_gpr_country -> security_exposure_no_nato_flexible\n")
cat("Pooled:\n")
F_pool <- first_stage_F(secexp, gpr, controls_stress, df_iv)
cat(sprintf("  Pooled first-stage F: %.2f %s\n",
            F_pool, ifelse(F_pool >= 10, "(above threshold)", "(BELOW threshold)")))

cat("\nFE (country dummies added):\n")
df_iv_fe <- df_iv %>% mutate(country_f = as.factor(country))
controls_fe <- c("country_f", trade, size, equity, stress, trump)
F_fe <- first_stage_F(secexp, gpr, controls_fe, df_iv_fe)
cat(sprintf("  FE first-stage F: %.2f %s\n\n",
            F_fe, ifelse(F_fe >= 10, "(above threshold)", "(BELOW threshold — FE absorbs GPR variation)")))

cat("IV3 — instruments: log_gpr_country + trump_x_gpr -> secexp + trump_x_secexp\n")
cat("Pooled first stage for secexp:\n")
F_pool_int1 <- first_stage_F(secexp, gpr,
                             c(controls_stress, "trump_x_gpr"), df_iv)
cat("Pooled first stage for trump_x_secexp:\n")
F_pool_int2 <- first_stage_F("trump_x_secexp", "trump_x_gpr",
                             c(controls_stress, gpr), df_iv)
cat("\n")

# 6. Model formulas

# OLS counterparts (on IV sample for direct comparison)
f_OLS1 <- as.formula(paste(dv, "~", trade, "+", size, "+", equity, "+",
                           nato, "+", secexp, "+", trump))

f_OLS2 <- as.formula(paste(dv, "~", trade, "+", size, "+", equity, "+",
                           stress, "+", nato, "+", secexp, "+", trump))

f_OLS3 <- as.formula(paste(dv, "~", trade, "+", size, "+", equity, "+",
                           stress, "+", nato, "+", secexp, "+",
                           "trump_x_secexp", "+", trump))

# IV formulas using ivreg two-part syntax:
# y ~ exogenous + endogenous | exogenous + instruments
f_IV1 <- as.formula(paste(
  dv, "~", trade, "+", size, "+", equity, "+", nato, "+", secexp, "+", trump,
  "|",
  trade, "+", size, "+", equity, "+", nato, "+", gpr, "+", trump
))

f_IV2 <- as.formula(paste(
  dv, "~", trade, "+", size, "+", equity, "+", stress, "+",
  nato, "+", secexp, "+", trump,
  "|",
  trade, "+", size, "+", equity, "+", stress, "+",
  nato, "+", gpr, "+", trump
))

f_IV3 <- as.formula(paste(
  dv, "~", trade, "+", size, "+", equity, "+", stress, "+",
  nato, "+", secexp, "+", "trump_x_secexp", "+", trump,
  "|",
  trade, "+", size, "+", equity, "+", stress, "+",
  nato, "+", gpr, "+", "trump_x_gpr", "+", trump
))

# FE IV
f_IV1_FE <- as.formula(paste(
  dv, "~", "country_f", "+", trade, "+", size, "+", equity, "+",
  stress, "+", secexp, "+", trump,
  "|",
  "country_f", "+", trade, "+", size, "+", equity, "+",
  stress, "+", gpr, "+", trump
))

# 7. Run models

cat("=== RUNNING MODELS ===\n")

OLS1    <- run_ols(f_OLS1, df_iv)
OLS2    <- run_ols(f_OLS2, df_iv)
OLS3    <- run_ols(f_OLS3, df_iv)

IV1     <- run_iv(f_IV1,    df_iv)
IV2     <- run_iv(f_IV2,    df_iv)
IV3     <- run_iv(f_IV3,    df_iv)
IV1_FE  <- run_iv(f_IV1_FE, df_iv_fe)

cat("All models done.\n\n")

# 8. Assemble output table

coef_map <- c(
  "log_us_trade_gdp"                            = "Log US trade/GDP",
  "size_ratio_us_to_holder"                     = "US-to-holder GDP size ratio",
  "log_equity_turnover_gdp"                     = "Log equity turnover/GDP",
  "fin_stress_3part_raw_flexible"               = "Financial stress (3-part raw flexible)",
  "nato_member_dummy"                           = "NATO member",
  "security_exposure_no_nato_flexible"          = "Security exposure (ex-NATO)",
  "trump_x_secexp"                              = "Trump 2025 x Security exposure",
  "trump_2025"                                  = "Trump 2025"
)

gof_map <- data.frame(
  raw   = c("nobs", "r.squared", "adj.r.squared"),
  clean = c("N", "R2", "Adj. R2"),
  fmt   = c(0, 3, 3)
)

# Build first-stage F and cluster count rows
fs_F <- c(round(F_pool, 2), round(F_pool, 2),
          round(F_pool_int1, 2), round(F_fe, 2))
nc   <- as.integer(sapply(list(OLS1, IV1, OLS2, IV2, OLS3, IV3, IV1_FE),
                          `[[`, "n_clusters"))

rows_iv <- data.frame(
  term         = c("Estimator", "First-stage F (sec. exp.)", "Country FE", "N clusters"),
  "OLS1"       = c("OLS",  "-",              "No",  as.character(nc[1])),
  "IV1"        = c("2SLS", as.character(round(F_pool, 2)), "No", as.character(nc[2])),
  "OLS2"       = c("OLS",  "-",              "No",  as.character(nc[3])),
  "IV2"        = c("2SLS", as.character(round(F_pool, 2)), "No", as.character(nc[4])),
  "OLS3"       = c("OLS",  "-",              "No",  as.character(nc[5])),
  "IV3"        = c("2SLS", as.character(round(F_pool_int1, 2)), "No", as.character(nc[6])),
  "IV1 FE"     = c("2SLS", as.character(round(F_fe, 2)), "Yes", as.character(nc[7])),
  check.names  = FALSE
)

fits_all <- list(
  "OLS1"   = OLS1$fit,
  "IV1"    = IV1$fit,
  "OLS2"   = OLS2$fit,
  "IV2"    = IV2$fit,
  "OLS3"   = OLS3$fit,
  "IV3"    = IV3$fit,
  "IV1 FE" = IV1_FE$fit
)

vcovs_all <- list(
  OLS1$cl_vcov, IV1$cl_vcov,
  OLS2$cl_vcov, IV2$cl_vcov,
  OLS3$cl_vcov, IV3$cl_vcov,
  IV1_FE$cl_vcov
)

modelsummary(
  fits_all,
  vcov      = vcovs_all,
  coef_map  = coef_map,
  coef_omit = "country_f",
  gof_map   = gof_map,
  add_rows  = rows_iv,
  stars     = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
  title     = "Table 12: IV robustness — Geopolitical Risk Index instrument",
  notes     = paste0(
    "Instrument: log Geopolitical Risk Index (Caldara & Iacoviello 2022). ",
    "Clustered SEs by country (HC1). Raw non-imputed data throughout. ",
    "IV1-FE first-stage F collapses due to near-zero within-country GPR variation — ",
    "FE IV results are uninformative and reported for illustration only. ",
    "Exclusion restriction is contestable: GPR may affect reserves through channels ",
    "beyond security exposure. Results treated as exploratory."
  ),
  output    = "Table12_IV_models.html"
)

cat("Done. Output written to: Table12_IV_models.html\n")