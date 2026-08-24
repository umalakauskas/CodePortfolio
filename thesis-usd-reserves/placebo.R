

library(readxl)
library(dplyr)
library(lmtest)
library(sandwich)
library(modelsummary)

# 0. Load data

df <- read_excel(
  "MASTER_PANEL_DATA_GLOBAL_RISK_corrected_finstress2_flagged.xlsx"
) %>%
  rename_with(tolower) %>%
  filter(!is.na(usd_share_pct))

# 1. Variable definitions

dv     <- "usd_share_pct"
trade  <- "log_us_trade_gdp"
size   <- "size_ratio_us_to_holder"
equity <- "log_equity_turnover_gdp"
stress <- "fin_stress_3part_raw_flexible"
nato   <- "nato_member_dummy"
secexp <- "security_exposure_no_nato_flexible"

# 2. Helper functions

run_ols <- function(formula, data, cluster_var = "country") {
  fit     <- lm(formula, data = data)
  cl_vcov <- vcovCL(fit, cluster = data[[cluster_var]], type = "HC1")
  list(fit = fit, cl_vcov = cl_vcov,
       n_clusters = length(unique(data[[cluster_var]])),
       n_obs = nobs(fit))
}

# 3. Build placebo dummies and run models

placebo_years <- c(2013, 2017, 2020, 2023, 2025)

fits_pooled <- list()
vcovs_pooled <- list()
fits_fe <- list()
vcovs_fe <- list()
nc_pooled <- c()
nc_fe     <- c()

for (yr in placebo_years) {
  
  # Build dummy: 1 in the placebo year, 0 otherwise
  df_p <- df %>% mutate(placebo_dummy = as.integer(year == yr))
  
  f_pooled <- as.formula(paste(
    dv, "~", trade, "+", size, "+", equity, "+", stress, "+",
    nato, "+", secexp, "+", "placebo_dummy"
  ))
  
  f_fe <- as.formula(paste(
    dv, "~", "as.factor(country)", "+",
    trade, "+", size, "+", equity, "+", stress, "+",
    secexp, "+", "placebo_dummy"
  ))
  
  m_pool <- run_ols(f_pooled, df_p)
  m_fe   <- run_ols(f_fe, df_p)
  
  label <- ifelse(yr == 2025, "Actual 2025", paste0("Placebo ", yr))
  fits_pooled[[label]]  <- m_pool$fit
  vcovs_pooled[[label]] <- m_pool$cl_vcov
  fits_fe[[label]]      <- m_fe$fit
  vcovs_fe[[label]]     <- m_fe$cl_vcov
  nc_pooled             <- c(nc_pooled, m_pool$n_clusters)
  nc_fe                 <- c(nc_fe, m_fe$n_clusters)
  
  # Print dummy coefficient for quick read
  coef_p <- coeftest(m_pool$fit, vcov = m_pool$cl_vcov)["placebo_dummy", ]
  coef_f <- coeftest(m_fe$fit,   vcov = m_fe$cl_vcov)["placebo_dummy", ]
  cat(sprintf("%-15s  Pooled: %+.3f (SE %.3f, p=%.3f)   FE: %+.3f (SE %.3f, p=%.3f)\n",
              label,
              coef_p["Estimate"], coef_p["Std. Error"], coef_p["Pr(>|t|)"],
              coef_f["Estimate"], coef_f["Std. Error"], coef_f["Pr(>|t|)"]))
}

cat("\n")

# 4. coef_map — rename placebo_dummy to year-specific label
# modelsummary will show it as "placebo_dummy" for all — that's correct
# we distinguish by column name

coef_map_placebo <- c(
  "log_us_trade_gdp"                   = "Log US trade/GDP",
  "size_ratio_us_to_holder"            = "US-to-holder GDP size ratio",
  "log_equity_turnover_gdp"            = "Log equity turnover/GDP",
  "fin_stress_3part_raw_flexible"      = "Financial stress (3-part raw flexible)",
  "nato_member_dummy"                  = "NATO member",
  "security_exposure_no_nato_flexible" = "Security exposure (ex-NATO)",
  "placebo_dummy"                      = "Shock year dummy"
)

gof_map <- data.frame(
  raw   = c("nobs", "r.squared", "adj.r.squared"),
  clean = c("N", "R2", "Adj. R2"),
  fmt   = c(0, 3, 3)
)

# Cluster count rows
rows_pooled <- data.frame(
  term            = "N clusters",
  "Placebo 2013"  = as.character(nc_pooled[1]),
  "Placebo 2017"  = as.character(nc_pooled[2]),
  "Placebo 2020"  = as.character(nc_pooled[3]),
  "Placebo 2023"  = as.character(nc_pooled[4]),
  "Actual 2025"   = as.character(nc_pooled[5]),
  check.names = FALSE
)

rows_fe <- data.frame(
  term            = c("Country FE", "N clusters"),
  "Placebo 2013"  = c("Yes", as.character(nc_fe[1])),
  "Placebo 2017"  = c("Yes", as.character(nc_fe[2])),
  "Placebo 2020"  = c("Yes", as.character(nc_fe[3])),
  "Placebo 2023"  = c("Yes", as.character(nc_fe[4])),
  "Actual 2025"   = c("Yes", as.character(nc_fe[5])),
  check.names = FALSE
)

# 5. Output tables

modelsummary(
  fits_pooled,
  vcov     = vcovs_pooled,
  coef_map = coef_map_placebo,
  gof_map  = gof_map,
  add_rows = rows_pooled,
  stars    = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
  title    = "Table A1: Placebo year tests — Pooled OLS",
  notes    = paste0(
    "Each column replaces trump_2025 with a placebo dummy = 1 in the stated year, 0 otherwise. ",
    "M3-equivalent specification throughout. Clustered SEs by country (HC1). ",
    "Raw non-imputed data. 'Actual 2025' is the real result for comparison."
  ),
  output   = "TableA1_placebo_pooled.html"
)

modelsummary(
  fits_fe,
  vcov      = vcovs_fe,
  coef_map  = coef_map_placebo,
  coef_omit = "as.factor",
  gof_map   = gof_map,
  add_rows  = rows_fe,
  stars     = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
  title     = "Table A2: Placebo year tests — Country FE",
  notes     = paste0(
    "Each column replaces trump_2025 with a placebo dummy = 1 in the stated year, 0 otherwise. ",
    "FE3-equivalent specification. Country dummies not shown. Clustered SEs by country (HC1). ",
    "Raw non-imputed data. 'Actual 2025' is the real result for comparison."
  ),
  output   = "TableA2_placebo_fe.html"
)

cat("Done. Files written:\n")
cat("  TableA1_placebo_pooled.html\n")
cat("  TableA2_placebo_fe.html\n")