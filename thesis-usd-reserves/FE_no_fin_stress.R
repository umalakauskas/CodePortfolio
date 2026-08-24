
# 0. Packages

packages <- c(
  "readxl",
  "dplyr",
  "tibble",
  "modelsummary",
  "fixest"
)

to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install)
}

library(readxl)
library(dplyr)
library(tibble)
library(modelsummary)
library(fixest)

# 1. Load flagged dataset

file_path <- "MASTER_PANEL_DATA_GLOBAL_RISK_corrected_finstress2_flagged.xlsx"

df_all <- read_excel(file_path, sheet = "Sheet1") %>%
  rename_with(tolower) %>%
  mutate(
    country = as.character(country),
    year = as.integer(year)
  ) %>%
  filter(!is.na(usd_share_pct))

# 2. Variable names

dv <- "usd_share_pct"
country_var <- "country"
year_var <- "year"

trade <- "log_us_trade_gdp"
size <- "size_ratio_us_to_holder"
equity <- "log_equity_turnover_gdp"
equity_imp <- "equity_turnover_to_gdp_pct_miss"

secexp <- "security_exposure_no_nato_flexible"
trump <- "trump_2025"

required_cols <- c(
  dv,
  country_var,
  year_var,
  trade,
  size,
  equity,
  equity_imp,
  secexp,
  trump
)

missing_cols <- setdiff(required_cols, names(df_all))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# 3. Build common FE no-financial-stress sample

common_vars <- c(
  dv,
  country_var,
  year_var,
  trade,
  size,
  equity,
  secexp,
  trump
)

df_no_finstress_fe <- df_all %>%
  filter(.data[[equity_imp]] == 0) %>%
  filter(if_all(all_of(common_vars), ~ !is.na(.x)))

cat("\nNo-financial-stress FE sample checks:\n")
cat("N:", nrow(df_no_finstress_fe), "\n")
cat("Clusters:", n_distinct(df_no_finstress_fe[[country_var]]), "\n")
cat("Imputed equity rows:", sum(df_no_finstress_fe[[equity_imp]] == 1, na.rm = TRUE), "\n\n")

# 4. Estimate FE models

f_fe_mercury_no_stress <- as.formula(
  paste(
    dv,
    "~",
    paste(
      c(
        trade,
        size,
        equity,
        trump
      ),
      collapse = " + "
    ),
    "| country"
  )
)

f_fe_mars_no_stress <- as.formula(
  paste(
    dv,
    "~",
    paste(
      c(
        trade,
        size,
        equity,
        secexp,
        trump
      ),
      collapse = " + "
    ),
    "| country"
  )
)

FE_mercury_no_stress <- feols(
  f_fe_mercury_no_stress,
  data = df_no_finstress_fe,
  vcov = ~ country
)

FE_mars_no_stress <- feols(
  f_fe_mars_no_stress,
  data = df_no_finstress_fe,
  vcov = ~ country
)

# 5. Output HTML table only

models_fe_no_finstress <- list(
  "FE Mercury no financial stress" = FE_mercury_no_stress,
  "FE Mars added, no financial stress" = FE_mars_no_stress
)

coef_map_fe_no_finstress <- c(
  "log_us_trade_gdp" = "Log US trade/GDP",
  "size_ratio_us_to_holder" = "US-to-holder GDP size ratio",
  "log_equity_turnover_gdp" = "Log equity turnover/GDP",
  "security_exposure_no_nato_flexible" = "Security exposure (ex-NATO)",
  "trump_2025" = "Trump 2025"
)

extra_rows_fe_no_finstress <- tibble(
  term = c(
    "Country FE",
    "Financial-stress index included",
    "Standalone NATO dummy",
    "Common sample",
    "Equity imputed rows",
    "N clusters"
  ),
  `FE Mercury no financial stress` = c(
    "Yes",
    "No",
    "Excluded",
    "Yes",
    "0",
    n_distinct(df_no_finstress_fe[[country_var]])
  ),
  `FE Mars added, no financial stress` = c(
    "Yes",
    "No",
    "Excluded",
    "Yes",
    "0",
    n_distinct(df_no_finstress_fe[[country_var]])
  )
)

modelsummary(
  models_fe_no_finstress,
  coef_map = coef_map_fe_no_finstress,
  statistic = "({std.error})",
  stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "AIC|BIC|Log.Lik.|RMSE|Std.Errors|FE",
  add_rows = extra_rows_fe_no_finstress,
  notes = "Country fixed effects included. Clustered standard errors by country in parentheses. Both models use a common FE sample and exclude observations where equity turnover relies on imputation. The financial-stress index is excluded from both specifications. The standalone NATO dummy is excluded because it is mostly time-invariant and absorbed by country fixed effects.",
  output = "Extra_FE_OLS_Mars_Mercury_no_finstress.html"
)
