
# 0. Packages

packages <- c(
  "readxl",
  "dplyr",
  "tibble",
  "modelsummary",
  "sandwich"
)

to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install)
}

library(readxl)
library(dplyr)
library(tibble)
library(modelsummary)
library(sandwich)

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

trade <- "log_us_trade_gdp"
size <- "size_ratio_us_to_holder"
equity <- "log_equity_turnover_gdp"
equity_imp <- "equity_turnover_to_gdp_pct_miss"

nato <- "nato_member_dummy"
secexp <- "security_exposure_no_nato_flexible"
trump <- "trump_2025"

required_cols <- c(
  dv,
  country_var,
  trade,
  size,
  equity,
  equity_imp,
  nato,
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

# 3. Build common no-financial-stress pooled sample

common_vars <- c(
  dv,
  country_var,
  trade,
  size,
  equity,
  nato,
  secexp,
  trump
)

df_no_finstress <- df_all %>%
  filter(.data[[equity_imp]] == 0) %>%
  filter(if_all(all_of(common_vars), ~ !is.na(.x)))

cat("\nNo-financial-stress pooled sample checks:\n")
cat("N:", nrow(df_no_finstress), "\n")
cat("Clusters:", n_distinct(df_no_finstress[[country_var]]), "\n")
cat("Imputed equity rows:", sum(df_no_finstress[[equity_imp]] == 1, na.rm = TRUE), "\n\n")

# 4. Clustered OLS helper

run_ols <- function(formula, data, cluster_var = "country") {
  needed_vars <- unique(c(all.vars(formula), cluster_var))
  
  model_data <- data %>%
    select(all_of(needed_vars)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  
  fit <- lm(formula, data = model_data)
  
  cl_vcov <- sandwich::vcovCL(
    fit,
    cluster = model_data[[cluster_var]],
    type = "HC1"
  )
  
  list(
    fit = fit,
    vcov = cl_vcov,
    n_obs = nobs(fit),
    n_clusters = n_distinct(model_data[[cluster_var]]),
    data = model_data
  )
}

# 5. Model formulas

f_mercury_no_stress <- as.formula(
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
    )
  )
)

f_mars_no_stress <- as.formula(
  paste(
    dv,
    "~",
    paste(
      c(
        trade,
        size,
        equity,
        nato,
        secexp,
        trump
      ),
      collapse = " + "
    )
  )
)

# 6. Estimate pooled OLS models

M_mercury_no_stress <- run_ols(f_mercury_no_stress, df_no_finstress)
M_mars_no_stress <- run_ols(f_mars_no_stress, df_no_finstress)

models_no_finstress <- list(
  "Mercury no financial stress" = M_mercury_no_stress$fit,
  "Mars no financial stress" = M_mars_no_stress$fit
)

vcovs_no_finstress <- list(
  M_mercury_no_stress$vcov,
  M_mars_no_stress$vcov
)

# 7. Output HTML table only

coef_map_no_finstress <- c(
  "log_us_trade_gdp" = "Log US trade/GDP",
  "size_ratio_us_to_holder" = "US-to-holder GDP size ratio",
  "log_equity_turnover_gdp" = "Log equity turnover/GDP",
  "nato_member_dummy" = "NATO member",
  "security_exposure_no_nato_flexible" = "Security exposure (ex-NATO)",
  "trump_2025" = "Trump 2025"
)

extra_rows_no_finstress <- tibble(
  term = c(
    "Financial-stress index included",
    "Common sample",
    "Equity imputed rows",
    "N clusters"
  ),
  `Mercury no financial stress` = c(
    "No",
    "Yes",
    "0",
    M_mercury_no_stress$n_clusters
  ),
  `Mars no financial stress` = c(
    "No",
    "Yes",
    "0",
    M_mars_no_stress$n_clusters
  )
)

modelsummary(
  models_no_finstress,
  vcov = vcovs_no_finstress,
  coef_map = coef_map_no_finstress,
  statistic = "({std.error})",
  stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "AIC|BIC|Log.Lik.|F|RMSE|Std.Errors",
  add_rows = extra_rows_no_finstress,
  notes = "Clustered standard errors by country (HC1) in parentheses. Both models use a common pooled OLS sample and exclude observations where equity turnover relies on imputation. The financial-stress index is excluded from both specifications.",
  output = "Extra_pooled_OLS_Mars_Mercury_no_finstress.html"
)

