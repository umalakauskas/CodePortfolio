

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

stress <- "fin_stress_3part_raw_flexible"
stress_imp <- "fin_stress_3part_raw_flexible_any_imp"

rel_stress <- "fin_stress_3part_rel_raw_flexible"
rel_stress_imp <- "fin_stress_3part_rel_raw_flexible_any_imp"

equity_imp <- "equity_turnover_to_gdp_pct_miss"

nato <- "nato_member_dummy"
secexp <- "security_exposure_no_nato_flexible"
trump <- "trump_2025"

required_cols <- c(
  dv,
  country_var,
  year_var,
  trade,
  size,
  equity,
  stress,
  stress_imp,
  rel_stress,
  rel_stress_imp,
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

# 3. Build strict Mundlak sample with financial stress

core_vars_strict <- c(
  dv,
  country_var,
  year_var,
  trade,
  size,
  equity,
  stress,
  nato,
  secexp,
  trump
)

df_mundlak_strict <- df_all %>%
  filter(
    .data[[stress_imp]] == 0,
    .data[[equity_imp]] == 0
  ) %>%
  filter(if_all(all_of(core_vars_strict), ~ !is.na(.x)))

cat("\nStrict Mundlak sample checks:\n")
cat("N:", nrow(df_mundlak_strict), "\n")
cat("Clusters:", n_distinct(df_mundlak_strict[[country_var]]), "\n")
cat("Imputed stress rows:", sum(df_mundlak_strict[[stress_imp]] == 1, na.rm = TRUE), "\n")
cat("Imputed equity rows:", sum(df_mundlak_strict[[equity_imp]] == 1, na.rm = TRUE), "\n\n")

sample_check_strict <- tibble(
  item = c(
    "Strict Mundlak sample",
    "N clusters",
    "Financial stress imputed rows",
    "Equity imputed rows"
  ),
  n = c(
    nrow(df_mundlak_strict),
    n_distinct(df_mundlak_strict[[country_var]]),
    sum(df_mundlak_strict[[stress_imp]] == 1, na.rm = TRUE),
    sum(df_mundlak_strict[[equity_imp]] == 1, na.rm = TRUE)
  )
)

write.csv(
  sample_check_strict,
  "Mundlak_STRICT_sample_check.csv",
  row.names = FALSE
)

# 4. Create Mundlak within and country-mean variables

make_mundlak_data <- function(data, vars, group_var = "country") {
  data %>%
    group_by(.data[[group_var]]) %>%
    mutate(
      across(
        all_of(vars),
        list(
          cm = ~ mean(.x, na.rm = TRUE),
          w = ~ .x - mean(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    ungroup()
}

mundlak_vars_strict <- c(
  trade,
  size,
  equity,
  stress,
  secexp,
  nato
)

df_mundlak_strict <- make_mundlak_data(
  data = df_mundlak_strict,
  vars = mundlak_vars_strict,
  group_var = country_var
)

# 5. Main strict Mundlak / CRE models

f_MU1 <- as.formula(
  paste(
    dv,
    "~",
    paste(
      c(
        paste0(trade, "_w"),
        paste0(trade, "_cm"),
        paste0(size, "_w"),
        paste0(size, "_cm"),
        paste0(equity, "_w"),
        paste0(equity, "_cm"),
        paste0(stress, "_w"),
        paste0(stress, "_cm"),
        paste0(secexp, "_w"),
        paste0(secexp, "_cm"),
        paste0(nato, "_w"),
        paste0(nato, "_cm"),
        trump
      ),
      collapse = " + "
    )
  )
)

f_MU2 <- as.formula(
  paste(
    dv,
    "~",
    paste(
      c(
        paste0(trade, "_w"),
        paste0(trade, "_cm"),
        paste0(size, "_w"),
        paste0(size, "_cm"),
        paste0(equity, "_w"),
        paste0(equity, "_cm"),
        paste0(stress, "_w"),
        paste0(stress, "_cm"),
        paste0(secexp, "_w"),
        paste0(secexp, "_cm"),
        paste0(nato, "_w"),
        paste0(nato, "_cm"),
        trump,
        paste0(trump, ":", nato)
      ),
      collapse = " + "
    )
  )
)

f_MU3 <- as.formula(
  paste(
    dv,
    "~",
    paste(
      c(
        paste0(trade, "_w"),
        paste0(trade, "_cm"),
        paste0(size, "_w"),
        paste0(size, "_cm"),
        paste0(equity, "_w"),
        paste0(equity, "_cm"),
        paste0(stress, "_w"),
        paste0(stress, "_cm"),
        paste0(secexp, "_w"),
        paste0(secexp, "_cm"),
        paste0(nato, "_w"),
        paste0(nato, "_cm"),
        trump,
        paste0(trump, ":", secexp)
      ),
      collapse = " + "
    )
  )
)

f_MU4 <- as.formula(
  paste(
    dv,
    "~",
    paste(
      c(
        paste0(trade, "_w"),
        paste0(trade, "_cm"),
        paste0(size, "_w"),
        paste0(size, "_cm"),
        paste0(equity, "_w"),
        paste0(equity, "_cm"),
        paste0(stress, "_w"),
        paste0(stress, "_cm"),
        paste0(secexp, "_w"),
        paste0(secexp, "_cm"),
        paste0(nato, "_w"),
        paste0(nato, "_cm"),
        trump,
        paste0(trump, ":", stress)
      ),
      collapse = " + "
    )
  )
)

MU1 <- feols(
  f_MU1,
  data = df_mundlak_strict,
  vcov = ~ country
)

MU2 <- feols(
  f_MU2,
  data = df_mundlak_strict,
  vcov = ~ country
)

MU3 <- feols(
  f_MU3,
  data = df_mundlak_strict,
  vcov = ~ country
)

MU4 <- feols(
  f_MU4,
  data = df_mundlak_strict,
  vcov = ~ country
)

models_mundlak_strict <- list(
  "MU1 Core" = MU1,
  "MU2 Trump x NATO" = MU2,
  "MU3 Trump x Security" = MU3,
  "MU4 Trump x Stress" = MU4
)

coef_map_mundlak_strict <- c(
  "log_us_trade_gdp_w" = "Log US trade/GDP (within)",
  "log_us_trade_gdp_cm" = "Log US trade/GDP (country mean)",
  "size_ratio_us_to_holder_w" = "US-to-holder GDP size ratio (within)",
  "size_ratio_us_to_holder_cm" = "US-to-holder GDP size ratio (country mean)",
  "log_equity_turnover_gdp_w" = "Log equity turnover/GDP (within)",
  "log_equity_turnover_gdp_cm" = "Log equity turnover/GDP (country mean)",
  "fin_stress_3part_raw_flexible_w" = "Financial stress (within)",
  "fin_stress_3part_raw_flexible_cm" = "Financial stress (country mean)",
  "security_exposure_no_nato_flexible_w" = "Security exposure ex-NATO (within)",
  "security_exposure_no_nato_flexible_cm" = "Security exposure ex-NATO (country mean)",
  "nato_member_dummy_w" = "NATO member (within)",
  "nato_member_dummy_cm" = "NATO member (country mean)",
  "trump_2025" = "Trump 2025",
  "trump_2025:nato_member_dummy" = "Trump 2025 x NATO",
  "trump_2025:security_exposure_no_nato_flexible" = "Trump 2025 x Security exposure",
  "trump_2025:fin_stress_3part_raw_flexible" = "Trump 2025 x Financial stress"
)

extra_rows_mundlak_strict <- tibble(
  term = c(
    "Estimator",
    "Country FE",
    "Financial-stress index included",
    "Strict observed sample",
    "Financial stress imputed rows",
    "Equity imputed rows",
    "N clusters"
  ),
  `MU1 Core` = c(
    "Mundlak / CRE",
    "No",
    "Yes",
    "Yes",
    "0",
    "0",
    n_distinct(df_mundlak_strict[[country_var]])
  ),
  `MU2 Trump x NATO` = c(
    "Mundlak / CRE",
    "No",
    "Yes",
    "Yes",
    "0",
    "0",
    n_distinct(df_mundlak_strict[[country_var]])
  ),
  `MU3 Trump x Security` = c(
    "Mundlak / CRE",
    "No",
    "Yes",
    "Yes",
    "0",
    "0",
    n_distinct(df_mundlak_strict[[country_var]])
  ),
  `MU4 Trump x Stress` = c(
    "Mundlak / CRE",
    "No",
    "Yes",
    "Yes",
    "0",
    "0",
    n_distinct(df_mundlak_strict[[country_var]])
  )
)

modelsummary(
  models_mundlak_strict,
  coef_map = coef_map_mundlak_strict,
  statistic = "({std.error})",
  stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "AIC|BIC|Log.Lik.|RMSE|Std.Errors|FE",
  add_rows = extra_rows_mundlak_strict,
  notes = "Clustered standard errors by country in parentheses. Within variables are country-demeaned values. Country-mean variables capture between-country differences. The strict observed-data sample excludes rows where the 3-part flexible financial-stress index or equity-turnover variable relies on imputation. Trump interactions use the observed/raw exposure variable in the interaction term.",
  output = "Table11_Mundlak_STRICT.html"
)

# 6. Build Mundlak no-financial-stress sample

core_vars_no_finstress <- c(
  dv,
  country_var,
  year_var,
  trade,
  size,
  equity,
  nato,
  secexp,
  trump
)

df_mundlak_no_finstress <- df_all %>%
  filter(.data[[equity_imp]] == 0) %>%
  filter(if_all(all_of(core_vars_no_finstress), ~ !is.na(.x)))

cat("\nNo-financial-stress Mundlak sample checks:\n")
cat("N:", nrow(df_mundlak_no_finstress), "\n")
cat("Clusters:", n_distinct(df_mundlak_no_finstress[[country_var]]), "\n")
cat("Imputed equity rows:", sum(df_mundlak_no_finstress[[equity_imp]] == 1, na.rm = TRUE), "\n\n")

sample_check_no_finstress <- tibble(
  item = c(
    "No-financial-stress Mundlak sample",
    "N clusters",
    "Equity imputed rows"
  ),
  n = c(
    nrow(df_mundlak_no_finstress),
    n_distinct(df_mundlak_no_finstress[[country_var]]),
    sum(df_mundlak_no_finstress[[equity_imp]] == 1, na.rm = TRUE)
  )
)

write.csv(
  sample_check_no_finstress,
  "Mundlak_no_finstress_sample_check.csv",
  row.names = FALSE
)

mundlak_vars_no_finstress <- c(
  trade,
  size,
  equity,
  secexp,
  nato
)

df_mundlak_no_finstress <- make_mundlak_data(
  data = df_mundlak_no_finstress,
  vars = mundlak_vars_no_finstress,
  group_var = country_var
)

# 7. Extra Mundlak models: Mercury and Mars without financial stress

f_MU_mercury_no_stress <- as.formula(
  paste(
    dv,
    "~",
    paste(
      c(
        paste0(trade, "_w"),
        paste0(trade, "_cm"),
        paste0(size, "_w"),
        paste0(size, "_cm"),
        paste0(equity, "_w"),
        paste0(equity, "_cm"),
        trump
      ),
      collapse = " + "
    )
  )
)

f_MU_mars_no_stress <- as.formula(
  paste(
    dv,
    "~",
    paste(
      c(
        paste0(trade, "_w"),
        paste0(trade, "_cm"),
        paste0(size, "_w"),
        paste0(size, "_cm"),
        paste0(equity, "_w"),
        paste0(equity, "_cm"),
        paste0(secexp, "_w"),
        paste0(secexp, "_cm"),
        paste0(nato, "_w"),
        paste0(nato, "_cm"),
        trump
      ),
      collapse = " + "
    )
  )
)

MU_mercury_no_stress <- feols(
  f_MU_mercury_no_stress,
  data = df_mundlak_no_finstress,
  vcov = ~ country
)

MU_mars_no_stress <- feols(
  f_MU_mars_no_stress,
  data = df_mundlak_no_finstress,
  vcov = ~ country
)

models_mundlak_no_finstress <- list(
  "MU Mercury no financial stress" = MU_mercury_no_stress,
  "MU Mars added, no financial stress" = MU_mars_no_stress
)

coef_map_mundlak_no_finstress <- c(
  "log_us_trade_gdp_w" = "Log US trade/GDP (within)",
  "log_us_trade_gdp_cm" = "Log US trade/GDP (country mean)",
  "size_ratio_us_to_holder_w" = "US-to-holder GDP size ratio (within)",
  "size_ratio_us_to_holder_cm" = "US-to-holder GDP size ratio (country mean)",
  "log_equity_turnover_gdp_w" = "Log equity turnover/GDP (within)",
  "log_equity_turnover_gdp_cm" = "Log equity turnover/GDP (country mean)",
  "security_exposure_no_nato_flexible_w" = "Security exposure ex-NATO (within)",
  "security_exposure_no_nato_flexible_cm" = "Security exposure ex-NATO (country mean)",
  "nato_member_dummy_w" = "NATO member (within)",
  "nato_member_dummy_cm" = "NATO member (country mean)",
  "trump_2025" = "Trump 2025"
)

extra_rows_mundlak_no_finstress <- tibble(
  term = c(
    "Estimator",
    "Country FE",
    "Financial-stress index included",
    "Common sample",
    "Equity imputed rows",
    "N clusters"
  ),
  `MU Mercury no financial stress` = c(
    "Mundlak / CRE",
    "No",
    "No",
    "Yes",
    "0",
    n_distinct(df_mundlak_no_finstress[[country_var]])
  ),
  `MU Mars added, no financial stress` = c(
    "Mundlak / CRE",
    "No",
    "No",
    "Yes",
    "0",
    n_distinct(df_mundlak_no_finstress[[country_var]])
  )
)

modelsummary(
  models_mundlak_no_finstress,
  coef_map = coef_map_mundlak_no_finstress,
  statistic = "({std.error})",
  stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "AIC|BIC|Log.Lik.|RMSE|Std.Errors|FE",
  add_rows = extra_rows_mundlak_no_finstress,
  notes = "Clustered standard errors by country in parentheses. Within variables are country-demeaned values. Country-mean variables capture between-country differences. Both models use a common no-financial-stress Mundlak sample and exclude observations where equity turnover relies on imputation. The financial-stress index is excluded from both specifications.",
  output = "Extra_Mundlak_Mars_Mercury_no_finstress.html"
)

# 8. Save model objects

saveRDS(
  list(
    sample_check_strict = sample_check_strict,
    sample_check_no_finstress = sample_check_no_finstress,
    df_mundlak_strict = df_mundlak_strict,
    df_mundlak_no_finstress = df_mundlak_no_finstress,
    models_mundlak_strict = list(
      MU1 = MU1,
      MU2 = MU2,
      MU3 = MU3,
      MU4 = MU4
    ),
    models_mundlak_no_finstress = list(
      MU_mercury_no_stress = MU_mercury_no_stress,
      MU_mars_no_stress = MU_mars_no_stress
    )
  ),
  "Mundlak_STRICT_results.rds"
)

cat("\nDONE: Mundlak / CRE analysis completed.\n")
cat("Key outputs:\n")
cat("- Table11_Mundlak_STRICT.html\n")
cat("- Extra_Mundlak_Mars_Mercury_no_finstress.html\n")
cat("- Mundlak_STRICT_sample_check.csv\n")
cat("- Mundlak_no_finstress_sample_check.csv\n")
cat("- Mundlak_STRICT_results.rds\n")