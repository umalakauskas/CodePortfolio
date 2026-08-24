
# 0. Packages

packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "modelsummary",
  "sandwich",
  "lmtest",
  "car"
)

to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install)
}

library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(modelsummary)
library(sandwich)
library(lmtest)
library(car)

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
  dv, country_var, year_var,
  trade, size, equity,
  stress, stress_imp,
  rel_stress, rel_stress_imp,
  equity_imp,
  nato, secexp, trump
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

# 3. Build strict observed-data samples

old_m3_vars <- c(
  dv, trade, size, equity, stress, nato, secexp, trump, country_var, year_var
)

df_old_m3 <- df_all %>%
  filter(if_all(all_of(old_m3_vars), ~ !is.na(.x)))

core_vars_abs <- c(
  dv, trade, size, equity, stress, nato, secexp, trump, country_var, year_var
)

df_core <- df_all %>%
  filter(
    .data[[stress_imp]] == 0,
    .data[[equity_imp]] == 0
  ) %>%
  filter(if_all(all_of(core_vars_abs), ~ !is.na(.x)))

core_vars_rel <- c(
  dv, trade, size, equity, rel_stress, nato, secexp, trump, country_var, year_var
)

df_rel_core <- df_all %>%
  filter(
    .data[[rel_stress_imp]] == 0,
    .data[[equity_imp]] == 0
  ) %>%
  filter(if_all(all_of(core_vars_rel), ~ !is.na(.x)))

sample_check <- tibble(
  item = c(
    "Dependent-variable sample",
    "Old M3 complete-case sample before flag exclusions",
    "Old M3 rows with imputed financial-stress components",
    "Old M3 rows with imputed equity turnover",
    "Strict absolute-stress pooled sample",
    "Strict relative-stress pooled sample"
  ),
  n = c(
    nrow(df_all),
    nrow(df_old_m3),
    sum(df_old_m3[[stress_imp]] == 1, na.rm = TRUE),
    sum(df_old_m3[[equity_imp]] == 1, na.rm = TRUE),
    nrow(df_core),
    nrow(df_rel_core)
  )
)

print(sample_check)
write.csv(sample_check, "pooled_OLS_strict_sample_check.csv", row.names = FALSE)

cat("\nStrict absolute-stress sample checks:\n")
cat("Imputed stress rows:", sum(df_core[[stress_imp]] == 1, na.rm = TRUE), "\n")
cat("Imputed equity rows:", sum(df_core[[equity_imp]] == 1, na.rm = TRUE), "\n")
cat("N:", nrow(df_core), "\n")
cat("Clusters:", n_distinct(df_core[[country_var]]), "\n\n")

# 4. Clustered OLS helper

run_ols <- function(formula, data, cluster_var = "country", id_vars = c("country", "year")) {
  needed_vars <- unique(c(all.vars(formula), cluster_var, id_vars))
  needed_vars <- needed_vars[needed_vars %in% names(data)]
  
  model_data <- data %>%
    select(all_of(needed_vars)) %>%
    filter(if_all(all_of(all.vars(formula)), ~ !is.na(.x))) %>%
    filter(!is.na(.data[[cluster_var]]))
  
  fit <- lm(formula, data = model_data)
  
  cl_vcov <- sandwich::vcovCL(
    fit,
    cluster = model_data[[cluster_var]],
    type = "HC1"
  )
  
  list(
    fit = fit,
    vcov = cl_vcov,
    n_clusters = n_distinct(model_data[[cluster_var]]),
    n_obs = nobs(fit),
    data = model_data
  )
}

make_formula <- function(rhs_terms) {
  as.formula(
    paste(dv, "~", paste(rhs_terms, collapse = " + "))
  )
}

# 5. Core pooled OLS specifications

f_M1 <- make_formula(c(
  trade,
  size,
  equity,
  stress,
  trump
))

f_M2 <- make_formula(c(
  trade,
  size,
  nato,
  secexp,
  trump
))

f_M3 <- make_formula(c(
  trade,
  size,
  equity,
  stress,
  nato,
  secexp,
  trump
))

f_M4 <- make_formula(c(
  trade,
  size,
  equity,
  stress,
  nato,
  secexp,
  trump,
  paste0(trump, ":", nato)
))

f_M5 <- make_formula(c(
  trade,
  size,
  equity,
  stress,
  nato,
  secexp,
  trump,
  paste0(trump, ":", secexp)
))

f_M6 <- make_formula(c(
  trade,
  size,
  equity,
  stress,
  nato,
  secexp,
  trump,
  paste0(trump, ":", stress)
))

M1 <- run_ols(f_M1, df_core)
M2 <- run_ols(f_M2, df_core)
M3 <- run_ols(f_M3, df_core)
M4 <- run_ols(f_M4, df_core)
M5 <- run_ols(f_M5, df_core)
M6 <- run_ols(f_M6, df_core)

models_core <- list(
  "M1 Mercury" = M1$fit,
  "M2 Mars" = M2$fit,
  "M3 Full" = M3$fit,
  "M4 Trump x NATO" = M4$fit,
  "M5 Trump x SecExp" = M5$fit,
  "M6 Trump x Stress" = M6$fit
)

vcovs_core <- list(
  M1$vcov,
  M2$vcov,
  M3$vcov,
  M4$vcov,
  M5$vcov,
  M6$vcov
)

coef_map_core <- c(
  "log_us_trade_gdp" = "Log US trade/GDP",
  "size_ratio_us_to_holder" = "US-to-holder GDP size ratio",
  "log_equity_turnover_gdp" = "Log equity turnover/GDP",
  "fin_stress_3part_raw_flexible" = "Financial stress (3-part raw flexible)",
  "nato_member_dummy" = "NATO member",
  "security_exposure_no_nato_flexible" = "Security exposure (ex-NATO)",
  "trump_2025" = "Trump 2025",
  "trump_2025:nato_member_dummy" = "Trump 2025 x NATO",
  "trump_2025:security_exposure_no_nato_flexible" = "Trump 2025 x Security exposure",
  "trump_2025:fin_stress_3part_raw_flexible" = "Trump 2025 x Financial stress"
)

extra_rows_core <- tibble(
  term = c("Strict observed sample", "Financial stress imputed rows", "Equity imputed rows", "N clusters"),
  `M1 Mercury` = c("Yes", "0", "0", M1$n_clusters),
  `M2 Mars` = c("Yes", "0", "0", M2$n_clusters),
  `M3 Full` = c("Yes", "0", "0", M3$n_clusters),
  `M4 Trump x NATO` = c("Yes", "0", "0", M4$n_clusters),
  `M5 Trump x SecExp` = c("Yes", "0", "0", M5$n_clusters),
  `M6 Trump x Stress` = c("Yes", "0", "0", M6$n_clusters)
)

modelsummary(
  models_core,
  vcov = vcovs_core,
  coef_map = coef_map_core,
  statistic = "({std.error})",
  stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "AIC|BIC|Log.Lik.|F|RMSE|Std.Errors",
  add_rows = extra_rows_core,
  notes = "Clustered standard errors by country (HC1) in parentheses. All core pooled models use a common strict observed-data sample excluding rows where the 3-part flexible financial-stress index or equity-turnover variable relies on imputation.",
  output = "Table7_pooled_OLS_core_STRICT.html"
)

modelsummary(
  models_core,
  vcov = vcovs_core,
  coef_map = coef_map_core,
  statistic = "({std.error})",
  stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "AIC|BIC|Log.Lik.|F|RMSE|Std.Errors",
  add_rows = extra_rows_core,
  notes = "Clustered standard errors by country (HC1) in parentheses. All core pooled models use a common strict observed-data sample excluding rows where the 3-part flexible financial-stress index or equity-turnover variable relies on imputation.",
  output = "Table7_pooled_OLS_core_STRICT.tex"
)

cat("\nTable 7: Core pooled OLS\n")
modelsummary(
  models_core,
  vcov = vcovs_core,
  coef_map = coef_map_core,
  statistic = "({std.error})",
  stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "AIC|BIC|Log.Lik.|F|RMSE|Std.Errors",
  add_rows = extra_rows_core,
  output = "markdown"
)

# 6. Pooled robustness specifications

df_noHK <- df_core %>%
  filter(country != "Hong Kong")

df_noHK_URU <- df_core %>%
  filter(!country %in% c("Hong Kong", "Uruguay"))

df_noMEX <- df_core %>%
  filter(country != "Mexico")

df_noHK_MEX <- df_core %>%
  filter(!country %in% c("Hong Kong", "Mexico"))

df_noSWE <- df_core %>%
  filter(country != "Sweden")

f_M3_rel <- make_formula(c(
  trade,
  size,
  equity,
  rel_stress,
  nato,
  secexp,
  trump
))

M3_full <- run_ols(f_M3, df_core)
M3_noHK <- run_ols(f_M3, df_noHK)
M3_noHK_URU <- run_ols(f_M3, df_noHK_URU)
M3_noMEX <- run_ols(f_M3, df_noMEX)
M3_noHK_MEX <- run_ols(f_M3, df_noHK_MEX)
M3_noSWE <- run_ols(f_M3, df_noSWE)
M3_rel <- run_ols(f_M3_rel, df_rel_core)

models_robust <- list(
  "M3 Full" = M3_full$fit,
  "M3 ex-HK" = M3_noHK$fit,
  "M3 ex-HK/URU" = M3_noHK_URU$fit,
  "M3 ex-MEX" = M3_noMEX$fit,
  "M3 ex-HK/MEX" = M3_noHK_MEX$fit,
  "M3 ex-SWE" = M3_noSWE$fit,
  "M3 Rel Stress" = M3_rel$fit
)

vcovs_robust <- list(
  M3_full$vcov,
  M3_noHK$vcov,
  M3_noHK_URU$vcov,
  M3_noMEX$vcov,
  M3_noHK_MEX$vcov,
  M3_noSWE$vcov,
  M3_rel$vcov
)

coef_map_robust <- c(
  "log_us_trade_gdp" = "Log US trade/GDP",
  "size_ratio_us_to_holder" = "US-to-holder GDP size ratio",
  "log_equity_turnover_gdp" = "Log equity turnover/GDP",
  "fin_stress_3part_raw_flexible" = "Financial stress (3-part raw flexible)",
  "fin_stress_3part_rel_raw_flexible" = "Rel. financial stress (3-part raw flexible)",
  "nato_member_dummy" = "NATO member",
  "security_exposure_no_nato_flexible" = "Security exposure (ex-NATO)",
  "trump_2025" = "Trump 2025"
)

extra_rows_robust <- tibble(
  term = c("Strict observed sample", "Excluded country/countries", "N clusters"),
  `M3 Full` = c("Yes", "None", M3_full$n_clusters),
  `M3 ex-HK` = c("Yes", "Hong Kong", M3_noHK$n_clusters),
  `M3 ex-HK/URU` = c("Yes", "Hong Kong; Uruguay", M3_noHK_URU$n_clusters),
  `M3 ex-MEX` = c("Yes", "Mexico", M3_noMEX$n_clusters),
  `M3 ex-HK/MEX` = c("Yes", "Hong Kong; Mexico", M3_noHK_MEX$n_clusters),
  `M3 ex-SWE` = c("Yes", "Sweden", M3_noSWE$n_clusters),
  `M3 Rel Stress` = c("Yes", "None", M3_rel$n_clusters)
)

modelsummary(
  models_robust,
  vcov = vcovs_robust,
  coef_map = coef_map_robust,
  statistic = "({std.error})",
  stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "AIC|BIC|Log.Lik.|F|RMSE|Std.Errors",
  add_rows = extra_rows_robust,
  notes = "Clustered standard errors by country (HC1) in parentheses. Robustness specifications use strict observed-data samples. The relative-stress model uses the relative financial-stress imputation flag.",
  output = "Table8_pooled_OLS_robustness_STRICT.html"
)

modelsummary(
  models_robust,
  vcov = vcovs_robust,
  coef_map = coef_map_robust,
  statistic = "({std.error})",
  stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "AIC|BIC|Log.Lik.|F|RMSE|Std.Errors",
  add_rows = extra_rows_robust,
  notes = "Clustered standard errors by country (HC1) in parentheses. Robustness specifications use strict observed-data samples. The relative-stress model uses the relative financial-stress imputation flag.",
  output = "Table8_pooled_OLS_robustness_STRICT.tex"
)

cat("\nTable 8: Pooled robustness\n")
modelsummary(
  models_robust,
  vcov = vcovs_robust,
  coef_map = coef_map_robust,
  statistic = "({std.error})",
  stars = c("+" = 0.1, "*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "AIC|BIC|Log.Lik.|F|RMSE|Std.Errors",
  add_rows = extra_rows_robust,
  output = "markdown"
)

# 7. Mars correlation matrix

mars_components <- c(
  "nato_member_dummy",
  "interpreted_security_score_norm",
  "mil_presence_index_flexible",
  "defence_dependence_norm",
  "security_exposure_no_nato_flexible",
  "security_exposure_full"
)

mars_labels <- c(
  "NATO member",
  "Alliance score (norm)",
  "Military presence index",
  "Defence dependence (norm)",
  "Security exposure (ex-NATO)",
  "Security exposure (full)"
)

mars_missing <- setdiff(mars_components, names(df_core))

if (length(mars_missing) > 0) {
  warning(
    paste0(
      "Skipping Mars correlation matrix. Missing columns: ",
      paste(mars_missing, collapse = ", ")
    )
  )
} else {
  cat("\nMars correlation matrix: strict pooled sample\n")
  
  mars_data <- df_core %>%
    select(all_of(mars_components)) %>%
    na.omit()
  
  cor_matrix <- cor(mars_data)
  
  rownames(cor_matrix) <- mars_labels
  colnames(cor_matrix) <- mars_labels
  
  print(round(cor_matrix, 3))
  
  cor_long <- as.data.frame(as.table(cor_matrix))
  colnames(cor_long) <- c("Var1", "Var2", "Correlation")
  
  cor_long$Var1 <- factor(cor_long$Var1, levels = mars_labels)
  cor_long$Var2 <- factor(cor_long$Var2, levels = rev(mars_labels))
  
  p_corr <- ggplot(cor_long, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlation, 2)), size = 3, color = "black") +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#d6604d",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Correlation"
    ) +
    labs(
      title = "Correlation matrix: Mars security components",
      subtitle = "Strict observed-data sample",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid = element_blank()
    )
  
  ggsave(
    "Appendix_corr_matrix_mars_STRICT.png",
    p_corr,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  write.csv(
    round(cor_matrix, 3),
    "Appendix_corr_matrix_mars_STRICT.csv"
  )
}

# 8. Manual VIF calculation for strict M3

cat("\nVIFs from M3: strict pooled sample\n")

vif_vars <- c(
  trade,
  size,
  equity,
  stress,
  nato,
  secexp,
  trump
)

vif_labels <- c(
  "Log US trade/GDP",
  "US-to-holder GDP size ratio",
  "Log equity turnover/GDP",
  "Financial stress (3-part raw flexible)",
  "NATO member",
  "Security exposure (ex-NATO)",
  "Trump 2025"
)

m3_data <- model.frame(M3_full$fit)

vif_results <- data.frame(
  Variable = vif_labels,
  VIF = NA_real_
)

for (i in seq_along(vif_vars)) {
  v <- vif_vars[i]
  others <- setdiff(vif_vars, v)
  
  f_aux <- as.formula(
    paste(v, "~", paste(others, collapse = " + "))
  )
  
  r2_aux <- summary(lm(f_aux, data = m3_data))$r.squared
  vif_results$VIF[i] <- round(1 / (1 - r2_aux), 3)
}

print(vif_results)

write.csv(
  vif_results,
  "Appendix_VIF_M3_STRICT.csv",
  row.names = FALSE
)

vif_results$Variable <- factor(
  vif_results$Variable,
  levels = vif_results$Variable[order(vif_results$VIF)]
)

p_vif <- ggplot(vif_results, aes(x = Variable, y = VIF, fill = VIF > 5)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "orange", linewidth = 0.8) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_text(aes(label = round(VIF, 2)), hjust = -0.15, size = 3.2) +
  scale_fill_manual(
    values = c("TRUE" = "#d6604d", "FALSE" = "#4393c3"),
    guide = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  coord_flip() +
  labs(
    title = "Variance Inflation Factors — M3 full pooled OLS",
    subtitle = "Strict observed-data sample; dashed lines at VIF = 5 and VIF = 10",
    x = NULL,
    y = "VIF"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  "Appendix_VIF_M3_STRICT.png",
  p_vif,
  width = 8,
  height = 5,
  dpi = 300
)

# 9. DFBETA diagnostics for trade coefficient

cat("\nDFBETA: influence on trade coefficient, strict M3 sample\n")

dfb <- dfbetas(M3_full$fit)

if (!(trade %in% colnames(dfb))) {
  stop(
    paste0(
      "The trade coefficient column was not found in dfbetas(). Available columns are: ",
      paste(colnames(dfb), collapse = ", ")
    )
  )
}

trade_dfb <- dfb[, trade]

df_influence <- M3_full$data %>%
  mutate(
    dfbeta_trade = as.numeric(trade_dfb),
    abs_dfbeta = abs(dfbeta_trade),
    flag = country %in% c("Hong Kong", "Mexico", "Canada", "Uruguay"),
    label = ifelse(
      flag | abs_dfbeta > 0.2,
      paste0(country, " ", year),
      ""
    )
  )

n <- nrow(df_influence)
thr <- 2 / sqrt(n)

cat(sprintf("Conventional DFBETA threshold: 2/sqrt(%d) = %.3f\n", n, thr))

top_infl <- df_influence %>%
  arrange(desc(abs_dfbeta)) %>%
  select(country, year, dfbeta_trade, abs_dfbeta) %>%
  head(15)

cat("Top 15 influential observations on trade coefficient:\n")
print(top_infl)

write.csv(
  df_influence %>%
    arrange(desc(abs_dfbeta)) %>%
    select(country, year, dfbeta_trade, abs_dfbeta),
  "Appendix_DFBETA_trade_observations_STRICT.csv",
  row.names = FALSE
)

df_influence <- df_influence %>%
  mutate(obs_index = row_number())

p_dfb <- ggplot(
  df_influence,
  aes(
    x = obs_index,
    y = dfbeta_trade,
    color = flag,
    label = label
  )
) +
  geom_hline(yintercept = thr, linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_hline(yintercept = -thr, linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_text(
    data = df_influence %>% filter(label != ""),
    aes(label = label),
    size = 2.6,
    vjust = -0.6,
    hjust = 0.5
  ) +
  scale_color_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "#d6604d"),
    guide = "none"
  ) +
  labs(
    title = "DFBETA: influence of each observation on the trade coefficient",
    subtitle = paste0(
      "Strict observed-data sample. Red dashed lines = conventional threshold ±",
      round(thr, 3),
      " (2/sqrt(N)). Red points = HK, Mexico, Canada, Uruguay."
    ),
    x = "Observation index",
    y = "DFBETA (trade)"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  "Appendix_DFBETA_trade_STRICT.png",
  p_dfb,
  width = 10,
  height = 5,
  dpi = 300
)

# 10. Country-level DFBETA summary for trade coefficient

dfbeta_country <- df_influence %>%
  group_by(country) %>%
  summarise(
    n_obs = n(),
    mean_dfbeta = round(mean(dfbeta_trade), 4),
    max_abs_dfbeta = round(max(abs_dfbeta), 4),
    n_exceed_thr = sum(abs_dfbeta > thr),
    net_direction = ifelse(mean(dfbeta_trade) > 0, "positive", "negative"),
    .groups = "drop"
  ) %>%
  arrange(desc(max_abs_dfbeta))

cat("\nCountry-level DFBETA summary:\n")
print(dfbeta_country, n = Inf)

write.csv(
  dfbeta_country,
  "Appendix_DFBETA_by_country_trade_STRICT.csv",
  row.names = FALSE
)

dfbeta_range <- df_influence %>%
  group_by(country) %>%
  summarise(
    mean_dfb = mean(dfbeta_trade),
    min_dfb = min(dfbeta_trade),
    max_dfb = max(dfbeta_trade),
    .groups = "drop"
  ) %>%
  arrange(mean_dfb) %>%
  mutate(country = factor(country, levels = country))

p_country <- ggplot(
  dfbeta_range,
  aes(
    x = country,
    y = mean_dfb,
    ymin = min_dfb,
    ymax = max_dfb,
    color = mean_dfb > 0
  )
) +
  geom_hline(yintercept = thr, linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_hline(yintercept = -thr, linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_linerange(linewidth = 0.7, alpha = 0.6) +
  geom_point(size = 2.2) +
  scale_color_manual(
    values = c("TRUE" = "#d6604d", "FALSE" = "#4393c3"),
    guide = "none"
  ) +
  coord_flip() +
  labs(
    title = "DFBETA by country: influence on trade coefficient",
    subtitle = paste0(
      "Strict observed-data sample. Dot = country mean DFBETA. Line = min/max range. Red dashed = ±",
      round(thr, 3),
      "."
    ),
    x = NULL,
    y = "DFBETA (trade)"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(
  "Appendix_DFBETA_by_country_STRICT.png",
  p_country,
  width = 9,
  height = 8,
  dpi = 300
)

exceed <- df_influence %>%
  filter(abs_dfbeta > thr) %>%
  select(country, year, dfbeta_trade, abs_dfbeta) %>%
  arrange(desc(abs_dfbeta))

cat("\nAll observations exceeding DFBETA threshold:\n")
print(exceed)

write.csv(
  exceed,
  "Appendix_DFBETA_exceed_threshold_STRICT.csv",
  row.names = FALSE
)

# 11. Save model objects

saveRDS(
  list(
    sample_check = sample_check,
    df_core = df_core,
    df_rel_core = df_rel_core,
    models_core = list(
      M1 = M1,
      M2 = M2,
      M3 = M3,
      M4 = M4,
      M5 = M5,
      M6 = M6
    ),
    models_robust = list(
      M3_full = M3_full,
      M3_noHK = M3_noHK,
      M3_noHK_URU = M3_noHK_URU,
      M3_noMEX = M3_noMEX,
      M3_noHK_MEX = M3_noHK_MEX,
      M3_noSWE = M3_noSWE,
      M3_rel = M3_rel
    ),
    vif_results = vif_results,
    dfbeta_country = dfbeta_country,
    df_influence = df_influence
  ),
  "pooled_OLS_STRICT_results.rds"
)

cat("\nDONE: Strict pooled OLS analysis completed.\n")
cat("Key outputs:\n")
cat("- pooled_OLS_strict_sample_check.csv\n")
cat("- Table7_pooled_OLS_core_STRICT.html\n")
cat("- Table7_pooled_OLS_core_STRICT.tex\n")
cat("- Table8_pooled_OLS_robustness_STRICT.html\n")
cat("- Table8_pooled_OLS_robustness_STRICT.tex\n")
cat("- Appendix_corr_matrix_mars_STRICT.csv / .png\n")
cat("- Appendix_VIF_M3_STRICT.csv / .png\n")
cat("- Appendix_DFBETA_trade_observations_STRICT.csv\n")
cat("- Appendix_DFBETA_trade_STRICT.png\n")
cat("- Appendix_DFBETA_by_country_trade_STRICT.csv\n")
cat("- Appendix_DFBETA_by_country_STRICT.png\n")
cat("- Appendix_DFBETA_exceed_threshold_STRICT.csv\n")
cat("- pooled_OLS_STRICT_results.rds\n")