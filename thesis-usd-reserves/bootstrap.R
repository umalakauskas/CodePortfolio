
# wild_bootstrap.R
#
# Wild cluster bootstrap p-values for key coefficients
# Addresses the few-cluster (N=22) problem with standard clustered SEs
# Uses Cameron, Gelbach & Miller (2008) wild bootstrap
#
# Models tested:
#   M3  - Full pooled OLS (trade, size, equity, stress, NATO, secexp, trump)
#   FE3 - Full country FE
#   MU1 - Mundlak full

library(readxl)
library(dplyr)
library(lmtest)
library(sandwich)
install.packages(
  "fwildclusterboot",
  repos = c("https://s3alfisc.r-universe.dev", "https://cloud.r-project.org")
)
library(fwildclusterboot)


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
trump  <- "trump_2025"

# 2. Add Mundlak country means

mundlak_vars <- c(trade, size, equity, stress, secexp)
for (v in mundlak_vars) {
  df[[paste0(v, "_cm")]] <- ave(df[[v]], df$country,
                                FUN = function(x) mean(x, na.rm = TRUE))
}

df$country_f <- as.factor(df$country)

# 3. Model formulas

f_M3 <- as.formula(paste(
  dv, "~", trade, "+", size, "+", equity, "+", stress, "+",
  nato, "+", secexp, "+", trump
))

f_FE3 <- as.formula(paste(
  dv, "~", "country_f +",
  trade, "+", size, "+", equity, "+", stress, "+", secexp, "+", trump
))

cm_terms <- paste(paste0(mundlak_vars, "_cm"), collapse = " + ")
f_MU1 <- as.formula(paste(
  dv, "~", trade, "+", size, "+", equity, "+", stress, "+",
  nato, "+", secexp, "+", trump, "+", cm_terms
))

# 4. Fit models on complete cases
# boottest requires no NAs in the model frame

df_m3  <- df %>% filter(!is.na(.data[[stress]]) & !is.na(.data[[secexp]]))
df_fe3 <- df_m3
df_mu1 <- df_m3

# country must be in the data passed to lm for boottest to find it
fit_M3  <- lm(f_M3,  data = df_m3,  model = TRUE)
fit_FE3 <- lm(f_FE3, data = df_fe3, model = TRUE)
fit_MU1 <- lm(f_MU1, data = df_mu1, model = TRUE)

# Attach country to model frame so boottest can find the cluster variable
fit_M3$model$country  <- df_m3$country[as.integer(rownames(fit_M3$model))]
fit_FE3$model$country <- df_fe3$country[as.integer(rownames(fit_FE3$model))]
fit_MU1$model$country <- df_mu1$country[as.integer(rownames(fit_MU1$model))]

# 5. Standard clustered SEs for comparison

cl_M3  <- vcovCL(fit_M3,  cluster = df_m3$country,  type = "HC1")
cl_FE3 <- vcovCL(fit_FE3, cluster = df_fe3$country, type = "HC1")
cl_MU1 <- vcovCL(fit_MU1, cluster = df_mu1$country, type = "HC1")

cat("=== STANDARD CLUSTERED SEs (HC1) — for comparison ===\n")
cat("\nM3 trade p-value: ",
    coeftest(fit_M3, cl_M3)[trade, "Pr(>|t|)"], "\n")
cat("FE3 trade p-value:",
    coeftest(fit_FE3, cl_FE3)[trade, "Pr(>|t|)"], "\n")
cat("MU1 stress_cm p-value:",
    coeftest(fit_MU1, cl_MU1)[paste0(stress,"_cm"), "Pr(>|t|)"], "\n\n")

# 6. Helper: run boottest
# clustid must be a character scalar (column name in the model data)
# The cluster column must exist in the dataframe used to fit the model

do_boot <- function(fit, param, B = 9999, seed = 42) {
  set.seed(seed)
  boottest(
    fit,
    clustid     = "country",
    param       = param,
    B           = B,
    type        = "rademacher",
    impose_null = TRUE
  )
}

cat("Running wild cluster bootstrap (B=9999, Rademacher)...\n")
cat("This will take a few minutes.\n\n")

# 7. M3 Pooled OLS

cat("--- M3 Pooled OLS ---\n")
b_m3_trade  <- do_boot(fit_M3, trade)
b_m3_stress <- do_boot(fit_M3, stress)
b_m3_secexp <- do_boot(fit_M3, secexp)
b_m3_trump  <- do_boot(fit_M3, trump)

cat(sprintf("  %-45s  coef=%+7.3f  cluster_p=%.3f  bootstrap_p=%.3f\n",
            trade,  coef(fit_M3)[trade],
            coeftest(fit_M3,cl_M3)[trade,"Pr(>|t|)"],  b_m3_trade$p_val))
cat(sprintf("  %-45s  coef=%+7.3f  cluster_p=%.3f  bootstrap_p=%.3f\n",
            stress, coef(fit_M3)[stress],
            coeftest(fit_M3,cl_M3)[stress,"Pr(>|t|)"], b_m3_stress$p_val))
cat(sprintf("  %-45s  coef=%+7.3f  cluster_p=%.3f  bootstrap_p=%.3f\n",
            secexp, coef(fit_M3)[secexp],
            coeftest(fit_M3,cl_M3)[secexp,"Pr(>|t|)"], b_m3_secexp$p_val))
cat(sprintf("  %-45s  coef=%+7.3f  cluster_p=%.3f  bootstrap_p=%.3f\n",
            trump,  coef(fit_M3)[trump],
            coeftest(fit_M3,cl_M3)[trump,"Pr(>|t|)"],  b_m3_trump$p_val))

# 8. FE3 Country FE

cat("\n--- FE3 Country FE ---\n")
b_fe3_trade  <- do_boot(fit_FE3, trade)
b_fe3_stress <- do_boot(fit_FE3, stress)
b_fe3_secexp <- do_boot(fit_FE3, secexp)
b_fe3_trump  <- do_boot(fit_FE3, trump)

cat(sprintf("  %-45s  coef=%+7.3f  cluster_p=%.3f  bootstrap_p=%.3f\n",
            trade,  coef(fit_FE3)[trade],
            coeftest(fit_FE3,cl_FE3)[trade,"Pr(>|t|)"],  b_fe3_trade$p_val))
cat(sprintf("  %-45s  coef=%+7.3f  cluster_p=%.3f  bootstrap_p=%.3f\n",
            stress, coef(fit_FE3)[stress],
            coeftest(fit_FE3,cl_FE3)[stress,"Pr(>|t|)"], b_fe3_stress$p_val))
cat(sprintf("  %-45s  coef=%+7.3f  cluster_p=%.3f  bootstrap_p=%.3f\n",
            secexp, coef(fit_FE3)[secexp],
            coeftest(fit_FE3,cl_FE3)[secexp,"Pr(>|t|)"], b_fe3_secexp$p_val))
cat(sprintf("  %-45s  coef=%+7.3f  cluster_p=%.3f  bootstrap_p=%.3f\n",
            trump,  coef(fit_FE3)[trump],
            coeftest(fit_FE3,cl_FE3)[trump,"Pr(>|t|)"],  b_fe3_trump$p_val))

# 9. MU1 Mundlak — within and between stress

stress_cm <- paste0(stress, "_cm")
trade_cm  <- paste0(trade, "_cm")

cat("\n--- MU1 Mundlak ---\n")
b_mu1_trade     <- do_boot(fit_MU1, trade)
b_mu1_stress    <- do_boot(fit_MU1, stress)
b_mu1_stress_cm <- do_boot(fit_MU1, stress_cm)
b_mu1_trade_cm  <- do_boot(fit_MU1, trade_cm)
b_mu1_trump     <- do_boot(fit_MU1, trump)

cat(sprintf("  %-55s  coef=%+7.3f  cluster_p=%.3f  bootstrap_p=%.3f\n",
            trade,      coef(fit_MU1)[trade],
            coeftest(fit_MU1,cl_MU1)[trade,"Pr(>|t|)"],      b_mu1_trade$p_val))
cat(sprintf("  %-55s  coef=%+7.3f  cluster_p=%.3f  bootstrap_p=%.3f\n",
            stress,     coef(fit_MU1)[stress],
            coeftest(fit_MU1,cl_MU1)[stress,"Pr(>|t|)"],     b_mu1_stress$p_val))
cat(sprintf("  %-55s  coef=%+7.3f  cluster_p=%.3f  bootstrap_p=%.3f\n",
            stress_cm,  coef(fit_MU1)[stress_cm],
            coeftest(fit_MU1,cl_MU1)[stress_cm,"Pr(>|t|)"],  b_mu1_stress_cm$p_val))
cat(sprintf("  %-55s  coef=%+7.3f  cluster_p=%.3f  bootstrap_p=%.3f\n",
            trade_cm,   coef(fit_MU1)[trade_cm],
            coeftest(fit_MU1,cl_MU1)[trade_cm,"Pr(>|t|)"],   b_mu1_trade_cm$p_val))
cat(sprintf("  %-55s  coef=%+7.3f  cluster_p=%.3f  bootstrap_p=%.3f\n",
            trump,      coef(fit_MU1)[trump],
            coeftest(fit_MU1,cl_MU1)[trump,"Pr(>|t|)"],      b_mu1_trump$p_val))

cat("\n=== INTERPRETATION ===\n")
cat("bootstrap_p > cluster_p  => few-cluster bias was inflating significance\n")
cat("bootstrap_p ≈ cluster_p  => standard clustered SEs were reliable\n")
cat("Key result: if trade bootstrap_p stays below 0.05 across M3 and MU1,\n")
cat("the headline Mercury finding is robust to the few-cluster correction.\n")