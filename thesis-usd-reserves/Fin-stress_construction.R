
# add_finstress_flags.R

library(readxl)
library(dplyr)
library(writexl)

# 0. Load data (adjust path as needed)
df <- read_excel(
  "MASTER_PANEL_DATA_GLOBAL_RISK_corrected_finstress2.xlsx",
  sheet = "panel_with_gpr"
)

us_data <- read_excel(
  "MASTER_PANEL_DATA_GLOBAL_RISK_corrected_finstress2.xlsx",
  sheet = "US mercury"
) %>%
  rename(year = Year)

# 1. Per-component imputation flags
#    flag = 1 if the value in the model came from imputation
#    (i.e. raw is NA but _imp is not NA)

df <- df %>%
  mutate(
    # Bank Z-spread: no existing miss flag, derive from raw vs imp
    bank_zspread_imp_flag = as.integer(
      is.na(bank_zspread_avg_bps) & !is.na(bank_zspread_avg_bps_imp)
    ),
    
    # Sovereign CDS: use existing miss flag if present, else derive
    sov_cds_imp_flag = as.integer(
      is.na(sovereign_cds_5y_bps) & !is.na(sovereign_cds_5y_bps_imp)
    ),
    
    # Inflation gap
    inflation_gap_imp_flag = as.integer(
      is.na(inflation_gap_pp) & !is.na(inflation_gap_pp_imp)
    )
  )

cat("Imputation flag counts:\n")
cat("  bank_zspread_imp_flag == 1:", sum(df$bank_zspread_imp_flag, na.rm = TRUE), "\n")
cat("  sov_cds_imp_flag == 1:     ", sum(df$sov_cds_imp_flag,     na.rm = TRUE), "\n")
cat("  inflation_gap_imp_flag == 1:", sum(df$inflation_gap_imp_flag, na.rm = TRUE), "\n\n")

# 2. Re-normalise RAW components using global min/max of
#    RAW values only (so the scale is consistent with the
#    existing indices but uses no imputed observations)

minmax_raw <- function(x) {
  mn <- min(x, na.rm = TRUE)
  mx <- max(x, na.rm = TRUE)
  (x - mn) / (mx - mn)
}

df <- df %>%
  mutate(
    bankcred_norm_raw = minmax_raw(bank_zspread_avg_bps),
    sovcred_norm_raw  = minmax_raw(sovereign_cds_5y_bps),
    moncred_norm_raw  = minmax_raw(inflation_gap_pp)
  )

# 3. Absolute financial stress indices — RAW only

df <- df %>%
  mutate(
    

    # Strict: both raw components must be present
    fin_stress_2part_raw_strict = ifelse(
      !is.na(bankcred_norm_raw) & !is.na(sovcred_norm_raw),
      (bankcred_norm_raw + sovcred_norm_raw) / 2,
      NA_real_
    ),
    # Flag: 1 if either component was imputed in the ORIGINAL index
    fin_stress_2part_raw_strict_any_imp = as.integer(
      bank_zspread_imp_flag == 1 | sov_cds_imp_flag == 1
    ),
    
    # Flexible: average whichever raw components are available
    fin_stress_2part_raw_flexible = case_when(
      !is.na(bankcred_norm_raw) & !is.na(sovcred_norm_raw) ~
        (bankcred_norm_raw + sovcred_norm_raw) / 2,
      !is.na(bankcred_norm_raw) ~ bankcred_norm_raw,
      !is.na(sovcred_norm_raw)  ~ sovcred_norm_raw,
      TRUE ~ NA_real_
    ),
    # Flag: 1 if ANY imputed value was used in corresponding flexible index
    fin_stress_2part_raw_flexible_any_imp = as.integer(
      bank_zspread_imp_flag == 1 | sov_cds_imp_flag == 1
    ),
    

    # Strict: all three raw components present
    fin_stress_3part_raw_strict = ifelse(
      !is.na(bankcred_norm_raw) & !is.na(sovcred_norm_raw) & !is.na(moncred_norm_raw),
      (bankcred_norm_raw + sovcred_norm_raw + moncred_norm_raw) / 3,
      NA_real_
    ),
    fin_stress_3part_raw_strict_any_imp = as.integer(
      bank_zspread_imp_flag == 1 | sov_cds_imp_flag == 1 | inflation_gap_imp_flag == 1
    ),
    
    # Flexible: average available raw components
    fin_stress_3part_raw_flexible = {
      n_avail <- (!is.na(bankcred_norm_raw)) +
        (!is.na(sovcred_norm_raw))  +
        (!is.na(moncred_norm_raw))
      total <- ifelse(!is.na(bankcred_norm_raw), bankcred_norm_raw, 0) +
        ifelse(!is.na(sovcred_norm_raw),  sovcred_norm_raw,  0) +
        ifelse(!is.na(moncred_norm_raw),  moncred_norm_raw,  0)
      ifelse(n_avail > 0, total / n_avail, NA_real_)
    },
    fin_stress_3part_raw_flexible_any_imp = as.integer(
      bank_zspread_imp_flag == 1 | sov_cds_imp_flag == 1 | inflation_gap_imp_flag == 1
    )
  )

cat("Raw absolute index coverage (non-missing obs):\n")
cat("  fin_stress_2part_raw_strict:   ",
    sum(!is.na(df$fin_stress_2part_raw_strict)), "\n")
cat("  fin_stress_2part_raw_flexible: ",
    sum(!is.na(df$fin_stress_2part_raw_flexible)), "\n")
cat("  fin_stress_3part_raw_strict:   ",
    sum(!is.na(df$fin_stress_3part_raw_strict)), "\n")
cat("  fin_stress_3part_raw_flexible: ",
    sum(!is.na(df$fin_stress_3part_raw_flexible)), "\n\n")

cat("Any-imputation flags (rows where index contains imputed value):\n")
cat("  2part strict:   ", sum(df$fin_stress_2part_raw_strict_any_imp,   na.rm=TRUE), "\n")
cat("  2part flexible: ", sum(df$fin_stress_2part_raw_flexible_any_imp, na.rm=TRUE), "\n")
cat("  3part strict:   ", sum(df$fin_stress_3part_raw_strict_any_imp,   na.rm=TRUE), "\n")
cat("  3part flexible: ", sum(df$fin_stress_3part_raw_flexible_any_imp, na.rm=TRUE), "\n\n")

# 4. Relative financial stress indices — RAW only
#    Relative = holder raw component minus US component,
#    re-normalised. Mirrors the logic used for the existing
#    fin_stress_*_rel columns.

# Merge US raw components by year
df <- df %>%
  left_join(
    us_data %>% select(
      year,
      us_bank_zspread  = `US Bank Z-spread Avg (bps)`,
      us_sov_cds       = `US Sovereign CDS 5Y (bps)`,
      us_inflation_gap = `US Inflation Gap |BEI−π*| (pp)`
    ),
    by = "year"
  )

# Raw relative differences (holder raw minus US)
df <- df %>%
  mutate(
    bankcred_rel_raw  = bank_zspread_avg_bps - us_bank_zspread,
    sovcred_rel_raw   = sovereign_cds_5y_bps - us_sov_cds,
    moncred_rel_raw   = inflation_gap_pp     - us_inflation_gap
  )

# Normalise relative raw components
df <- df %>%
  mutate(
    bankcred_norm_rel_raw = minmax_raw(bankcred_rel_raw),
    sovcred_norm_rel_raw  = minmax_raw(sovcred_rel_raw),
    moncred_norm_rel_raw  = minmax_raw(moncred_rel_raw)
  )

# Build relative indices
df <- df %>%
  mutate(
    
    # TWO-PART RELATIVE
    fin_stress_2part_rel_raw_strict = ifelse(
      !is.na(bankcred_norm_rel_raw) & !is.na(sovcred_norm_rel_raw),
      (bankcred_norm_rel_raw + sovcred_norm_rel_raw) / 2,
      NA_real_
    ),
    fin_stress_2part_rel_raw_strict_any_imp = as.integer(
      bank_zspread_imp_flag == 1 | sov_cds_imp_flag == 1
    ),
    
    fin_stress_2part_rel_raw_flexible = case_when(
      !is.na(bankcred_norm_rel_raw) & !is.na(sovcred_norm_rel_raw) ~
        (bankcred_norm_rel_raw + sovcred_norm_rel_raw) / 2,
      !is.na(bankcred_norm_rel_raw) ~ bankcred_norm_rel_raw,
      !is.na(sovcred_norm_rel_raw)  ~ sovcred_norm_rel_raw,
      TRUE ~ NA_real_
    ),
    fin_stress_2part_rel_raw_flexible_any_imp = as.integer(
      bank_zspread_imp_flag == 1 | sov_cds_imp_flag == 1
    ),
    
    # THREE-PART RELATIVE
    fin_stress_3part_rel_raw_strict = ifelse(
      !is.na(bankcred_norm_rel_raw) & !is.na(sovcred_norm_rel_raw) & !is.na(moncred_norm_rel_raw),
      (bankcred_norm_rel_raw + sovcred_norm_rel_raw + moncred_norm_rel_raw) / 3,
      NA_real_
    ),
    fin_stress_3part_rel_raw_strict_any_imp = as.integer(
      bank_zspread_imp_flag == 1 | sov_cds_imp_flag == 1 | inflation_gap_imp_flag == 1
    ),
    
    fin_stress_3part_rel_raw_flexible = {
      n_avail <- (!is.na(bankcred_norm_rel_raw)) +
        (!is.na(sovcred_norm_rel_raw))  +
        (!is.na(moncred_norm_rel_raw))
      total <- ifelse(!is.na(bankcred_norm_rel_raw), bankcred_norm_rel_raw, 0) +
        ifelse(!is.na(sovcred_norm_rel_raw),  sovcred_norm_rel_raw,  0) +
        ifelse(!is.na(moncred_norm_rel_raw),  moncred_norm_rel_raw,  0)
      ifelse(n_avail > 0, total / n_avail, NA_real_)
    },
    fin_stress_3part_rel_raw_flexible_any_imp = as.integer(
      bank_zspread_imp_flag == 1 | sov_cds_imp_flag == 1 | inflation_gap_imp_flag == 1
    )
  )

cat("Raw relative index coverage (non-missing obs):\n")
cat("  fin_stress_2part_rel_raw_strict:   ",
    sum(!is.na(df$fin_stress_2part_rel_raw_strict)), "\n")
cat("  fin_stress_2part_rel_raw_flexible: ",
    sum(!is.na(df$fin_stress_2part_rel_raw_flexible)), "\n")
cat("  fin_stress_3part_rel_raw_strict:   ",
    sum(!is.na(df$fin_stress_3part_rel_raw_strict)), "\n")
cat("  fin_stress_3part_rel_raw_flexible: ",
    sum(!is.na(df$fin_stress_3part_rel_raw_flexible)), "\n\n")

# 5. 2025 audit: which rows in 2025 carry any imputed value
df_2025 <- df %>%
  filter(year == 2025) %>%
  select(country,
         bank_zspread_imp_flag,
         sov_cds_imp_flag,
         inflation_gap_imp_flag,
         fin_stress_2part_raw_flexible_any_imp,
         fin_stress_3part_raw_flexible_any_imp)

cat("2025 imputation audit:\n")
print(df_2025)
cat("\n")

# 6. New columns added — summary
new_cols <- c(
  # flags
  "bank_zspread_imp_flag", "sov_cds_imp_flag", "inflation_gap_imp_flag",
  # raw absolute
  "fin_stress_2part_raw_strict",          "fin_stress_2part_raw_strict_any_imp",
  "fin_stress_2part_raw_flexible",        "fin_stress_2part_raw_flexible_any_imp",
  "fin_stress_3part_raw_strict",          "fin_stress_3part_raw_strict_any_imp",
  "fin_stress_3part_raw_flexible",        "fin_stress_3part_raw_flexible_any_imp",
  # raw relative
  "fin_stress_2part_rel_raw_strict",      "fin_stress_2part_rel_raw_strict_any_imp",
  "fin_stress_2part_rel_raw_flexible",    "fin_stress_2part_rel_raw_flexible_any_imp",
  "fin_stress_3part_rel_raw_strict",      "fin_stress_3part_rel_raw_strict_any_imp",
  "fin_stress_3part_rel_raw_flexible",    "fin_stress_3part_rel_raw_flexible_any_imp"
)

cat("New columns added to df:\n")
cat(paste(" ", new_cols, collapse = "\n"), "\n\n")
cat("df now has", ncol(df), "columns and", nrow(df), "rows.\n")

# 7. Save
write_xlsx(df, "MASTER_PANEL_DATA_GLOBAL_RISK_corrected_finstress2_flagged.xlsx")