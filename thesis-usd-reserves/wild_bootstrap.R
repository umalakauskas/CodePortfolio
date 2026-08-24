

packages <- c("readxl", "dplyr", "tibble", "sandwich")
to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(readxl); library(dplyr); library(tibble); library(sandwich)

set.seed(12345)
B           <- 9999   # switch to 999 for test
cluster_var <- "country"

# Load data

df_all <- read_excel(
  "MASTER_PANEL_DATA_GLOBAL_RISK_corrected_finstress2_flagged.xlsx",
  sheet = "Sheet1"
) %>%
  rename_with(tolower) %>%
  mutate(country = as.character(country), year = as.integer(year)) %>%
  filter(!is.na(usd_share_pct))

dv         <- "usd_share_pct"
trade      <- "log_us_trade_gdp"
size       <- "size_ratio_us_to_holder"
equity     <- "log_equity_turnover_gdp"
stress     <- "fin_stress_3part_raw_flexible"
stress_imp <- "fin_stress_3part_raw_flexible_any_imp"
equity_imp <- "equity_turnover_to_gdp_pct_miss"
nato       <- "nato_member_dummy"
secexp     <- "security_exposure_no_nato_flexible"
trump      <- "trump_2025"

df_strict <- df_all %>%
  filter(.data[[stress_imp]] == 0, .data[[equity_imp]] == 0) %>%
  filter(if_all(all_of(c(dv,"country","year",trade,size,equity,
                         stress,nato,secexp,trump)), ~ !is.na(.x))) %>%
  mutate(country = as.character(country))

df_nfs <- df_all %>%
  filter(.data[[equity_imp]] == 0) %>%
  filter(if_all(all_of(c(dv,"country","year",trade,size,equity,
                         nato,secexp,trump)), ~ !is.na(.x))) %>%
  mutate(country = as.character(country))

cat("Strict N:", nrow(df_strict),
    "| clusters:", n_distinct(df_strict$country), "\n")
cat("No-stress N:", nrow(df_nfs),
    "| clusters:", n_distinct(df_nfs$country), "\n\n")

# Mundlak prep

make_mundlak_data <- function(data, vars) {
  data %>%
    group_by(country) %>%
    mutate(across(all_of(vars),
                  list(cm = ~ mean(.x, na.rm=TRUE), w = ~ .x - mean(.x, na.rm=TRUE)),
                  .names = "{.col}_{.fn}")) %>%
    ungroup()
}

ols_f <- function(rhs) as.formula(paste(dv,"~",paste(rhs,collapse="+")))
fe_f  <- function(rhs) as.formula(paste(dv,"~",paste(rhs,collapse="+"),
                                        "+ factor(country) - 1"))

# Bootstrap function 

wild_p_fast <- function(model, data, term, B = 9999, cluster = "country",
                        verbose = FALSE) {
  
  X <- model.matrix(model)
  y <- model.response(model.frame(model))
  n <- nrow(X); k <- ncol(X)
  j <- which(colnames(X) == term)
  
  if (length(j) == 0) stop(paste("Term not in model matrix:", term))
  
  # (X'X)^{-1} — use pseudoinverse if singular
  XtX <- crossprod(X)
  XtX_inv <- tryCatch(
    solve(XtX),
    error = function(e) {
      if (verbose) cat("  [INFO] X'X singular, using MASS::ginv\n")
      MASS::ginv(XtX)
    }
  )
  
  # Restricted residuals: remove column j, refit with lm.fit (handles rank def)
  X_r    <- X[, -j, drop = FALSE]
  fit_r  <- lm.fit(X_r, y)          # lm.fit not .lm.fit — safer
  e_r    <- fit_r$residuals
  
  if (verbose) {
    cat("  Restricted residuals: mean =", round(mean(e_r), 6),
        "| any NA?", any(is.na(e_r)), "\n")
  }
  
  if (any(is.na(e_r))) {
    stop("Restricted model produced NA residuals — model matrix may be rank deficient")
  }
  
  # Cluster info
  clusters    <- data[[cluster]]
  cluster_ids <- unique(clusters)
  G           <- length(cluster_ids)
  hc1         <- (G / (G - 1)) * (n / (n - k))
  
  # Observed t-statistic
  beta_hat <- coef(model)[term]
  cl_vcov  <- sandwich::vcovCL(model, cluster = clusters, type = "HC1")
  se_hat   <- sqrt(diag(cl_vcov)[term])
  t_obs    <- beta_hat / se_hat
  
  if (verbose) cat("  t_obs =", round(t_obs, 4), "\n")
  
  # Cluster-level score contributions for beta_j
  xj_row   <- XtX_inv[j, , drop = FALSE]
  scores_i <- as.numeric(X %*% t(xj_row)) * e_r
  S        <- tapply(scores_i, clusters, sum)
  
  denom_const <- sqrt(hc1 * sum(S^2))
  
  if (verbose) {
    cat("  sum(S^2) =", round(sum(S^2), 8),
        "| denom_const =", round(denom_const, 8), "\n")
  }
  
  if (denom_const == 0 || is.na(denom_const)) {
    stop(paste("Zero/NA denominator — score vector S is all zeros for term:", term))
  }
  
  # Bootstrap loop
  t_boot <- numeric(B)
  for (b in seq_len(B)) {
    w         <- sample(c(-1, 1), G, replace = TRUE)
    t_boot[b] <- sum(w * S) / denom_const
  }
  
  mean(abs(t_boot) >= abs(t_obs))
}

# DIAGNOSTIC: test on one coefficient first 

cat("=== DIAGNOSTIC: testing trade coefficient in simple pooled OLS ===\n")

test_model <- lm(as.formula(paste(dv, "~", trade, "+", size, "+", equity,
                                  "+", stress, "+", nato, "+", secexp, "+", trump)),
                 data = df_strict)

cat("Model fitted. Coefficients found:\n")
print(names(coef(test_model)))

diag_result <- tryCatch({
  p <- wild_p_fast(test_model, df_strict, trade, B = 99, verbose = TRUE)
  cat("\nDiagnostic PASSED. Wild p-value (B=99):", round(p, 4), "\n\n")
  TRUE
}, error = function(e) {
  cat("\nDiagnostic FAILED. Error message:\n")
  cat(conditionMessage(e), "\n\n")
  FALSE
})

if (!diag_result) {
  cat("STOPPING: fix the error above before running the full bootstrap.\n")
  cat("Most likely causes:\n")
  cat("1. MASS package not installed — run: install.packages('MASS')\n")
  cat("2. Rank-deficient model matrix — check for collinear variables\n")
  cat("3. NA values in model data — check df_strict for missing values\n")
  stop("Bootstrap diagnostic failed — see messages above.")
}

cat("Diagnostic passed. Proceeding with full run (B =", B, ")...\n\n")

# Runner 

run_boot <- function(model, data, term, family, name, label) {
  cat(sprintf("%-30s | %-22s | %s\n", family, name, term))
  
  cn <- names(coef(model))
  if (!(term %in% cn)) {
    cat("  [SKIP] term not found in model\n")
    return(tibble(model_family=family, model_name=name, coefficient=label,
                  term=term, estimate=NA_real_, cluster_se=NA_real_,
                  cluster_t=NA_real_, cluster_p=NA_real_, wild_p=NA_real_,
                  B=B, status="term not found"))
  }
  
  cl_vcov <- sandwich::vcovCL(model, cluster=data[[cluster_var]], type="HC1")
  est     <- coef(model)[term]
  se      <- sqrt(diag(cl_vcov)[term])
  tval    <- est / se
  G       <- n_distinct(data[[cluster_var]])
  pval    <- 2 * pt(abs(tval), df=G-1, lower.tail=FALSE)
  
  wild_p <- tryCatch(
    wild_p_fast(model, data, term, B=B, cluster=cluster_var, verbose=FALSE),
    error = function(e) {
      # Print the actual error — this is what was missing before
      cat("  [ERROR]", conditionMessage(e), "\n")
      NA_real_
    }
  )
  
  tibble(model_family=family, model_name=name, coefficient=label, term=term,
         estimate=round(est,4), cluster_se=round(se,4),
         cluster_t=round(tval,4), cluster_p=round(pval,4),
         wild_p=round(wild_p,4), B=B,
         status=if(!is.na(wild_p)) "ok" else "bootstrap failed")
}

# Fit all models 

cat("Fitting all models...\n")

mu_vars     <- c(trade,size,equity,stress,secexp,nato)
mu_vars_nfs <- c(trade,size,equity,secexp,nato)
df_mu       <- make_mundlak_data(df_strict, mu_vars)
df_mu_nfs   <- make_mundlak_data(df_nfs, mu_vars_nfs)

mu_base <- c(paste0(trade,"_w"),paste0(trade,"_cm"),
             paste0(size,"_w"),paste0(size,"_cm"),
             paste0(equity,"_w"),paste0(equity,"_cm"),
             paste0(stress,"_w"),paste0(stress,"_cm"),
             paste0(secexp,"_w"),paste0(secexp,"_cm"),
             paste0(nato,"_w"),paste0(nato,"_cm"), trump)

mu_base_nfs <- c(paste0(trade,"_w"),paste0(trade,"_cm"),
                 paste0(size,"_w"),paste0(size,"_cm"),
                 paste0(equity,"_w"),paste0(equity,"_cm"), trump)

mu_mars_nfs <- c(paste0(trade,"_w"),paste0(trade,"_cm"),
                 paste0(size,"_w"),paste0(size,"_cm"),
                 paste0(equity,"_w"),paste0(equity,"_cm"),
                 paste0(secexp,"_w"),paste0(secexp,"_cm"),
                 paste0(nato,"_w"),paste0(nato,"_cm"), trump)

P_M3 <- lm(ols_f(c(trade,size,equity,stress,nato,secexp,trump)), df_strict)
P_M4 <- lm(ols_f(c(trade,size,equity,stress,nato,secexp,trump,paste0(nato,":",trump))), df_strict)
P_M5 <- lm(ols_f(c(trade,size,equity,stress,nato,secexp,trump,paste0(secexp,":",trump))), df_strict)
P_M6 <- lm(ols_f(c(trade,size,equity,stress,nato,secexp,trump,paste0(stress,":",trump))), df_strict)

FE3 <- lm(fe_f(c(trade,size,equity,stress,secexp,trump)), df_strict)
FE4 <- lm(fe_f(c(trade,size,equity,stress,secexp,trump,paste0(secexp,":",trump))), df_strict)
FE5 <- lm(fe_f(c(trade,size,equity,stress,secexp,trump,paste0(stress,":",trump))), df_strict)
FE6 <- lm(fe_f(c(trade,size,equity,stress,secexp,trump,paste0(nato,":",trump))), df_strict)

MU1 <- lm(ols_f(mu_base), df_mu)
MU2 <- lm(ols_f(c(mu_base,paste0(nato,":",trump))), df_mu)
MU3 <- lm(ols_f(c(mu_base,paste0(secexp,":",trump))), df_mu)
MU4 <- lm(ols_f(c(mu_base,paste0(stress,":",trump))), df_mu)

P_merc_nfs  <- lm(ols_f(c(trade,size,equity,trump)), df_nfs)
P_mars_nfs  <- lm(ols_f(c(trade,size,equity,nato,secexp,trump)), df_nfs)
FE_merc_nfs <- lm(fe_f(c(trade,size,equity,trump)), df_nfs)
FE_mars_nfs <- lm(fe_f(c(trade,size,equity,secexp,trump)), df_nfs)
MU_merc_nfs <- lm(ols_f(mu_base_nfs), df_mu_nfs)
MU_mars_nfs <- lm(ols_f(mu_mars_nfs), df_mu_nfs)

cat("Models fitted.\n\n")

get_int <- function(m, a, b) {
  cn <- names(coef(m))
  c1 <- paste0(a,":",b); c2 <- paste0(b,":",a)
  if(c1 %in% cn) return(c1); if(c2 %in% cn) return(c2); NA_character_
}

int_P_nato    <- get_int(P_M4,nato,trump)
int_P_secexp  <- get_int(P_M5,secexp,trump)
int_P_stress  <- get_int(P_M6,stress,trump)
int_FE_secexp <- get_int(FE4,secexp,trump)
int_FE_stress <- get_int(FE5,stress,trump)
int_FE_nato   <- get_int(FE6,nato,trump)
int_MU_nato   <- get_int(MU2,nato,trump)
int_MU_secexp <- get_int(MU3,secexp,trump)
int_MU_stress <- get_int(MU4,stress,trump)

# Run all bootstrap tests 

cat(sprintf("=== Wild cluster bootstrap B=%d ===\n\n", B))

results <- tibble()
add <- function(m,d,term,fam,nm,lab)
  results <<- bind_rows(results, run_boot(m,d,term,fam,nm,lab))

add(P_M3,df_strict,trade,  "Pooled OLS strict","M3 Full","Log US trade/GDP")
add(P_M3,df_strict,size,   "Pooled OLS strict","M3 Full","US-to-holder GDP size ratio")
add(P_M3,df_strict,equity, "Pooled OLS strict","M3 Full","Log equity turnover/GDP")
add(P_M3,df_strict,stress, "Pooled OLS strict","M3 Full","Financial stress")
add(P_M3,df_strict,nato,   "Pooled OLS strict","M3 Full","NATO member")
add(P_M3,df_strict,secexp, "Pooled OLS strict","M3 Full","Security exposure ex-NATO")
add(P_M3,df_strict,trump,  "Pooled OLS strict","M3 Full","Trump 2025")

if(!is.na(int_P_nato))   add(P_M4,df_strict,int_P_nato,  "Pooled OLS strict","M4 Trump x NATO",    "Trump 2025 x NATO")
if(!is.na(int_P_secexp)) add(P_M5,df_strict,int_P_secexp,"Pooled OLS strict","M5 Trump x Security","Trump 2025 x Security exposure")
if(!is.na(int_P_stress)) add(P_M6,df_strict,int_P_stress,"Pooled OLS strict","M6 Trump x Stress",  "Trump 2025 x Financial stress")

add(FE3,df_strict,trade,  "Country FE strict","FE3 Full","Log US trade/GDP")
add(FE3,df_strict,size,   "Country FE strict","FE3 Full","US-to-holder GDP size ratio")
add(FE3,df_strict,equity, "Country FE strict","FE3 Full","Log equity turnover/GDP")
add(FE3,df_strict,stress, "Country FE strict","FE3 Full","Financial stress")
add(FE3,df_strict,secexp, "Country FE strict","FE3 Full","Security exposure ex-NATO")
add(FE3,df_strict,trump,  "Country FE strict","FE3 Full","Trump 2025")

if(!is.na(int_FE_secexp)) add(FE4,df_strict,int_FE_secexp,"Country FE strict","FE4 Trump x Security","Trump 2025 x Security exposure")
if(!is.na(int_FE_stress)) add(FE5,df_strict,int_FE_stress,"Country FE strict","FE5 Trump x Stress",  "Trump 2025 x Financial stress")
if(!is.na(int_FE_nato))   add(FE6,df_strict,int_FE_nato,  "Country FE strict","FE6 Trump x NATO",    "Trump 2025 x NATO")

add(MU1,df_mu,paste0(trade,"_w"),  "Mundlak strict","MU1 Core","Log US trade/GDP within")
add(MU1,df_mu,paste0(trade,"_cm"), "Mundlak strict","MU1 Core","Log US trade/GDP country mean")
add(MU1,df_mu,paste0(size,"_w"),   "Mundlak strict","MU1 Core","US size ratio within")
add(MU1,df_mu,paste0(size,"_cm"),  "Mundlak strict","MU1 Core","US size ratio country mean")
add(MU1,df_mu,paste0(equity,"_w"), "Mundlak strict","MU1 Core","Log equity turnover/GDP within")
add(MU1,df_mu,paste0(equity,"_cm"),"Mundlak strict","MU1 Core","Log equity turnover/GDP country mean")
add(MU1,df_mu,paste0(stress,"_w"), "Mundlak strict","MU1 Core","Financial stress within")
add(MU1,df_mu,paste0(stress,"_cm"),"Mundlak strict","MU1 Core","Financial stress country mean")
add(MU1,df_mu,paste0(secexp,"_w"), "Mundlak strict","MU1 Core","Security exposure within")
add(MU1,df_mu,paste0(secexp,"_cm"),"Mundlak strict","MU1 Core","Security exposure country mean")
add(MU1,df_mu,paste0(nato,"_w"),   "Mundlak strict","MU1 Core","NATO within")
add(MU1,df_mu,paste0(nato,"_cm"),  "Mundlak strict","MU1 Core","NATO country mean")
add(MU1,df_mu,trump,               "Mundlak strict","MU1 Core","Trump 2025")

if(!is.na(int_MU_nato))   add(MU2,df_mu,int_MU_nato,  "Mundlak strict","MU2 Trump x NATO",    "Trump 2025 x NATO")
if(!is.na(int_MU_secexp)) add(MU3,df_mu,int_MU_secexp,"Mundlak strict","MU3 Trump x Security","Trump 2025 x Security exposure")
if(!is.na(int_MU_stress)) add(MU4,df_mu,int_MU_stress,"Mundlak strict","MU4 Trump x Stress",  "Trump 2025 x Financial stress")

add(P_merc_nfs,df_nfs,trade, "Pooled no financial stress","Mercury",   "Log US trade/GDP")
add(P_merc_nfs,df_nfs,equity,"Pooled no financial stress","Mercury",   "Log equity turnover/GDP")
add(P_merc_nfs,df_nfs,trump, "Pooled no financial stress","Mercury",   "Trump 2025")
add(P_mars_nfs,df_nfs,trade, "Pooled no financial stress","Mars added","Log US trade/GDP")
add(P_mars_nfs,df_nfs,equity,"Pooled no financial stress","Mars added","Log equity turnover/GDP")
add(P_mars_nfs,df_nfs,secexp,"Pooled no financial stress","Mars added","Security exposure ex-NATO")
add(P_mars_nfs,df_nfs,trump, "Pooled no financial stress","Mars added","Trump 2025")

add(FE_merc_nfs,df_nfs,trade, "FE no financial stress","Mercury",   "Log US trade/GDP")
add(FE_merc_nfs,df_nfs,size,  "FE no financial stress","Mercury",   "US-to-holder GDP size ratio")
add(FE_merc_nfs,df_nfs,equity,"FE no financial stress","Mercury",   "Log equity turnover/GDP")
add(FE_merc_nfs,df_nfs,trump, "FE no financial stress","Mercury",   "Trump 2025")
add(FE_mars_nfs,df_nfs,trade, "FE no financial stress","Mars added","Log US trade/GDP")
add(FE_mars_nfs,df_nfs,size,  "FE no financial stress","Mars added","US-to-holder GDP size ratio")
add(FE_mars_nfs,df_nfs,equity,"FE no financial stress","Mars added","Log equity turnover/GDP")
add(FE_mars_nfs,df_nfs,secexp,"FE no financial stress","Mars added","Security exposure ex-NATO")
add(FE_mars_nfs,df_nfs,trump, "FE no financial stress","Mars added","Trump 2025")

add(MU_merc_nfs,df_mu_nfs,paste0(trade,"_w"),  "Mundlak no financial stress","Mercury",   "Log US trade/GDP within")
add(MU_merc_nfs,df_mu_nfs,paste0(trade,"_cm"), "Mundlak no financial stress","Mercury",   "Log US trade/GDP country mean")
add(MU_merc_nfs,df_mu_nfs,paste0(equity,"_w"), "Mundlak no financial stress","Mercury",   "Log equity turnover/GDP within")
add(MU_merc_nfs,df_mu_nfs,paste0(equity,"_cm"),"Mundlak no financial stress","Mercury",   "Log equity turnover/GDP country mean")
add(MU_merc_nfs,df_mu_nfs,trump,               "Mundlak no financial stress","Mercury",   "Trump 2025")
add(MU_mars_nfs,df_mu_nfs,paste0(trade,"_w"),  "Mundlak no financial stress","Mars added","Log US trade/GDP within")
add(MU_mars_nfs,df_mu_nfs,paste0(trade,"_cm"), "Mundlak no financial stress","Mars added","Log US trade/GDP country mean")
add(MU_mars_nfs,df_mu_nfs,paste0(equity,"_w"), "Mundlak no financial stress","Mars added","Log equity turnover/GDP within")
add(MU_mars_nfs,df_mu_nfs,paste0(equity,"_cm"),"Mundlak no financial stress","Mars added","Log equity turnover/GDP country mean")
add(MU_mars_nfs,df_mu_nfs,paste0(secexp,"_w"), "Mundlak no financial stress","Mars added","Security exposure within")
add(MU_mars_nfs,df_mu_nfs,paste0(secexp,"_cm"),"Mundlak no financial stress","Mars added","Security exposure country mean")
add(MU_mars_nfs,df_mu_nfs,paste0(nato,"_w"),   "Mundlak no financial stress","Mars added","NATO within")
add(MU_mars_nfs,df_mu_nfs,paste0(nato,"_cm"),  "Mundlak no financial stress","Mars added","NATO country mean")
add(MU_mars_nfs,df_mu_nfs,trump,               "Mundlak no financial stress","Mars added","Trump 2025")

# Save and summarise 

results_clean <- results %>%
  mutate(
    sig_cl_10    = ifelse(!is.na(cluster_p) & cluster_p < 0.10, "Yes", "No"),
    sig_wl_10    = ifelse(!is.na(wild_p)    & wild_p    < 0.10, "Yes", "No"),
    sig_cl_05    = ifelse(!is.na(cluster_p) & cluster_p < 0.05, "Yes", "No"),
    sig_wl_05    = ifelse(!is.na(wild_p)    & wild_p    < 0.05, "Yes", "No"),
    cl_not_wl_05 = sig_cl_05 == "Yes" & sig_wl_05 == "No"
  )

write.csv(results_clean, "WildClusterBootstrap_results_DIAG.csv",  row.names = FALSE)
write.csv(filter(results_clean, status != "ok"),
          "WildClusterBootstrap_failed_DIAG.csv", row.names = FALSE)

# HTML table
he <- function(x) { x <- as.character(x)
gsub("<","&lt;", gsub("&","&amp;",x,fixed=TRUE), fixed=TRUE) }
df_d <- results_clean %>% mutate(across(
  c(estimate,cluster_se,cluster_t,cluster_p,wild_p),
  ~ ifelse(is.na(.x),"",sprintf("%.4f",.x))))
hdr  <- names(df_d)
rows <- apply(df_d,1,function(r) paste0("<tr>",paste0("<td>",he(r),"</td>",collapse=""),"</tr>"))
writeLines(paste0(
  "<!DOCTYPE html><html><head><meta charset='UTF-8'>",
  "<style>body{font-family:Arial;margin:24px}",
  "table{border-collapse:collapse;width:100%;font-size:13px}",
  "th,td{border:1px solid #ddd;padding:5px}th{background:#f2f2f2}",
  "tr:nth-child(even){background:#fafafa}</style></head><body>",
  "<h2>Wild cluster bootstrap</h2><table><thead><tr>",
  paste0("<th>",he(hdr),"</th>",collapse=""),
  "</tr></thead><tbody>",paste0(rows,collapse="\n"),
  "</tbody></table><p style='font-size:12px;color:#555'>",
  "CGM (2008), Rademacher, B=",B,
  ". Score-based. HC1. FE via country dummies.</p>",
  "</body></html>"),
  "WildClusterBootstrap_results_DIAG.html")

cat("\n=== STATUS ===\n"); print(table(results_clean$status))
cat("\n=== Cluster sig but wild not (p<0.05) ===\n")
print(filter(results_clean,cl_not_wl_05) %>%
        select(model_family,model_name,coefficient,estimate,cluster_p,wild_p), n=Inf)
cat("\n=== Pooled M3 ===\n")
print(filter(results_clean,model_family=="Pooled OLS strict",model_name=="M3 Full") %>%
        select(coefficient,estimate,cluster_p,wild_p,status), n=Inf)
cat("\n=== FE3 ===\n")
print(filter(results_clean,model_family=="Country FE strict",model_name=="FE3 Full") %>%
        select(coefficient,estimate,cluster_p,wild_p,status), n=Inf)
cat("\n=== MU1 ===\n")
print(filter(results_clean,model_family=="Mundlak strict",model_name=="MU1 Core") %>%
        select(coefficient,estimate,cluster_p,wild_p,status), n=Inf)

saveRDS(list(results=results_clean,B=B),
        "WildClusterBootstrap_results_DIAG.rds")

cat("\nDONE.\n")
cat("Outputs: WildClusterBootstrap_results_DIAG.csv/html/rds\n")