# Economic and Security Determinants of US Dollar Reserve Holdings

MSc Economics thesis, Vrije Universiteit Amsterdam (2025–2026).

Panel dataset covering 22 reserve-holding economies from 2010 to 2025, combining central-bank, financial-market, and geopolitical data to study what drives a currency's international reserve status. Applies an Eichengreen "Mars vs Mercury" framework, testing whether security alignment (NATO membership, military presence, defence dependence) and financial stability/credibility explain USD reserve share variation.

Key finding: a country's trade exposure to the US was the most consistent predictor of its dollar reserve share across specifications.

## Scripts

| File | Purpose |
|---|---|
| `Fin-stress_construction.R` | Constructs the financial-stress index variable, merging panel and GPR (geopolitical risk) data with imputation flags |
| `OLS_pooled.R` | Main pooled OLS specification |
| `OLS_pooled_and_correlation.R` | Pooled OLS with correlation diagnostics |
| `OLS_pooled_no_fin_stress.R` | Robustness check excluding the financial-stress regressor |
| `FE_OLS.R` | Country fixed-effects model (main specification) |
| `FE_no_fin_stress.R` | Fixed-effects robustness check excluding financial stress |
| `Strict_Mundlak.R` | Mundlak correlated-random-effects (CRE) specification |
| `IV.R` | Instrumental-variables estimation |
| `placebo.R` | Placebo / falsification tests |
| `bootstrap.R`, `wild_bootstrap.R` | Wild cluster bootstrap inference, addressing the small-cluster (N=22) problem in standard clustered standard errors, following Cameron, Gelbach & Miller (2008) |

Data file referenced throughout: `MASTER_PANEL_DATA_GLOBAL_RISK_corrected_finstress2_flagged.xlsx` (not included — hand-assembled from central-bank reports, LSEG Refinitiv, and other sources described in the thesis).
