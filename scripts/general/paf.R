#!/usr/bin/env Rscript

# ==============================
# Attributable Fraction (GLM)
# ==============================

suppressPackageStartupMessages({
  library(dplyr)
})

# ------------------------------
# Parse command-line arguments
# ------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript attributable_fraction.R <metadata.tsv> <ordered_vars.txt>")
}

metadata_file <- args[1]
ordered_vars_file <- args[2]
cancer_type <- args[3]

# Read ordered variables (one per line)
#ordered_vars <- read_lines(ordered_vars_file) %>%
#  trimws() %>%
#  .[. != ""]
ordered_vars <- readLines(ordered_vars_file)
ordered_vars <- trimws(ordered_vars)
ordered_vars <- ordered_vars[ordered_vars != ""]

# Baseline covariates
covariates <- c("PC1", "PC2", "PC3", "PC4")

# Conditionally add passed_sv_quality_control
if (cancer_type != "cervix") {
  covariates <- c(covariates, "passed_sv_quality_control")
}

# ------------------------------
# Read data
# ------------------------------
#df <- read_tsv(metadata_file, show_col_types = FALSE) %>%
#  mutate(case = ifelse(tolower(original_dx) == "control", 0, 1))
df <- read.delim(metadata_file, sep='\t', stringsAsFactors = FALSE)
df$case <- ifelse(tolower(df$original_dx) == "control", 0, 1)

# Add sex_binary if applicable
if ("inferred_sex" %in% colnames(df) && n_distinct(df$inferred_sex) > 1) {
  df <- df %>%
    mutate(
      sex_binary = case_when(
        tolower(inferred_sex) == "male" ~ 1,
        tolower(inferred_sex) == "female" ~ 0,
        TRUE ~ NA_real_
      )
    )
  covariates <- c(covariates, "sex_binary")
}

# ------------------------------
# Transform R2 to liability scale
# ------------------------------
# Transform R2 to liability scale per Equation 23 in Lee et al., AJHG, 2011
lee.transform <- function(r2, K, P){
  r2 * ( (K * (1-K)) / (dnorm(norminv(K))^2) ) * ( (K * (1-K)) / (P * (1-P)) )
}

# ------------------------------
# Nagelkerke R² (manual)
# ------------------------------
nagelkerke_r2 <- function(model, null_model) {
  llf <- as.numeric(logLik(model))
  llnull <- as.numeric(logLik(null_model))
  n <- nobs(model)
  
  (1 - exp((2 / n) * (llnull - llf))) /
    (1 - exp((2 / n) * llnull))
}

# ------------------------------
# Attributable fraction function
# ------------------------------
compute_attributable_fraction <- function(df, outcome_col,
                                          predictor_cols, covariates) {
  
  full_formula <- as.formula(
    paste(outcome_col, "~", paste(c(predictor_cols, covariates), collapse = "+"))
  )
  
  reduced_formula <- if (length(covariates) > 0) {
    as.formula(paste(outcome_col, "~", paste(covariates, collapse = "+")))
  } else {
    as.formula(paste(outcome_col, "~ 1"))
  }
  
  null_formula <- as.formula(paste(outcome_col, "~ 1"))
  
  full_model <- glm(full_formula, data = df, family = binomial())
  reduced_model <- glm(reduced_formula, data = df, family = binomial())
  null_model <- glm(null_formula, data = df, family = binomial())
  
  r2_full <- nagelkerke_r2(full_model, null_model)
  r2_reduced <- nagelkerke_r2(reduced_model, null_model)
  
  data.frame(
    R2_full = r2_full,
    R2_reduced = r2_reduced,
    fraction_explained = r2_full - r2_reduced
  )
}

# ------------------------------
# Iterate over predictors
# ------------------------------
results <- list()

for (i in seq_along(ordered_vars)) {
  current_predictors <- ordered_vars[1:i]
  newest_predictor <- ordered_vars[i]
  
  cat("Processing:", newest_predictor, "\n")
  
  res <- compute_attributable_fraction(
    df,
    outcome_col = "case",
    predictor_cols = current_predictors,
    covariates = covariates
  )
  
  res$added_predictor <- newest_predictor
  results[[i]] <- res
}

# ------------------------------
# Save output
# ------------------------------
results_df <- bind_rows(results)

write.table(results_df,paste0(cancer_type, ".attributable_fraction_results.tsv"),sep = "\t",row.names = FALSE, quote = FALSE)

cat("Finished! Results written to attributable_fraction_results.tsv\n")
