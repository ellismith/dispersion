suppressPackageStartupMessages({
  library(dglm)
  library(edgeR)
  library(parallel)
})

counts_file <- "/scratch/easmit31/dispersion/dglm/disp_age__raw_counts/GABAergic_neurons/filtered/GABAergic_neurons_0_ACC_filtered_cutoff0.5.csv"
meta_file <- "/scratch/easmit31/dispersion/dglm/disp_age__raw_counts/GABAergic_neurons/GABAergic_neurons_0_metadata.csv"

counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

meta <- meta[meta$region == "ACC", ]
meta <- meta[match(colnames(counts), meta$animal_id), ]
rownames(meta) <- meta$animal_id

stopifnot(identical(colnames(counts), rownames(meta)))

meta$sex <- factor(meta$sex)
meta$age <- as.numeric(meta$age)
meta$mean_n_umi <- as.numeric(meta$mean_n_umi)
meta$n_cells <- as.numeric(meta$n_cells)

y <- DGEList(counts = as.matrix(counts))
y <- calcNormFactors(y, method = "TMM")
logCPM <- cpm(y, log = TRUE, prior.count = 0.5)

fit_one_gene <- function(i, age_vec) {
  d <- data.frame(
    age = age_vec,
    sex = meta$sex,
    mean_n_umi = meta$mean_n_umi,
    n_cells = meta$n_cells,
    e = as.numeric(logCPM[i, ])
  )

  warn <- character(0)

  f <- withCallingHandlers(
    try(
      dglm(
        e ~ age + sex + mean_n_umi + n_cells,
        ~ age,
        family = gaussian(),
        dlink = "log",
        data = d
      ),
      silent = TRUE
    ),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (inherits(f, "try-error")) {
    return(c(beta = NA, se = NA, z = NA, mean_converged = FALSE,
             disp_converged = FALSE, had_warning = TRUE))
  }

  tb <- try(coef(summary(f)$dispersion.summary), silent = TRUE)

  if (inherits(tb, "try-error") ||
      !("age" %in% rownames(tb)) ||
      !("Std. Error" %in% colnames(tb))) {
    return(c(beta = NA, se = NA, z = NA,
             mean_converged = isTRUE(f$fit$converged),
             disp_converged = isTRUE(f$dispersion.fit$converged),
             had_warning = length(warn) > 0))
  }

  beta <- unname(tb["age", "Estimate"])
  se <- unname(tb["age", "Std. Error"])

  c(
    beta = beta,
    se = se,
    z = beta / se,
    mean_converged = isTRUE(f$fit$converged),
    disp_converged = isTRUE(f$dispersion.fit$converged),
    had_warning = length(warn) > 0
  )
}

run_one <- function(age_vec, label) {
  cat("\n=== ", label, " ===\n", sep = "")

  ans <- t(vapply(
    seq_len(nrow(logCPM)),
    function(i) fit_one_gene(i, age_vec),
    numeric(6)
  ))

  ans <- as.data.frame(ans)
  ans$gene <- rownames(logCPM)

  ans$mean_converged <- as.logical(ans$mean_converged)
  ans$disp_converged <- as.logical(ans$disp_converged)
  ans$had_warning <- as.logical(ans$had_warning)

  valid <- ans[
    ans$mean_converged &
    ans$disp_converged &
    !ans$had_warning &
    is.finite(ans$beta) &
    is.finite(ans$se) &
    ans$se > 0 &
    is.finite(ans$z),
    ,
    drop = FALSE
  ]

  cat("Total genes: ", nrow(ans), "\n", sep = "")
  cat("Both fits converged: ",
      sum(ans$mean_converged & ans$disp_converged, na.rm = TRUE),
      "\n", sep = "")
  cat("No-warning, converged, valid genes: ", nrow(valid), "\n", sep = "")
  cat("Percent |z| > 1.96 in valid genes: ",
      round(100 * mean(abs(valid$z) > 1.96), 2),
      "%\n", sep = "")
  cat("Percent p < 0.05 in valid genes: ",
      round(100 * mean(2 * pnorm(-abs(valid$z)) < 0.05), 2),
      "%\n", sep = "")

  ans
}

real_results <- run_one(meta$age, "REAL AGE")

set.seed(1)
perm_results <- run_one(sample(meta$age), "PERMUTED AGE")

write.csv(
  real_results,
  "/scratch/easmit31/dispersion/dglm/scripts/NORM_FIX_TEST/converged_only_real.csv",
  row.names = FALSE
)

write.csv(
  perm_results,
  "/scratch/easmit31/dispersion/dglm/scripts/NORM_FIX_TEST/converged_only_permuted.csv",
  row.names = FALSE
)
