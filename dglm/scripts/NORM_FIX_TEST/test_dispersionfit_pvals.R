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

meta$age <- as.numeric(meta$age)
meta$sex <- factor(meta$sex)
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

  f <- try(
    suppressWarnings(
      dglm(
        e ~ age + sex + mean_n_umi + n_cells,
        ~ age,
        family = gaussian(),
        dlink = "log",
        data = d
      )
    ),
    silent = TRUE
  )

  if (inherits(f, "try-error")) {
    return(c(beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_))
  }

  tb <- try(coef(summary(f$dispersion.fit)), silent = TRUE)

  if (inherits(tb, "try-error") ||
      !("age" %in% rownames(tb)) ||
      !all(c("Estimate", "Std. Error", "Pr(>|t|)") %in% colnames(tb))) {
    return(c(beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_))
  }

  c(
    beta = unname(tb["age", "Estimate"]),
    se = unname(tb["age", "Std. Error"]),
    t = unname(tb["age", "t value"]),
    p = unname(tb["age", "Pr(>|t|)"])
  )
}

run_all <- function(age_vec, label) {
  cat("\n=== ", label, " ===\n", sep = "")

  ans <- t(mclapply(
    seq_len(nrow(logCPM)),
    fit_one_gene,
    age_vec = age_vec,
    mc.cores = 16
  ))

  ans <- as.data.frame(ans)
  ans$gene <- rownames(logCPM)

  valid <- ans[
    is.finite(ans$beta) &
    is.finite(ans$se) &
    ans$se > 0 &
    is.finite(ans$t) &
    is.finite(ans$p),
    ,
    drop = FALSE
  ]

  cat("Genes with valid dispersion-fit statistics: ", nrow(valid), "\n", sep = "")
  cat("Percent p < 0.05: ", round(100 * mean(valid$p < 0.05), 2), "%\n", sep = "")
  cat("Median p: ", round(median(valid$p), 4), "\n", sep = "")

  ans
}

real_results <- run_all(meta$age, "REAL AGE")

set.seed(1)
perm_results <- run_all(sample(meta$age), "PERMUTED AGE")

write.csv(
  real_results,
  "/scratch/easmit31/dispersion/dglm/scripts/NORM_FIX_TEST/dispersionfit_real.csv",
  row.names = FALSE
)

write.csv(
  perm_results,
  "/scratch/easmit31/dispersion/dglm/scripts/NORM_FIX_TEST/dispersionfit_permuted.csv",
  row.names = FALSE
)
