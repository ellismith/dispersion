suppressPackageStartupMessages({
  library(dglm)
  library(edgeR)
})

counts <- read.csv(
  "/scratch/easmit31/dispersion/dglm/disp_age__raw_counts/GABAergic_neurons/filtered/GABAergic_neurons_0_ACC_filtered_cutoff0.5.csv",
  row.names = 1,
  check.names = FALSE
)

meta <- read.csv(
  "/scratch/easmit31/dispersion/dglm/disp_age__raw_counts/GABAergic_neurons/GABAergic_neurons_0_metadata.csv",
  stringsAsFactors = FALSE
)

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

fit_one <- function(i, age_vec) {
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

  if (inherits(f, "try-error")) return(NA_real_)

  tb <- try(coef(summary(f$dispersion.fit)), silent = TRUE)

  if (inherits(tb, "try-error") ||
      !("age" %in% rownames(tb)) ||
      !("Pr(>|t|)" %in% colnames(tb))) {
    return(NA_real_)
  }

  unname(tb["age", "Pr(>|t|)"])
}

idx <- seq_len(min(1000, nrow(logCPM)))

real_p <- vapply(idx, fit_one, numeric(1), age_vec = meta$age)

set.seed(1)
perm_p <- vapply(idx, fit_one, numeric(1), age_vec = sample(meta$age))

cat("Genes tested: ", length(idx), "\n", sep = "")
cat(
  "Real age: valid = ", sum(is.finite(real_p)),
  "; p < 0.05 = ",
  round(100 * mean(real_p < 0.05, na.rm = TRUE), 2),
  "%\n",
  sep = ""
)
cat(
  "Fake age: valid = ", sum(is.finite(perm_p)),
  "; p < 0.05 = ",
  round(100 * mean(perm_p < 0.05, na.rm = TRUE), 2),
  "%\n",
  sep = ""
)
