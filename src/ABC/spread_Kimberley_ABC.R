# spread_Kimberley_ABC.R
#
# Rejection-ABC runs for lambda on the Kimberley invasion 2009-2023, with the same model structure as
# the TCZ runs (docs/roadmap.md, Phase 1): SILO standard and extreme kernels, the real 2010-2023
# extreme-year sequence, c >= 2, no drying (so rollup is on). Needs out/setup_kimberley_silo.RData
# from src/ABC/kimberley-setup.R.
#
# Runs are done in batches, in parallel (forked, so the large pairs matrices are shared). Each batch
# has its own seed and is saved to <out.dir>/ABC_<batch>.RData as a table of lambda and scores
# (src/ABC/abc-score.R), not the full population matrices (~0.7 MB a run). Batches already on disk
# are skipped, so a stopped run resumes.
#
# Usage (from the repo root):
#   Rscript src/ABC/spread_Kimberley_ABC.R 'prior.mean=5.5' 'prior.sd=0.5' 'n.batches=20' \
#     'runs.per.batch=10' 'out.dir="out/ABC-silo/pilot"' 'cores=4'
# The prior is log10(lambda) ~ N(prior.mean, prior.sd), drawn by run_sims().
#
# Large point set: set R_MAX_VSIZE high (usethis::edit_r_environ()).

prior.mean <- 5.5
prior.sd <- 0.5
n.batches <- 20
runs.per.batch <- 10
out.dir <- "out/ABC-silo/pilot"
cores <- 4
seed.base <- 20261127
for (a in commandArgs(TRUE)) eval(parse(text = a))

source("src/pprocess_functions.R")
source("src/ABC/abc-score.R")
library(parallel)

load("out/setup_kimberley_silo.RData") # spread.table, pairs, pairs.extreme, extreme.years, kim.start.year, kim.gens
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

sample.lambda <- prior.mean # run_sims() draws lambda from these globals
sample.lambda.sd <- prior.sd

run_batch <- function(b) {
  f <- file.path(out.dir, sprintf("ABC_%03d.RData", b))
  if (file.exists(f)) return(f)
  set.seed(seed.base + b)
  t0 <- Sys.time()
  sims <- run_sims(n.sims = runs.per.batch, gens = kim.gens, plot = FALSE, rollup = TRUE,
                   pairs.extreme = pairs.extreme, extreme.years = extreme.years)
  abc <- do.call(rbind, lapply(sims, function(s) c(lambda = unname(s$pars[, "lambda"]), abc_score(s$popmatrix))))
  abc <- data.frame(abc, batch = b, prior.mean = prior.mean, prior.sd = prior.sd)
  secs.per.run <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / runs.per.batch
  save(abc, secs.per.run, file = f)
  f
}

done <- mclapply(seq_len(n.batches), run_batch, mc.cores = cores, mc.preschedule = FALSE)
failed <- vapply(done, function(x) inherits(x, "try-error") || is.null(x), logical(1))
if (any(failed)) stop(sum(failed), " batch(es) failed: ", paste(which(failed), collapse = ", "))
cat("Done:", n.batches, "batches in", out.dir, "\n")
