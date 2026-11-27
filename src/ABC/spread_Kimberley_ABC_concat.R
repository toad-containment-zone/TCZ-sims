# spread_Kimberley_ABC_concat.R
#
# Collects the Kimberley ABC batches from spread_Kimberley_ABC.R, keeps the top 5% of runs by score
# (agreement with the 2022 and 2023 presence data, src/ABC/abc-score.R) and saves the log10(lambda)
# posterior as sample.lambda / sample.lambda.sd, the objects run_sims() reads.
#
# Usage: Rscript src/ABC/spread_Kimberley_ABC_concat.R 'in.dir="out/ABC-silo/main"' \
#          'post.file="dat/Posteriors_2026-12.RData"'
# post.file = NULL (the default) only reports and plots, without saving a posterior.

in.dir <- "out/ABC-silo/main"
post.file <- NULL
accept.q <- 0.95 # keep the top 5%
for (a in commandArgs(TRUE)) eval(parse(text = a))

flist <- list.files(in.dir, pattern = "^ABC_.*\\.RData$", full.names = TRUE)
stopifnot(length(flist) > 0)
abc <- do.call(rbind, lapply(flist, function(f) { load(f); abc }))
abc$log.lambda <- log10(abc$lambda)
cat(nrow(abc), "runs from", length(flist), "batches\n")

cut.q <- quantile(abc$score, probs = accept.q)
acc <- abc[abc$score >= cut.q, ]
sample.lambda <- mean(acc$log.lambda)
sample.lambda.sd <- sd(acc$log.lambda)
cat(sprintf("Accepted %d runs (score >= %.4f) | log10 lambda mean %.3f, sd %.3f\n",
            nrow(acc), cut.q, sample.lambda, sample.lambda.sd))
cat(sprintf("Prior log10 lambda N(%.3f, %.3f); accepted range %.3f - %.3f; all runs %.3f - %.3f\n",
            abc$prior.mean[1], abc$prior.sd[1], min(acc$log.lambda), max(acc$log.lambda),
            min(abc$log.lambda), max(abc$log.lambda)))

# the score against lambda: the accepted runs should sit well inside the prior, and the peak should
# not be flat-topped (a flat top means the data can't pin lambda down at this resolution)
png(file.path(in.dir, "score-vs-lambda.png"), width = 18, height = 12, units = "cm", res = 200)
plot(score ~ log.lambda, data = abc, pch = 16, cex = 0.5, col = "grey50",
     xlab = expression(log[10](lambda)), ylab = "Agreement with 2022-23 presence")
points(score ~ log.lambda, data = acc, pch = 16, cex = 0.6, col = "firebrick")
abline(h = cut.q, lty = 2)
lines(lowess(abc$log.lambda, abc$score, f = 0.2), lwd = 2)
dev.off()

if (!is.null(post.file)) {
  abc.runs <- nrow(abc)
  save(sample.lambda, sample.lambda.sd, abc.runs, cut.q, file = post.file)
  cat("Saved", post.file, "\n")
}
