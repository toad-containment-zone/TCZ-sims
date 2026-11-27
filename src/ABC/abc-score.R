# abc-score.R
#
# Scores one Kimberley ABC run against the observed presence/absence. The run starts in 2009 and
# runs 14 generations to 2023, so a point present in 2023 has age > 0 and one present in 2022 has
# age - 1 > 0 (as in the original spread_Kimberley_ABC_concat.R). The score is the fraction of
# observed points (2022 and 2023 together) whose simulated state matches.

abc_score <- function(popmatrix) {
  obs.2023 <- popmatrix[, "X2023_pres"]
  obs.2022 <- popmatrix[, "X2022_pres"]
  sim.2023 <- as.numeric(popmatrix[, "age"] > 0)
  sim.2022 <- as.numeric((popmatrix[, "age"] - 1) > 0)
  m.2023 <- sum(sim.2023 == obs.2023, na.rm = TRUE)
  m.2022 <- sum(sim.2022 == obs.2022, na.rm = TRUE)
  n.2023 <- sum(is.finite(obs.2023))
  n.2022 <- sum(is.finite(obs.2022))
  c(score = (m.2023 + m.2022) / (n.2023 + n.2022),
    score.2022 = m.2022 / n.2022,
    score.2023 = m.2023 / n.2023,
    n.colonised = sum(popmatrix[, "Pres"] == 1))
}
