# Synthetic demo/regression script for the control-step (management-driven extinction) mechanics
# added to src/pprocess_functions.R. Uses made-up point data (no DATA_PATH dependency) so it can be
# run standalone to sanity-check the new setup()/spread.pilb() arguments before wiring real TCZ
# infrastructure-audit data through load-data.R.
#
# Three checks:
#   A. Regression: control.step=FALSE reproduces today's permanent-occupation behavior exactly,
#      even though the new columns are present in spread.table.
#   B. Attrition: with control.step=TRUE, tank/fence points go extinct at their expected annual
#      rates (1%/5%) and natural points (default surv.prob=1) never do.
#   D. Eradication responses: each occupied fence found at the annual audit is counted once.
#   E. Extreme-wet-year kernel: extreme.years all FALSE reproduces a standard run exactly (same seed);
#      all TRUE, with a wetter extreme kernel, colonises distant targets more often.
#   F. Natural drying: a drying matrix of 1s empties every pre-colonised natural point, 0s none;
#      draw_drying_matrix() gives p.avg at zero rainfall anomaly and 0 for "never" points.
#   C. Fence breach-thinning: fenced target points (r.eff == r, isolating p.breach's effect) get
#      colonised less often than otherwise-identical unfenced targets at the same distance, by
#      about the amount Binom(c, p.breach) thinning predicts. Uses a realistic, small p.breach
#      (~0.1) so that thinning in the wrong direction (by 1 - p.breach) fails the check.

source("src/pprocess_functions.R")
set.seed(20260730)

r <- R_DETECT # detection radius (m), the model's shared value

######## Dataset 1: attrition test (tank / fence / natural points, all pre-colonised) ########
n.per.type <- 800
make_grid <- function(n, spacing = 1000, offset = 0) {
  side <- ceiling(sqrt(n))
  g <- expand.grid(gx = seq_len(side), gy = seq_len(side))[seq_len(n), ]
  data.frame(X = g$gx * spacing + offset, Y = g$gy * spacing)
}

tank.pts    <- cbind(make_grid(n.per.type, offset = 0),       control_type = "tank",  art_nat = 1)
fence.pts   <- cbind(make_grid(n.per.type, offset = 1e6),     control_type = "fence", art_nat = 1)
natural.pts <- cbind(make_grid(n.per.type, offset = 2e6),     control_type = "none",  art_nat = 0)

attrition.data <- rbind(tank.pts, fence.pts, natural.pts)
attrition.data$present     <- 1 # all pre-colonised
attrition.data$fence_area  <- ifelse(attrition.data$control_type=="fence", pi*r^2, NA) # r.eff == r
attrition.data$fence_units <- ifelse(attrition.data$control_type=="fence", 54, 0)      # p.breach ~ 0.5 (unused here, no new colonisation)

setup(point.data = attrition.data,
      X.id = "X", Y.id = "Y", present.id = "present", artificial.natural.id = "art_nat",
      remove_duplicates = FALSE, constant.rain = 100, TCZ = FALSE,
      control.type.id = "control_type", fence.area.id = "fence_area", fence.units.id = "fence_units")
attrition.table <- spread.table
attrition.pairs <- pairs

cat("\n== A. Regression check (control.step=FALSE) ==\n")
out.A <- spread.pilb(pop = attrition.table, gens = 10, pairs = attrition.pairs,
                      delta = 0, r = r, rollup = FALSE, control.step = FALSE)
cat("All points still colonised:", all(out.A$popmatrix[, "Pres"]==1), "\n")
cat("occ.years == age for every point:", all(out.A$popmatrix[, "occ.years"]==out.A$popmatrix[, "age"]), "\n")
cat("age advanced by 10 (from its initial value of 1) for every point:", all(out.A$popmatrix[, "age"]==11), "\n")

cat("\n== B. Attrition rates (control.step=TRUE, delta=0, 1 generation) ==\n")
# control.fail.prob is P(control fails, i.e. toads persist), so extinction probability is
# 1-control.fail.prob: ~99% for tanks (1% failure rate per the qmd), ~95% for fences (5% eradication
# failure rate) -- both are deliberately HIGH per-generation extinction probabilities; it's the fence's
# breach-thinning at the colonisation step (tested separately below) that limits how often a fence
# point gets colonised in the first place.
out.B <- spread.pilb(pop = attrition.table, gens = 1, pairs = attrition.pairs,
                      delta = 0, r = r, rollup = FALSE, control.step = TRUE)
pm <- out.B$popmatrix
rate <- function(type) 1 - mean(pm[pm[,"control.type"]==type, "Pres"])
cat(sprintf("Tank extinction rate:    %.4f (expected ~0.99)\n", rate(CONTROL_TANK)))
cat(sprintf("Fence extinction rate:   %.4f (expected ~0.95)\n", rate(CONTROL_FENCE)))
cat(sprintf("Natural extinction rate: %.4f (expected 0, surv.prob defaults to 1)\n", rate(CONTROL_NONE)))

cat("\n== D. Eradication responses (same run as B) ==\n")
# every occupied fence is found at the audit, so each pre-colonised fence has exactly one response
# after one generation; tanks and natural points never trigger one
stopifnot(all(pm[pm[, "control.type"]==CONTROL_FENCE, "erad.events"] == 1),
          all(pm[pm[, "control.type"]!=CONTROL_FENCE, "erad.events"] == 0))
cat("One eradication response per occupied fence, none elsewhere: TRUE\n")

cat("\n== occ.years / first.col.gen sanity (control.step=TRUE, 10 generations) ==\n")
out.C <- spread.pilb(pop = attrition.table, gens = 10, pairs = attrition.pairs,
                      delta = 0, r = r, rollup = FALSE, control.step = TRUE)
pm <- out.C$popmatrix
cat("occ.years <= age for every point:", all(pm[, "occ.years"] <= pm[, "age"]), "\n")
cat("first.col.gen is 0 for all (all pre-colonised at t=0):", all(pm[, "first.col.gen"]==0), "\n")

cat("\n== rollup + control.step incompatibility is rejected ==\n")
rejected <- tryCatch({
  spread.pilb(pop = attrition.table, gens = 1, pairs = attrition.pairs, delta = 0, r = r,
              rollup = TRUE, control.step = TRUE)
  FALSE
}, error = function(e) TRUE)
cat("stop() raised when rollup=TRUE & control.step=TRUE:", rejected, "\n")

######## Dataset 2: fence breach-thinning test ########
n.targets <- 30
angles <- seq(0, 2*pi, length.out = n.targets*2 + 1)[seq_len(n.targets*2)]
dist.to.seed <- 1000

# fence_units = 8.2 (km + gates) gives p.breach = 1-(1-0.05)^(8.2/4) ~ 0.10. Fenced points get
# control_fail_prob = 1 (a breached fence is never eradicated), so the end-of-generation control
# step can't knock them back and the comparison measures breach-thinning alone.
fence.units.C <- 8.2
seed.pt      <- data.frame(X = 0, Y = 0, control_type = "none", art_nat = 1, present = 1,
                           fence_area = NA, fence_units = 0, control_fail_prob = NA)
unfenced.pts <- data.frame(X = dist.to.seed*cos(angles[1:n.targets]),
                            Y = dist.to.seed*sin(angles[1:n.targets]),
                            control_type = "none", art_nat = 1, present = 0,
                            fence_area = NA, fence_units = 0, control_fail_prob = NA)
fenced.pts   <- data.frame(X = dist.to.seed*cos(angles[(n.targets+1):(2*n.targets)]),
                            Y = dist.to.seed*sin(angles[(n.targets+1):(2*n.targets)]),
                            control_type = "fence", art_nat = 1, present = 0,
                            fence_area = pi*r^2, fence_units = fence.units.C, # r.eff == r
                            control_fail_prob = 1)

breach.data <- rbind(seed.pt, unfenced.pts, fenced.pts)

setup(point.data = breach.data,
      X.id = "X", Y.id = "Y", present.id = "present", artificial.natural.id = "art_nat",
      remove_duplicates = FALSE, constant.rain = 100, TCZ = FALSE,
      control.type.id = "control_type", fence.area.id = "fence_area", fence.units.id = "fence_units",
      control.fail.prob.id = "control_fail_prob")
breach.table <- spread.table
breach.pairs <- pairs
p.breach.val <- unique(breach.table[breach.table[,"control.type"]==CONTROL_FENCE, "p.breach"])
cat(sprintf("\n== C. Fence breach-thinning test (p.breach = %.3f, r.eff == r == %d) ==\n", p.breach.val, r))

# rows are in input order (remove_duplicates=FALSE): 1 = seed, 2:(n.targets+1) = unfenced targets,
# (n.targets+2):(2*n.targets+1) = fenced targets
unfenced.idx <- breach.table[, "ID"] %in% 2:(n.targets+1)
fenced.idx   <- breach.table[, "ID"] %in% (n.targets+2):(2*n.targets+1)
stopifnot(all(breach.table[unfenced.idx, "control.type"]==CONTROL_NONE),
          all(breach.table[fenced.idx, "control.type"]==CONTROL_FENCE))

n.reps <- 300
unfenced.col <- numeric(n.reps); fenced.col <- numeric(n.reps)
for (rr in seq_len(n.reps)) {
  out <- spread.pilb(pop = breach.table, gens = 1, pairs = breach.pairs,
                      delta = 200, r = r, rollup = FALSE, control.step = TRUE)
  pm <- out$popmatrix
  unfenced.col[rr] <- mean(pm[unfenced.idx, "Pres"])
  fenced.col[rr]   <- mean(pm[fenced.idx, "Pres"])
}
cat(sprintf("Unfenced target colonisation rate: %.3f\n", mean(unfenced.col)))
cat(sprintf("Fenced target colonisation rate:   %.3f\n", mean(fenced.col)))
# Expected rates, from the same multinomial-then-threshold logic: each target's arrivals are
# ~Poisson(m) for some m, unfenced colonised when >= COLONISE_MIN arrive, fenced when >= COLONISE_MIN
# of Binom(n, p.breach) get through, i.e. ~Poisson(m * p.breach). Back out m from the observed unfenced rate.
pois.col <- function(m) 1 - ppois(COLONISE_MIN - 1, m)
m.hat <- uniroot(function(m) pois.col(m) - mean(unfenced.col), c(1e-6, 1e4))$root
exp.right <- pois.col(m.hat * p.breach.val)       # thinning by p.breach (intended)
exp.wrong <- pois.col(m.hat * (1 - p.breach.val)) # thinning by 1 - p.breach (the old bug)
cat(sprintf("Expected fenced rate if thinned by p.breach:     %.3f\n", exp.right))
cat(sprintf("Expected fenced rate if thinned by 1 - p.breach: %.3f\n", exp.wrong))
stopifnot(abs(mean(fenced.col) - exp.right) < abs(mean(fenced.col) - exp.wrong),
          mean(fenced.col) < 0.5 * mean(unfenced.col))
cat("Breach-thinning direction check passed\n")

# age==0 iff first.col.gen is NA holds for any point never yet colonised (using the last rep's result,
# which always has some never-colonised target points since colonisation isn't certain)
cat("age==0 exactly where first.col.gen is NA (last rep):",
    all((pm[, "age"]==0) == is.na(pm[, "first.col.gen"])), "\n")


######## Dataset 3: extreme-wet-year kernel ########
cat("\n== E. Extreme-wet-year kernel ==\n")
n.far <- 40
far.angles <- seq(0, 2*pi, length.out = n.far + 1)[seq_len(n.far)]
far.dist <- 1500 # the extreme kernel is flatter: fewer arrivals within ~1 km, more beyond (expected ~0.12 vs ~0.18 here)
kernel.data <- rbind(data.frame(X = 0, Y = 0, present = 1),
                     data.frame(X = far.dist*cos(far.angles), Y = far.dist*sin(far.angles), present = 0))
kernel.data$art_nat <- 1
kernel.data$control_type <- "none"
kernel.data$rain_extreme <- 80 # rainy days in an extreme year -> rain_to_days(80) ~ 181 movement days
setup(point.data = kernel.data,
      X.id = "X", Y.id = "Y", present.id = "present", artificial.natural.id = "art_nat",
      remove_duplicates = FALSE, constant.rain = 100, TCZ = FALSE, # standard: 100 movement days
      control.type.id = "control_type", rain.extreme.id = "rain_extreme")
far.idx <- spread.table[, "ID"] > 1
run_E <- function(extreme, seed.E) {
  set.seed(seed.E)
  spread.pilb(pop = spread.table, gens = 1, pairs = pairs, delta = 200, r = r, rollup = FALSE,
              control.step = TRUE, pairs.extreme = if (is.null(extreme)) NULL else pairs.extreme,
              extreme.years = extreme)$popmatrix
}
same <- all(sapply(1:20, function(k) identical(run_E(NULL, k), run_E(FALSE, k))))
stopifnot(same)
cat("extreme.years = FALSE reproduces a standard run exactly:", same, "\n")
std.rate <- mean(sapply(1:200, function(k) mean(run_E(FALSE, k)[far.idx, "Pres"])))
ext.rate <- mean(sapply(1:200, function(k) mean(run_E(TRUE, k)[far.idx, "Pres"])))
cat(sprintf("Colonisation rate of targets %d m away: standard %.3f, extreme %.3f\n", far.dist, std.rate, ext.rate))
stopifnot(ext.rate > std.rate)
cat("Extreme kernel check passed\n")

######## F. Natural-point drying ########
cat("\n== F. Natural-point drying ==\n")
source("src/tcz-sims/drying-model.R")
n.att <- nrow(attrition.table)
nat.rows <- attrition.table[, "nats"] == 1
run_F <- function(p) spread.pilb(pop = attrition.table, gens = 1, pairs = attrition.pairs, delta = 0, r = r,
                                 rollup = FALSE, control.step = TRUE,
                                 drying.matrix = matrix(p, nrow = 1, ncol = n.att))$popmatrix
stopifnot(all(run_F(1)[nat.rows, "Pres"] == 0), all(run_F(0)[nat.rows, "Pres"] == 1))
cat("Drying matrix of 1s empties natural points, 0s leaves them: TRUE\n")
fake.spec <- list(class = c("never", "sometimes", "sometimes"), sometimes = c(FALSE, TRUE, TRUE),
                  x1 = matrix(c(0, 0, 0, 0, 1, -1), nrow = 2, byrow = TRUE),
                  x2 = matrix(0, nrow = 2, ncol = 3), p.avg = 0.25, beta = c(-1.3, -0.8))
dm <- draw_drying_matrix(fake.spec, idx = c(1, 2))
stopifnot(all(dm[, 1] == 0), isTRUE(all.equal(dm[1, 2:3], c(0.25, 0.25))),
          dm[2, 2] < 0.25, dm[2, 3] > 0.25)
cat("draw_drying_matrix(): p.avg at zero anomaly, wetter -> less drying, never -> 0: TRUE\n")
