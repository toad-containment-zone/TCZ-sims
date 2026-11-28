# Null model for invasion speed: for each year with an empirical front-position
# estimate, start the point-process spread model from that front and run it
# forward one generation (~one wet season = one calendar year -- see rain_to_days()
# in pprocess_functions.R), many times, then compare the simulated front
# displacement against the empirical year-over-year displacement estimated in the
# invasion-front-monitoring repo (time-to-detection-multi-year_JAGS.R).
#
# Requires invasion-front-parameters-multi-year.Rdata (mod.multi, post.samples,
# bbox.fit, mean.coord, scale, years.all) to already exist in $DATA_PATH --
# produced by (re-)running that script in invasion-front-monitoring.
#
# STATUS (2026-11-28, Phase 2 of docs/roadmap.md): on the shared calibration and model structure.
# Rain comes from SILO via rain-years.R (Nov-May rainy days, season.days = 212), with two scenarios:
#   - static: each point's mean and extreme kernels over the 2002-2026 reference years, switching
#     to the extreme kernel in years that were extreme by the top-20% rule;
#   - actual rainfall: each season's own Nov-May rainy days per cell, one kernel.
# Results are an out-of-sample check of the Kimberley calibration against observed front movement.
#
# Notes:
# - The starting front for year k is fixed at its posterior mean (a, b_k) -- not
#   resampled per rep -- so the resulting spread of b_sim across reps reflects only
#   the stochastic spread process (propagule draws, colonisation draws, lambda
#   resampling in run_sims()), not front-location uncertainty. That's compared
#   against the empirical delta[k] posterior, which instead reflects detection-data
#   uncertainty -- these are two different kinds of uncertainty.
# - b_sim is recovered per rep via a deterministic misclassification-minimising
#   threshold on the simulated Pres pattern (fit_b_threshold()), not a smooth
#   logistic fit -- this matches the deterministic step-function front definition
#   used by score_colonised()/score_colonised_multi_year(), and avoids
#   (quasi-)separation problems a logistic fit would hit when only a narrow ring of
#   points near the front changes status in a single generation.
# - The threshold fit is restricted to points inside bbox.fit (the Albers-metres
#   bounding box of the data the empirical front line was actually fit to), even
#   though the simulation itself runs on the full point network -- the front is
#   only approximately straight over that shorter span.

######## load functions, libraries, and point data ########
source("src/tcz-sims/load-data.R") # produces all.points (colonisation_year, Albers X/Y already attached)

######## load empirical front parameters (produced by invasion-front-monitoring) ########
# Moved up from its original position (after all.points/survey.point.bbox construction) --
# survey.point.bbox below needs a/mean.coord_m/b_metres() to build a front-relative cutoff.
load(file.path(Sys.getenv("DATA_PATH"), "invasion-front-parameters-multi-year.Rdata"))
# objects brought into scope: mod.multi, post.samples, bbox.fit, mean.coord, scale, years.all, run.info

a <- mod.multi[[1]]["a", 1]        # posterior-mean slope (shared across years, fixed across reps)
mean.coord_m <- mean.coord * scale # mean.coord (km, model space) -> Albers metres

# posterior-mean b (Albers metres) for years.all[tt] -- same transform score_colonised_multi_year() uses
b_metres <- function(tt) {
  b_raw <- mod.multi[[1]][paste0("b[", tt, "]"), 1]
  (b_raw + mean.coord["Y"]) * scale
}
b1 <- b_metres(1) # earliest modelled year's intercept, used to build survey.point.bbox below

# replace old.lagrange.points with points from a different geographic area, to cover the area
# covered by toad monitoring surveys plus an 80km buffer -- except on the already-colonised
# (east) side, where the boundary follows a line parallel to the invasion front's axis, 80km
# (perpendicular) behind the earliest modelled front (years.all[1]). A plain axis-aligned bbox
# buffer can't guarantee uniform perpendicular clearance from a diagonal front (slope `a`); this
# does. West/north/south sides aren't reported as problematic and keep the original 80km buffer.
# Note, invasion front analysis needs to be run first, and its repo needs to be sitting in same
# parent directory as this repo.
load(file = "../invasion-front-monitoring/out/merged-visual-surveys.RData")
survey.bb <- df |>
  st_transform(crs = 3577) |>
  st_bbox()

# west/north/south: original 80km axis-aligned buffer. East side left generously unconstrained
# (1000km) -- irrelevant once intersected with front.cutoff.poly below (that's always the
# binding edge in practice), so its exact value doesn't need tuning.
survey.box <- st_bbox(c(xmin = survey.bb[["xmin"]] - 80000,
                         ymin = survey.bb[["ymin"]] - 80000,
                         xmax = survey.bb[["xmax"]] + 1e6,
                         ymax = survey.bb[["ymax"]] + 80000),
                       crs = 3577) |>
  st_as_sfc()

# front-relative cutoff: half-plane D <= b1.offset, D = Y - a*(X - meanX_m) (score_colonised()'s
# convention: colonised iff D > b). b1.offset shifts b1 80km into colonised (increasing-D)
# territory: a perpendicular offset `dist` in D-space is `dist * sqrt(1+a^2)` (same conversion
# used later in this script for disp.sim.km).
b1.offset <- b1 + 80000 * sqrt(1 + a^2)

# Offset line D = b1.offset, built generously far along the front's axis (u = X - meanX_m).
# singleSide buffers give the LEFT side of the line's direction of travel for positive dist,
# RIGHT for negative (?st_buffer). For a line traversed in the direction of increasing u (tangent
# (1,a)), the left-hand normal is (-a,1) == grad(D) (D's increasing direction) -- so the *right*
# side (negative dist) is the decreasing-D/uncolonised-ish side we want, for any sign of `a`.
u.vals <- c(-2e6, 2e6) # +/-2000km along the front's axis -- generously covers WA
front.offset.line <- st_linestring(cbind(mean.coord_m[["X"]] + u.vals,
                                          a * u.vals + b1.offset))
front.cutoff.poly <- st_sfc(front.offset.line, crs = 3577) |>
  st_buffer(dist = -2e6, singleSide = TRUE, endCapStyle = "FLAT")

survey.point.bbox <- st_intersection(survey.box, front.cutoff.poly)
stopifnot(length(survey.point.bbox) > 0, !st_is_empty(survey.point.bbox))
# essentially all df survey points should fall inside the final box -- if not, something is
# badly wrong with the cutoff direction/parameters, not just a few points on a tight edge
n.survey.outside <- sum(!st_within(st_transform(df, crs = 3577), survey.point.bbox, sparse = FALSE))
if (n.survey.outside > 0) {
  warning(n.survey.outside, " of ", nrow(df), " survey points fall outside survey.point.bbox")
}

# replace all.points with the relevant set of waterpoints from combined waterpoint data
all.points <- st_read(file.path(spatial.dir, "Edited-layers/merged-points.shp")) |> 
  select(fcsubtype_, full_name, perennia_1, origin_des, watercou_1, area_m, st_perimet) |> 
  st_transform(crs = 4326) |>
  bind_rows(d_extra) |> # add additional points from aerial imagery
  bind_rows(df |> select(full_name = display_label, geometry)) |>  # add additional points from on-ground surveys
  st_transform(crs = 3577) |>
  st_filter(survey.point.bbox, .predicate = st_within)
  

####### Bring in estimated invasion front(s) and scored colonised points #######
colnsd <- score_colonised_multi_year(all.points)
all.points <- colnsd$scored.points # includes colonised (0/1) and colonisation_year
inv.front <- colnsd$fronts # one front line per modelled year
fp <- colnsd$fp; rm(colnsd)

# get coordinates in albers
albers <- all.points |> 
  st_transform(crs = 3577) |> 
  st_coordinates()
all.points <- bind_cols(all.points, albers)



######## SILO rain: the shared two-kernel inputs and per-year rainy days (rain-years.R) ########
# The front region is on the Kimberley side, so its cells share the Kimberley ABC's SILO cache
# (keyed per cell; any new cells are fetched into it). This replaces the script's own SILO fetch
# and dat/silo-actual-rain-days-cache.RData (Phase 2, 2026-11-28).
source("src/tcz-sims/rain-years.R")
all.points$art_nat <- as.integer(!coalesce(all.points$origin_des == "Natural", FALSE)) # 1 = artificial, 0 = natural (no control step here)
silo <- attach_silo_kernels(all.points, cache = "dat/silo-kimberley-rain-cache.RData")
all.points <- silo$points
rain.years <- silo$rain.years
rm(silo)

######## load simulation posteriors ########
load_posteriors() # the shared calibration, POSTERIORS_FILE in pprocess_functions.R

## Make a static figure here to show selected points, and the invasion front lines on a map.
dir.create("out/null-model", showWarnings = FALSE)
bbox <- st_bbox(all.points)
selected.points.fig <- ggplot() +
  geom_sf(data = wa.coast, fill = NA, color = "grey30") +
  geom_sf(data = tcz.boundary, fill = NA, color = "grey60", linetype = "dashed") +
  geom_sf(data = inv.front, aes(group = factor(year)), color = "red", lty = 2) +
  geom_sf(data = all.points, aes(color = factor(colonisation_year)), size = 0.7) +
  geom_sf(data = fp, color = "red") +
  labs(color = "Colonisation year") +
  # all.points/inv.front/fp are Albers (3577) here, unlike wa.coast/tcz.boundary (4326) -- pin the
  # plot CRS explicitly so the Albers-metre bbox isn't misread as lon/lat by coord_sf's default (first-layer) CRS
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]),
           crs = st_crs(all.points))
ggsave(filename = "out/null-model/selected-points-and-fronts.pdf", plot = selected.points.fig,
       width = 180, height = 160, units = "mm")


######## helper: recover b from a colonisation pattern ########
# Given perpendicular-axis position D = Y - a*(X - meanX_m) and colonisation status Pres
# (0/1) for a set of points, find the b that best separates colonised from uncolonised
# points under the rule "colonised if D > b" (score_colonised()'s own rule) -- i.e. the
# threshold minimising misclassifications against the observed Pres pattern. Ties broken
# by taking the middle of the tied set. Returns NA if every point falls on one side (no
# separating threshold exists within the data).
fit_b_threshold <- function(D, Pres) {
  ord <- order(D)
  D <- D[ord]; Pres <- Pres[ord]
  n <- length(D)
  if (n < 2) return(NA_real_)
  total1 <- sum(Pres)
  cum1 <- c(0, cumsum(Pres))          # cum1[i+1] = colonised count among first i points, i = 0..n
  i.seq <- 0:n
  errors <- 2 * cum1 - i.seq - total1 + n  # errors[i+1] = misclassifications if 1..i called uncolonised
  best <- which(errors == min(errors)) - 1L # back to i-indexing (0..n)
  i.star <- best[ceiling(length(best) / 2)]
  if (i.star <= 0L || i.star >= n) return(NA_real_) # degenerate: whole window on one side
  mean(D[i.star:(i.star + 1L)])
}

######## quick self-test of fit_b_threshold() on synthetic data ########
# (independent of real data/geometry -- just checks the fitting logic itself)
{
  set.seed(1)
  D.test <- sort(runif(500, 0, 1000))
  b.true <- 550
  Pres.test <- as.integer(D.test > b.true)
  b.rec <- fit_b_threshold(D.test, Pres.test)
  stopifnot(abs(b.rec - b.true) < diff(range(D.test)) / length(D.test) * 2)
}

######## years to simulate forward from (every year except the last has a "next year") ########
# A front fitted to year yr's surveys moves in the following wet season: season year yr + 1
# (Jul yr to Jun yr + 1, rain-years.R), so that season's rain and extreme flag drive the step.
sim.years <- years.all[-length(years.all)]
season.of <- sim.years + 1
if (!all(season.of %in% rain.years$ref$season.year))
  stop("season years ", paste(setdiff(season.of, rain.years$ref$season.year), collapse = ", "),
       " are outside the SILO reference years (rain-years.R, rain.ref.years)")
extreme.of <- rain.years$ref$extreme[match(season.of, rain.years$ref$season.year)]
cat("Season years", paste(season.of, collapse = ", "), "| extreme:",
    paste(season.of[extreme.of], collapse = ", "), "\n")
n.sims <- 100 # sim reps per year -- raise for a smoother null distribution (runtime scales linearly)

######## one year: start from that year's front, run one generation, recover the new front ########
# Uses the global spread.table/pairs (and pairs.extreme) from the last setup() call. The
# colonisation-dependent columns are overwritten in place: exactly what setup() would compute for
# them given col.year.id = NULL (see pprocess_functions.R), without rebuilding the dispersal matrix.
simulate_year <- function(tt, extreme = NULL, pairs.extreme = NULL) {
  yr <- sim.years[tt]
  kept.ids <- spread.table[, "ID"] # all.points row for each thinned spread.table row
  present_k_full <- ifelse(!is.na(all.points$colonisation_year) &
                             all.points$colonisation_year <= yr, 1L, 0L)
  present_k <- present_k_full[kept.ids]
  spread.table[, "Pres"] <<- present_k
  spread.table[, "age"] <<- present_k
  spread.table[, "occ.years"] <<- present_k
  spread.table[, "first.col.gen"] <<- ifelse(present_k == 1, 0, NA_real_)

  ######## sanity check: threshold-fit should recover ~b_k from the (deterministic) starting state ########
  b_k <- b_metres(tt)
  in.bbox0 <- spread.table[, "X"] >= bbox.fit["xmin"] & spread.table[, "X"] <= bbox.fit["xmax"] &
    spread.table[, "Y"] >= bbox.fit["ymin"] & spread.table[, "Y"] <= bbox.fit["ymax"]
  D0 <- spread.table[, "Y"] - a * (spread.table[, "X"] - mean.coord_m["X"])
  b0.check <- fit_b_threshold(D0[in.bbox0], spread.table[in.bbox0, "Pres"])
  cat("  Sanity check: b recovered from starting state =", round(b0.check),
      "m vs fitted b_k =", round(b_k), "m (diff =", round(b0.check - b_k), "m)\n")

  ######## run one generation (= one wet season) forward, many reps ########
  sim_out <- run_sims(n.sims = n.sims, gens = 1, plot = FALSE, rollup = FALSE,
                      stop.on.target = FALSE, pairs.extreme = pairs.extreme, extreme.years = extreme)

  ######## fit b_sim per rep, restricted to the empirical fitting bbox ########
  disp.sim.km <- vapply(sim_out, function(rep) {
    pm <- rep$popmatrix
    in.bbox <- pm[, "X"] >= bbox.fit["xmin"] & pm[, "X"] <= bbox.fit["xmax"] &
      pm[, "Y"] >= bbox.fit["ymin"] & pm[, "Y"] <= bbox.fit["ymax"]
    pm <- pm[in.bbox, , drop = FALSE]
    D <- pm[, "Y"] - a * (pm[, "X"] - mean.coord_m["X"])
    b_sim <- fit_b_threshold(D, pm[, "Pres"])
    (b_sim - b_k) / sqrt(1 + a^2) / 1000 # metres -> km; signed, same convention as empirical delta
  }, numeric(1))

  n.degenerate <- sum(is.na(disp.sim.km))
  if (n.degenerate > 0) {
    cat("  Warning:", n.degenerate, "of", n.sims, "reps had no separating threshold in bbox.fit (dropped)\n")
  }

  ######## empirical comparison: delta[tt] = displacement from years.all[tt] to years.all[tt+1] ########
  list(year = yr,
       season.year = season.of[tt],
       sim = disp.sim.km,
       empirical = post.samples[, paste0("delta[", tt, "]")],
       n.degenerate = n.degenerate,
       n.sims = n.sims)
}

######## STATIC SCENARIO: SILO mean and extreme kernels, real extreme years ########
# Built once: every point's standard kernel from its mean Nov-May rainy days over the reference
# years and its extreme kernel from the maximum (attach_silo_kernels(), as in the TCZ runs and the
# forecast); each year uses the extreme kernel if its season year was extreme by the top-20% rule.
# present_k is a placeholder (overwritten per year); no targets and no control step.
all.points$present_k <- 0L
setup(point.data = all.points,
      X.id = "X",
      Y.id = "Y",
      present.id = "present_k",
      artificial.natural.id = "art_nat",
      rain.id = "rain.mean",
      rain.extreme.id = "rain.extreme",
      constant.rain = NULL,
      trunc.dist = TRUE,
      TCZ = FALSE)

null.results <- vector("list", length(sim.years))
for (tt in seq_along(sim.years)) {
  cat("\n==== [static] Year", sim.years[tt], "( season", season.of[tt], if (extreme.of[tt]) "extreme" else "", ") ====\n")
  null.results[[tt]] <- simulate_year(tt, extreme = extreme.of[tt], pairs.extreme = pairs.extreme)
}
names(null.results) <- sim.years

######## PER-YEAR SCENARIO: each season's actual Nov-May rainy days ########
# How much of the empirical/null discrepancy is real inter-annual rainfall? The kernel now varies
# by year, so setup() rebuilds spread.table/pairs each year from that season's rainy days per cell
# (rain.years$cell.year), with one kernel (the season's own rain, not an extreme switch).
null.results.actual.rain <- vector("list", length(sim.years))
for (tt in seq_along(sim.years)) {
  cat("\n==== [actual rainfall] Year", sim.years[tt], "( season", season.of[tt], ") ====\n")
  cy <- rain.years$cell.year |> filter(season.year == season.of[tt])
  all.points$rain.season <- cy$rain.days[match(rain.years$point.cell, cy$silo_cell)]
  stopifnot(!anyNA(all.points$rain.season))
  cat("  Mean Nov-May rainy days:", round(mean(all.points$rain.season), 1),
      "(reference mean", round(mean(all.points$rain.mean), 1), ")\n")
  setup(point.data = all.points,
        X.id = "X",
        Y.id = "Y",
        present.id = "present_k",
        artificial.natural.id = "art_nat",
        rain.id = "rain.season",
        constant.rain = NULL,
        trunc.dist = TRUE,
        TCZ = FALSE)
  null.results.actual.rain[[tt]] <- simulate_year(tt)
  null.results.actual.rain[[tt]]$mean.rain.days <- mean(all.points$rain.season)
}
names(null.results.actual.rain) <- sim.years

######## combine both scenarios with empirical into one three-way comparison ########
comparison.df.all <- do.call(rbind, c(
  lapply(null.results, function(x) rbind(
    data.frame(year = x$year, source = "empirical", distance_km = abs(x$empirical)),
    data.frame(year = x$year, source = "simulated (null)", distance_km = abs(x$sim[!is.na(x$sim)])))),
  lapply(null.results.actual.rain, function(x)
    data.frame(year = x$year, source = "simulated (actual rainfall)", distance_km = abs(x$sim[!is.na(x$sim)])))))
comparison.df.all$source <- factor(comparison.df.all$source,
                                   levels = c("empirical", "simulated (null)", "simulated (actual rainfall)"))

# per year: medians, and where the empirical median falls in each simulated distribution
summary.df <- do.call(rbind, lapply(seq_along(sim.years), function(tt) {
  emp <- median(abs(null.results[[tt]]$empirical))
  st <- abs(na.omit(null.results[[tt]]$sim)); ac <- abs(na.omit(null.results.actual.rain[[tt]]$sim))
  data.frame(year = sim.years[tt], season.year = season.of[tt], extreme = extreme.of[tt],
             empirical.median.km = emp,
             static.median.km = median(st), static.emp.pctile = mean(st <= emp),
             actual.median.km = median(ac), actual.emp.pctile = mean(ac <= emp),
             actual.mean.rain.days = null.results.actual.rain[[tt]]$mean.rain.days)
}))
print(summary.df, digits = 3)

dir.create("out/null-model", showWarnings = FALSE)
write.csv(summary.df, "out/null-model/front-displacement-summary.csv", row.names = FALSE)
save(null.results, null.results.actual.rain, comparison.df.all, summary.df,
     file = "out/null-model/front-displacement-null-vs-empirical-vs-actual-rainfall.RData")

print(
  ggplot(comparison.df.all, aes(x = factor(year), y = distance_km, fill = source)) +
    geom_violin(alpha = 0.6, colour = NA, position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    labs(x = "Starting year", y = "Perpendicular front movement (km)", fill = NULL) +
    theme_bw()
)
ggsave("out/null-model/front-displacement-three-way.pdf", width = 180, height = 120, units = "mm")

######## degenerate-rep summary across all years ########
degen <- rbind(static = vapply(null.results, `[[`, numeric(1), "n.degenerate"),
               actual = vapply(null.results.actual.rain, `[[`, numeric(1), "n.degenerate"))
if (any(degen > 0)) {
  cat("\nDegenerate (whole-bbox-one-sided) reps by year:\n")
  print(degen)
}
