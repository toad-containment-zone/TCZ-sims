# drying-model.R
#
# Natural-waterpoint drying for the TCZ simulations (ms/tcz-optimisation-methods.qmd, "Auditing
# natural waterpoints"). Supersedes dea-drying-probs.R and rainfall-drying-model.R.
#
# Traditional Owner knowledge decides WHETHER a point dries; DEA decides only HOW drying responds
# to rainfall (project lead's decision, 2026-11-26). In DEA Waterbodies, the Living Waters points
# TOs class as "Always Wet" reach 0% open water in most dry seasons -- DEA sees open surface water,
# not the persistent water toads use -- so DEA can't calibrate whether a point dries.
#
# Drying classes (per simulated point):
#   never     -- Living Waters "Always Wet"; natural points with no TO assessment (LaGrange and
#                aerial-imagery points, all outside the TCZ); artificial points (not used)
#   sometimes -- Living Waters "Sometimes Dry", "Don't Know" or blank: in a season of average
#                rainfall at that point the drying probability is p.avg (0.25); it moves with the
#                wet-season rainfall anomaly and the previous season's, using the DEA slopes:
#                logit p = logit(p.avg) + b1 * x1 + b2 * x2
#
# Season year Y: wet season Dec (Y-1) to May (Y), dry season Jun to Dec (Y) (rain-years.R).
# DEA queries and time series are cached in dat/dea-drying-cache.RData (gitignored).

library(sf)
library(dplyr)
library(httr2)
library(lubridate)
library(lme4)

DEA_MAX_DIST_M     <- 100  # Living Waters point must be within 100 m of the DEA polygon
DEA_DRY_MONTHS     <- 6:12 # June-December
DEA_MIN_OBS        <- 3    # valid dry-season observations needed to score a year
dea.cache          <- "dat/dea-drying-cache.RData"

# Match Living Waters points to DEA waterbodies and score each matched site's dry seasons
# (dry = any valid observation with 0% wet). Returns one row per matched point and year.
dea_drying <- function(lw) {
  if (file.exists(dea.cache)) {
    load(dea.cache) # dea.wb, dea.ts
  } else {
    bb <- st_bbox(st_transform(lw, 4326)); buf <- 0.05
    url <- paste0("https://geoserver.dea.ga.gov.au/geoserver/wfs?service=WFS&version=1.1.0",
                  "&request=GetFeature&typeName=dea:DigitalEarthAustraliaWaterbodies_v3&bbox=",
                  paste(bb["xmin"] - buf, bb["ymin"] - buf, bb["xmax"] + buf, bb["ymax"] + buf,
                        "EPSG:4326", sep = ","),
                  "&outputFormat=application/json&maxFeatures=5000")
    dea.wb <- request(url) |> req_perform() |> resp_body_string() |> st_read(quiet = TRUE)
    if (nrow(dea.wb) >= 5000) stop("DEA WFS hit maxFeatures; raise it")
    dea.ts <- list()
  }
  lwa <- st_transform(lw, 3577)
  near <- st_nearest_feature(lwa, st_transform(dea.wb, 3577))
  dist <- as.numeric(st_distance(lwa, st_transform(dea.wb, 3577)[near, ], by_element = TRUE))
  m <- which(dist <= DEA_MAX_DIST_M)
  uids <- unique(dea.wb$uid[near[m]])
  for (u in setdiff(uids, names(dea.ts))) {
    ts <- read.csv(url(dea.wb$timeseries[dea.wb$uid == u][1]))
    ts$date <- as.Date(substr(ts$date, 1, 10))
    dea.ts[[u]] <- ts
  }
  save(dea.wb, dea.ts, file = dea.cache)

  annual <- lapply(m, function(k) {
    ts <- dea.ts[[dea.wb$uid[near[k]]]]
    ts <- ts[!is.na(ts$pc_wet) & month(ts$date) %in% DEA_DRY_MONTHS, ]
    ts |>
      group_by(season.year = year(date)) |>
      summarise(n.obs = n(), dry = any(pc_wet == 0), .groups = "drop") |>
      filter(n.obs >= DEA_MIN_OBS) |>
      mutate(site = dea.wb$uid[near[k]], lw_record_id = lw$lw_record_id[k],
             group_assessment = lw$group_assessment[k], dist.m = dist[k])
  })
  out <- bind_rows(annual)
  attr(out, "points") <- lw[m, ]
  out
}

# glmer(dry ~ wet_sc + lag_sc + (1 | site)) on the DEA-matched sites, with each site's wet-season
# rainfall from SILO (rain.years). Only the slopes b1, b2 (and their covariance) are used.
fit_drying_response <- function(annual, rain.years.lw) {
  cells <- data.frame(lw_record_id = attr(annual, "points")$lw_record_id,
                      silo_cell = rain.years.lw$point.cell)
  wet <- rain.years.lw$cell.year
  # previous season's wet mm, by explicit year match (robust to gaps)
  wet$wet.lag <- wet$wet.mm[match(paste(wet$silo_cell, wet$season.year - 1),
                                  paste(wet$silo_cell, wet$season.year))]
  d <- annual |>
    left_join(cells, by = "lw_record_id") |>
    inner_join(wet |> select(silo_cell, season.year, wet.mm, wet.lag), by = c("silo_cell", "season.year")) |>
    filter(!is.na(wet.lag))
  sc <- list(wet = c(center = mean(d$wet.mm), scale = sd(d$wet.mm)),
             lag = c(center = mean(d$wet.lag), scale = sd(d$wet.lag)))
  d$wet_sc <- (d$wet.mm - sc$wet["center"]) / sc$wet["scale"]
  d$lag_sc <- (d$wet.lag - sc$lag["center"]) / sc$lag["scale"]
  fit <- glmer(dry ~ wet_sc + lag_sc + (1 | site), data = d, family = binomial)
  beta <- fixef(fit)[c("wet_sc", "lag_sc")]
  list(fit = fit, beta = beta, vcov = as.matrix(vcov(fit))[c("wet_sc", "lag_sc"), c("wet_sc", "lag_sc")],
       scale = sc, n.sites = n_distinct(d$site), n.site.years = nrow(d),
       years = range(d$season.year), singular = isSingular(fit), data = d)
}

# Per-point drying class and rainfall anomalies over the reference years (rows = reference years).
# `points` and `cells` (SILO cell per point) must be aligned with spread.table rows.
drying_spec <- function(points, cells, rain.years, response, p.avg = 0.25) {
  nat <- points$art_nat == 0
  is.lw <- !is.na(points$lw_record_id)
  cls <- ifelse(nat & is.lw & !(points$group_assessment %in% "Always Wet"), "sometimes", "never")
  ref <- rain.years$ref$season.year
  cy <- rain.years$cell.year
  wet.of <- function(years) {
    sapply(seq_len(nrow(points)), function(j)
      cy$wet.mm[match(paste(cells[j], years), paste(cy$silo_cell, cy$season.year))])
  }
  W  <- wet.of(ref)     # reference years x points
  WL <- wet.of(ref - 1) # previous season
  stopifnot(!anyNA(W[, cls == "sometimes"]), !anyNA(WL[, cls == "sometimes"]))
  # anomalies from each point's own mean over the reference years, so p = p.avg in an average season
  x1 <- sweep(W, 2, colMeans(W)) / response$scale$wet["scale"]
  x2 <- sweep(WL, 2, colMeans(WL)) / response$scale$lag["scale"]
  list(class = cls, sometimes = cls == "sometimes", x1 = x1, x2 = x2, p.avg = p.avg,
       beta = response$beta, vcov = response$vcov, ref.years = ref)
}

# Drying probabilities for one run: rows = generations (reference years in order idx), cols = points.
draw_drying_matrix <- function(spec, idx, beta = spec$beta) {
  out <- matrix(0, nrow = length(idx), ncol = length(spec$class))
  s <- spec$sometimes
  if (any(s)) {
    out[, s] <- plogis(qlogis(spec$p.avg) + beta[1] * spec$x1[idx, s, drop = FALSE] +
                         beta[2] * spec$x2[idx, s, drop = FALSE])
  }
  out
}

# one run's coefficients: (b1, b2) ~ MVN(estimate, covariance), carrying the slope uncertainty
draw_drying_beta <- function(spec) {
  as.numeric(spec$beta + t(chol(spec$vcov)) %*% rnorm(2))
}
