# rain-years.R
#
# Year-by-year rainfall for the TCZ simulation region, from SILO DataDrill (0.05-degree cells):
#   - the 25-year reference sequence of simulated years, and which of them are "extreme" wet years
#     (top 20% by Dec-May rainfall averaged over the region's cells); the extreme years switch the
#     dispersal kernel to the extreme kernel (ms/tcz-optimisation-methods.qmd, "Modelled scenario")
#   - per point, the mean and maximum rainy days (>= 1 mm, Nov-May only: the movable season in
#     rain_to_days()) over the reference years, which set the standard and extreme kernels
#   - per cell and year, wet-season rainfall, kept for the step 5 drying matrix (which needs the two
#     years before the reference for its lag terms)
#
# Season year Y runs Jul (Y-1) to Jun (Y); its wet season is Dec (Y-1) to May (Y).
# Daily rain is cached in dat/silo-tcz-rain-cache.RData (gitignored; rain_years(cache=) picks
# another file for another point set); needs SILO_API_KEY for cells
# not yet cached. Cell keying matches src/null-model/predict-front-displacement.R.

library(sf)
library(dplyr)
library(weatherOz)

rain.ref.years  <- 2002:2026 # 25 most recent complete season years
rain.lag.years  <- 2000:2001 # extra years for the drying model's lag terms (step 5)
rain.extreme.q  <- 0.8       # top 20% of reference years by wet-season rainfall are "extreme"
silo.cell.res   <- 0.05
silo.tcz.cache  <- "dat/silo-tcz-rain-cache.RData"

silo_cells <- function(points) {
  ll <- st_coordinates(st_transform(points, 4326))
  lon <- round(ll[, 1] / silo.cell.res) * silo.cell.res
  lat <- round(ll[, 2] / silo.cell.res) * silo.cell.res
  # fixed-width key: round()*0.05 can leave floating-point crumbs (122.05000000000001)
  data.frame(silo_cell = sprintf("%.2f_%.2f", lon, lat), silo_lon = lon, silo_lat = lat)
}

fetch_silo_cell <- function(cell, lon, lat, start, end, tries = 3) {
  cat("  Fetching", cell, "\n")
  # SILO occasionally fails or returns no rows under load (239 of 4,325 Kimberley cells on
  # 2026-11-27, all of which returned data on a retry), so retry with a pause before giving up
  for (attempt in seq_len(tries)) {
    out <- tryCatch(
      get_data_drill(longitude = lon, latitude = lat, start_date = start, end_date = end,
                     values = "rain", api_key = Sys.getenv("SILO_API_KEY")) |>
        as.data.frame() |>
        select(year, month, day, rainfall) |>
        mutate(silo_cell = cell),
      error = function(e) {
        warning("SILO fetch failed for ", cell, " (attempt ", attempt, "): ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(out) && nrow(out) > 0) return(out)
    if (attempt < tries) Sys.sleep(5 * attempt)
  }
  NULL
}

# daily SILO rain for every cell under `points`, fetching only cells not already in the cache.
# `cache` lets another point set (e.g. the Kimberley ABC, src/ABC/kimberley-setup.R) keep its own
# file; the cache is saved every `save.every` fetched cells so a long fetch can be resumed
silo_daily_for <- function(points, cache = silo.tcz.cache, save.every = 100) {
  fetch.start <- as.Date(sprintf("%d-07-01", min(rain.lag.years) - 1))
  fetch.end   <- as.Date(sprintf("%d-06-30", max(rain.ref.years)))
  cells <- distinct(silo_cells(points))
  daily <- data.frame()
  silo.tcz.nodata <- character(0) # cells SILO has no data for, even on retry (not refetched)
  if (file.exists(cache)) {
    load(cache) # silo.tcz.daily, silo.tcz.range, and silo.tcz.nodata if any
    if (silo.tcz.range["start"] <= fetch.start && silo.tcz.range["end"] >= fetch.end) daily <- silo.tcz.daily
  }
  to.fetch <- cells |> filter(!silo_cell %in% unique(daily$silo_cell), !silo_cell %in% silo.tcz.nodata)
  cat(nrow(to.fetch), "of", nrow(cells), "SILO cells need fetching\n")
  if (nrow(to.fetch) > 0) {
    chunks <- split(seq_len(nrow(to.fetch)), ceiling(seq_len(nrow(to.fetch)) / save.every))
    for (ch in chunks) {
      new <- Map(fetch_silo_cell, to.fetch$silo_cell[ch], to.fetch$silo_lon[ch], to.fetch$silo_lat[ch],
                 MoreArgs = list(start = fetch.start, end = fetch.end)) |> bind_rows()
      daily <- bind_rows(daily, new)
      silo.tcz.daily <- daily
      silo.tcz.range <- c(start = fetch.start, end = fetch.end)
      save(silo.tcz.daily, silo.tcz.range, silo.tcz.nodata, file = cache)
      cat("  saved", length(unique(daily$silo_cell)), "cells to", cache, "\n")
    }
    silo.tcz.nodata <- union(silo.tcz.nodata, setdiff(to.fetch$silo_cell, daily$silo_cell))
    save(silo.tcz.daily, silo.tcz.range, silo.tcz.nodata, file = cache)
  }
  # a few cells have no SILO data at all, even on retry: small islands and coastal cells outside its
  # land mask (5 around King Sound and the Buccaneer and Bonaparte archipelagos, 21 Kimberley points,
  # 2026-11-27).
  # They take the nearest cell with data, within max.sub.cells cells; the substitution is redone each
  # run, not cached
  missing <- setdiff(cells$silo_cell, daily$silo_cell)
  if (length(missing) > 0) {
    have <- distinct(cells |> filter(!silo_cell %in% missing))
    max.sub.cells <- 6 # 0.3 deg: one offshore-island point (124.20_-15.35) is 0.255 deg from data
    for (m in missing) {
      mm <- cells[match(m, cells$silo_cell), ]
      dd <- sqrt((have$silo_lon - mm$silo_lon)^2 + (have$silo_lat - mm$silo_lat)^2)
      if (min(dd) > max.sub.cells * silo.cell.res + 1e-9)
        stop("no SILO data for ", m, " and no cell with data within ", max.sub.cells, " cells")
      src <- have$silo_cell[which.min(dd)]
      cat("  no SILO data for", m, "-- using", src, sprintf("(%.2f deg away)\n", min(dd)))
      daily <- bind_rows(daily, daily |> filter(silo_cell == src) |> mutate(silo_cell = m))
    }
  }
  daily |> filter(silo_cell %in% cells$silo_cell)
}

rain_years <- function(points, cache = silo.tcz.cache) {
  daily <- silo_daily_for(points, cache) |>
    mutate(season.year = ifelse(month >= 7, year + 1, year))

  cell.year <- daily |>
    filter(season.year %in% c(rain.lag.years, rain.ref.years)) |>
    group_by(silo_cell, season.year) |>
    summarise(rain.days = sum(rainfall >= 1 & month %in% c(11, 12, 1:5)), # Nov-May: the movable season in rain_to_days()
              wet.mm = sum(rainfall[month %in% c(12, 1:5)]),
              n.days = n(), .groups = "drop")
  stopifnot(all(cell.year$n.days >= 365))

  ref <- cell.year |>
    filter(season.year %in% rain.ref.years) |>
    group_by(season.year) |>
    summarise(wet.mm = mean(wet.mm), rain.days = mean(rain.days), .groups = "drop") |>
    arrange(season.year)
  ref$extreme <- ref$wet.mm > quantile(ref$wet.mm, rain.extreme.q)
  stopifnot(nrow(ref) == length(rain.ref.years), sum(ref$extreme) == round(length(rain.ref.years) * (1 - rain.extreme.q)))

  cell.max <- cell.year |>
    filter(season.year %in% rain.ref.years) |>
    group_by(silo_cell) |>
    summarise(rain.days.max = max(rain.days), rain.days.mean = mean(rain.days), .groups = "drop")
  pts <- silo_cells(points) |> left_join(cell.max, by = "silo_cell")

  cat("Reference years", min(ref$season.year), "-", max(ref$season.year), "| extreme:",
      paste(ref$season.year[ref$extreme], collapse = ", "), "\n")
  list(ref = ref, cell.year = cell.year,
       point.cell = pts$silo_cell,
       point.rain.days.max = pts$rain.days.max,
       point.rain.days.mean = pts$rain.days.mean)
}

# The shared SILO kernel inputs (docs/roadmap.md, "One calibration and model structure"): adds each
# point's mean (`rain.mean`, the standard kernel) and maximum (`rain.extreme`, the extreme kernel)
# Nov-May rainy days over the reference years. Used by tcz-setup.R, forecast.R and the null model,
# which then pass rain.id = "rain.mean", rain.extreme.id = "rain.extreme" to setup().
attach_silo_kernels <- function(points, cache = silo.tcz.cache) {
  ry <- rain_years(points, cache = cache)
  points$rain.mean <- ry$point.rain.days.mean
  points$rain.extreme <- ry$point.rain.days.max
  stopifnot(!anyNA(points$rain.mean), !anyNA(points$rain.extreme))
  list(points = points, rain.years = ry)
}
