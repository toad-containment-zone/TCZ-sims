# rain-years.R
#
# Year-by-year rainfall for the TCZ simulation region, from SILO DataDrill (0.05-degree cells):
#   - the 25-year reference sequence of simulated years, and which of them are "extreme" wet years
#     (top 20% by Dec-May rainfall averaged over the region's cells); the extreme years switch the
#     dispersal kernel to the extreme kernel (ms/tcz-optimisation-methods.qmd, "Modelled scenario")
#   - per point, the maximum rainy days (>= 1 mm) in any reference year, which sets that extreme kernel
#   - per cell and year, wet-season rainfall, kept for the step 5 drying matrix (which needs the two
#     years before the reference for its lag terms)
#
# Season year Y runs Jul (Y-1) to Jun (Y); its wet season is Dec (Y-1) to May (Y).
# Daily rain is cached in dat/silo-tcz-rain-cache.RData (gitignored); needs SILO_API_KEY for cells
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

fetch_silo_cell <- function(cell, lon, lat, start, end) {
  cat("  Fetching", cell, "\n")
  tryCatch(
    get_data_drill(longitude = lon, latitude = lat, start_date = start, end_date = end,
                   values = "rain", api_key = Sys.getenv("SILO_API_KEY")) |>
      as.data.frame() |>
      select(year, month, day, rainfall) |>
      mutate(silo_cell = cell),
    error = function(e) {
      warning("SILO fetch failed for ", cell, ": ", conditionMessage(e))
      NULL
    }
  )
}

# daily SILO rain for every cell under `points`, fetching only cells not already in the cache
silo_daily_for <- function(points) {
  fetch.start <- as.Date(sprintf("%d-07-01", min(rain.lag.years) - 1))
  fetch.end   <- as.Date(sprintf("%d-06-30", max(rain.ref.years)))
  cells <- distinct(silo_cells(points))
  daily <- data.frame()
  if (file.exists(silo.tcz.cache)) {
    load(silo.tcz.cache) # silo.tcz.daily, silo.tcz.range
    if (silo.tcz.range["start"] <= fetch.start && silo.tcz.range["end"] >= fetch.end) daily <- silo.tcz.daily
  }
  to.fetch <- cells |> filter(!silo_cell %in% unique(daily$silo_cell))
  cat(nrow(to.fetch), "of", nrow(cells), "SILO cells need fetching\n")
  if (nrow(to.fetch) > 0) {
    new <- Map(fetch_silo_cell, to.fetch$silo_cell, to.fetch$silo_lon, to.fetch$silo_lat,
               MoreArgs = list(start = fetch.start, end = fetch.end)) |> bind_rows()
    daily <- bind_rows(daily, new)
    silo.tcz.daily <- daily
    silo.tcz.range <- c(start = fetch.start, end = fetch.end)
    save(silo.tcz.daily, silo.tcz.range, file = silo.tcz.cache)
  }
  missing <- setdiff(cells$silo_cell, daily$silo_cell)
  if (length(missing) > 0) stop("no SILO data for ", length(missing), " cell(s), e.g. ", missing[1])
  daily |> filter(silo_cell %in% cells$silo_cell)
}

rain_years <- function(points) {
  daily <- silo_daily_for(points) |>
    mutate(season.year = ifelse(month >= 7, year + 1, year))

  cell.year <- daily |>
    filter(season.year %in% c(rain.lag.years, rain.ref.years)) |>
    group_by(silo_cell, season.year) |>
    summarise(rain.days = sum(rainfall >= 1),
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
