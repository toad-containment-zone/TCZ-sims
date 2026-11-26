# export_web_data.R
#
# Convert shiny/data/app_data.rds (written by prep_app_data.R) into the static
# files the JavaScript forecast site reads (shiny/web/data/):
#
#   forecast.json      points, per-year colonisation probabilities, metadata
#   tcz.geojson        TCZ boundary polygon
#   pastoral.geojson   pastoral lease boundaries (geometry only)
#
# prep_app_data.R sources this at the end, so a normal forecast run refreshes
# both the .rds (used by the legacy Shiny app, app.R) and the web data. It can
# also be run on its own from the repo root:
#   Rscript shiny/export_web_data.R

suppressPackageStartupMessages({
  library(sf)
  library(jsonlite)
})

export_web_data <- function(rds     = "shiny/data/app_data.rds",
                            out_dir = "shiny/web/data") {
  stopifnot("app_data.rds not found -- run from the repo root" = file.exists(rds))
  d <- readRDS(rds)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- Points: one row per waterpoint, in all_points order ----
  pts <- d$all_points
  pd  <- d$prob_data
  arrival <- pd$mean_arrival[match(pts$ID, pd$ID)]

  # ---- Probabilities: one vector per year, aligned to pts ----
  # A point missing from prob_data for a year means probability 0.
  prob <- lapply(d$years, function(yr) {
    rows <- pd[pd$year == yr, , drop = FALSE]
    p <- rows$prob[match(pts$ID, rows$ID)]
    p[is.na(p)] <- 0
    round(p, 4)
  })

  b <- d$initial_bounds
  start_year <- if (!is.null(d$start_year)) d$start_year else min(d$years)

  forecast <- list(
    updated       = format(d$updated, "%Y-%m-%d"),
    updated_label = format(d$updated, "%B %Y"),
    years         = d$years,
    start_year    = start_year,
    n_sims        = d$n_sims,
    bounds        = list(c(b$lat1, b$lng1), c(b$lat2, b$lng2)),
    points        = list(lon = round(pts$lon, 5), lat = round(pts$lat, 5),
                         arrival = arrival),
    prob          = prob
  )
  write_json(forecast, file.path(out_dir, "forecast.json"),
             auto_unbox = TRUE, digits = NA)

  # ---- Boundaries as GeoJSON (WGS84, ~10 m coordinate precision) ----
  write_geojson <- function(x, path) {
    x <- st_sf(geometry = st_geometry(st_transform(x, 4326)))
    if (file.exists(path)) unlink(path)
    st_write(x, path, driver = "GeoJSON", quiet = TRUE,
             layer_options = c("COORDINATE_PRECISION=4", "RFC7946=YES"))
  }
  tcz <- d$tcz_boundary
  write_geojson(tcz, file.path(out_dir, "tcz.geojson"))
  write_geojson(d$pastoral_boundaries, file.path(out_dir, "pastoral.geojson"))

  sizes <- file.size(file.path(out_dir, c("forecast.json", "tcz.geojson", "pastoral.geojson")))
  cat(
    "Wrote web data to ", out_dir, "\n",
    "  points: ", nrow(pts), ", years: ", length(d$years), "\n",
    "  sizes (KB): forecast ", round(sizes[1] / 1024), ", tcz ", round(sizes[2] / 1024),
    ", pastoral ", round(sizes[3] / 1024), "\n",
    sep = ""
  )
  invisible(out_dir)
}

export_web_data()
