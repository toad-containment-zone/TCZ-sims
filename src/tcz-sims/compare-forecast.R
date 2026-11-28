# compare-forecast.R
#
# Compares the current forecast (out/forecast.RData, from forecast.R) with a baseline copy, by default
# the September 2026 forecast kept in out/forecast-baseline-2026-09/ (old λ, raster rain, one
# kernel, season.days = 365, r = 100 m, 100 m de-duplication). Phase 2 of docs/roadmap.md.
#
#   - TCZ arrival: per run, the year the first TCZ point is colonised (start year + `gen`, the
#     generation at which forecast.R's run stopped on its first target)
#   - per point: mean arrival year (conditional on colonisation in the run), for points within
#     MATCH_M of a point in the other set (the point sets differ: thinning replaced de-duplication)
#
# Writes out/forecast-comparison/ (tcz-arrival.csv, point-arrival.csv, two figures).

library(dplyr)
library(ggplot2)

start.year     <- 2026
baseline.file  <- "out/forecast-baseline-2026-09/forecast.RData"
current.file   <- "out/forecast.RData"
MATCH_M        <- 100
out.dir        <- "out/forecast-comparison"
dir.create(out.dir, showWarnings = FALSE)

load_forecast <- function(file) { e <- new.env(); load(file, envir = e); e$output }

tcz_arrival <- function(output, label) {
  gen <- vapply(output, function(x) if (is.null(x$gen)) NA_real_ else as.numeric(x$gen), numeric(1))
  data.frame(forecast = label, run = seq_along(gen), arrival = start.year + gen)
}

point_arrival <- function(output) {
  do.call(rbind, lapply(output, function(x) {
    m <- x$popmatrix[x$popmatrix[, "age"] > 0, c("ID", "X", "Y", "first.col.gen"), drop = FALSE]
    as.data.frame(m)
  })) |>
    group_by(ID) |>
    summarise(X = mean(X), Y = mean(Y), arrival = mean(start.year + first.col.gen), n.col = n(), .groups = "drop")
}

base <- load_forecast(baseline.file)
curr <- load_forecast(current.file)

######## TCZ arrival year ########
arr <- rbind(tcz_arrival(base, "baseline (2026-09)"), tcz_arrival(curr, "current"))
arr.summary <- arr |>
  group_by(forecast) |>
  summarise(runs = n(), never = sum(is.na(arrival)),
            median = median(arrival, na.rm = TRUE),
            q10 = quantile(arrival, 0.1, na.rm = TRUE), q90 = quantile(arrival, 0.9, na.rm = TRUE),
            p.by.2030 = mean(arrival <= 2030, na.rm = FALSE),
            p.by.2035 = mean(arrival <= 2035, na.rm = FALSE),
            .groups = "drop")
print(as.data.frame(arr.summary))
write.csv(arr.summary, file.path(out.dir, "tcz-arrival.csv"), row.names = FALSE)
print(table(arr$forecast, arr$arrival, useNA = "ifany"))

ggplot(arr |> filter(!is.na(arrival)), aes(x = arrival, fill = forecast)) +
  geom_bar(aes(y = after_stat(prop), group = forecast), position = position_dodge(preserve = "single")) +
  labs(x = "Year toads first reach the TCZ", y = "Share of runs", fill = NULL) +
  theme_bw()
ggsave(file.path(out.dir, "tcz-arrival.pdf"), width = 160, height = 100, units = "mm")

######## per-point mean arrival year, matched by location ########
pb <- point_arrival(base); pc <- point_arrival(curr)
nn <- vapply(seq_len(nrow(pc)), function(i) {
  d2 <- (pb$X - pc$X[i])^2 + (pb$Y - pc$Y[i])^2
  j <- which.min(d2)
  if (d2[j] <= MATCH_M^2) j else NA_integer_
}, integer(1))
pts <- pc |>
  mutate(arrival.base = pb$arrival[nn], n.col.base = pb$n.col[nn]) |>
  filter(!is.na(arrival.base)) |>
  mutate(shift = arrival - arrival.base)
cat(nrow(pts), "of", nrow(pc), "current points matched a baseline point within", MATCH_M, "m\n")
# points colonised at forecast start carry their real colonisation year in both: only the rest move
moving <- pts |> filter(arrival.base >= start.year | arrival >= start.year)
cat("Mean arrival-year shift (current - baseline) at points colonised after", start.year, ":",
    round(mean(moving$shift), 2), "years; median", median(moving$shift), "\n")
write.csv(pts, file.path(out.dir, "point-arrival.csv"), row.names = FALSE)

ggplot(moving, aes(X, Y, colour = shift)) +
  geom_point(size = 0.8) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b", name = "Arrival shift (y)") +
  coord_equal() +
  labs(title = "Mean arrival year, current - baseline", x = NULL, y = NULL) +
  theme_bw()
ggsave(file.path(out.dir, "point-arrival-shift.pdf"), width = 180, height = 160, units = "mm")
