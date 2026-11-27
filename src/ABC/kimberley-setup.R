# kimberley-setup.R
#
# Builds the Kimberley ABC setup with the same model structure as the TCZ runs
# (docs/roadmap.md, Phase 1):
#   - both kernels from SILO, by rain_years(): the standard kernel from each cell's mean rainy days
#     over the 2002-2026 reference years, the extreme kernel from the maximum;
#   - a fixed year sequence following the real invasion: generation i is season year 2009 + i, so the
#     14 generations from the 2009 start run 2010-2023, using the extreme kernel in the years that
#     were extreme by the top-20% rule;
#   - no natural-point drying (no TO assessment in the Kimberley), so the control step is off and
#     the ABC can use rollup.
#
# Writes out/setup_kimberley_silo.RData (spread.table, pairs, pairs.extreme, extreme.years).
# The SILO fetch for the ~4,300 cells is cached in dat/silo-kimberley-rain-cache.RData.
# Large point set: set R_MAX_VSIZE high (usethis::edit_r_environ()) before running.

source("src/pprocess_functions.R")
source("src/tcz-sims/rain-years.R")

kim.start.year <- 2009
kim.gens <- 14
kim.season.years <- kim.start.year + seq_len(kim.gens) # 2010-2023

k.points <- read.csv("dat/waterpoint-data_Kimberley_trimmed.csv")
k.points$origin_des <- as.numeric(as.factor(k.points$origin_des)) - 1 # Manmade = 0, Natural = 1 (unused: no control step)
k.points$colonised[is.na(k.points$colonised)] <- 0
k.points$colonised[k.points$colonised == 2009] <- 1

rain.kim <- rain_years(st_as_sf(k.points, coords = c("X", "Y"), crs = 3577, remove = FALSE),
                       cache = "dat/silo-kimberley-rain-cache.RData")
k.points$rain.mean <- rain.kim$point.rain.days.mean
k.points$rain.extreme <- rain.kim$point.rain.days.max
stopifnot(!anyNA(k.points$rain.mean), !anyNA(k.points$rain.extreme))

extreme.years <- rain.kim$ref$extreme[match(kim.season.years, rain.kim$ref$season.year)]
stopifnot(length(extreme.years) == kim.gens, !anyNA(extreme.years))
cat("Extreme season years in", min(kim.season.years), "-", max(kim.season.years), ":",
    paste(kim.season.years[extreme.years], collapse = ", "), "\n")

# step 1 diagnostic: the old calibration's ndays_1 against SILO mean rainy days at the same points
cat(sprintf("ndays_1 mean %.1f | SILO mean %.1f | SILO max %.1f | ratio SILO/ndays_1 %.2f | r = %.2f\n",
            mean(k.points$ndays_1), mean(k.points$rain.mean), mean(k.points$rain.extreme),
            mean(k.points$rain.mean) / mean(k.points$ndays_1), cor(k.points$ndays_1, k.points$rain.mean)))
cat("Movement days, standard kernel:\n"); print(summary(rain_to_days(k.points$rain.mean)))
cat("Movement days, extreme kernel:\n");  print(summary(rain_to_days(k.points$rain.extreme)))

setup(point.data = k.points,
      X.id = "X",
      Y.id = "Y",
      present.id = "colonised",
      artificial.natural.id = "origin_des",
      rain.id = "rain.mean",
      rain.extreme.id = "rain.extreme",
      observations = c("X2021_pres", "X2022_pres", "X2023_pres"),
      constant.rain = NULL,
      trunc.dist = TRUE,
      TCZ = FALSE)
rm(k.points); gc()

cat("pairs:", format(object.size(pairs), units = "MB"),
    "| pairs.extreme:", format(object.size(pairs.extreme), units = "MB"), "\n")
save(spread.table, pairs, pairs.extreme, extreme.years, kim.start.year, kim.gens,
     file = "out/setup_kimberley_silo.RData")
