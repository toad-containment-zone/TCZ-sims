# Script to run sims from the current invasion front, until the tcz is reached (the "do nothing"
# forecast behind the public site).
#
# Same model structure as the TCZ runs (docs/roadmap.md, Phase 2): both kernels from SILO
# (attach_silo_kernels()), season.days = 212, a random reference-year start per run, the
# extreme-wet-year kernel in extreme years, and natural-point drying (build_drying_spec()). The
# control step is on only for drying: nothing here is a tank or fence.
#
# The public site (shiny/prep_app_data.R -> shiny/web/data/) is regenerated in Phase 5, with the
# final text, so this script no longer writes it unless UPDATE_SITE is TRUE.

######## load functions and libraries ########
source("src/tcz-sims/load-data.R")
source("src/tcz-sims/rain-years.R")
source("src/tcz-sims/drying-model.R")

if (!exists("UPDATE_SITE")) UPDATE_SITE <- FALSE
if (!exists("N_SIMS")) N_SIMS <- 1000

######## define target points and point flags ########
all.points <- all.points |>
  mutate(colonised = ifelse(st_within(geometry, tcz.boundary, sparse = FALSE)[, 1], 2L, colonised),
         natural = coalesce(origin_des == "Natural", FALSE),
         art_nat = as.integer(!natural), # setup() wants 1 = artificial, 0 = natural (natural points get the drying test)
         control.type = "none",          # do nothing: no tanks or fences
         keep.thin = inside_tcz & natural) # natural water inside the TCZ is never thinned away, as in tcz-setup.R

######## SILO kernels (shared with tcz-setup.R) ########
silo <- attach_silo_kernels(all.points)
all.points <- silo$points
rain.years <- silo$rain.years
rm(silo)

######## load simulation posteriors ########
load_posteriors() # the shared calibration, POSTERIORS_FILE in pprocess_functions.R

start_year <- 2026

######## tcz-arrival scenario ########
scen.name <- "forecast"
setup(point.data = all.points,
      X.id = "X",
      Y.id = "Y",
      present.id = "colonised",
      artificial.natural.id = "art_nat",
      rain.id = "rain.mean",            # SILO mean Nov-May rainy days
      rain.extreme.id = "rain.extreme", # builds pairs.extreme
      constant.rain = NULL,
      trunc.dist = TRUE,
      TCZ = FALSE,
      control.type.id = "control.type",
      keep.id = "keep.thin",
      col.year.id = "colonisation_year",
      start.year = start_year)
stopifnot(sum(spread.table[, "target"]) > 0, !anyNA(spread.table[, c("u", "v", "max.dist", "u.x", "v.x", "max.dist.x")]))

######## natural-point drying (shared with tcz-setup.R) ########
drying.spec <- build_drying_spec(all.points, spread.table[, "ID"], rain.years, living.waters, p.avg = 0.25)
cat("Drying classes:", paste(names(table(drying.spec$class)), table(drying.spec$class), collapse = ", "),
    "| sometimes-dry outside the TCZ:",
    sum(drying.spec$sometimes & !all.points$inside_tcz[spread.table[, "ID"]]), "\n")

# write out points for basemapping
write.csv(spread.table, file = "out/basemap_points.csv", row.names = FALSE)
# run sims..
sim_out <- run_sims(n.sims = N_SIMS, gens = 20, plot = FALSE, rollup = FALSE,
                    control.step = TRUE, # for drying only
                    pairs.extreme = pairs.extreme,
                    extreme.ref = rain.years$ref$extreme,
                    drying.spec = drying.spec)
save_outputs(output = sim_out,
             path = "out",
             scenario.name = scen.name,
             start.year=start_year,
             plot.time = TRUE,
             extinction.aware = TRUE)
make_plots(scen.name, plot.year = TRUE, tcz.boundary = tcz.boundary)

## Output data to the public site (Phase 5 only)
if (UPDATE_SITE) source("shiny/prep_app_data.R")
