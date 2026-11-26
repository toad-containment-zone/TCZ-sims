# Shared setup for the TCZ scenarios (tcz-sims.R, tcz-optimise.R): loads data, defines the
# simulation region and colonised/target points, runs setup(), and checks the intervention units.
# Leaves spread.table, pairs, pairs.extreme, rain.years and tcz.units.sim (unit metadata per
# spread.table row) in the session.

######## load functions and libraries ########
source("src/tcz-sims/load-data.R") # provides tcz.points: TCZ audit records as intervention units (build-tcz-units.R)

######## filter to TCZ region ########
# buffer TCZ by 80km (project to Albers for metric distances, then back)
tcz_buffer_80 <- tcz.boundary |>
  st_transform(3577) |>
  st_buffer(80000) |>
  st_union() |>
  st_transform(st_crs(tcz.points))

tcz.points <- tcz.points[st_within(tcz.points, tcz_buffer_80, sparse = FALSE)[, 1], ]

######## define colonised and target points ########
# Use TCZ bounding box to define north/south: 
#   north of TCZ → already colonised (colonised = 1)
#   south of TCZ → target points    (colonised = 2)
tcz_bbox   <- st_bbox(st_transform(tcz.boundary, st_crs(tcz.points)))
pt_y       <- st_coordinates(tcz.points)[, "Y"]
within_tcz <- st_within(tcz.points, tcz.boundary, sparse = FALSE)[, 1]

tcz.points <- tcz.points |>
  mutate(
    colonised = case_when(
      within_tcz                 ~ 0L,  # inside TCZ → uncolonised simulation space
      pt_y > tcz_bbox[["ymax"]] ~ 1L,  # north of TCZ → already colonised
      pt_y < tcz_bbox[["ymin"]] ~ 2L,  # south of TCZ → target points
      TRUE ~ colonised
    ),
    origin_des = (origin_des == "Natural"), # to match logical required of setup()
    control.type = case_when(
      within_tcz & !origin_des ~ control.type, # TCZ intervention unit: "fence" or "tank" from build_tcz_units()
      TRUE                     ~ "none"        # natural, or outside TCZ
    ),
    fence.units = coalesce(distance, 0) / 1000 + coalesce(n.gates, 0), # distance is metres -> km, + gate count
    fence.area = if_else(control.type == "fence", (distance / 4)^2, NA_real_), # distance = full perimeter (m); assume square -> side = distance/4, area = side^2
    art_nat = as.integer(!origin_des) # setup() wants 1 = artificial, 0 = natural (natural points get the drying test)
  ) |>
  rename(TCZ = inside_tcz) # to match required name in setup()

# Plot to make sure all is as it should be
bbox <- st_bbox(tcz_buffer_80)

ggplot() +
  geom_sf(data = wa.coast, fill = NA, color = "grey30") +
  geom_sf(data = tcz_buffer_80, fill = NA, color = "blue", lty = 2) +
  geom_sf(data = tcz.boundary, fill = NA, color = "red") +
  geom_sf(data = tcz.points,
          aes(color = factor(colonised,
                             levels = c(0, 1, 2),
                             labels = c("Uncolonised (TCZ)", "Colonised (N)", "Target (S)")),
              shape = origin_des),
          size = 0.8) +
  scale_color_manual(values = c("Uncolonised (TCZ)" = "grey50",
                                "Colonised (N)"     = "firebrick",
                                "Target (S)"        = "steelblue"),
                     name = "Status") +
  scale_shape_manual(values = c(`TRUE` = 8, `FALSE` = 16),
                     labels = c("Natural", "Manmade"),
                     name = "Origin") +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]))


######## year-by-year rainfall: extreme wet years and extreme-year rainy days (rain-years.R) ########
source("src/tcz-sims/rain-years.R")
rain.years <- rain_years(tcz.points)
# both kernels come from SILO: the standard kernel from each point's mean rainy days over the
# reference years (replacing the rainfall raster used by forecast.R), the extreme kernel from its
# maximum (project lead's decision, 2026-11-26)
tcz.points$rain.mean <- rain.years$point.rain.days.mean
tcz.points$rain.extreme <- rain.years$point.rain.days.max
stopifnot(!anyNA(tcz.points$rain.mean), !anyNA(tcz.points$rain.extreme))

######## load simulation posteriors ########
load("dat/Posteriors_2026.RData")


######## full-implementation point set ########
setup(point.data = tcz.points, 
      X.id = "X", 
      Y.id = "Y",
      present.id = "colonised",
      artificial.natural.id = "art_nat",
      rain.id = "rain.mean", # SILO mean rainy days (not the raster's `rainfall`)
      threshold = 100,
      constant.rain = NULL,
      trunc.dist = TRUE,
      TCZ = FALSE, # keep TCZ artificial points in: they're modelled as tanks/fences via the control step, not removed
      control.type.id = "control.type",
      fence.area.id = "fence.area",
      fence.units.id = "fence.units",
      rain.extreme.id = "rain.extreme") # builds pairs.extreme

# every TCZ intervention unit must survive setup()'s spatial de-duplication, and nothing outside
# the TCZ may carry a control
stopifnot(all(!is.na(tcz.points$unit_id[!is.na(tcz.points$control.type) & tcz.points$control.type != "none"])),
          sum(tcz.points$control.type != "none") == nrow(tcz.units))
sim.unit <- tcz.points$unit_id[spread.table[, "ID"]]
if (!setequal(na.omit(sim.unit), tcz.units$unit_id)) {
  stop("TCZ units dropped by setup() de-duplication: ",
       paste(setdiff(tcz.units$unit_id, sim.unit), collapse = ", "))
}
stopifnot(all(spread.table[is.na(sim.unit), "control.type"] == CONTROL_NONE))
# unit metadata aligned row-for-row with spread.table (NA rows = not a TCZ unit), for the optimiser
tcz.units.sim <- st_drop_geometry(tcz.units)[match(sim.unit, tcz.units$unit_id),
                                             c("unit_id", "cost", "candidate", "control.type")]
