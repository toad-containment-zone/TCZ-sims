# Script to load all background data required for TCZ sims

######## load functions and libraries ########
source("src/pprocess_functions.R") # functions for point process model
source("src/tcz-sims/get-invasion-front.R") # function for bringing in current invasion front
source("src/tcz-sims/build-tcz-units.R") # audit records -> intervention units for tcz-sims.R

######## load posteriors ########
load_posteriors() # the shared calibration, POSTERIORS_FILE in pprocess_functions.R

######## define the data directories ########
data.dir <- file.path(Sys.getenv("DATA_PATH"), "Toads/TCZ/infrastructure")
spatial.dir <- file.path(Sys.getenv("DATA_PATH"), "GIS - General", "GIS_layers_read_only/")

######## Load point data ########
# load TCZ site data (downloaded as geojson)
data.dump.id <- "089e7bef-6ea8-45d7-a819-a1c964d552fb"
tcz.sites <- st_read(file.path(data.dir, data.dump.id, "water_point_audit.geojson")) |> 
  mutate(origin_des = "Manmade", inside_tcz = TRUE) |> 
  select(-observer, -station_name, -traditional_owner_country, -whole_site_photo)

# load costing data (9 Jul 2026 activity-steps workbook, 4 crews)
cost.file <- file.path(data.dir, "optimisation-data/Activity Steps Toad project - 4 crew - 9Jul26.xlsx")
cost.sheet <- readxl::read_xlsx(path = cost.file, sheet = "Analysis", .name_repair = "unique_quiet")
## Project totals, for reporting. Cost/Site spreads only the direct infrastructure cost across
## sites, so site costs sum to the direct construction cost, not the total project cost.
cost.totals <- c(total.project = as.numeric(cost.sheet$`Total Project Cost`[1]),
                 direct.construction = as.numeric(cost.sheet$`Direct Construction Cost`[1]))
## Site-level costs (one row per audit record, keyed on ParentID = X_record_id)
tcz.costs <- cost.sheet |>
  filter(!is.na(ParentID)) |>
  select(1:6)
names(tcz.costs) <- tolower(make.names(names(tcz.costs)))

## Fence lengths for long fencing areas
fence.lengths <- readxl::read_xlsx(path = file.path(data.dir, "optimisation-data/Draft QRA data - Fencing unit rates.xlsx"),
                               sheet = "fence-lengths")
# fences with no Site.ID aren't in the audit yet (e.g. "Shamrock Stn Pivots", 2026-11-27): no location
# or cost, so they're left out until they have one
if (any(is.na(fence.lengths$Site.ID)))
  message("Skipping fence(s) with no Site.ID: ", paste(fence.lengths$Fence[is.na(fence.lengths$Site.ID)], collapse = ", "))
fence.lengths <- fence.lengths |> filter(!is.na(Site.ID))
# eradication: a fence's erad.units (crews of two + quad bike needed to cover it in a night) times the
# total cost per erad-unit per response (the sheet's day rate x n-days, from its
# Eradication-days-calculations sheet)
erad.sheet <- readxl::read_xlsx(path = file.path(data.dir, "optimisation-data/Draft QRA data - Fencing unit rates.xlsx"),
                                sheet = "eradication-costs-per-erad-unit")
erad.cost.per.unit <- erad.sheet$`cost-per erad-unit`[erad.sheet$Item %in% "Total-cost-per-erad-unit"]
stopifnot(length(erad.cost.per.unit) == 1, erad.cost.per.unit > 0, !anyNA(fence.lengths$erad.units))

# merge costs to infrastructure
tcz.sites <- tcz.sites |>
  left_join(tcz.costs, by = join_by(X_record_id == parentid)) |> 
  left_join(fence.lengths |> select(Site.ID, distance, n.gates), by = join_by(X_record_id == Site.ID))

# load background points (from Southwell et al)
old.lagrange.points <- st_read(file.path(spatial.dir, "Edited-layers/merged-points_rainfall_LaGrange.shp")) |> 
  select(fcsubtype_, full_name, perennia_1, origin_des, watercou_1, area_m, st_perimet) |> 
  st_transform(4326) # switch to WGS84 to match other data sources

## Note the manmade points in the TCZ here ^ have been replaced by the audit in tcz.sites
# these old.lagrange.points also do a bad job of identifying natural waterpoints in TCZ.  Need to replace with something else
# Living waters points make the most sense
# (Fulcrum export, Nov 2026: ground-truthing status and the drying survey's group_assessment)
living.waters <- st_read(file.path(Sys.getenv("DATA_PATH"), "Toads/TCZ/living-waters/living_waters_points.geojson")) |>
  select(lw_record_id = X_record_id, lw_status = X_status, group_assessment) |> # drop free-text and people fields
  mutate(origin_des = "Natural") |>
  st_transform(4326) # switch to WGS84 to match other data sources


# read in additional points from aerial imagery ()
d_extra <- st_read(file.path(spatial.dir, "additional-waterpoints/additional-waterpoints.kml")) |> 
  mutate(
    origin_des = case_when(
      # "point" in name → Manmade (catches watering point, waterpoint, etc.)
      str_detect(Name, regex("point", ignore_case = TRUE)) ~ "Manmade",
      # Natural water features
      str_detect(Name, regex("lagoon|lgoon|lake|geegully|creek|seep|spring|soak|salt", ignore_case = TRUE)) ~ "Natural",
      # Remaining manmade infrastructure
      str_detect(Name, regex("water|dam|tank|trough|irrigat|borrow|pond|ditch|resort|homestead|campsite|station|park", ignore_case = TRUE)) ~ "Manmade",
      TRUE ~ NA_character_
    )
  )
  # add additional records onto old.lagrange.points
old.lagrange.points <- bind_rows(old.lagrange.points, d_extra)


######## Load rasters ########
# load rainfall data (Number of days where at least 1mm of rain falls)
rainfall.raster <- terra::rast(file.path(spatial.dir, "annual-rainfall/rdann-1.asc"))

######## Load polygons ########
# load TCZ poly
tcz.boundary <- st_read(file.path(spatial.dir, "TCZ_boundary/Toad_Containment_Zone_Boundary_July 26.shp"))
# load wa coastline poly
wa.coast <- st_read(file.path(spatial.dir, "coast-poly_wa.shp")) |> 
  st_transform(crs = 4326) |>
  st_make_valid() |>
  st_intersection(old.lagrange.points |> st_transform(4326) |> st_bbox() |> st_as_sfc())
# load pastoral boundaries
pastoral.boundaries <- st_read(file.path(spatial.dir, "Pastoral_Stations_DPLH_083.gdb")) |> 
  st_transform(crs = 4326) |>
  st_make_valid() |>
  st_intersection(old.lagrange.points |> st_transform(4326) |> st_bbox() |> st_as_sfc())

# write out trimmed pastoral boundaries for use elsewhere
st_write(pastoral.boundaries, "out/pastoral-boundaries.shp", append = FALSE)

####### Merge and filter datasets #######
# 1. Flag which lagrange points fall inside the TCZ boundary
inside_tcz <- st_within(old.lagrange.points, tcz.boundary, sparse = FALSE)[, 1]
old.lagrange.points <- old.lagrange.points |>
  mutate(inside_tcz = inside_tcz)

# 2. Remove Manmade points that are inside the TCZ
lagrange.filtered <- old.lagrange.points |>
  filter(!inside_tcz) # remove all these points, natural and manmade

# 3. Intervention units for tcz-sims.R: audit records merged within 100 m and into fenced
#    enclosures, with summed costs (see build-tcz-units.R)
tcz.units <- build_tcz_units(tcz.sites, tcz.costs, fence.lengths) |>
  mutate(origin_des = "Manmade", inside_tcz = TRUE)
write.csv(st_drop_geometry(tcz.units), "out/tcz-units.csv", row.names = FALSE)
tcz.units <- tcz.units |> select(-X, -Y) |> st_transform(4326) # X/Y are re-added below

# 4. Rainfall, colonisation status and Albers coordinates for a merged point set
finish_points <- function(pts) {
  pts$inside_tcz <- st_within(pts, tcz.boundary, sparse = FALSE)[, 1]
  rainfall.vals <- terra::extract(rainfall.raster, terra::vect(pts))
  pts <- pts |> mutate(rainfall = rainfall.vals[[2]])
  colnsd <- score_colonised_multi_year(pts)
  pts <- colnsd$scored.points # includes colonised (0/1) and colonisation_year
  albers <- pts |> st_transform(crs = 3577) |> st_coordinates()
  list(points = bind_cols(pts, albers), fronts = colnsd$fronts, fp = colnsd$fp)
}

# all.points: every audit record as its own point (forecast.R, the do-nothing scenario)
fin <- finish_points(bind_rows(lagrange.filtered, tcz.sites, living.waters))
all.points <- fin$points
inv.front <- fin$fronts # one front line per modelled year
fp <- fin$fp
# tcz.points: TCZ audit records replaced by intervention units (tcz-sims.R)
tcz.points <- finish_points(bind_rows(lagrange.filtered, tcz.units, living.waters))$points
rm(fin)

####### Plot it #######
bbox <- st_bbox(all.points)

time.to.tcz.fig <-ggplot() +
  geom_sf(data = wa.coast, fill = NA, color = "grey30") +
  geom_sf(data = tcz.boundary, fill = NA, color = "red") +
  geom_sf(data = inv.front, aes(group = factor(year)), color = "red", lty = 2) +
  geom_sf(data = all.points, aes(color = factor(colonised))) +
  scale_color_manual(values = c("0" = "#132B43", "1" = "#56B1F7")) +
  geom_sf(data = fp, color = "red") +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]))
ggsave(filename = "out/time-to-tcz.png", plot = time.to.tcz.fig)

bbox <- st_bbox(tcz.boundary)
tcz.waterpoints.fig <- ggplot() +
  geom_sf(data = wa.coast, fill = NA, color = "grey30") +
  geom_sf(data = tcz.boundary, fill = NA, color = "red") +
  geom_sf(data = all.points, aes(color = origin_des)) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]))
ggsave(filename = "out/tcz-waterpoints.pdf", plot = tcz.waterpoints.fig)
