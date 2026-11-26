# build-tcz-units.R
#
# Turns TCZ infrastructure-audit records into "intervention units": one simulated waterpoint per
# unit, carrying the unit's total cost and control type. The greedy optimiser removes interventions
# one unit at a time, so every intervention has to be its own point with the right cost.
#
# Why units rather than raw audit records:
#   - setup()'s 100 m spatial de-duplication keeps only the first of any records closer than 100 m,
#     silently dropping the rest along with their costs. Merging those records into one unit first
#     leaves setup() nothing to drop.
#   - A fenced area can enclose several audit records (bores, tanks, the homestead) plus the
#     fence-line records that carry the fence's cost. They all form one fenced unit: toads that get
#     past the fence reach all of them, and the fence is the one intervention.
#
# Rules:
#   1. Records listed in `fence.enclosures` go into their fenced unit, placed at the centroid of
#      the enclosed water points (fence-line records carry cost but not location).
#   2. Remaining records closer than `merge.dist` (single linkage) are merged, placed at their
#      centroid.
#   3. Unit cost = sum of members' Cost/Site (direct infrastructure cost). Control type = fence if
#      any member is on the fence-lengths sheet, else tank. "Nothing"-status records (no works)
#      are controlled points with $0 cost.
#   4. candidate = cost > 0: only units with a cost are ever removed by the greedy optimiser.
#
# Inputs:  tcz.sites (sf, the audit, with X_record_id/name_of_site/X_status),
#          costs (ParentID + Cost/Site), fence.lengths (Site.ID/Fence/distance/n.gates/erad.units)
# Output:  sf POINT (Albers, EPSG:3577) with one row per unit.

library(sf)
library(dplyr)

# Fenced areas whose fence encloses other audit records. Keyed by X_record_id; `name` is only a
# check that the IDs still point at the expected records. `role`: "anchor" = the record the
# fence-lengths sheet names; "enclosed" = a water point inside the fence; "fence-line" = a record
# for the fence itself (holds cost, no water). Membership confirmed by the project lead 2026-11-26.
fence.enclosures <- tribble(
  ~fence,             ~X_record_id,                           ~name,                               ~role,
  "Shelemar Station", "4b523d19-0507-4f30-a6fa-57198c64b1db", "Homestead",                         "anchor",
  "Shelemar Station", "b0494d15-3276-4a99-bfe8-5062215a205f", "P6 bore",                           "enclosed",
  "Shelemar Station", "1726d2c6-8086-4a9e-815b-29ee7fbcf26a", "P3 bore",                           "enclosed",
  "Shelemar Station", "84f4db18-592b-4465-ae3d-aad472f707a6", "Roper 1 bore",                      "enclosed",
  "Shelemar Station", "58ebddd8-4b98-4c82-8d44-c28ccacae11d", "Roper 2 bore",                      "enclosed",
  "Shelemar Station", "c03563db-a23d-42f2-889e-d656ab5402fe", "Typical valve installation",        "enclosed",
  "Shelemar Station", "62bcba07-9446-461b-aba0-644ab6b24354", "Perimeter fence",                   "fence-line",
  # assumed part of the same fence: no cost of its own, ~6.8 km from "Perimeter fence"
  # (about one side of the 27.6 km fence's square)
  "Shelemar Station", "8f4d571c-d68e-4ada-98f4-ba8b55445225", "Shelemar Perimeter Fence West - 3", "fence-line"
)

build_tcz_units <- function(tcz.sites, costs, fence.lengths, enclosures = fence.enclosures,
                            merge.dist = 100) {
  sites <- tcz.sites |>
    st_transform(3577) |>
    select(X_record_id, name_of_site, X_status) |>
    mutate(cost = costs$cost.site[match(X_record_id, costs$parentid)],
           works = X_status != "Nothing")

  # -- checks on the inputs --
  stopifnot(!anyDuplicated(sites$X_record_id), !anyDuplicated(costs$parentid))
  enc <- enclosures |> mutate(found = name_of_site_lookup(sites, X_record_id))
  bad <- enc |> filter(is.na(found) | found != name)
  if (nrow(bad) > 0) stop("fence.enclosures IDs no longer match the audit: ", paste(bad$name, collapse = ", "))
  stopifnot(all(enc$X_record_id[enc$role == "anchor"] %in% fence.lengths$Site.ID))

  # -- 1. enclosure units --
  sites$unit_key <- NA_character_
  sites$unit_key[match(enc$X_record_id, sites$X_record_id)] <- paste0("fence:", enc$fence)
  sites$locates <- TRUE # whether this record contributes to its unit's location
  sites$locates[match(enc$X_record_id[enc$role == "fence-line"], sites$X_record_id)] <- FALSE

  # -- 2. single-linkage merge of everything else within merge.dist --
  rest <- which(is.na(sites$unit_key))
  if (length(rest) > 1) {
    d <- as.matrix(st_distance(sites[rest, ])) |> units::drop_units()
    grp <- cutree(hclust(as.dist(d), method = "single"), h = merge.dist - 1e-9)
  } else grp <- seq_along(rest)
  sites$unit_key[rest] <- paste0("near:", grp)

  # -- 3. one row per unit --
  xy <- st_coordinates(sites)
  sites$x <- xy[, 1]; sites$y <- xy[, 2]
  fence.row <- match(sites$X_record_id, fence.lengths$Site.ID)
  sites$fence.distance <- fence.lengths$distance[fence.row]
  sites$fence.gates <- fence.lengths$n.gates[fence.row]
  sites$fence.erad.units <- fence.lengths$erad.units[fence.row]

  units <- sites |>
    st_drop_geometry() |>
    group_by(unit_key) |>
    summarise(
      record_ids = paste(X_record_id, collapse = ";"),
      names = paste(name_of_site, collapse = "; "),
      n_records = n(),
      n_works = sum(works),
      n_uncosted_works = sum(works & is.na(cost)),
      cost = sum(cost, na.rm = TRUE),
      control.type = if (any(!is.na(fence.distance))) "fence" else "tank",
      distance = if (any(!is.na(fence.distance))) sum(fence.distance, na.rm = TRUE) else NA_real_,
      n.gates = if (any(!is.na(fence.gates))) sum(fence.gates, na.rm = TRUE) else NA_real_,
      erad.units = if (any(!is.na(fence.erad.units))) sum(fence.erad.units, na.rm = TRUE) else NA_real_,
      X = mean(x[locates]),
      Y = mean(y[locates]),
      .groups = "drop"
    ) |>
    arrange(desc(cost), unit_key) |>
    mutate(unit_id = sprintf("U%03d", row_number()),
           candidate = cost > 0) |>
    relocate(unit_id)

  # -- checks on the output --
  member.ids <- unlist(strsplit(units$record_ids, ";"))
  stopifnot(setequal(member.ids, sites$X_record_id), !anyDuplicated(member.ids))
  stopifnot(isTRUE(all.equal(sum(units$cost), sum(sites$cost, na.rm = TRUE))))
  ud <- as.matrix(dist(units[, c("X", "Y")])); diag(ud) <- Inf
  if (any(ud < merge.dist)) stop("two units are closer than merge.dist (", merge.dist, " m)")
  stopifnot(all(units$control.type[grepl("^fence:", units$unit_key)] == "fence"))

  # -- report --
  orphan.costs <- costs |> filter(!parentid %in% sites$X_record_id)
  cat(sprintf("TCZ units: %d from %d audit records (%d merged within %g m, %d into fenced enclosures)\n",
              nrow(units), nrow(sites),
              sum(units$n_records[grepl("^near:", units$unit_key)] - 1), merge.dist,
              sum(units$n_records[grepl("^fence:", units$unit_key)] - 1)))
  cat(sprintf("  %d fence, %d tank | %d candidates (cost > 0), total candidate cost $%s\n",
              sum(units$control.type == "fence"), sum(units$control.type == "tank"),
              sum(units$candidate), format(round(sum(units$cost)), big.mark = ",")))
  uncosted <- units |> filter(n_uncosted_works > 0)
  if (nrow(uncosted) > 0) {
    cat("  Units containing works with no cost:\n")
    print(as.data.frame(uncosted |> select(unit_id, names, n_uncosted_works, cost, candidate)), row.names = FALSE)
  }
  if (nrow(orphan.costs) > 0) {
    cat("  Cost rows matching no audit record (excluded):\n")
    print(as.data.frame(orphan.costs |> select(station, site, cost.site)), row.names = FALSE)
  }

  units |>
    select(-unit_key) |>
    st_as_sf(coords = c("X", "Y"), crs = 3577, remove = FALSE)
}

name_of_site_lookup <- function(sites, ids) sites$name_of_site[match(ids, sites$X_record_id)]
