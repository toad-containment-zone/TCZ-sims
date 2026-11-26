# Script to run sims over the TCZ: the full-implementation scenario (every intervention unit in place)
# The design optimisation is in tcz-optimise.R; both share tcz-setup.R.

source("src/tcz-sims/tcz-setup.R")

scen.name <- "tcz-full"
# write out points for basemapping
write.csv(spread.table, file = "out/basemap_points.csv", row.names = FALSE)
# run sims..
# stop.on.target = FALSE runs the full 50 years after any breach, so point-years occupied are counted over the whole horizon
sim_out <- run_sims(n.sims = 100, gens = 50, plot = FALSE, rollup = FALSE, control.step = TRUE, stop.on.target = FALSE,
                    pairs.extreme = pairs.extreme, extreme.ref = rain.years$ref$extreme) # extreme-wet-year kernel, random start per rep
save_outputs(output = sim_out,
             path = "out",
             scenario.name = scen.name,
             start.year=2032,
             plot.time = TRUE,
             extinction.aware = TRUE)
make_plots(scen.name, plot.year = TRUE, tcz.boundary = tcz.boundary)
