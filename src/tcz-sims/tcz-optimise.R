# tcz-optimise.R
#
# Optimises the TCZ design (ms/tcz-optimisation-methods.qmd, "Cost optimisation"):
# With the extreme-wet-year kernel (rain-years.R) and fence eradication responses counted.
#   1. do-nothing counterfactual (every intervention removed) and full design, 100 full-horizon runs each
#   2. backward greedy search (most expensive first, 20 runs per test) + 100-run confirmation,
#      restoring removed units in reverse order while the confirmation finds breaches
#   3. optimised design, 100 full-horizon runs; breach probability, point-years avoided, cost and CE
# Outputs go to out/optimisation/.
#
# Set SMOKE <- TRUE before sourcing to check the plumbing quickly (5 candidates, tiny run counts).

if (!exists("SMOKE")) SMOKE <- FALSE

source("src/tcz-sims/tcz-setup.R")
source("src/tcz-sims/optimise-functions.R")

out.dir <- if (SMOKE) "out/optimisation/smoke" else "out/optimisation"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

n.eval    <- if (SMOKE) 4 else 100 # full-horizon evaluation runs per design
n.screen  <- if (SMOKE) 2 else 20
n.confirm <- if (SMOKE) 4 else 100
max.cand  <- if (SMOKE) 5 else Inf
seed      <- 20261126
eradication.cost <- 0 # $ per fence eradication response; no figure yet

st <- spread.table
tcz.rows <- !is.na(tcz.units.sim$unit_id)
none <- rep(FALSE, nrow(st))

# sanity checks on the design machinery
stopifnot(identical(apply_design(st, none), st),
          all(apply_design(st, tcz.rows)[, "control.type"] == CONTROL_NONE))

summarise_design <- function(name, ev, remove) {
  data.frame(design = name,
             units.retained = sum(tcz.rows & !remove),
             cost = design_cost(tcz.units.sim, remove, eradication.cost = eradication.cost,
                                eradication.events = mean(ev$n.erad)),
             erad.events = mean(ev$n.erad),
             runs = nrow(ev),
             breaches = sum(ev$breach),
             p.breach = mean(ev$breach),
             p.breach.upper95 = binom.test(sum(ev$breach), nrow(ev))$conf.int[2],
             py.total = mean(ev$py.total),
             py.tcz = mean(ev$py.tcz))
}
eval_full <- function(remove, seed.offset) {
  eval_design(apply_design(st, remove), pairs, n.eval, stop.on.target = FALSE,
              seed = seed + seed.offset, tcz.rows = tcz.rows,
              pairs.extreme = pairs.extreme, extreme.ref = rain.years$ref$extreme)
}

######## 1. counterfactual and full design ########
t0 <- Sys.time()
cat("Do-nothing counterfactual...\n")
ev.nothing <- eval_full(tcz.rows, 2e6)
cat("Full design...\n")
ev.full <- eval_full(none, 2e6)

######## 2. greedy search ########
cat("Greedy search...\n")
greedy <- greedy_backward(st, pairs, tcz.units.sim, n.screen = n.screen, n.confirm = n.confirm,
                          seed = seed, max.candidates = max.cand,
                          pairs.extreme = pairs.extreme, extreme.ref = rain.years$ref$extreme)

######## 3. optimised design ########
cat("Optimised design...\n")
ev.opt <- eval_full(greedy$remove, 2e6)

######## summaries ########
summary.tab <- rbind(summarise_design("do-nothing", ev.nothing, tcz.rows),
                     summarise_design("full", ev.full, none),
                     summarise_design("optimised", ev.opt, greedy$remove))
summary.tab$py.avoided <- summary.tab$py.total[1] - summary.tab$py.total
summary.tab$ce.per.million <- ifelse(summary.tab$cost > 0, summary.tab$py.avoided / (summary.tab$cost / 1e6), NA)
print(summary.tab)

design <- data.frame(tcz.units.sim[tcz.rows, c("unit_id", "control.type", "cost", "candidate")],
                     retained = !greedy$remove[tcz.rows])
design$names <- st_drop_geometry(tcz.units)$names[match(design$unit_id, tcz.units$unit_id)]

# check: each candidate tested once; capital cost removed adds up
stopifnot(!anyDuplicated(greedy$log$unit_id),
          isTRUE(all.equal(design_cost(tcz.units.sim, none) - design_cost(tcz.units.sim, greedy$remove),
                           sum(design$cost[!design$retained], na.rm = TRUE))))

write.csv(greedy$log, file.path(out.dir, "greedy-log.csv"), row.names = FALSE)
write.csv(greedy$restorations, file.path(out.dir, "restorations.csv"), row.names = FALSE)
write.csv(design, file.path(out.dir, "design.csv"), row.names = FALSE)
write.csv(summary.tab, file.path(out.dir, "summary.csv"), row.names = FALSE)
save(greedy, ev.nothing, ev.full, ev.opt, summary.tab, design, file = file.path(out.dir, "optimisation.RData"))
cat(sprintf("Done in %.1f min. Outputs in %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), out.dir))
