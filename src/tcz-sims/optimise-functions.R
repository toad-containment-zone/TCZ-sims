# optimise-functions.R
#
# Design evaluation and the backward greedy optimiser for the TCZ (ms/tcz-optimisation-methods.qmd,
# "Cost optimisation"). Works on the spread.table/pairs left by tcz-setup.R, plus tcz.units.sim
# (intervention-unit metadata aligned row-for-row with spread.table; NA unit_id = not a TCZ unit).
#
# A "design" is a logical vector over spread.table rows: TRUE = intervention removed. Removing an
# intervention leaves the water in place, uncontrolled: once colonised the point persists and
# emits propagules like any other artificial waterpoint. A removed fence keeps its effective
# radius (the wet area is still there); it loses only the breach thinning and the control test.

library(parallel)

# copy of spread.table with the interventions at rows `remove` taken out
apply_design <- function(st, remove) {
  stopifnot(is.logical(remove), length(remove) == nrow(st))
  if (any(remove)) {
    st[remove, "control.type"] <- CONTROL_NONE
    # r.eff is kept: an unfenced wet area is still as big a target as the fenced one was (project
    # lead's decision, 2026-11-26), it just loses the fence's thinning and the control test
    st[remove, "p.breach"] <- 0
    st[remove, "control.fail.prob"] <- 0
  }
  st
}

# n stochastic 50-year runs of one design, split across cores.
# Common random numbers: run r always uses set.seed(seed + r), so every design is tested against the
# same lambda draws and random streams, and differences between designs come from the design.
# Returns one row per run. With stop.on.target = TRUE a run stops at its first breach (fast; only
# `breach` is then meaningful); use FALSE for point-years over the full horizon.
# `runs` picks which run numbers to do (default 1:n), e.g. to replay only the runs that breached.
eval_design <- function(st, pairs, n, gens = 50, stop.on.target = TRUE, seed = 1,
                        tcz.rows = rep(FALSE, nrow(st)), drying.matrix = NULL,
                        pairs.extreme = NULL, extreme.ref = NULL, drying.spec = NULL,
                        cores = max(1, detectCores() - 2), runs = seq_len(n),
                        erad.cost.row = rep(0, nrow(st))) {
  one_run <- function(r) {
    set.seed(seed + r)
    lambda <- 10^rnorm(1, mean = sample.lambda, sd = sample.lambda.sd) # as in run_sims()
    start <- NA; extreme.years <- NULL; drying.run <- drying.matrix
    n.ref <- if (!is.null(extreme.ref)) length(extreme.ref) else if (!is.null(drying.spec)) nrow(drying.spec$x1) else NA
    if (!is.na(n.ref)) { # random start in the reference years, as in run_sims()
      start <- sample.int(n.ref, 1)
      idx <- year_index(n.ref, gens, start)
      if (!is.null(extreme.ref)) extreme.years <- extreme.ref[idx]
      if (!is.null(drying.spec)) drying.run <- draw_drying_matrix(drying.spec, idx, draw_drying_beta(drying.spec))
    }
    invisible(capture.output( # silence spread.pilb()'s progress bar
      out <- spread.pilb(pop = st, gens = gens, pairs = pairs, delta = lambda, r = 10^2,
                         plot = FALSE, rollup = FALSE, control.step = TRUE,
                         stop.on.target = stop.on.target, drying.matrix = drying.run,
                         pairs.extreme = pairs.extreme, extreme.years = extreme.years)
    ))
    pm <- out$popmatrix
    data.frame(run = r, lambda = lambda, start = start, breach = !is.na(out$gen), gen = out$gen,
               py.total = sum(pm[, "occ.years"]), py.tcz = sum(pm[tcz.rows, "occ.years"]),
               n.erad = sum(pm[, "erad.events"]),
               erad.cost = sum(pm[, "erad.events"] * erad.cost.row))
  }
  res <- mclapply(runs, one_run, mc.cores = cores, mc.preschedule = FALSE)
  failed <- vapply(res, inherits, logical(1), "try-error")
  if (any(failed)) stop("eval_design: ", sum(failed), " run(s) failed: ", res[[which(failed)[1]]])
  do.call(rbind, res)
}

# total design cost over the horizon: capital for retained units, plus recurrent costs (0 until there
# are figures) and the run-mean eradication cost (eval_design()'s erad.cost, fence by fence)
design_cost <- function(units.sim, remove, years = 50, maintenance.per.year = 0, eradication = 0) {
  kept <- !is.na(units.sim$unit_id) & !remove
  capital <- sum(units.sim$cost[kept], na.rm = TRUE)
  capital + maintenance.per.year * sum(kept) * years + eradication
}

# backward greedy search (single pass, most expensive first). A candidate whose removal produces any
# breach in n.screen runs is kept and never retested. The final design is confirmed with n.confirm
# runs on a separate seed block; while that finds breaches, one removed unit is restored and the
# design re-confirmed:
#   restore = "best"    (default; project lead's decision, 2026-11-27) replay the runs that breached
#                       with each still-removed unit restored in turn, and restore the unit that fixes
#                       the most breaches per dollar (ties: cheaper, then unit_id). If no single unit
#                       fixes any, restore the most recently removed one ("fallback").
#   restore = "reverse" restore the most recently removed unit (the steps 3-5 rule).
greedy_backward <- function(st, pairs, units.sim, n.screen = 20, n.confirm = 100, gens = 50,
                            seed = 20261126, max.candidates = Inf,
                            cores = max(1, detectCores() - 2), drying.matrix = NULL,
                            pairs.extreme = NULL, extreme.ref = NULL, drying.spec = NULL,
                            restore = c("best", "reverse")) {
  restore <- match.arg(restore)
  ev_args <- list(pairs = pairs, gens = gens, stop.on.target = TRUE, cores = cores,
                  drying.matrix = drying.matrix, pairs.extreme = pairs.extreme,
                  extreme.ref = extreme.ref, drying.spec = drying.spec)
  ev <- function(remove, n, seed, runs = seq_len(n))
    do.call(eval_design, c(list(st = apply_design(st, remove), n = n, seed = seed, runs = runs), ev_args))
  money <- function(x) format(round(x), big.mark = ",")

  cand <- which(units.sim$candidate %in% TRUE)
  ord <- cand[order(-units.sim$cost[cand], units.sim$unit_id[cand])]
  ord <- ord[seq_len(min(length(ord), max.candidates))]
  remove <- rep(FALSE, nrow(st))
  log <- vector("list", length(ord))
  t0 <- Sys.time()

  for (k in seq_along(ord)) {
    row <- ord[k]
    trial <- remove; trial[row] <- TRUE
    n.breach <- sum(ev(trial, n.screen, seed)$breach)
    if (n.breach == 0) remove <- trial
    log[[k]] <- data.frame(step = k, unit_id = units.sim$unit_id[row], row = row,
                           control.type = units.sim$control.type[row], cost = units.sim$cost[row],
                           n.runs = n.screen, n.breach = n.breach,
                           decision = if (n.breach == 0) "removed" else "kept",
                           cost.removed.cum = sum(units.sim$cost[remove], na.rm = TRUE))
    cat(sprintf("[%3d/%d] %s %-5s $%10s  breaches %2d/%d -> %-7s (removed so far: %d, $%s; %.1f min)\n",
                k, length(ord), units.sim$unit_id[row], units.sim$control.type[row],
                money(units.sim$cost[row]), n.breach, n.screen,
                log[[k]]$decision, sum(remove), money(log[[k]]$cost.removed.cum),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  log <- do.call(rbind, log)

  # confirmation, restoring one unit at a time until it finds no breaches
  removal.order <- log$row[log$decision == "removed"]
  restorations <- list()
  conf.seed <- seed + 1e6
  repeat {
    conf <- ev(remove, n.confirm, conf.seed)
    n.breach <- sum(conf$breach)
    cat(sprintf("Confirmation: %d/%d breaches with %d units removed\n", n.breach, n.confirm, sum(remove)))
    if (n.breach == 0) break
    still.removed <- removal.order[remove[removal.order]]
    if (length(still.removed) == 0) {
      warning("confirmation still finds breaches with every unit restored")
      break
    }
    rule <- "reverse"; fixed <- NA_integer_
    back <- tail(still.removed, 1)
    if (restore == "best") {
      # replay only the breaching runs (same run numbers = same random streams) with each unit back
      bad <- conf$run[conf$breach]
      fixes <- vapply(still.removed, function(u) {
        trial <- remove; trial[u] <- FALSE
        length(bad) - sum(ev(trial, n.confirm, conf.seed, runs = bad)$breach)
      }, integer(1))
      if (any(fixes > 0)) {
        score <- fixes / units.sim$cost[still.removed]
        pick <- order(-score, units.sim$cost[still.removed], units.sim$unit_id[still.removed])[1]
        back <- still.removed[pick]; rule <- "best"
      } else {
        rule <- "fallback"
      }
      fixed <- fixes[match(back, still.removed)]
    }
    remove[back] <- FALSE
    restorations[[length(restorations) + 1]] <- data.frame(
      unit_id = units.sim$unit_id[back], row = back, cost = units.sim$cost[back],
      n.breach.before = n.breach, n.runs = n.confirm, rule = rule,
      breaches.fixed = fixed, n.tested = if (restore == "best") length(still.removed) else NA_integer_)
    cat(sprintf("  restored %s ($%s; %s%s)\n", units.sim$unit_id[back], money(units.sim$cost[back]), rule,
                if (is.na(fixed)) "" else sprintf(", fixes %d of %d breaching runs", fixed, n.breach)))
  }
  restorations <- if (length(restorations)) do.call(rbind, restorations) else
    data.frame(unit_id = character(), row = integer(), cost = numeric(), n.breach.before = integer(),
               n.runs = integer(), rule = character(), breaches.fixed = integer(), n.tested = integer())

  list(log = log, remove = remove, restorations = restorations, confirmation = conf,
       restore = restore, minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
