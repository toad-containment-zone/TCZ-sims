# TCZ design optimisation — plan

Goal: find the cheapest set of TCZ interventions with zero breaches of the zone over 50 years, and
report cost-effectiveness (toad point-years avoided per dollar), as described in
`ms/tcz-optimisation-methods.qmd` ("Optimising the design of the TCZ").

This is the overall plan. Each step gets a detailed plan (a section under "Detailed step plans")
when we start it, and this file is updated as steps change or finish.

## Where we are (2026-11-27)

- **Built:** the control step in `spread.pilb()` (`src/pprocess_functions.R`):
  - natural points dry out via `surv.prob`, or a per-year `drying.matrix`;
  - tanks and fences persist via `control.fail.prob`;
  - fences get an effective radius `r.eff` and arrivals thinned by `p.breach`.

  Also `stop.on.target=FALSE`, the `occ.years` point-years count, `extinction.aware` outputs, and
  the synthetic check `src/tcz-sims/control-step-demo.R`.
- **Built:** `load-data.R` joins the post-PCR001 per-site cost and the fence lengths and gate counts
  (8 fences) to the audit records. `tcz-sims.R` gives each artificial TCZ point a control type
  (fence if it has a fence-length record, otherwise tank).
- **Built, but separate from the simulation:** `dea-drying-probs.R` → `rainfall-drying-model.R`
  (`fit_drying`, `predict_drying()`).
- **Data:** the drying survey's per-point answers are in the updated Living Waters layer (`group_assessment`), loaded in step 2 but not yet used.
- **Baseline:** `out/tcz-full.*` was last regenerated in step 6 (2026-11-27) with the final code.
- **All six steps are done.** The final design and its numbers are under "Step 6" below. The open
  questions that could still change them are listed under "Open questions".

## Steps

1. **Fix fence thinning; re-run the full-implementation baseline.** ✅ Done 2026-11-26 (see below).
   - `pprocess_functions.R:557,563` let arrivals through a fence with probability `1 - p.breach`.
     The methods say `c* ~ Binom(p_b, c)`. Fix the code, and extend `control-step-demo.R` so it
     would catch this.
   - Make the colonisation threshold the same in the code (`realised_n > 2`) and the methods
     ("failed when c ≥ 2").
   - Re-run `tcz-sims.R` with `stop.on.target=FALSE`.
2. **Tidy the infrastructure data and costs (`load-data.R`).** ✅ Done 2026-11-26 (see below).
   - Inside the TCZ: one row per artificial site (`X_record_id` in `tcz.sites`) and one per natural
     point (Living Waters, plus soaks).
   - Outside the TCZ: every point gets no control, so colonised artificial points outside the TCZ
     polygon stay colonised through the control step.
   - Sum costs per site (conversion and maintenance). Check that every TCZ site has a cost and a
     control type.
3. **Cost and objective accounting, then the greedy optimiser.** ✅ Pipeline built 2026-11-26 (see below). Results are provisional until steps 4–5.
   - Total 50-year cost: capital, recurrent, and eradication events.
   - Point-years avoided against a do-nothing run on the same point set; a breach flag; CE.
   - Backward greedy search:
     - try removing interventions, most expensive first, with 20 runs per candidate;
     - keep any intervention whose removal causes a breach;
     - confirm the final design with 100 runs.
   - Work out the run budget and whether it needs the cluster (`src/cluster/`).
   - Build this first with a static `surv.prob`, so it works end to end before the drying data
     goes in.
4. **Fence eradication branch and the extreme-wet-year kernel.** ✅ Done 2026-11-26 (see below).
   - When a fence is breached, the annual audit triggers eradication, which has its own cost and a
     5% chance of failing. The fence only emits toads (λ) the next year if eradication fails.
     At the moment the code treats fences like tanks, with a flat `control.fail.prob`.
   - In the top 20% of wet years, use a dispersal kernel built from the most rainy days on record.
5. **Drying data into the model.** ✅ Done 2026-11-26 (see below). The approach changed from the original plan: Traditional Owner assessment decides whether a point dries, DEA only how drying responds to rainfall.
   - Load the survey data.
   - Calibrate the survey's `group_assessment` (Always Wet / Sometimes Dry / Don't Know) against DEA
     drying, using the long-run drying probability from `fit_drying`. Rewrite the methods text to
     match these categories.
   - Predict drying for points without a DEA match.
   - Build the waterpoint × year drying matrix (a 25-year rainfall reference, repeated to 50
     years) and pass it to `run_sims(drying.matrix=)`.
   - Carry the uncertainty through by drawing a calibration sample for each run.
6. **Final runs and write-up.** ✅ Done 2026-11-27 (see below).
   - Merge `main` into the branch. It will conflict with main's removal of the old Shiny app in
     `shiny/`.
   - Run the final optimisation.
   - Fill the XX placeholders and missing references in the methods, then render.

## Open questions

- **Costs:** eradication responses are costed from step 6 at $37,469 per erad-unit per response.
- **To do: maintenance costs for fences and tanks** (project lead, 2026-11-27). Calculate a
  50-year maintenance cost for each unit and add it to that unit's capital cost, so it counts in
  the cost totals and in the screening and restoration ranking.
  - Fences would likely scale with fence length and gates (`distance`, `n.gates` are already on
    the units). Tanks would be a per-unit figure.
  - `design_cost(maintenance.per.year=)` is a flat per-unit hook that defaults to 0. Replace it
    with a per-unit `maintenance` column on `tcz-units`, built in `build_tcz_units()` from the
    cost workbook.
  - Do this together with the ranking item below, since both change the cost each unit is ranked
    on.
- **Screening and restoration rank units on capital cost only.** Eradication is now most of a
  fence's cost: the full design's $18.3M of eradication against $16.5M of capital. Ranking
  fences by capital + expected eradication cost could change which units are tested first and
  which are restored. Worth doing before the design is final.
- **"Shamrock Stn Pivots" fence** (on `fence-lengths` from 2026-11-27: 19.5 km, 3 gates, 7
  erad-units) has no `Site.ID`, audit record or cost, so it's left out (project lead,
  2026-11-27). Add it once it's in the audit. At 7 erad-units it would cost about $262k per
  response.
- **Uncosted works:** "Erin B yards" and "Refuse tip" have works planned but no cost, so they are
  fixed points. Five more uncosted records sit in units that do have a cost. Two cost rows ("Two
  tree well", "Sea Ranger base", $20k each) match no audit record and are excluded.
- **Gates:** `fence-lengths` already has `n.gates` for all 8 fenced areas. Confirm that's complete.
- **Follow-up (now Future work 1): is the Kimberley calibration still valid with SILO rainy days?** λ
  (`dat/Posteriors_2026.RData`) was fitted by ABC (`src/ABC/`) using rainy days from the
  Kimberley point data's `ndays_1` column, not SILO (its source isn't checked yet). The TCZ runs
  now take rainy days from SILO, which counts about 50% more than the rainfall raster at the TCZ
  points (42.6 against 28.3, r = 0.90). If `ndays_1` behaves like the raster, the TCZ model gives
  longer movement windows than the calibration assumed, and toads may spread faster than the
  calibrated λ implies. First check: compare `ndays_1` with SILO at the Kimberley points. Possible fixes: refit the ABC with SILO rainy days for the
  Kimberley points, or rescale SILO to the calibration source. (Flagged by the project lead,
  2026-11-26.)

## Future work (flagged by the project lead, 2026-11-27)

These come after the optimisation and are in this order, because each depends on the one before.

1. **Refit the Kimberley calibration on SILO rainy days.** SILO counts about 50% more rainy days
   than the BoM raster at the TCZ points (42.6 against 28.3, r = 0.90). λ
   (`dat/Posteriors_2026.RData`) was fitted by ABC (`src/ABC/spread_Kimberley_ABC*.R`) with the
   Kimberley points' `ndays_1`, so the TCZ runs, which now use SILO, move toads further each year
   than the calibration assumed.
   - Refit with SILO rainy days at the Kimberley points, using the same source as `rain-years.R`,
     so that the parameters carry over to the TCZ.
   - This replaces the "is the Kimberley calibration still valid?" open question above.
2. **Switch the colonisation threshold back to c ≥ 2 in the same refit.** Two arrivals are enough
   to found a population, so ≥ 2 makes more biological sense.
   - The code uses `realised_n > 2` (`pprocess_functions.R`, `spread.pilb()`). The methods
     describe it as c > 2 for now, to match the current results.
   - `spread.pilb()` is shared by the ABC fit, the TCZ runs and `forecast.R`. Change it only
     together with the refit, because λ was calibrated under the current threshold.
   - Then update the methods (simulation model and fence-failure text) and re-run the
     optimisation.
3. **Thread the refit into the annual forecast (`forecast.R`).**
   - Use the new posteriors, and SILO mean rainy days in place of the raster
     (`rain.id = "rainfall"`).
   - Add inter-annual rainfall variation as the optimisation runs do (`rain-years.R`): a random
     start year in the reference period per run, the extreme-wet-year kernel in the top 20% of
     years, and the natural-point drying matrix where it applies.
   - Then regenerate the public forecast site (`shiny/prep_app_data.R` → `shiny/web/`).

## Detailed step plans

### Step 1 — fence thinning fix and baseline re-run

Findings:
- `spread.pilb()` thins arrivals at a not-yet-colonised fence with `prob = 1 - p.breach` (in both
  the multinomial branch and the expected-value branch). `p.breach` is the chance that one toad
  gets through, so the fence currently stops only a fraction `p.breach` of arrivals. For a typical
  fence (2.3 km + 3 gates, g = 5.3) `p.breach` ≈ 0.07, so 93% of toads get through instead of 7%.
- `control-step-demo.R` test C uses `fence_units = 54`, which gives `p.breach` ≈ 0.5. At 0.5 the
  two directions give the same result, so the test can't catch the bug. Its printed explanation
  also assumes `1 - p.breach`.
- The `realised_n > 2` threshold has been in the code since the matrix rewrite (`fdec393`,
  2024-05). It applies to every point and was in place for the λ calibration, so it stays. The
  methods text changes to match.

Changes:
1. `src/pprocess_functions.R`: thin with `prob = p.breach` (and `expected_n * p.breach`). Correct
   the comments in `setup()`/`spread.pilb()` that describe the thinning.
2. `src/tcz-sims/control-step-demo.R` test C: use a realistic, asymmetric fence (`p.breach` ≈ 0.1).
   Add a `stopifnot()` that fails when the thinning goes the wrong way: the fenced colonisation
   rate must be well below the unfenced rate. Fix the printed explanation.
3. `ms/tcz-optimisation-methods.qmd`: "failed when $c_i \geq 2$" → failed when
   $c^*_i > 2$ (at least three colonists, the same threshold as every other point).
4. `src/tcz-sims/tcz-sims.R`: `run_sims(..., stop.on.target = FALSE)`, so the baseline runs the
   full 50 years and records point-years.
5. Run `control-step-demo.R`, then `tcz-sims.R` (100 runs × 50 years). Record in this file how
   often the full design is breached, and the mean point-years occupied.

Done when: the demo passes and would fail on the old code; `out/tcz-full.*` is regenerated; and
the baseline numbers are written here.

**Done (2026-11-26).**
- Thinning fixed. Demo test C now uses `p.breach` = 0.10 and switches off the fence
  persistence test (`control_fail_prob = 1`), so it measures thinning alone.
  - With the fix: unfenced 0.547, fenced 0.003 (predicted 0.003).
  - On the old code: fenced 0.475, and the check fails as it should.
- **Second bug found and fixed in `tcz-sims.R`.** `setup()` reads the natural/artificial column as
  1 = artificial, 0 = natural. `tcz-sims.R` passed `origin_des == "Natural"`, which is the
  opposite, together with `TCZ = TRUE`. As a result:
  - every artificial point inside the TCZ was removed, so tanks and fences were never simulated
    (the 31 Jul `tcz-full` run was really the old "remove all artificial water" scenario);
  - natural points were flagged as artificial, so they would never have dried in step 5.

  Now `tcz-sims.R` passes `art_nat = !natural` with `TCZ = FALSE`. `forecast.R` is unaffected,
  because the flag is only used by the control step. The convention is corrected in `CLAUDE.md`.
- **Baseline: full implementation** (100 runs × 50 years, `stop.on.target = FALSE`,
  natural points never dry):

  | | |
  |---|---|
  | Points in run | 727: 286 inside the TCZ (259 tanks, 8 fences, 19 natural), 82 targets |
  | Breaches | **0 / 100** |
  | Point-years occupied, all points | mean 16,298 (almost all from points already colonised north of the TCZ) |
  | Point-years occupied, inside the TCZ | mean 441 (range 54–758): natural 410, tank 25, fence 6 |
  | TCZ points ever colonised | mean 38 of 286 |

  Inside the TCZ, the point-years are almost all on natural points. Once colonised they never
  dry, because `surv.prob` = 1 until step 5. Expect this to fall once the drying data is in.
- Run time: about 3 min per 100 runs × 50 years on a laptop (about 2 s per run). For step 3:
  one greedy pass over about 267 interventions × 20 runs is about 3 h.

### Step 2 — infrastructure units and costs

Findings (2026-11-26):
- The cost workbook moved. The new one (`Activity Steps Toad project - 4 crew - 9Jul26.xlsx`,
  sheet "Analysis") has the same 202 sites. Each `Cost/Site` is 0.583 × the old value, because it
  now spreads only the direct infrastructure cost ($16.56M), not the whole project ($28.41M).
- `setup()`'s 100 m de-duplication was silently dropping 10 of the 277 audit records, 5 of them
  with costs (about $0.76M).
- Shelemar's fence was split across records. The fence-lengths row points at "Homestead" (no
  cost); the $1.83M cost sits on "Perimeter fence", which was modelled as a tank. The project lead
  confirmed that the homestead, P6 bore, P3 bore, Roper 1 and 2 bores and Typical valve
  installation are inside the fence. The points near the Nita Downs Pivot 2 and Shamrock Gardens
  fences are outside them.
- Decision: the 69 "Nothing" records are fixed controlled points (tank, 1% failure, $0, never
  removed).
- The updated Living Waters layer (`DATA_PATH/Toads/TCZ/living-waters/living_waters_points.geojson`)
  has 35 points (the 33 old ones plus 2 new), 20 inside the TCZ. It carries the survey's
  `group_assessment` and ground-truthing status.

Changes:
- New `src/tcz-sims/build-tcz-units.R` (`build_tcz_units()`) builds intervention units:
  - records in a fenced enclosure (table `fence.enclosures`) become one fenced unit, placed at
    the centroid of its water points;
  - other records within 100 m are merged;
  - each unit gets the summed cost, a control type, and `candidate = cost > 0`;
  - it checks that every record is in exactly one unit, costs add up, and no units are within 100 m.
- `load-data.R` reads the new workbook and Living Waters layer, and builds two point sets:
  `all.points` (every record; `forecast.R`) and `tcz.points` (units; `tcz-sims.R`). It writes
  `out/tcz-units.csv`.
- `tcz-sims.R` runs on `tcz.points`. It checks that every unit survives `setup()` and that nothing
  outside the TCZ has a control, and builds `tcz.units.sim` (unit ID, cost, candidate), aligned
  row-for-row with `spread.table`, for step 3.
- "Shelemar Perimeter Fence West - 3" is part of the same fence (confirmed by the project lead).

Decisions (project lead, 2026-11-26):
- **CE uses the direct cost** (`Cost/Site`, $16.56M in total), not the whole-project total.
- **Survey consent is recorded outside this layer.** The blank `model_use` field doesn't block
  modelling `group_assessment`.
- **The methods text changes to match the survey's categories** (Always Wet / Sometimes Dry /
  Don't Know), rather than the calibration being made to fit never / some years / most years /
  every year. Do this in step 5, with the calibration.

Results:
- **261 units** from 277 records (9 merged within 100 m, 7 folded into the Shelemar fence):
  8 fences and 253 tanks. **197 candidates**, $16,521,766 in total (the direct cost less the two
  orphan rows). 64 units are fixed.
- **Baseline: full implementation** (100 runs × 50 years, natural points never dry):

  | | Step 1 | Step 2 |
  |---|---|---|
  | Points in run (inside TCZ: tank / fence / natural) | 727 (259 / 8 / 19) | 723 (253 / 8 / 20) |
  | Breaches | 0 / 100 | **0 / 100** |
  | Point-years inside the TCZ, mean | 441 | 504 |
  | of which natural / tank / fence | 410 / 25 / 6 | 472 / 26 / 6 |
  | TCZ points ever colonised, mean | 38 | 40 |

  The rise in point-years is almost all on natural points: one more natural point inside the
  TCZ, and natural points still never dry.

### Step 3 — cost accounting and greedy optimiser

Plan:
- Split `tcz-sims.R` into `tcz-setup.R` (shared setup and unit checks) and the baseline run.
- New `optimise-functions.R`:
  - `apply_design()` removes interventions from a copy of `spread.table`. The water stays,
    uncontrolled: `control.type` none, `r.eff` NA, `p.breach` 0, `control.fail.prob` 0.
  - `eval_design()` runs n runs in parallel with common random numbers (run r always uses
    `set.seed(seed + r)`).
  - `design_cost()`: capital, plus maintenance and eradication hooks that default to 0.
  - `greedy_backward()`: one pass, most expensive first, 20 runs per test, with a run stopping
    at its first breach. Then a 100-run confirmation on a separate seed block, restoring removed
    units in reverse order of removal while it finds breaches (project lead's decision).
- New `tcz-optimise.R`: do-nothing counterfactual, full design, greedy search, then the optimised
  design (100 full-horizon runs each). Writes `out/optimisation/`. `SMOKE <- TRUE` gives a quick
  plumbing check.
- Methods text: costs are direct costs; the reverse-order restoration rule; common random numbers.

First full run (2026-11-26, 2 min 11 s on 8 cores; natural points never dry, no eradication
branch, no extreme-wet-year kernel):

| Design | Units kept | Cost | Breaches /100 (95% upper) | Point-years (all points) | Point-years avoided | CE (avoided per $M) |
|---|---|---|---|---|---|---|
| Do nothing | 0 | $0 | 100 | 33,327 | — | — |
| Full | 261 | $16.52M | 0 (3.6%) | 16,343 | 16,984 | 1,028 |
| Optimised | 249 | $11.72M | 0 (3.6%) | 16,765 | 16,561 | 1,413 |

- The screening removed 54 of the 197 candidates. Every removal had passed 0/20 on its own, but
  together they gave 17/100 breaches at confirmation. Restoring in reverse order then put 42
  units back ($2.89M) before reaching 0/100. Most of those restorations (cheap units, removed
  late) barely changed the breach rate (17 → 16/100). The rate only fell once the restoring
  reached units U014–U023. **Reverse-order restoration is inefficient here.** One alternative is
  to restore by each unit's effect on breaches.
- The final design removes 12 units: the Shelemar, Shamrock Gardens, Nita Downs Homestead and
  Pivot 1/2 fences (5 of the 8 fences), and 7 tanks.
- **Problem: in the model, removing a fence makes its site *safer*.** A fenced unit is a big
  target: its effective radius `r.eff` comes from the fenced area (Shelemar 3.9 km), thinned by
  `p.breach` (Shelemar 0.31). `apply_design()` turns a removed fence into one ordinary point with
  r = 100 m, so its catchment area falls about 1,500-fold. Removing a fence therefore cuts
  arrivals rather than raising them. This is why the most expensive fences go first. Fence
  removal has to be modelled differently before the fence results mean anything.
- Many single removals breach on their own (143 of 197 candidates kept). Natural points never
  dry yet, so they act as permanent stepping stones. Expect this to change a lot in step 5.

Decisions (project lead, 2026-11-26):
- **A removed fence keeps its effective radius**, everywhere (Shelemar included). The wet area is
  still as big a target; only the breach thinning and the control test go. `apply_design()` no
  longer clears `r.eff`.
- **Restoration rule: decide after step 5.** Keep reverse order until then.

Second run, with removed fences keeping their radius (2026-11-26, 2 min 12 s):

| Design | Units kept | Cost | Breaches /100 (95% upper) | Point-years (all points) | Point-years avoided | CE (avoided per $M) |
|---|---|---|---|---|---|---|
| Do nothing | 0 | $0 | 100 | 33,775 | — | — |
| Full | 261 | $16.52M | 0 (3.6%) | 16,343 | 17,432 | 1,055 |
| Optimised | 252 | $12.27M | **2** (7.0%) | 16,973 | 16,802 | 1,370 |

- The screening removed 47 units (4 of 8 fences, 43 tanks). The confirmation found 27/100
  breaches, and restoring in reverse order put back 38 units ($2.65M) before it reached 0/100.
- The final design removes 9 units, saving $4.26M: the Shelemar, Shamrock Gardens and Nita Downs
  Homestead fences, plus 6 tanks (U003–U006, U008, U009). The Shelemar and Shamrock Gardens
  fences can still go even with their radius kept, which suggests these sites are off the main
  invasion path.
- The optimised design passed its confirmation (0/100), but the independent evaluation (a
  separate seed block) found **2/100 breaches**. "Zero in 100" only bounds the breach risk at
  about 3.6%, and this design's true risk looks to be around 2%. Worth weighing when the
  restoration rule is decided.

### Step 4 — fence eradication and extreme-wet-year kernel

Plan and decisions (project lead, 2026-11-26):
- Count an eradication response for every occupied fence at each annual audit (`erad.events` in
  `spread.pilb()`). The fence biology (5% eradication failure, then emitting λ) was already in the
  control step. Eradication cost is a parameter: `eradication.cost <- 0` in `tcz-optimise.R`
  until there's a figure.
- Extreme wet years are the top 20% of the 25-year reference (2002–2026) by Dec–May SILO rainfall
  averaged over the cells under the simulated points. In those years toads disperse through
  `pairs.extreme`, built from each point's cell's maximum rainy days in any reference year.
- Each run starts at a random reference year and wraps round it. Run r gets the same start in
  every design. Step 5's drying matrix will use the same year sequence.

Changes:
- New `src/tcz-sims/rain-years.R` (SILO fetch for 306 cells, cached in
  `dat/silo-tcz-rain-cache.RData`).
- `setup(rain.extreme.id=)` builds `pairs.extreme`, and `erad.events` is always present.
- `spread.pilb(pairs.extreme=, extreme.years=)`.
- `run_sims(extreme.ref=)` and `eval_design()` draw a random start per run.
- Demo checks D (eradication counting) and E (`extreme.years = FALSE` reproduces a standard run
  exactly; the extreme kernel colonises targets 1.5 km away at 0.177 vs 0.116).
- Methods text: reference years, the random start, the wet-year definition, the 250-day kernel
  cap, and eradication responses.

Results (2026-11-26; natural points still never dry):
- **Extreme years:** 2004, 2011, 2017, 2018, 2023. Per point, the wettest-year rainy days range
  from 45 to 103 (mean 78), against a SILO mean of 43. Points with more than about 91 rainy days
  go past the 250-day kernel table and are capped: 147 of 848 input points.
- The extreme kernel is flatter: fewer arrivals within about 1 km, more beyond (at 1.5 km,
  expected arrivals rise from 1.19 to 1.46).
- **Baseline: full implementation:** 0/100 breaches. Mean point-years inside the TCZ 610 (step 2:
  504): natural 572, tank 31, fence 7.
- **Eradication responses:** the full design triggers a mean of **145 per 50-year run** (range
  70–192), roughly 3 a year across the 8 fences. The big fenced areas are colonised often. Once a
  cost per response exists this will weigh heavily on fences.
- **Optimiser** (1.3 min):

  | Design | Units kept | Cost | Breaches /100 (95% upper) | Point-years avoided | CE (avoided per $M) | Eradications / run |
  |---|---|---|---|---|---|---|
  | Do nothing | 0 | $0 | 100 | — | — | 0 |
  | Full | 261 | $16.52M | 0 (3.6%) | 17,649 | 1,068 | 145 |
  | Optimised | 256 | $13.15M | 0 (3.6%) | 17,245 | 1,312 | 150 |

  - The screening removed 25 units. The confirmation found 30/100 breaches, and reverse-order
    restoration put back 20 ($1.91M).
  - The final design removes 5 units, saving $3.37M: the Shelemar and Shamrock Gardens fences,
    plus tanks U003–U005.
  - With the extreme kernel, fewer removals survive than in step 3 (9 units removed then).

Decision (project lead, 2026-11-26): **take all TCZ rainfall from SILO.** The standard kernel now
uses each point's mean SILO rainy days over the 25 reference years (`rain.id = "rain.mean"` in
`tcz-setup.R`), not the rainfall raster. The extreme kernel uses the SILO maximum, as before.
The raster averaged 28.3 rainy days at these points against SILO's 42.6 (r = 0.90), so mixing
the two sources had exaggerated the jump to the extreme kernel. `forecast.R` still uses the raster.

Re-run with SILO for both kernels:
- Standard movement days: mean 140 (range 81–190), up from 101 with the raster. Dispersal is
  faster in ordinary years too.
- **Baseline: full implementation:** 0/100 breaches. Mean point-years inside the TCZ 647: natural
  604, tank 35, fence 7.

  | Design | Units kept | Cost | Breaches /100 (95% upper) | Point-years avoided | CE (avoided per $M) | Eradications / run |
  |---|---|---|---|---|---|---|
  | Do nothing | 0 | $0 | 100 | — | — | 0 |
  | Full | 261 | $16.52M | 0 (3.6%) | 17,862 | 1,081 | 152 |
  | Optimised | 256 | $13.15M | 0 (3.6%) | 17,443 | 1,327 | 166 |

- The screening removed 27 units. The confirmation found 20/100 breaches, and restoring put back
  22 ($2.19M).
- The final design is the same as with the raster: the Shelemar and Shamrock Gardens fences, plus
  tanks U003–U005, saving $3.37M.

Decision (project lead, 2026-11-26): **keep `days.per.rain = 4`.** It's what both Kimberley ABC
calibrations used, and it matches the original Tingley et al. formula (`r + 3r(1 − overlap)`:
movement on the rainy day plus the three following days, a 4-day window). The methods text said
"three days of movement" and now describes the 4-day window instead. No code change.

### Step 5 — natural-point drying

Findings (2026-11-26):
- 10 of the 35 Living Waters points are within 100 m of a DEA waterbody (463 in the region): 7
  Always Wet, 1 Sometimes Dry, 1 Don't Know, 1 blank.
- **DEA contradicts the TO assessment.** The waterbodies at the "Always Wet" points reach 0% open
  water in 75–100% of dry seasons (1987–2026), and the "Sometimes Dry" one in every season. DEA
  sees open surface water across the polygon, not the persistent water toads use.
- `dea-drying-probs.R` didn't match the methods: it matched within 1,000 m, and the comment on its
  months didn't match the code. All 141 LaGrange natural points outside the TCZ are flagged
  "Perennial".

Decisions (project lead, 2026-11-26):
- **TO-primary.** "Always Wet" never dries. "Sometimes Dry", "Don't Know" and blank Living Waters
  points dry with probability logit p = logit(**0.25**) + b₁·x₁ + b₂·x₂, where x₁ and x₂ are the
  current and previous wet-season rainfall anomalies from the point's own reference mean. b₁ and
  b₂ come from the DEA glmer, so DEA sets only the rainfall response.
- Natural points with no TO assessment (LaGrange and aerial-imagery points, all outside the TCZ)
  never dry. The Living Waters rule takes precedence, so blank Living Waters records outside the
  TCZ are treated as Sometimes Dry.

Changes:
- New `src/tcz-sims/drying-model.R`:
  - `dea_drying()`: 100 m match, Jun–Dec, cached in `dat/dea-drying-cache.RData`;
  - `fit_drying_response()`: the glmer, with SILO rainfall from `rain-years.R`;
  - `drying_spec()`, `draw_drying_matrix()` and `draw_drying_beta()`.
- `run_sims(drying.spec=)` and `eval_design(drying.spec=)` draw (b₁, b₂) from their sampling
  distribution for each run, and use that run's random start. Drying and the extreme kernel follow
  the same year sequence.
- `tcz-setup.R` builds `drying.spec` aligned with `spread.table`, and stops if the glmer is
  singular or a slope isn't negative.
- `dea-drying-probs.R` and `rainfall-drying-model.R` are marked superseded.
- Demo check F.
- Methods: "Auditing natural waterpoints" is rewritten (**a draft for the project lead to review**).

Results:
- **glmer:** 10 sites, 258 site-years (2001–2026), not singular.
  - b₁ (current wet season) = −1.30 per SD (213 mm), SE 0.29;
  - b₂ (previous wet season) = −0.78 per SD (253 mm), SE 0.26;
  - site SD 1.48.
- **Drying classes:**
  - inside the TCZ, 10 natural points never dry (Always Wet) and 10 sometimes dry (7 Sometimes
    Dry, 3 Don't Know);
  - outside the TCZ, 8 sometimes dry (Living Waters) and 155 never dry.
- **Drying probability by reference year** (mean over "sometimes" points): it ranges from 0.01
  (2018) to 0.72 (2016), and the mean over years is 0.31. The extreme years dry least (0.01–0.18).
- **Baseline: full implementation:** 0/100 breaches. Mean point-years inside the TCZ 522 (step 4:
  647): natural 483 (was 604), tank 32, fence 7.
- **Optimiser** (1.6 min):

  | Design | Units kept | Cost | Breaches /100 (95% upper) | Point-years avoided | CE (avoided per $M) | Eradications / run |
  |---|---|---|---|---|---|---|
  | Do nothing | 0 | $0 | 100 | — | — | 0 |
  | Full | 261 | $16.52M | 0 (3.6%) | 17,756 | 1,075 | 142 |
  | Optimised | 255 | $12.91M | **1** (5.4%) | 17,358 | 1,345 | 153 |

  - The screening removed 38 units (5 fences, 33 tanks). The confirmation found 18/100 breaches,
    and reverse-order restoration put back 32 ($2.75M).
  - The final design removes 6 units, saving $3.61M: the Shelemar and Shamrock Gardens fences,
    plus tanks U003–U006.
  - The independent 100-run evaluation found 1 breach.
- Drying cuts natural point-years by about 20%, but the 10 "Always Wet" points inside the TCZ
  stay permanent stepping stones. The design changes little from step 4 (one more tank, U006,
  removed).

### Step 6 — final runs and write-up

Decisions (project lead, 2026-11-27):
- **Restoration rule: best value.** While the confirmation finds breaches:
  - replay the breaching runs (same run numbers, so the same random streams) once for each unit
    still removed, with that unit restored;
  - restore the unit that fixes the most breaches per capital dollar (ties: cheaper, then
    `unit_id`);
  - if no single unit fixes any breach, restore the most recently removed one ("fallback");
  - re-confirm.

  `greedy_backward(restore = "reverse")` keeps the old rule.
- **Run sizes:** 100 runs per screening test, a 1,000-run confirmation, and 1,000-run evaluations
  of each design. "0 in 1,000" bounds the breach risk at about 0.3%.
- **Eradication costs** (added by the project lead to "Draft QRA data - Fencing unit rates.xlsx"):
  - each fence has `erad.units` (two-person crews with a quad bike needed to cover it in a
    night): Shelemar 8, Shamrock Gardens 2, Shamrock Stn Homestead 1, the others 0.5;
  - an erad-unit costs $1,800 a day (personnel $1,400, quad bike $150, vehicle $250);
  - a response takes 20.8 days per erad-unit (the workbook's "Eradication-days-calculations":
    N₀ = 700, detection probability 0.27), so **$37,469 per erad-unit per response** (updated by
    the project lead, 2026-11-27, replacing a first run that assumed one night per response).

Changes:
- Merged `main` (`59da87c`): the static forecast site replaces the Shiny app. The conflicts were
  in `CLAUDE.md` and `.gitignore`.
- `optimise-functions.R`:
  - `eval_design(runs=)` replays chosen run numbers;
  - `eval_design(erad.cost.row=)` returns each run's eradication cost, fence by fence;
  - `greedy_backward(restore=)`, with `restorations.csv` gaining `rule`, `breaches.fixed` and
    `n.tested`;
  - `design_cost(eradication=)` takes the run-mean eradication cost.
- `load-data.R` reads the eradication sheet. `build-tcz-units.R` carries `erad.units` onto the
  fence units, and `tcz-setup.R` builds `erad.cost.row`.
- `tcz-optimise.R`: runs at 100 / 1,000 / 1,000, with the summary split into capital, eradication
  and total cost.
- Methods:
  - best-value restoration, the run sizes, the independent final evaluation and eradication
    costing;
  - filled the telemetry counts (67 toads with at least 5 turning angles, so 670,000 d*ₙ per n)
    and the Tingley 2013, DEA (Krause et al. 2021) and SILO (Jeffrey et al. 2001) citations;
  - corrected the colonisation threshold to c > 2, as in the code. Line 16 still said ≥ 2.
  - The PDF renders (`ms/output/ms/`, gitignored).

Checks:
- Demo checks A–F pass.
- With `restore = "reverse"` at step 5's sizes and seed, the code gives exactly step 5's design
  (U001–U006 removed, 32 restorations).
- Forcing step 5's 38 screened removals and applying best-value restoration put back 16 units
  ($1.87M), against reverse order's 32 ($2.75M).

Results (2026-11-27, 6.5 min on 8 cores, at $37,469 per erad-unit per response):

| Design | Units kept | Capital | Eradication | Total cost | Breaches /1000 (95% upper) | Point-years avoided | CE (avoided per $M) | Eradications / run |
|---|---|---|---|---|---|---|---|---|
| Do nothing | 0 | $0 | $0 | $0 | 1000 | — | — | 0 |
| Full | 261 | $16.52M | $18.32M | $34.85M | 0 (0.37%) | 17,748 | 509 | 141 |
| Optimised | 244 | $12.41M | $2.33M | $14.74M | **1** (0.56%) | 17,408 | 1,181 | 103 |

- The screening removed 23 units (2 fences, 21 tanks). The confirmation found 38/1000 breaches.
  Best-value restoration put back 6 units ($0.84M): U042, U159, U010, U009, U060 and U006. It
  reached 0/1000 after six rounds.
- **The final design removes 17 units, saving $4.12M capital and $20.11M in total:**
  - the Shelemar (U001) and Shamrock Gardens (U002) fences;
  - 15 tanks: U003, U004, U008, U022, U026, U039, U069, U071, U120, U121, U123, U136, U158, U161,
    U190.
- The independent 1,000-run evaluation found **1 breach** (p ≈ 0.001, 95% upper 0.56%). That is
  consistent with the confirmation's 0/1000, and the residual risk is far below step 5's
  (1/100).
- **Eradication is now the largest cost of the full design.** Shelemar (8 erad-units, about
  $300k per response) and Shamrock Gardens (2) are colonised often. Removing them cuts the
  eradication cost from $18.3M to $2.3M. Most of the gain in CE (509 → 1,181) comes from that,
  not from the capital saved.
- The design is the same as in the first step 6 run (one night per response). Screening and
  restoration rank units on capital cost, and eradication cost only enters the reported totals.
  See "Open questions".

Still open for the write-up (methods):
- the Geoscience Australia dataset names ("XX and YY") and the "step out details";
- the soak count (XX);
- full references for Dunlop et al., Brown et al. and Phillips, Dempster et al.;

Fixed since: Southwell et al. is now 2017 throughout, and the telemetry dates read "wet seasons
(January–March) of 2005–2010", to match `gregsData.txt` (project lead, 2026-11-27). The remaining
gaps the project lead will fill.
