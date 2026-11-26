# TCZ design optimisation — plan

Goal: find the cheapest set of TCZ interventions with zero breaches of the zone over 50 years, and
report cost-effectiveness (toad point-years avoided per dollar), as described in
`ms/tcz-optimisation-methods.qmd` ("Optimising the design of the TCZ").

This is the overall plan. Each step gets a detailed plan (a section under "Detailed step plans")
when we start it, and this file is updated as steps change or finish.

## Where we are (2026-11-26)

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
- **Baseline:** `out/tcz-full.*` was regenerated in step 1 (2026-11-26) with tanks and fences in the run.

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
4. **Fence eradication branch and the extreme-wet-year kernel.**
   - When a fence is breached, the annual audit triggers eradication, which has its own cost and a
     5% chance of failing. The fence only emits toads (λ) the next year if eradication fails.
     At the moment the code treats fences like tanks, with a flat `control.fail.prob`.
   - In the top 20% of wet years, use a dispersal kernel built from the most rainy days on record.
5. **Drying data into the model.**
   - Load the survey data.
   - Calibrate the survey's `group_assessment` (Always Wet / Sometimes Dry / Don't Know) against DEA
     drying, using the long-run drying probability from `fit_drying`. Rewrite the methods text to
     match these categories.
   - Predict drying for points without a DEA match.
   - Build the waterpoint × year drying matrix (a 25-year rainfall reference, repeated to 50
     years) and pass it to `run_sims(drying.matrix=)`.
   - Carry the uncertainty through by drawing a calibration sample for each run.
6. **Final runs and write-up.**
   - Merge `main` into the branch. It will conflict with main's removal of the old Shiny app in
     `shiny/`.
   - Run the final optimisation.
   - Fill the XX placeholders and missing references in the methods, then render.

## Open questions

- **Costs:** are there maintenance costs, or a cost per eradication response? The cost sheet has
  only `Cost/Site`. Until there are, these are parameters that default to 0.
- **Uncosted works:** "Erin B yards" and "Refuse tip" have works planned but no cost, so they are
  fixed points. Five more uncosted records sit in units that do have a cost. Two cost rows ("Two
  tree well", "Sea Ranger base", $20k each) match no audit record and are excluded.
- **Gates:** `fence-lengths` already has `n.gates` for all 8 fenced areas. Confirm that's complete.

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
