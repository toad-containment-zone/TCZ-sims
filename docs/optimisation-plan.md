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
- **Data:** the Traditional Owner drying survey data has arrived (2026-11) but isn't loaded yet.
- **Baseline:** `out/tcz-full.*` was regenerated in step 1 (2026-11-26) with tanks and fences in the run.

## Steps

1. **Fix fence thinning; re-run the full-implementation baseline.** ✅ Done 2026-11-26 (see below).
   - `pprocess_functions.R:557,563` let arrivals through a fence with probability `1 - p.breach`.
     The methods say `c* ~ Binom(p_b, c)`. Fix the code, and extend `control-step-demo.R` so it
     would catch this.
   - Make the colonisation threshold the same in the code (`realised_n > 2`) and the methods
     ("failed when c ≥ 2").
   - Re-run `tcz-sims.R` with `stop.on.target=FALSE`.
2. **Tidy the infrastructure data and costs (`load-data.R`).**
   - Inside the TCZ: one row per artificial site (`X_record_id` in `tcz.sites`) and one per natural
     point (Living Waters, plus soaks).
   - Outside the TCZ: every point gets no control, so colonised artificial points outside the TCZ
     polygon stay colonised through the control step.
   - Sum costs per site (conversion and maintenance). Check that every TCZ site has a cost and a
     control type.
3. **Cost and objective accounting, then the greedy optimiser.**
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
   - Calibrate the survey categories (never / some / most / every year) against DEA drying, using
     an ordered logistic model on the long-run drying probability from `fit_drying`.
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

- **Survey data:** where it lives under `DATA_PATH`, how it's laid out, and whether it adds soaks
  beyond the 33 Living Waters points.
- **Costs:** are there maintenance costs, or a cost per eradication response? The PCR001 sheet has
  only `Cost/Site`. Until there are, these are parameters that default to 0.
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
