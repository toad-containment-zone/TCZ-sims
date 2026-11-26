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
- **Stale:** `out/tcz-full.*` (31 Jul) predates the fence work.

## Steps

1. **Fix fence thinning; re-run the full-implementation baseline.**
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

_(added as each step starts)_
