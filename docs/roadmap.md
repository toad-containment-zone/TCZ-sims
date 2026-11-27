# TCZ work — roadmap to finalise

Goal: finish the TCZ modelling, so that its results can be reported and acted on. That means:
- a Kimberley calibration that transfers to the TCZ;
- an annual forecast built on it;
- an optimised TCZ design built on checked inputs;
- a finished methods and results paper, a program summary report and an updated public forecast;
- everything merged into `main`.

This is the overall roadmap. Each phase gets a detailed plan (under "Detailed phase plans") when
we start it, and this file is updated as phases change or finish. It replaces
`docs/optimisation-plan.md`, which is archived as
[`docs/archive/optimisation-plan-2026-11.md`](archive/optimisation-plan-2026-11.md): the record of
optimisation steps 1–6, with every decision, result table and check.

## Where we are (2026-11-27)

- **Optimisation steps 1–6 are done** on branch `tcz-optimisation` (see the archive).
  - The model has the control step, intervention units, eradication responses, the
    extreme-wet-year kernel, natural-point drying and best-value restoration.
  - **Final design (provisional):** remove 17 of 261 units (the Shelemar and Shamrock Gardens
    fences and 15 tanks). That costs $14.74M ($12.41M capital and $2.33M eradication), against
    $34.85M for the full design. It found 0/1000 breaches at confirmation and 1/1000 in the
    independent evaluation, and CE is 1,181 against 509.
  - The design is provisional until Phases 1, 3 and 4: the calibration, the inputs and the cost
    model can all still change it.
- **The annual forecast (`forecast.R`) is unchanged.** It still uses the rainfall raster, the old
  posteriors and constant rainfall each year.

## Standing decisions (carried forward)

These are settled unless the project lead reopens them. Dates are when each was decided; the
archive has the reasons.

- Costs are direct infrastructure costs (`Cost/Site`), with no overheads (2026-11-26).
- A removed fence keeps its effective radius `r.eff` and loses only its breach thinning and
  control test (2026-11-26).
- `days.per.rain = 4`: movement on the rainy day plus the 3 days after it (2026-11-26).
- Movement days come from `rain_to_days()` with `season.days = 212`. That is the Nov–May movable
  season, since toads don't move from June to October. Only Nov–May rainy days count. The
  mapping is a rough, saturating guess: about 4 days per rainy day early in the season, levelling
  off at the season length. It isn't worth refining without data on movement against water
  balance (2026-11-27).
- All TCZ rainfall comes from SILO. The reference years are 2002–2026, extreme years are the top
  20% by wet-season rainfall, and each run starts at a random reference year (2026-11-26).
- Drying: the TO assessment decides whether a point dries, and DEA only how drying responds to
  rainfall. "Don't Know" and blank count as Sometimes Dry, p = 0.25 in an average year, and
  natural points with no TO assessment never dry (2026-11-26).
- Restoration is best value, with 100 runs per screen, a 1,000-run confirmation and 1,000-run
  evaluations (2026-11-27).
- Eradication costs $37,469 per erad-unit per response (2026-11-27).
- The Kimberley refit uses the same model structure as the TCZ runs (2026-11-27).
- **One calibration and model structure for every live application (2026-11-27).**
  - All live applications use the same λ posterior, kernel table, movement mapping, colonisation
    threshold and rain source.
  - The calibration is set in one place, `POSTERIORS_FILE` in `pprocess_functions.R`, and loaded
    with `load_posteriors()`.
  - When the calibration or the structure changes, update every live application in the same
    piece of work, and flag each change as it's made.
  - The live applications are listed under "Applications" below. Everything else is in
    `src/archive/` and is not maintained.
- The Kimberley fit has no natural-point drying. There's no TO assessment there, so under the
  current rule its natural points never dry. That keeps the control step off and rollup on
  (2026-11-27).
- The Kimberley SILO fetch covers 4,325 cells for 1999–2026. The fetch and the ABC both run
  locally (2026-11-27).
- **r = 300 m** (`R_DETECT` in `pprocess_functions.R`), changed from 100 m (2026-11-27).
  - 300 m is closer to biology.
  - The calibration only fixes λ·r². A fence's target (`r.eff`) is set by its area, not by r, so
    r sets how hard fences are pressed compared with open water.
  - A fence is never a smaller target than open water: the model uses `max(r, r.eff)`. This
    matters for 3 small fences, with `r.eff` of 179–219 m.
  - Points are thinned so that no detection discs overlap (`thin_overlapping()`, which replaced
    the 100 m dedup). Overlapping discs double-count arrivals and break the multinomial: at
    r = 300 with the 100 m dedup, the median disc was 44% double-counted, and 53 Kimberley
    sources had more than 100% of their propagules "landing" (2026-11-27).
    - Survivors are chosen by priority: units, then colonised points, then observed points,
      then targets, then the rest.
    - Units and natural waterpoints inside the TCZ are always kept, even where they overlap
      (44 pairs at the TCZ). Without this, a spring beside a tank would be hidden and the place
      would be modelled as controlled (`setup(keep.id=)`).
    - In the Kimberley this keeps 15,698 of 31,521 points; at the TCZ it keeps 619 of 848,
      including all 22 Living Waters points.

## Applications

The applications that must stay on the shared calibration and model structure. "Aligned" means
everything except λ, which changes for all of them when Phase 1 sets `POSTERIORS_FILE`.

| Application | Rain / kernels | Status |
|---|---|---|
| `src/ABC/` (the calibration itself) | SILO mean + extreme, real 2010–23 sequence | Phase 1 |
| `tcz-setup.R` → `tcz-sims.R`, `tcz-optimise.R` | SILO mean + extreme, random start, drying | Aligned; final runs in Phase 4 |
| `forecast.R` → public site | raster, one kernel, `season.days = 365` | Phase 2 |
| `null-model/predict-front-displacement.R` | raster / all-year SILO per year, `season.days = 365` | Phase 2 |
| `control-step-demo.R` (synthetic checks) | constant movement days | Aligned |

Archived on 2026-11-27 (`src/archive/README.md`): `tingley-sims/`, `pilbara-impact-sims/` (Dunlop
*et al.*, submitted 2025), `paruku-impact-sims/` and the original VRD ABC scripts.

## Phases

The phases are in dependency order: each uses only the outputs of earlier phases.

1. **Refit the Kimberley calibration** (`src/ABC/`, `src/cluster/`).
   - **Why:** λ (`dat/Posteriors_2026.RData`) was fitted with the Kimberley points' `ndays_1`
     rainy days. SILO, which the TCZ runs now use, counts about 50% more (42.6 against 28.3 at
     the TCZ points, r = 0.90). So the TCZ runs move toads further each year than the calibration
     assumed.
   - First check: compare `ndays_1` with SILO at the Kimberley points
     (`dat/waterpoint-data_Kimberley_trimmed.csv`), and find out where `ndays_1` came from.
   - Get daily SILO for the Kimberley cells. Reuse `rain_years()` (`src/tcz-sims/rain-years.R`,
     cached per cell). `dat/silo-actual-rain-days-cache.RData` (from `src/null-model/`) may
     already hold some of them.
   - Fit with the **same model structure as the TCZ runs**:
     - year-by-year SILO rainy days over the actual invasion years (following real history, not
       a random start);
     - the extreme-wet-year kernel by the same top-20% rule;
     - SILO for both kernels.
   - **Switch the colonisation threshold to c ≥ 2** (`realised_n >= 2` in `spread.pilb()`,
     `src/pprocess_functions.R`). Two arrivals are enough to found a population. The ABC, the TCZ
     runs and the forecast share this function, so it changes once, here, together with the
     refit.
   - Re-run the ABC (probably on the cluster) and score it against the 2022 and 2023 presence
     data, as `spread_Kimberley_ABC_concat.R` does.
   - Save the posteriors as `dat/Posteriors_<date>.RData`, compare them with
     `Posteriors_2026.RData`, and update the methods' calibration paragraph.
2. **Annual forecast** (`src/tcz-sims/forecast.R`).
   - Use the Phase 1 posteriors and SILO mean rainy days in place of the raster
     (`rain.id = "rainfall"`).
   - Add inter-annual variation as the optimisation runs have it: a random reference-year start
     per run, the extreme kernel and the drying matrix (via `run_sims(extreme.ref=, drying.spec=)`).
   - Factor the setup shared with `tcz-setup.R` into one path for both pipelines.
   - Move the null model (`src/null-model/predict-front-displacement.R`) onto the same path. Its
     static scenario takes SILO mean and extreme kernels, with real extreme years. Its per-year
     scenario counts Nov–May rainy days with `season.days = 212`. Then rerun it as an
     out-of-sample check of the Phase 1 calibration against the observed front movement.
   - Compare the new arrival-time forecast with the current one.
   - Regenerate the public site in Phase 5, not here, so it changes once, with the final text.
3. **Input-data review for the optimisation.** A checklist to work through with the project
   lead. Each item ends as confirmed, corrected, or an accepted assumption stated in the methods.
   - **Audit and units:**
     - uncosted works: "Erin B yards" and "Refuse tip" (fixed points now), plus 5 uncosted
       records inside costed units;
     - the orphan cost rows: "Two tree well" and "Sea Ranger base", $20k each, which match no
       audit record;
     - the fence enclosure memberships: Shelemar's 8 records, and nothing else inside Pivot 2 or
       Shamrock Gardens;
     - the 100 m merge rule.
   - **Fences:**
     - lengths, gates (is `n.gates` complete?) and `erad.units`;
     - **"Shamrock Stn Pivots"** (19.5 km, 3 gates, 7 erad-units), which needs a `Site.ID`, its
       enclosed records and a cost. It's left out for now and costs about $262k per response;
     - the square-fence assumption behind the fence area and `r.eff`.
   - **Costs:**
     - which workbook version is in use;
     - the eradication days calculation (N₀ = 700, detection probability 0.27, 20.8 days per
       erad-unit).
   - **Control parameters:** tank failure 1%, fence failure 5%, eradication failure 5%, and the
     timing of the annual audit.
   - **Natural points:**
     - the drying classes;
     - blank Living Waters records outside the TCZ, treated as Sometimes Dry;
     - the LaGrange points, all "Perennial" and never drying;
     - p = 0.25 in an average year.
   - **Setup:**
     - the invasion-front parameters (are they the latest fit?);
     - the 80 km buffer and the target points;
     - SILO cache freshness and the reference years.
4. **Remaining optimisation updates and final runs** (`src/tcz-sims/`).
   - **Maintenance costs:** a 50-year maintenance cost for each fence and tank, added to its
     capital cost.
     - Fences probably scale with length and gates (`distance` and `n.gates` are already on the
       units); tanks get a per-unit figure.
     - Replace the flat `design_cost(maintenance.per.year=)` hook with a per-unit column built in
       `build_tcz_units()` from the cost workbook.
   - **Rank on whole-life cost.** Rank on capital + maintenance + expected eradication cost, not
     capital alone. Eradication is now most of a fence's cost ($18.3M against $16.5M of capital
     in the full design).
     - `eval_design()` must return eradication events per fence row, not just their sum.
     - Use the new cost in the screening order and in best-value restoration.
   - Add Shamrock Stn Pivots once Phase 3 has its data.
   - Final runs with the Phase 1 posteriors and the Phase 3 corrections, at 100 / 1,000 / 1,000.
   - **Sensitivity analysis:** how far the design and its cost move as each assumption changes.
     The set is to be agreed in this phase's plan; candidates are the λ posterior spread,
     p.avg, the control failure rates and the eradication cost.
5. **Methods, results and reporting.**
   - **Methods** (`ms/tcz-optimisation-methods.qmd`):
     - the project lead fills the Geoscience Australia dataset names ("XX and YY", "step out
       details"), the soak count, and the Dunlop et al., Brown et al. and Phillips, Dempster et
       al. references;
     - the project lead reviews the "Auditing natural waterpoints" draft;
     - c ≥ 2 throughout (the simulation model and fence-failure text);
     - the refit and the forecast changes.
   - **Results:**
     - design, cost and CE tables;
     - breach risk with its bounds;
     - the eradication share of cost;
     - sensitivity;
     - figures: a design map and the arrival-time forecast.
   - **Program summary report:** short and plain-language. What to build, what it costs, how
     safe it is, and the main assumptions.
   - **Forecast site:** regenerate `shiny/web/data/` (`shiny/prep_app_data.R` →
     `shiny/export_web_data.R`), and update the About text for the new model.
   - **Merge:** a PR from `tcz-optimisation` into `main`.

## Open questions

- **Sensitivity set:** which parameters to vary, and over what ranges (Phase 4).
- **Summary report:** its audience and format (Phase 5).

## Detailed phase plans

(Filled in as each phase starts.)
