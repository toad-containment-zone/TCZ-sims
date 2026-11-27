# Archived model applications

These scripts applied earlier versions of the model and calibration. They are kept for the record but
are **not maintained against the current engine** (`src/pprocess_functions.R`). Since 2026-11-27 the
engine colonises at c ≥ 2 and maps rainy days to movement over a 212-day (Nov–May) season. These
scripts were written for c ≥ 3, a 365-day season, the `ndays_1` raster and the old VRD
posteriors (`dat/Posteriors.RData`), so running them now mixes old and new.

To reproduce one as it was, check out commit `b6432e9`, the last commit before the Phase 1 engine
changes, and run it from the repo root at its old path (`src/<folder>/...`).

- `pilbara-impact-sims/`: Dunlop *et al.* (submitted 2025). `timing-to-pilbara.R` is Method 3
  (the point-process model) and `spread-rate-rainfall.R` is Method 2. `timing-to-tcz.R` was
  already superseded by `src/tcz-sims/forecast.R`.
- `tingley-sims/`: the original Pilbara barrier scenarios.
- `paruku-impact-sims/`: arrival-time simulation at Paruku.
- `ABC/`: the original VRD ABC (`spread_model.R`, `spread_table_generator.R`, `SpreadABC_concat.R`)
  and a Kimberley demo run.

Bringing any of these back means moving it onto the shared calibration and model structure first
(`docs/roadmap.md`, "One calibration and model structure").
