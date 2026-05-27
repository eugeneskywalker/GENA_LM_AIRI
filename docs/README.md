# docs/

Interactive landing page that accompanies the Research Proposal — a 5-minute scrollytelling
companion to `../proposal/proposal.pdf` with interactive charts of the experiment results.

Once built, served at: <https://eugeneskywalker.github.io/GENA_LM_AIRI/>

## Data sources

The landing reads pre-aggregated JSON in `data/` (derived from `../results/*.json`):

| File | Contents |
|---|---|
| `data/e1_layerwise.json` | Layer-wise probing F1 across 13 layers × 4 tasks (E1) |
| `data/e3_comparison.json` | GENA-LM vs HyenaDNA, two views: last-layer (E3) and best-per-task-layer (E3-v2) |
| `data/e5_ctcf_saturation.json` | CTCF 19×4 saturation map + PWM reference + negative-control comparison |
| `data/experiments_summary.json` | Aggregate metadata across all 8 experiments, `weaknesses_closed` map, `limitations_scope` list |

Numbers are rounded to 3–4 significant figures for readability; the original JSONs in `../results/`
are authoritative for re-runs and citations.

## Build & deploy

The landing is a single static HTML file served from this folder via
GitHub Pages (`Settings → Pages → main /docs`). No build step is required.
