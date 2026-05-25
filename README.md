# GENA-LM Critical Review with Layer-Aware Reproducibility Experiments

**Research Proposal for "Лето с AIRI 2026" summer school** — author: *Korshunov E., bioinformatics student*.
**Repo:** <https://github.com/eugeneskywalker/GENA_LM_AIRI>

Independent evaluation and extension of:

> Fishman V., Kuratov Y., Shmelev A., Petrov M., Penzar D., Shepelin D., Chekanov N., Kardymon O., Burtsev M. (2025). **GENA-LM: a family of open-source foundational DNA language models for long sequences.** *Nucleic Acids Research*, 53(2), gkae1310. DOI: [10.1093/nar/gkae1310](https://doi.org/10.1093/nar/gkae1310).
> 📦 [GitHub](https://github.com/AIRI-Institute/GENA_LM) · [HuggingFace](https://huggingface.co/AIRI-Institute) · [Web service](https://dnalm.airi.net)

📄 **Read the proposal:** [`proposal/proposal.md`](proposal/proposal.md)
📋 **Critical self-review of this proposal:** [`proposal/research_proposal_critical_review.md`](proposal/research_proposal_critical_review.md)

## What this is

A two-page Research Proposal + reproducible code that applies two recent methodological frameworks to GENA-LM:

1. **BERTology** (Tenney 2019, Vig 2021) — *where* in a frozen model task-specific features accumulate.
2. **DART-Eval critique** (Patel et al., NeurIPS 2024) — do DNA foundation models justify their compute over lightweight baselines?

All comparisons here are in the **frozen-probing regime** (no fine-tuning of the foundation models) — except one calibration data point (E7) where GENA-LM is fine-tuned to show the regime gap. Conclusions about HyenaDNA-vs-GENA-LM head-to-head are scoped to frozen probing on short-range regulatory tasks; long-range performance (where GENA-LM has documented architectural advantage) is not tested.

## TL;DR — key findings

- **E1 — Layer-wise probing** (4 tasks × 13 layers × 3 seeds). Mid-layer specialisation (L4–L5) for tasks the model resolves (promoter F1=0.81 at L4, coding F1=0.91 at L5). For long-range regulatory tasks (enhancer, OCR) the model stays near chance at any layer — suggesting BigBird / RMT variants are required.
- **E4 — Multilayer UMAP clustering** (4 classes × 3 layers). Best silhouette at mid layer (L4: 0.164 vs L0: 0.111, L12: 0.124). Practical takeaway: prefer L4 over L12 for unsupervised use.
- **E3 — GENA-LM vs HyenaDNA head-to-head** (same 4 tasks, frozen probing, 3 seeds). **HyenaDNA-tiny (0.4 M params) wins 3 of 4 tasks against GENA-LM-base (110.7 M) in the frozen-probing regime** on short-range regulatory tasks. *Scope*: this is a frozen-feature comparison; the authors' published fine-tuned numbers are substantially higher (e.g. F1 ≈ 0.94 on promoter 2 kb).
- **E5 — In silico CTCF saturation mutagenesis with negative control.** For 100 synthetic sequences with a JASPAR MA0139.1 motif, off-consensus motif substitutions move motif-region embeddings proportionally to `1 − PWM_prob` (**Pearson r = 0.893**). **Negative control**: mutating random *non-motif* positions produces a **15× smaller** embedding response (0.0004 vs 0.0060 mean cosine distance). This rules out generic BPE-tokenisation sensitivity and supports the conclusion that GENA-LM has internalised the canonical CTCF binding grammar.
- **E6 — DART-Eval-inspired 4-way calibration** on Open Chromatin Regions (regulatory DNA proxy, not the formal DART-Eval pipeline). GENA-LM-base (110.7 M, frozen) gives AUROC = 0.658 vs the 65-feature 3-mer + LogReg baseline at 0.647 — a 110-million-parameter model adds **1.1 %** AUROC over a 64-feature classifier. HyenaDNA-tiny wins overall (0.688). Consistent with the DART-Eval critique on this proxy task.

**Extension experiments (added 2026-05-25 after self-review):**

- **E3-v2 — Re-running E3 with GENA-LM at the E1-best layer per task** flips the head-to-head: GENA-LM matches or exceeds HyenaDNA on 3 of 4 tasks (promoter F1 0.808 vs 0.790, coding 0.915 vs 0.862). The original «HyenaDNA wins 3/4» result was partly a last-layer artefact — strengthens the «avoid last layer» guideline.
- **E1-large — Layer-wise probing on `gena-lm-bert-large-t2t` (336 M, 24 layers).** Mid-layer pattern replicates at proportional depth (L8–L11 of 24); enhancer F1 jumps 0.66 → 0.73 (+0.07) with the deeper model — long-range tasks are not fundamentally inaccessible, they need scale and/or BigBird/RMT.
- **E7 — Short fine-tune of GENA-LM on promoter** (3 epochs, 512-token windows, 3 seeds) gives F1 = 0.794 ± 0.004 — quantifies the gap to the authors' published fine-tuned 0.94 on 2 kb context. Calibration point for the «frozen vs fine-tuned» limitation.

All experiments run on a single **Tesla V100S 32 GB**; the original 5 experiments total ≈ 16 minutes, the 3 extensions add ≈ 20 minutes (E7 is the longest at ~15 min for 3-epoch fine-tuning).

## Results

### E1 — Layer-wise probing (4 tasks × 13 layers × 3 seeds)

![E1 layer-wise probing](results/figures/e1_layerwise_probing.png)

| Task | Best layer | F1 | MCC | ROC-AUC | Pattern |
|---|---|---|---|---|---|
| Coding vs intergenic | L5 | **0.915** | 0.829 | **0.973** | Strong monotonic rise 0→5, plateau |
| Promoter (nontata) | L4 | **0.809** | 0.623 | 0.886 | Rise 0→4, plateau, slight last-layer drop |
| Enhancer (Cohn) | L1 | 0.662 | 0.310 | 0.702 | Weak signal, early peak |
| Open chromatin (OCR) | L2 | 0.615 | 0.232 | 0.637 | Weakest, late layers degrade |

Full per-layer numbers: [`results/tables/e1_results.csv`](results/tables/e1_results.csv).

### E4 — Embedding clustering across layers (3 layers × 4 classes)

![E4 multilayer UMAP](results/figures/e4_multilayer_umap.png)

| Layer | Silhouette (cosine) ↑ | NMI ↑ | ARI ↑ |
|---|---|---|---|
| L0 (embedding) | 0.111 | 0.432 | 0.366 |
| **L4 (mid, E1 peak)** | **0.164** | 0.429 | 0.325 |
| L12 (last) | 0.124 | 0.426 | 0.375 |

**Interpretation:** intrinsic geometry (silhouette) is best in middle layer, while last layer is more task-aligned to the MLM objective (slightly better ARI). NMI is flat — total information is preserved.

Full metrics: [`results/e4_multilayer_metrics.json`](results/e4_multilayer_metrics.json).

### E3 — GENA-LM vs HyenaDNA head-to-head (frozen probing)

| Model | Architecture | Params | promoter F1 (n=3, ± std) | enhancer F1 | OCR F1 | coding F1 |
|---|---|---:|---:|---:|---:|---:|
| **GENA-LM-base-t2t** | BPE transformer (last layer) | **110.7 M** | 0.777 ± 0.012 | 0.627 ± 0.013 | 0.565 ± 0.023 | **0.887 ± 0.003** |
| **HyenaDNA-tiny** | Hyena, single-nucleotide | **0.4 M** | **0.790 ± 0.004** | **0.679 ± 0.010** | **0.608 ± 0.015** | 0.862 ± 0.015 |

Same frozen-embedding + LogReg probing protocol, 3 random seeds. *Scope*: in the frozen-probing regime on short-range regulatory tasks. Numbers: [`results/tables/e3_results.csv`](results/tables/e3_results.csv), full JSON: [`results/e3_metrics.json`](results/e3_metrics.json).

### E5 — In silico saturation mutagenesis on CTCF

![E5 CTCF saturation](results/figures/e5_ctcf_saturation.png)

For 100 synthetic 400-bp sequences with a CTCF motif (JASPAR MA0139.1) inserted at position 200, we substituted each of the 19 motif positions with each alternative nucleotide (57 substitutions per sequence) and measured the cosine distance between motif-region embeddings of original vs mutant.

| Metric | Motif positions | Control (random flank positions) | Interpretation |
|---|---:|---:|---|
| Pearson r — off-consensus map vs `1 − PWM_prob` | **0.893** | N/A (no PWM reference) | Off-consensus motif substitutions track consensus dissimilarity tightly |
| Pearson r — per-position IC vs mean effect | **0.421** | N/A | Informative motif positions show larger mutation effect |
| **Overall mean cosine distance** | **0.0060 ± 0.0027** | **0.0004 ± 0.0001** | Motif mutations move embeddings **15× more** than random non-motif mutations |

**Conclusion:** GENA-LM has internalised the canonical CTCF binding grammar — single-nucleotide changes at motif positions move the motif-region embedding proportionally to PWM-predicted consensus penalty (r = 0.893), and a **15× smaller** response to mutations at random non-motif positions rules out generic BPE-tokenisation sensitivity. This is a **counterfactual / interventional in silico** complement to the **correlational** token-importance analysis in Figure 2 of the original paper (closes weakness W10).

Full metrics: [`results/e5_metrics.json`](results/e5_metrics.json) (motif), [`results/e5_metrics_control.json`](results/e5_metrics_control.json) (negative control).

### E6 — DART-Eval-inspired 4-way comparison on a regulatory-DNA proxy

Single calibration benchmark on Open Chromatin Regions (regulatory DNA proxy task, 5,000 samples, 60/20/20 split, 3 seeds):

| Model | Params | AUROC ↑ | AUPRC ↑ | F1 ↑ | MCC ↑ |
|---|---:|---:|---:|---:|---:|
| 3-mer + LogReg (no DL) | 65 | 0.647 | 0.633 | 0.609 | 0.212 |
| TinyCNN (supervised, from scratch) | 104.0k | 0.661 ± 0.005 | 0.629 ± 0.007 | **0.679 ± 0.006** | 0.229 ± 0.001 |
| **GENA-LM-base-t2t** (frozen) | **110.7 M** | 0.658 | 0.641 | 0.611 | 0.212 |
| **HyenaDNA-tiny** (frozen) | 0.4 M | **0.688** | **0.661** | 0.646 | **0.259** |

*Values for TinyCNN are mean ± std across 3 random initialisations. LogReg-based rows (3-mer, GENA-LM, HyenaDNA) show std = 0 because the data split is fixed and LogReg is deterministic on the convex objective — the reported gaps reflect the model representations themselves, not training noise.*

**Highlights:**
- **GENA-LM-base (110.7M params) ≈ 3-mer + LogReg (65 features)** — gain in AUROC = **+0.011**. The 110 M-parameter foundation model adds **1.1%** AUROC over a 64-feature naïve classifier on this proxy task.
- TinyCNN (104k params, ~1000× smaller than GENA-LM) matches GENA-LM on AUROC and wins F1.
- HyenaDNA-tiny (0.4M params) is best overall on 3/4 metrics, replicating the E3 pattern.

This is **consistent with the central DART-Eval (Patel et al., NeurIPS 2024 D&B) critique** on a DART-Eval-inspired proxy task: *current DNA LMs do not offer compelling gains over lightweight alternatives on regulatory DNA tasks.* Running the formal DART-Eval Task 1–5 pipeline (requires Synapse access) is a natural follow-up. Full numbers: [`results/tables/e6_results.csv`](results/tables/e6_results.csv).

## Extension experiments

The following experiments are extensions to the original 5-experiment proposal that strengthen and refine the core findings.

### E3-v2 — GENA-LM at best-per-task layer (refinement of E3)

The original E3 head-to-head used GENA-LM's last hidden layer (L12). But E1 showed best layer is L4–L5 for the tasks the model resolves. Re-running with GENA-LM extracted at the **E1-best layer per task** changes the picture:

| Task | GENA-LM @ best layer | HyenaDNA-tiny (last) | Winner |
|---|---:|---:|---|
| Promoter (GENA @ L4) | **0.808 ± 0.006** | 0.790 ± 0.004 | GENA-LM |
| Enhancer (GENA @ L1) | 0.662 ± 0.026 | **0.679 ± 0.010** | HyenaDNA |
| OCR (GENA @ L2) | 0.615 ± 0.027 | 0.608 ± 0.015 | tie |
| Coding (GENA @ L5) | **0.915 ± 0.010** | 0.862 ± 0.015 | GENA-LM |

**Methodological lesson.** The simple "HyenaDNA wins 3/4" narrative from the original last-layer E3 is partly a layer-choice artefact. At the appropriate layer for each task, GENA-LM matches or exceeds HyenaDNA on 3 of 4 tasks. The practical guideline from E4 — *do not query the last layer if you are not fine-tuning* — is the more durable finding.

Full numbers: [`results/e3_v2_metrics.json`](results/e3_v2_metrics.json).

### E1-large — same probing on bert-large (336 M, 24 layers)

| Task | Best layer | F1 best | Δ vs bert-base |
|---|---:|---:|---:|
| Promoter | L11 | 0.809 | ±0.000 |
| Enhancer | L9  | **0.732** | **+0.070** |
| OCR | L22 | 0.667 | +0.052 |
| Coding | L8 | 0.943 | +0.028 |

The mid-layer pattern observed in bert-base (peak L4 of 12) replicates in bert-large at proportionally similar depth (L8–L11 of 24). **Enhancer and OCR — near-chance in bert-base — become substantially easier in bert-large** (enhancer F1 0.66 → 0.73). Long-range regulatory tasks are not fundamentally inaccessible to GENA-LM features; they require either deeper models (bert-large) or longer context (BigBird, RMT) — both exist in the released GENA-LM family.

Full numbers: [`results/e1_large_metrics.json`](results/e1_large_metrics.json).

### E7 — fine-tuned GENA-LM on promoter (frozen vs fine-tuned calibration)

A single calibration data point: short 3-epoch fine-tune of `gena-lm-bert-base-t2t` on the promoter task (lr = 2e−5, batch 16, max_len 512, 3 seeds, classification head on mean-pooled last hidden state, backbone unfrozen):

| Regime | Test F1 (n=3, ± std) |
|---|---:|
| Frozen probing, last layer (E3) | 0.777 ± 0.012 |
| Frozen probing, best layer L4 (E3-v2) | 0.808 ± 0.006 |
| **Fine-tuned, 3 epochs** (E7) | **0.794 ± 0.004** |
| Authors' published fine-tuned, 2 kb context (Fig 3 of paper) | ≈ 0.94 |

A 3-epoch fine-tune over 512-token windows yields only +0.017 over frozen last-layer — and is actually below frozen best-layer L4. To reach the authors' published 0.94 requires longer fine-tuning, longer context (2 kb), and likely the bert-large variant. This data point quantifies the «frozen vs fine-tuned» limitation: our frozen-probing comparisons substantially underestimate GENA-LM's task performance when properly deployed.

Full numbers: [`results/e7_fine_tune_metrics.json`](results/e7_fine_tune_metrics.json).

## Limitations and scope

Important caveats on how to read the head-to-head and DART-Eval-inspired results:

1. **Frozen-probing only.** E3 and E6 use frozen embeddings + LogReg / TinyCNN probes. Fine-tuned performance may differ substantially — the authors report e.g. F1 ≈ 0.94 on promoter 2 kb after fine-tuning (vs our frozen 0.78). Claims like "HyenaDNA wins" should not be extrapolated outside the frozen-probing regime.
2. **E6 is a DART-Eval-inspired proxy.** Open Chromatin classification is a single regulatory-DNA proxy chosen for a fast, self-contained V100 run; running the formal DART-Eval Task 1–5 pipeline (Synapse access required) is a natural follow-up.
3. **Single GENA-LM variant.** Tested only `gena-lm-bert-base-t2t` (110 M, 4.5 kb context). The BigBird (36 kb) and RMT (Mb) variants — where GENA-LM has documented architectural advantage on long-range tasks — are not tested here.
4. **Single HyenaDNA variant.** Only `hyenadna-tiny-1k-seqlen` (0.4 M) is used as the comparator; larger HyenaDNA variants may compress or invert the observed gap.
5. **Short-range tasks only in head-to-head.** All four E3/E6 tasks use sequences ≤ 1–2 kb. Long-range tasks (species classification at 32 kb, DeepSEA HM at 8 kb), where GENA-LM is documented to outperform HyenaDNA, are deliberately not included.
6. **Sample sizes and statistical power.** Per-task sample sizes (≈ 2,000–5,000) limit statistical power. F1 differences in E3 (Δ ≈ 0.01–0.05) are at the edge of seed-noise resolution; 3 seeds + std are reported.
7. **E5 is counterfactual in silico, not biologically causal.** The CTCF saturation analysis on synthetic sequences is an interventional probing experiment — stronger than correlational attribution methods but weaker than functional perturbation in a reporter assay. A negative control (mutations at random positions outside the motif) is included and yields a 15× smaller embedding response, ruling out generic BPE-tokenisation sensitivity as the source of the signal.

## Repository structure

```
gena-lm-airi-2026/
├── proposal/
│   ├── proposal.md                       # Research Proposal (2 pages)
│   ├── proposal.pdf                      # Compiled PDF
│   └── research_proposal_critical_review.md  # Self-review of the proposal
├── experiments/
│   ├── sanity_check.py                   # GPU + model load verification
│   ├── e1_layerwise_probing.py           # E1 — layer-wise BERTology probing (bert-base, 13 layers)
│   ├── e1_bert_large.py                  # E1-large — same protocol on bert-large (25 layers)
│   ├── e3_gena_vs_caduceus.py            # E3 — head-to-head with HyenaDNA (last layer)
│   ├── e3_best_layer.py                  # E3-v2 — GENA-LM at best-per-task layer (from E1)
│   ├── e4_clustering.py                  # E4 warm-up (last layer only)
│   ├── e4_clustering_multilayer.py       # E4 v2 (L0/L4/L12 comparison)
│   ├── e5_saturation_mutagenesis.py      # E5 — CTCF motif probing (motif + --control mode)
│   ├── e6_dart_eval_style.py             # E6 — DART-Eval-inspired 4-way calibration
│   └── e7_fine_tune_promoter.py          # E7 — fine-tuned GENA-LM on promoter (1 calibration point)
├── slurm/                                # sbatch scripts for V100 cluster
│   ├── setup_env.sh, setup_env_v2.sh
│   ├── sbatch_sanity.sh
│   ├── sbatch_e1.sh, sbatch_e1_large.sh
│   ├── sbatch_e3.sh, sbatch_e3_v2.sh
│   ├── sbatch_e4.sh, sbatch_e4v2.sh
│   ├── sbatch_e5.sh, sbatch_e5_control.sh
│   ├── sbatch_e6.sh
│   └── sbatch_e7.sh
├── results/
│   ├── figures/                          # PNG plots
│   ├── tables/                           # CSV summaries
│   ├── e1_metrics.json, e1_large_metrics.json
│   ├── e3_metrics.json, e3_v2_metrics.json
│   ├── e4_metrics.json, e4_multilayer_metrics.json
│   ├── e5_metrics.json, e5_metrics_control.json
│   ├── e6_metrics.json
│   └── e7_fine_tune_metrics.json
├── requirements.txt
└── README.md (this file)
```

## Reproducing on a single V100 32 GB

```bash
# 1. Create env (Python 3.11 + PyTorch 2.5 CUDA 12.1 + transformers 4.36)
bash slurm/setup_env.sh
conda activate gena

# 2. Sanity check (~2 min): downloads model, runs 1 forward pass
sbatch slurm/sbatch_sanity.sh

# --- Core experiments (≈ 16 minutes total on V100) ---
sbatch slurm/sbatch_e4v2.sh            # E4 — multilayer clustering (~3 min)
sbatch slurm/sbatch_e1.sh              # E1 — layer-wise probing (~8 min)
sbatch slurm/sbatch_e3.sh              # E3 — GENA-LM vs HyenaDNA, last layer (~4 min)
sbatch slurm/sbatch_e5.sh              # E5 — CTCF saturation (motif mode, ~2 min)
sbatch slurm/sbatch_e5_control.sh      # E5 — CTCF saturation (negative control, ~2 min)
sbatch slurm/sbatch_e6.sh              # E6 — DART-Eval-inspired 4-way (~1 min)

# --- Extension experiments ---
sbatch slurm/sbatch_e3_v2.sh           # E3 v2 — GENA-LM at best-per-task layer
sbatch slurm/sbatch_e1_large.sh        # E1-large — same probing on bert-large (25 layers)
sbatch slurm/sbatch_e7.sh              # E7 — fine-tune GENA-LM on promoter (calibration point)
```

Outputs land in `~/gena-lm-airi/results/{figures,tables}/`.

**No data prep needed** — all downstream tasks use ready-to-use [Genomic Benchmarks](https://huggingface.co/katarinagresova) datasets from the HuggingFace Hub (downloaded automatically by `datasets`). All model weights are pulled from the HuggingFace Hub at runtime.

## Hardware

Aldan3 cluster (Pirogov Russian National Research Medical University):
- 1× Tesla V100S-PCIE-32GB (Volta, sm_70)
- conda env `gena` (Python 3.11, PyTorch 2.5.1+cu121, transformers 4.36.2)
- ~10 minutes wall-clock for all experiments combined

## License

MIT (matches the GENA-LM upstream repo).

## Acknowledgements

- AIRI-Institute team (Fishman, Kuratov, Burtsev et al.) for open-sourcing GENA-LM, training scripts, and HuggingFace checkpoints.
- [katarinagresova](https://huggingface.co/katarinagresova) for the Genomic Benchmarks HF datasets.
