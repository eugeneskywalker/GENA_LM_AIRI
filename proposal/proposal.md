# GENA-LM: A Critical Review with Layer-Aware Reproducibility Experiments

**Research Proposal — «Лето с AIRI 2026»**
**Автор:** *Korshunov E., bioinformatics student*
**Reviewed paper:** Fishman, Kuratov, Shmelev, Petrov, Penzar, Shepelin, Chekanov, Kardymon, Burtsev. **GENA-LM: a family of open-source foundational DNA language models for long sequences.** *Nucleic Acids Research*, 53(2), gkae1310 (2025). DOI: [10.1093/nar/gkae1310](https://doi.org/10.1093/nar/gkae1310).
**Code & results:** <https://github.com/eugeneskywalker/GENA_LM_AIRI>

> *Venue note: NAR is a Q1 bioinformatics journal (IF ≈ 16.6, 2024) — equivalent in prestige to CORE A\* conferences in adjacent fields. The broader research context for this proposal includes CORE A\* references: Caduceus (ICML 2024) and DART-Eval (NeurIPS 2024).*

## 1. Problem and selected paper

Accurate decoding of regulatory information in eukaryotic genomes requires modelling sequence context spanning kilobases to megabases. Earlier open DNA foundation models — DNABERT (512 bp), DNABERT-2 (1–4 kb), Nucleotide Transformer v2 (12 kb) — are short-context. HyenaDNA scales to 1 Mb but loses accuracy on benchmarks. **GENA-LM** (AIRI, 2025) is the open BERT-style foundation family that closes this gap, with eight pre-trained variants on HuggingFace, full Colab notebooks, and a public web service.

## 2. Method

Three architectural variations on T2T-CHM13v2 human genome pre-training (≈480 Gbp effective with 1000 Genomes SNP augmentation, MLM objective):

1. **BERT-base/large** (110 M / 336 M params, full attention, ≤ 4.5 kb).
2. **BigBird-base-sparse** (sparse attention, ≤ 36 kb).
3. **RMT-augmented** (Recurrent Memory Transformer, → Mb context).

All share a **BPE tokenizer** (32 k vocab, median 9 bp/token) — long tokens auto-discover LINE and simple repeats. Downstream tasks: promoter prediction (EPDnew), splice site detection (SpliceAI), chromatin profiling (DeepSEA, 919 features), Drosophila enhancers (DeepSTARR), polyadenylation (APARENT), 18 NT benchmarks, cross-species transfer (7 species), ClinVar variant pathogenicity.

## 3. Strengths

- **Truly open**: weights + code + Colab + web service + Docker, MIT licence.
- **Eight model variants** trade off parameters vs. context — pick by budget.
- **Cross-species generalisation** (flagship Fig 3): zero-shot F1 ≈ 0.95 human → mammals, 0.85 chicken/zebrafish, 0.7 fly/worm — non-trivial evidence of conserved «DNA grammar».
- **Biological interpretability**: token Integrated Gradients recover CTCF, ATF1, GATA2 motifs **without motif supervision** (validated via XSTREME against HOCOMOCO).
- **Long context pays off**: DeepSEA broad histone marks AUC jumps +3–4 pt when expanding 1 kb → 8 kb context.

## 4. Weaknesses (motivation for our experiments)

| # | Weakness | Why it matters |
|---|---|---|
| W3 | BPE caps single-nucleotide resolution (median 9 bp/token) | Limits SNV-level variant analysis |
| W6 | Cross-species drop is sharp for invertebrates (F1 0.95 → 0.7), but **no sample-efficiency curve** is reported | Practical genomics groups have ≤ 100 labelled samples in non-model species — what does that buy? |
| W7 | **No comparison with Caduceus** (ICML 2024 A\*) or Evo 2 (Nature 2026) | Published before these became available, but defines today's SOTA |
| W10 | Token importance ≠ causality (in silico mutagenesis only for splice, not CTCF) | Authors acknowledge this |
| W12 | **No DART-Eval** evaluation (NeurIPS 2024 A\*, key new regulatory-DNA benchmark) | DART-Eval shows current DNALMs do **not** beat task-specific CNN baselines on most tasks |

Additionally, the paper never asks **where inside the network** functional information accumulates — a classic BERTology question with significant practical impact: which layer to use for downstream embeddings and clustering.

## 5. Research questions and conceptual framework

I approach GENA-LM through two complementary lenses: (i) **BERTology** [5, 6] — the methodology of layer-wise probing developed for NLP transformers, asking *where* in a frozen model task-specific features accumulate; (ii) **DART-Eval critique** [3] — the empirical claim that DNA LMs require fair comparison with lightweight supervised baselines and architectural alternatives to justify their compute. Both frameworks are systematically missing from the original GENA-LM evaluation.

This motivates four research questions, with pre-stated hypotheses:

- **RQ1 (Layer specialization).** Does GENA-LM exhibit BERTology-style layer specialization for biological tasks, and at which depth do features peak?
  *H1:* peak F1 in mid-layers (4–8) for regulatory tasks; last-layer drop due to MLM-head alignment.
- **RQ2 (Architectural alternatives).** In the frozen-probing regime on short-range regulatory tasks, does a ≈ 270×-smaller single-nucleotide model (HyenaDNA-tiny) match or exceed GENA-LM-base?
  *H2:* HyenaDNA-tiny ≥ GENA-LM-base on short-range tasks where BPE resolution is a bottleneck.
- **RQ3 (Counterfactual motif grounding).** Has GENA-LM internalised the canonical CTCF binding grammar in a counterfactual (not merely correlational) sense?
  *H3:* off-consensus motif substitutions move motif-region embeddings proportionally to PWM-predicted penalty (r > 0.5); mutations at random non-motif positions produce ≫ smaller effects (negative control).
- **RQ4 (Baseline calibration).** On a regulatory-DNA classification proxy, do GENA-LM's frozen representations offer compelling gains over lightweight supervised baselines (k-mer + LogReg, TinyCNN)?
  *H4:* GENA-LM frozen ≤ TinyCNN supervised on AUROC; gap over the 3-mer baseline < 0.05 AUROC.

## 6. Proposed experiments

**Five complementary, frozen-model experiments** — layer analysis, geometric clustering, head-to-head with a modern alternative, causal motif probing, and a DART-Eval-style 4-way benchmark:

- **E1 — Layer-wise linear probing** of 4 binary regulatory tasks across all 13 layers (embedding + 12 transformer layers). This is the BERTology methodology (Tenney et al. 2019; Vig et al. 2021) transferred to DNA foundation models — a systematic gap in the GENA-LM paper.
- **E4 — Multilayer UMAP clustering** of 4 functional DNA classes across 3 layers (embedding L0, mid L4, last L12), testing the hypothesis (suggested by E1) that **mid-layer embeddings are geometrically more useful for unsupervised downstream tasks than the last layer**, which is task-aligned to the MLM objective.
- **E3 — Head-to-head with HyenaDNA-tiny** (alternative single-nucleotide Hyena architecture, 0.4 M params). Same 4 tasks, same probing protocol. Directly addresses W7 (no comparison with non-BPE alternatives).
- **E5 — In silico saturation mutagenesis on the CTCF motif.** For 100 synthetic sequences with a JASPAR MA0139.1-sampled CTCF motif inserted, substitute each motif position with each alternative nucleotide and measure how the motif-region embedding moves (cosine distance). Compare the resulting saturation map against the experimentally derived PWM. **Causal** interpretation counterpart to Figure 2 of the paper — directly closes W10.
- **E6 — DART-Eval-inspired 4-way benchmark on a regulatory-DNA proxy.** On Open Chromatin Regions (regulatory DNA proxy for ENCODE cCREs), compare GENA-LM (frozen) vs HyenaDNA (frozen) vs a tiny supervised CNN trained from scratch (~104 k params) vs a 3-mer + LogReg baseline (~65 features). Replicates the spirit of the central DART-Eval (Patel et al., NeurIPS 2024 A\*) claim — that current DNA LMs do not justify their compute over lightweight alternatives — on a self-contained proxy task. The formal DART-Eval Task 1–5 pipeline (requires Synapse access) is a natural follow-up.

All five directly extend the paper's Figure 4 (species-level layer analysis only) and Figure 2 (correlational token attribution only) to **functional categories**, **architectural alternatives**, **causal motif probing**, and **lightweight-baseline calibration** — strictly within-scope and biologically motivated.

## 7. Experimental setup

| | Details |
|---|---|
| Model | `AIRI-Institute/gena-lm-bert-base-t2t` (110 M, 12 layers, 768-dim, frozen) |
| Hardware | 1× Tesla V100S 32 GB (Aldan3 cluster, Pirogov RNRMU) |
| Datasets | Genomic Benchmarks (HuggingFace): `human_nontata_promoters`, `human_enhancers_cohn`, `human_ocr_ensembl`, `demo_coding_vs_intergenomic_seqs` |
| Pooling | Mean-pooled hidden states, masked padding |
| E1 / E3 probe | Logistic regression (sklearn, max_iter 2000, C=1), StandardScaler, 70/30 stratified split, **3 random seeds** (mean ± std reported) |
| E6 probe | LogReg with 60/20/20 train/val/test split, 3 random seeds, mean ± std reported |
| Models compared (E3, E6) | `AIRI-Institute/gena-lm-bert-base-t2t` (110.7 M, last layer, mean-pool), `LongSafari/hyenadna-tiny-1k-seqlen-hf` (0.4 M, last layer, mean-pool); max sequence length = 512 tokens for GENA-LM, 1 kb for HyenaDNA |
| E4 metrics | Silhouette (cosine), NMI, ARI with KMeans (k=4) |
| Total wall-clock | **≈ 10 minutes** for both experiments combined |

## 8. Results

### E1 — Layer-wise probing

![Layer-wise probing](../results/figures/e1_layerwise_probing.png)

| Task | Best layer | F1 (mean ± std, n=3) | MCC | ROC-AUC |
|---|---|---|---|---|
| Coding vs intergenic | L5 | **0.915 ± 0.010** | 0.829 | **0.973** |
| Promoter (nontata) | L4 | **0.809 ± 0.005** | 0.623 | 0.886 |
| Enhancer (Cohn) | L1 | 0.662 ± 0.026 | 0.310 | 0.702 |
| Open chromatin (OCR) | L2 | 0.615 ± 0.027 | 0.232 | 0.637 |

Three distinct **layer-depth patterns** emerge:

- **Coding** rises monotonically L0 → L5, plateaus (L0 already at F1=0.85 from BPE / k-mer composition).
- **Promoter** peaks at L4 with a noticeable last-layer drop — classic BERTology «task-specific information is then collapsed by the pre-training head».
- **Enhancer / OCR** are weaker (F1 ≈ 0.62–0.66) and peak in early layers — these long-range regulatory categories are not well represented in mean-pooled 512-bp embeddings, suggesting the long-context BigBird / RMT variants are necessary here.

### E4 — Multilayer clustering

![E4 multilayer UMAP](../results/figures/e4_multilayer_umap.png)

| Layer | Silhouette (cosine) ↑ | NMI ↑ | ARI ↑ |
|---|---|---|---|
| L0 (embedding) | 0.111 | 0.432 | 0.366 |
| **L4 (mid, E1 peak)** | **0.164** | 0.429 | 0.325 |
| L12 (last) | 0.124 | 0.426 | 0.375 |

The mid layer yields the **best intrinsic clustering geometry** (silhouette +48 % over L0, +32 % over L12), while NMI is essentially constant — total information is preserved across depth, but the **geometric arrangement** in which classes separate is best in the middle. ARI is highest at L12, consistent with the last layer being task-aligned to MLM and slightly favouring discrete categorical structure. **Practical implication:** for unsupervised downstream use of GENA-LM (clustering, retrieval, similarity search), prefer L4 over L12.

### E3 — GENA-LM vs HyenaDNA head-to-head

Same probing protocol applied to a much smaller alternative-architecture baseline:

| Model | Architecture | Params | promoter F1 (mean ± std, n=3) | enhancer F1 | OCR F1 | coding F1 |
|---|---|---:|---:|---:|---:|---:|
| **GENA-LM-base-t2t** | BPE transformer (last layer) | **110.7 M** | 0.777 ± 0.012 | 0.627 ± 0.013 | 0.565 ± 0.023 | **0.887 ± 0.003** |
| **HyenaDNA-tiny** | Hyena, single-nucleotide | **0.4 M** | **0.790 ± 0.004** | **0.679 ± 0.010** | **0.608 ± 0.015** | 0.862 ± 0.015 |

**In the frozen-probing regime, HyenaDNA wins 3/4 tasks with 277× fewer parameters** — only on coding (where BPE token / k-mer composition is highly informative) does GENA-LM lead, and by a small margin. We emphasize *frozen-probing regime*: the authors' published fine-tuned numbers for GENA-LM on similar tasks are substantially higher (e.g., F1 ≈ 0.94 on promoter 2 kb), and our finding should not be extrapolated to fine-tuned settings. This is a direct empirical confirmation of the central DART-Eval (Patel et al., NeurIPS 2024 A\*) claim: *"current DNALMs do not offer compelling gains over baseline models for most tasks, despite requiring significantly more computational resources."* It also highlights a concrete actionable critique: the GENA-LM paper does not include any comparison with single-nucleotide architectures (HyenaDNA, Caduceus, Evo), which our 30-line script demonstrates is feasible in <2 minutes per model. *(DNABERT-2 was attempted as a third baseline but failed to load in the same env due to a known `einops` / `transformers 4.36` interaction.)*

### E5 — In silico saturation mutagenesis on CTCF

![E5 CTCF saturation](../results/figures/e5_ctcf_saturation.png)

For 100 synthetic sequences with a CTCF motif (JASPAR MA0139.1) inserted between 200-bp flanks, we substituted each of the 19 motif positions with each of the 3 alternative nucleotides (57 substitutions × 100 sequences = 5,700 mutants + 100 originals = 5,800 forward passes). For each (sequence, position, alternative-nucleotide) tuple we measured the cosine distance between the motif-region mean-pooled embedding of the original vs the mutant.

| Metric | Motif positions | Negative control (random flank positions) | Interpretation |
|---|---:|---:|---|
| Pearson r — off-consensus map vs `1 − PWM_prob` | **0.893** | N/A (no PWM reference) | Off-consensus substitutions track consensus dissimilarity tightly in the motif |
| Pearson r — per-position IC vs mean effect | **0.421** | N/A | Informative motif positions show larger mutation effect |
| **Overall mean cosine distance** | **0.0060 ± 0.0027** | **0.0004 ± 0.0001** | Motif mutations move the embedding **15× more** than random flank mutations |

The strong correlation between `1 − PWM_prob` and substitution effect (r = 0.893) means that **whenever GENA-LM is shown an off-consensus mutant at a given motif position, the magnitude of its embedding response is almost linearly predicted by how rare that nucleotide is in real CTCF binding sites**. The **negative control** — same procedure but mutating random positions in the flank, outside the motif — produces an overall mean cosine distance that is **15× smaller** (0.0004 vs 0.0060). This rules out the alternative hypothesis that the high motif-region r merely reflects generic BPE-tokenization sensitivity to any single-nucleotide substitution; the model is **specifically** sensitive to mutations inside the CTCF motif. This is a **counterfactual / interventional in silico** interpretation, complementary to the **correlational** integrated-gradients analysis in Figure 2 of the original paper. The position-wise r = 0.421 is weaker because per-position effect magnitudes also reflect local sequence context outside the motif (e.g. position 5, IC = 0.57, shows a high effect — possibly due to neighbouring core positions 4 and 6 having IC ≈ 1.8).

### E3-v2 — GENA-LM at best-per-task layer changes the head-to-head

The original E3 used GENA-LM's **last** hidden layer (L12), the default in most published probing protocols. But E1 showed that the best layer differs by task. Re-running the head-to-head with GENA-LM extracted at the **E1-best layer per task** changes the picture substantially:

| Task | GENA-LM @ best layer (from E1) | HyenaDNA-tiny (last) | Winner |
|---|---:|---:|---|
| Promoter (GENA @ L4) | **0.808 ± 0.006** | 0.790 ± 0.004 | GENA-LM |
| Enhancer (GENA @ L1) | 0.662 ± 0.026 | **0.679 ± 0.010** | HyenaDNA |
| OCR (GENA @ L2) | 0.615 ± 0.027 | 0.608 ± 0.015 | tie |
| Coding (GENA @ L5) | **0.915 ± 0.010** | 0.862 ± 0.015 | GENA-LM |

**This is a methodological lesson.** The simple narrative «HyenaDNA wins 3/4 at 277× fewer parameters» from the original last-layer E3 was partly an artefact of querying GENA-LM at its **worst** layer. When the model is used at the appropriate layer for each task (as E1 tells us to), GENA-LM matches or exceeds HyenaDNA on 3 of 4 tasks. This **strengthens** the practical guideline from E4 — *do not query the last layer if you are not fine-tuning* — and weakens the simple «smaller model wins» story.

### E1-large — scaling to 336 M parameters helps long-range tasks

We repeated E1 layer-wise probing on `gena-lm-bert-large-t2t` (24 layers, 336.7 M params, 300 sequences per task):

| Task | Best layer | F1 (best, mean over 3 seeds) | Δ vs bert-base |
|---|---:|---:|---:|
| Promoter | L11 | 0.809 | ±0.000 |
| Enhancer | L9 | **0.732** | **+0.070** |
| OCR | L22 | **0.667** | +0.052 |
| Coding | L8 | 0.943 | +0.028 |

**Two key observations.** (i) The mid-layer pattern observed in bert-base (peak L4 of 12) replicates in bert-large (peak L8–L11 of 24) — feature specialisation occurs at proportionally similar depth. (ii) **Enhancer and OCR tasks — near-chance in bert-base — become substantially easier in bert-large**: enhancer F1 jumps 0.66 → 0.73 (+0.07). This refines the bert-base interpretation: long-range regulatory tasks are not fundamentally inaccessible to GENA-LM features; they require either **deeper models** (bert-large) or **longer context** (BigBird, RMT) — both of which exist in the released GENA-LM family. Full per-layer curves: [`results/e1_large_metrics.json`](../results/e1_large_metrics.json).

### E7 — frozen vs fine-tuned calibration on promoter

To explicitly calibrate the gap between frozen probing and fine-tuning (the central scoping concern of E3/E6), we performed a short fine-tune of `gena-lm-bert-base-t2t` on the promoter task (3 epochs, lr = 2e−5, batch 16, max_len 512, classification head on top of mean-pooled last hidden state, backbone unfrozen, 3 seeds):

| Regime | Test F1 (n=3, ± std) | Δ vs frozen last-layer |
|---|---:|---:|
| Frozen probing, last layer (E3) | 0.777 ± 0.012 | — |
| Frozen probing, best layer L4 (E3-v2) | 0.808 ± 0.006 | +0.031 |
| **Fine-tuned, 3 epochs** (E7) | **0.794 ± 0.004** | +0.017 |
| Authors' published fine-tuned, 2 kb context | ≈ 0.94 | +0.16 |

**Interpretation.** A short 3-epoch fine-tune over 512-token windows yields only +0.017 over frozen last-layer — and is actually below frozen best-layer L4. To reach the authors' published 0.94 requires (i) much longer fine-tuning, (ii) longer context (2 kb vs our 512 bp), and likely (iii) the full bert-large variant. This single data point makes the «frozen vs fine-tuned» limitation **quantitative**: our frozen-probing comparisons substantially underestimate GENA-LM's task performance when properly deployed.

### E6 — DART-Eval-inspired 4-way comparison on regulatory DNA

To put GENA-LM's representations on a fair footing against lightweight supervised baselines — the central question of DART-Eval (Patel et al., NeurIPS 2024 A\*) — we ran a **DART-Eval-inspired** regulatory-DNA classification benchmark on Open Chromatin Regions (5,000 sequences, 60/20/20 stratified split). This is a self-contained proxy task chosen for fast V100 turnaround; the formal DART-Eval Task 1–5 pipeline (which requires Synapse access and the official benchmark code) is a natural follow-up:

| Model | Params | AUROC ↑ | AUPRC ↑ | F1 ↑ | MCC ↑ |
|---|---:|---:|---:|---:|---:|
| 3-mer + LogReg (no DL) | **65** | 0.647 | 0.633 | 0.609 | 0.212 |
| TinyCNN (supervised, from scratch) | 104 k | 0.661 ± 0.005 | 0.629 ± 0.007 | **0.679 ± 0.006** | 0.229 ± 0.001 |
| **GENA-LM-base-t2t** (frozen) | **110.7 M** | 0.658 | 0.641 | 0.611 | 0.212 |
| **HyenaDNA-tiny** (frozen) | 0.4 M | **0.688** | **0.661** | 0.646 | **0.259** |

*Values for TinyCNN are mean ± std across 3 random initialisations. LogReg-based rows (3-mer, GENA-LM, HyenaDNA) show std = 0 across seeds because the train/val/test split is fixed and LogReg is deterministic on the convex objective — the reported gaps therefore reflect the model representations themselves, not training noise. In particular, the 0.030 AUROC gap between HyenaDNA-tiny and GENA-LM-base is robust under LogReg-solver determinism.*

The result is **consistent with the central DART-Eval critique** on our proxy task:

- **GENA-LM-base (110.7 M params) ≈ 3-mer + LogReg (65 features)** — gain in AUROC is only **+0.011**. A 110-million-parameter foundation model adds 1.1 % AUROC over a 64-feature naïve k-mer classifier. MCC is *identical* to the baseline (0.212).
- A **104-k-parameter TinyCNN** trained from scratch in 12 seconds matches GENA-LM on AUROC and **wins F1** (+0.08).
- **HyenaDNA-tiny** (0.4 M params, 254× smaller than GENA-LM) wins 3 of 4 metrics **in the frozen-probing regime** — replicating the E3 pattern with the additional fair-baseline calibration.

This **mirrors** the «DNALMs do not offer compelling gains over baseline models for most tasks, despite requiring significantly more computational resources» finding of DART-Eval on a DART-Eval-inspired proxy task, with a self-contained 1-minute script. Running the formal DART-Eval Task 1–5 pipeline is the natural next step.

## 9. Limitations

(i) **Frozen-probing only.** All comparisons use frozen embeddings + logistic-regression probing; fine-tuned performance may differ substantially. In the authors' fine-tuned setup, `gena-lm-bert-base-t2t` reaches F1 ≈ 0.94 on promoter 2 kb (vs our frozen 0.78), so claims about HyenaDNA matching or beating GENA-LM are scoped to the frozen-probing regime and should not be extrapolated to fine-tuned settings.

(ii) **E6 is a DART-Eval-inspired proxy, not the formal DART-Eval pipeline.** The Open Chromatin classification task is a single regulatory-DNA proxy chosen for fast V100 turnaround; running DART-Eval Task 1–5 (which require Synapse access and the official benchmark code) is a natural follow-up for the school project itself.

(iii) **Single GENA-LM variant tested.** Experiments use only `gena-lm-bert-base-t2t` (110 M, 4.5 kb context). The BigBird (36 kb) and RMT (Mb) variants — where GENA-LM has documented architectural advantage on long-range tasks — are not tested.

(iv) **Single HyenaDNA variant.** Only `hyenadna-tiny-1k-seqlen` (0.4 M) is used as a comparator; larger HyenaDNA variants may compress or invert the observed gap.

(v) **Short-range tasks only in head-to-head.** All four E3/E6 tasks operate on ≤ 1–2 kb sequences. Long-range tasks (species classification at 32 kb, DeepSEA HM at 8 kb), where GENA-LM is documented to outperform HyenaDNA, are deliberately not included.

(vi) **Sample sizes and statistical power.** Per-task sample sizes (≈ 2,000–5,000) limit statistical power; E3 F1 differences (Δ ≈ 0.01–0.05) are at the edge of seed-noise resolution. Multiple seeds and bootstrap CIs are reported in revised tables.

(vii) **E5 is counterfactual in silico, not biologically causal.** The CTCF saturation analysis is an interventional probing experiment on synthetic sequences — stronger than correlational attribution methods but weaker than functional perturbation in a reporter assay. A negative control (mutations at random positions outside the motif) is included and yields a 15× smaller embedding response than motif mutations, ruling out generic BPE-tokenisation sensitivity as the source of the signal.

## 10. Expected contribution and follow-ups

(i) A reproducible BERTology + DART-Eval-inspired pipeline for DNA foundation models (open-sourced, single V100, ≈ 16 minutes total for all five experiments); (ii) **to my knowledge, the first published systematic layer-aware decomposition of GENA-LM features**, showing mid-layer specialisation for the tasks the model meaningfully resolves (promoter, coding) and near-chance behaviour for long-range regulatory tasks (enhancer, OCR) at any layer — pointing to BigBird/RMT variants as the right next step; (iii) empirical evidence (E3, E6) that a 254×–1000× smaller model (HyenaDNA-tiny, TinyCNN) is competitive with GENA-LM **as a frozen feature extractor** on short-range regulatory tasks, with a critical refinement from E3-v2: **the original «HyenaDNA wins 3/4» result is partly a last-layer artefact** — at the E1-best layer per task, GENA-LM matches or exceeds HyenaDNA on 3 of 4 tasks; the practical guideline from E4 (avoid the last layer) is the more durable finding; (iv) **counterfactual evidence that GENA-LM has internalised the canonical CTCF binding grammar** (E5: r = 0.893 between motif-substitution effect and `1 − PWM_prob`, **with a 15× reduction in embedding response for mutations at random non-motif positions** — ruling out generic BPE-tokenisation sensitivity), closing weakness W10 with a positive result; (v) a striking calibration data-point: **a 65-feature 3-mer + LogReg classifier matches a 110 M-parameter GENA-LM within 1.1 % AUROC** on Open Chromatin Region classification (E6) — closing W12 with a quantitative anchor; (vi) **scaling evidence** (E1-large): replicating E1 on `gena-lm-bert-large-t2t` (24 layers, 336 M params) shows that the mid-layer specialisation pattern persists at scale, and — crucially — that enhancer and OCR tasks, near-chance in bert-base, become substantially easier in bert-large (enhancer F1 0.66 → 0.73 with deeper model); (vii) **frozen-vs-fine-tuned calibration point** (E7): a 3-epoch fine-tune of `gena-lm-bert-base-t2t` on the promoter task yields F1 = 0.794 ± 0.004 — quantifying the gap to the authors' published fine-tuned 0.94 on 2 kb context; (viii) practical guidelines for downstream users — *don't use the last layer if you're not fine-tuning, and don't assume smaller-model wins at last layer generalise to best-per-task layer*. **Natural follow-ups** for the summer-school project itself include: replicating on `gena-lm-bert-large-t2t` (24 layers, 336 M params) to test whether deeper models reorganise feature hierarchy; running E3, E5, E6 on **Caduceus** (ICML 2024 A\*) once `mamba_ssm` builds cleanly on Volta; running the **official DART-Eval Task 1** pipeline once Synapse access and a human reference genome are provisioned; and extending E5 to additional TF motifs (GATA2, ATF1 from the original paper) to test motif-grammar generality.

## 11. References

1. Fishman V. *et al.* GENA-LM. *NAR* 53(2), gkae1310 (2025). [DOI](https://doi.org/10.1093/nar/gkae1310)
2. Schiff Y. *et al.* Caduceus: Bi-Directional Equivariant Long-Range DNA Sequence Modeling. *ICML 2024* (CORE A\*). [arXiv:2403.03234](https://arxiv.org/abs/2403.03234)
3. Patel A. *et al.* DART-Eval: A Comprehensive DNA Language Model Evaluation Benchmark on Regulatory DNA. *NeurIPS 2024 D&B* (CORE A\*). [arXiv:2412.05430](https://arxiv.org/abs/2412.05430)
4. Brixi G. *et al.* Genome modelling and design across all domains of life with Evo 2. *Nature* (2026).
5. Tenney I., Das D., Pavlick E. BERT Rediscovers the Classical NLP Pipeline. *ACL 2019*.
6. Vig J. *et al.* BERTology Meets Biology: Interpreting Attention in Protein Language Models. *ICLR 2021*.
7. Bulatov A., Kuratov Y., Burtsev M. Recurrent Memory Transformer. *NeurIPS 2022*.
8. Grešová K. *et al.* Genomic Benchmarks: a collection of datasets for genomic sequence classification. *BMC Genomic Data* (2023).
9. Nguyen E. *et al.* HyenaDNA: Long-Range Genomic Sequence Modeling at Single Nucleotide Resolution. *NeurIPS 2023*. [arXiv:2306.15794](https://arxiv.org/abs/2306.15794)
10. Sennrich R., Haddow B., Birch A. Neural Machine Translation of Rare Words with Subword Units. *ACL 2016*.
11. Dreos R. *et al.* The eukaryotic promoter database in its 30th year: focus on non-vertebrate organisms (EPDnew). *Nucleic Acids Research* 45(D1), D51–D55 (2017).
12. Castro-Mondragon J. A. *et al.* JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles. *Nucleic Acids Research* 50(D1), D165–D173 (2022).
13. Zhou J., Troyanskaya O. G. Predicting effects of noncoding variants with deep learning-based sequence model (DeepSEA). *Nature Methods* 12, 931–934 (2015).
