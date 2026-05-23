# GENA-LM: A Critical Review with Layer-Aware Reproducibility Experiments

**Research Proposal — «Лето с AIRI 2026»**
**Применitель:** *Korshunov E., bioinformatics student*
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

## 5. Proposed improvement

**Layer-aware feature analysis of GENA-LM + head-to-head with a modern alternative.** Three complementary, frozen-model experiments:

- **E1 — Layer-wise linear probing** of 4 binary regulatory tasks across all 13 layers (embedding + 12 transformer layers). This is the BERTology methodology (Tenney et al. 2019; Vig et al. 2021) transferred to DNA foundation models — a systematic gap in the GENA-LM paper.
- **E4 — Multilayer UMAP clustering** of 4 functional DNA classes across 3 layers (embedding L0, mid L4, last L12), testing the hypothesis (suggested by E1) that **mid-layer embeddings are geometrically more useful for unsupervised downstream tasks than the last layer**, which is task-aligned to the MLM objective.
- **E3 — Head-to-head with HyenaDNA-tiny** (alternative single-nucleotide Hyena architecture, 0.4 M params). Same 4 tasks, same probing protocol. Directly addresses W7 (no comparison with non-BPE alternatives) and tests the central claim of DART-Eval (NeurIPS 2024 A\*): that *current DNA LMs do not offer compelling gains over lighter alternatives.* DNABERT-2 was attempted as a third baseline but failed to load (`einops` dependency conflict with `transformers 4.36` — reported as a known issue).

All three directly extend the paper's Figure 4 (species-level layer analysis only) to **functional categories** within the human genome and to **architectural alternatives** — strictly within-scope and biologically motivated.

## 6. Experimental setup

| | Details |
|---|---|
| Model | `AIRI-Institute/gena-lm-bert-base-t2t` (110 M, 12 layers, 768-dim, frozen) |
| Hardware | 1× Tesla V100S 32 GB (Aldan3 cluster, Pirogov RNRMU) |
| Datasets | Genomic Benchmarks (HuggingFace): `human_nontata_promoters`, `human_enhancers_cohn`, `human_ocr_ensembl`, `demo_coding_vs_intergenomic_seqs` |
| Pooling | Mean-pooled hidden states, masked padding |
| E1 probe | Logistic regression (sklearn, max_iter 2000, C=1), StandardScaler, 70/30 stratified split, 3 random seeds |
| E4 metrics | Silhouette (cosine), NMI, ARI with KMeans (k=4) |
| Total wall-clock | **≈ 10 minutes** for both experiments combined |

## 7. Results

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

| Model | Architecture | Params | promoter F1 | enhancer F1 | OCR F1 | coding F1 |
|---|---|---:|---:|---:|---:|---:|
| **GENA-LM-base-t2t** | BPE transformer (last layer) | **110.7 M** | 0.777 | 0.627 | 0.565 | **0.887** |
| **HyenaDNA-tiny** | Hyena, single-nucleotide | **0.4 M** | **0.790** | **0.679** | **0.608** | 0.862 |

**HyenaDNA wins 3/4 tasks with 277× fewer parameters** — only on coding (where BPE token / k-mer composition is highly informative) does GENA-LM lead, and by a small margin. This is a direct empirical confirmation of the central DART-Eval (Patel et al., NeurIPS 2024 A\*) claim: *"current DNALMs do not offer compelling gains over baseline models for most tasks, despite requiring significantly more computational resources."* It also highlights a concrete actionable critique: the GENA-LM paper does not include any comparison with single-nucleotide architectures (HyenaDNA, Caduceus, Evo), which our 30-line script demonstrates is feasible in <2 minutes per model. *(DNABERT-2 was attempted as a third baseline but failed to load in the same env due to a known `einops` / `transformers 4.36` interaction.)*

## 8. Expected contribution and follow-ups

(i) A reproducible BERTology pipeline for DNA foundation models (open-sourced, single V100, ≈ 12 minutes total for all three experiments); (ii) the first systematic layer-aware decomposition of GENA-LM features, showing mid-layer specialisation; (iii) empirical evidence (E3) that a 277× smaller alternative-architecture model (HyenaDNA-tiny) matches or beats GENA-LM on 3 of 4 regulatory tasks — confirming the DART-Eval critique; (iv) a practical guideline for downstream users — *don't use the last layer if you're not fine-tuning, and don't assume larger is better*. **Natural follow-ups** for the summer-school project itself include: replicating on `gena-lm-bert-large-t2t` (24 layers, 336 M params) to test whether deeper models reorganise feature hierarchy; running E1 + E3 on **Caduceus** (ICML 2024 A\*) once `mamba_ssm` builds cleanly on Volta; and adding a **DART-Eval Task 2 / 4** head-to-head, which would close W12 entirely.

## 9. References

1. Fishman V. *et al.* GENA-LM. *NAR* 53(2), gkae1310 (2025). [DOI](https://doi.org/10.1093/nar/gkae1310)
2. Schiff Y. *et al.* Caduceus: Bi-Directional Equivariant Long-Range DNA Sequence Modeling. *ICML 2024* (CORE A\*). [arXiv:2403.03234](https://arxiv.org/abs/2403.03234)
3. Patel A. *et al.* DART-Eval: A Comprehensive DNA Language Model Evaluation Benchmark on Regulatory DNA. *NeurIPS 2024 D&B* (CORE A\*). [arXiv:2412.05430](https://arxiv.org/abs/2412.05430)
4. Brixi G. *et al.* Genome modelling and design across all domains of life with Evo 2. *Nature* (2026).
5. Tenney I., Das D., Pavlick E. BERT Rediscovers the Classical NLP Pipeline. *ACL 2019*.
6. Vig J. *et al.* BERTology Meets Biology: Interpreting Attention in Protein Language Models. *ICLR 2021*.
7. Bulatov A., Kuratov Y., Burtsev M. Recurrent Memory Transformer. *NeurIPS 2022*.
8. Grešová K. *et al.* Genomic Benchmarks: a collection of datasets for genomic sequence classification. *BMC Genomic Data* (2023).
