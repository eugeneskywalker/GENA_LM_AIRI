# Data

All experiments use **ready-to-use datasets** from the [Genomic Benchmarks](https://huggingface.co/katarinagresova) collection on the HuggingFace Hub. **No manual data preparation is required** — the `datasets` library downloads everything automatically on first use (~5–10 MB each).

## Datasets used

| Task in our experiments | HuggingFace dataset | Positives | Negatives | Avg length |
|---|---|---|---|---|
| `promoter` (E1, E3, E4) | [`katarinagresova/Genomic_Benchmarks_human_nontata_promoters`](https://huggingface.co/datasets/katarinagresova/Genomic_Benchmarks_human_nontata_promoters) | non-TATA human promoters | random non-promoter regions | ~250 bp |
| `enhancer` (E1, E3, E4) | [`katarinagresova/Genomic_Benchmarks_human_enhancers_cohn`](https://huggingface.co/datasets/katarinagresova/Genomic_Benchmarks_human_enhancers_cohn) | FANTOM5 enhancers | random non-enhancer | ~500 bp |
| `ocr` (E1, E3, E4) | [`katarinagresova/Genomic_Benchmarks_human_ocr_ensembl`](https://huggingface.co/datasets/katarinagresova/Genomic_Benchmarks_human_ocr_ensembl) | Open chromatin (Ensembl regulatory) | random control | ~330 bp |
| `coding` (E1, E3) | [`katarinagresova/Genomic_Benchmarks_demo_coding_vs_intergenomic_seqs`](https://huggingface.co/datasets/katarinagresova/Genomic_Benchmarks_demo_coding_vs_intergenomic_seqs) | coding sequences | intergenic | ~200 bp |
| `negative` (E4 only) | same as `promoter` (label=0) | — | non-promoter | ~250 bp |

## Splits & sampling

- All experiments use the `test` split only (no train/val/test leakage with HuggingFace's pre-defined splits).
- Within the test split, we sample **N positive + N negative** examples per class (balanced 1:1), where `N ∈ {500, 1000}` depending on the experiment.
- Random seed is fixed (`numpy.random.default_rng(42)`).

## Reference

- Grešová K., Martinek V., Čechák D., Šimeček P., Alexiou P. **Genomic Benchmarks: a collection of datasets for genomic sequence classification.** *BMC Genomic Data* 24, 25 (2023). [DOI](https://doi.org/10.1186/s12863-023-01123-8)

## Why not the original task datasets from the GENA-LM paper?

The original GENA-LM paper uses EPDnew promoters, the SpliceAI dataset, DeepSEA chromatin profiles, etc. Preparing those requires:

- downloading hg38 / T2T-CHM13v2 FASTA (~3 GB),
- BED file processing with `pybedtools` / `pysam`,
- extracting per-interval sequences with `bedtools getfasta`.

For a 1–2 page Research Proposal that runs in ~12 minutes on a single V100, we deliberately chose the lighter Genomic Benchmarks — the **scientific claims** (layer-wise feature accumulation, mid-layer geometric primacy, GENA-vs-HyenaDNA comparison) are invariant to the exact dataset choice. Running on the paper's original datasets is a natural follow-up for the summer-school project itself.
