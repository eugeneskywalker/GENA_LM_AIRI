"""E5: In silico saturation mutagenesis on the CTCF motif (embedding-distance version).

Tests whether GENA-LM embeddings are sensitive to mutations in the canonical
CTCF binding motif by measuring how much the mean-pooled hidden state of the
motif region changes when each position is substituted.

WHY EMBEDDING-DISTANCE INSTEAD OF MLM-LIKELIHOOD:
The GENA-LM custom architecture uses pre/post-attention LayerNorm that does
not map cleanly to the standard `BertForMaskedLM` MLM head weights —
`AutoModelForMaskedLM.from_pretrained(..., trust_remote_code=True)` silently
initialises some LayerNorm weights randomly, breaking MLM scoring. Loading
through `AutoModel` keeps the custom architecture intact and gives us valid
hidden states. We use the distance between motif-region embeddings of
original vs mutant sequences as the importance signal — this is the standard
interpretability approach used in ESM-style protein language models.

PROTOCOL:
1. Load JASPAR MA0139.1 CTCF PWM (19 positions x 4 nucleotides).
2. Generate 100 synthetic 400-bp sequences with a motif sampled from the
   PWM inserted in the middle.
3. For each sequence and each (position, alternative_nucleotide):
     - Substitute position p with alt_nuc.
     - Forward pass; mean-pool last hidden state over tokens covering motif.
     - effect[i, p, n] = 1 - cosine_similarity(emb_orig, emb_mut)
                       (large = mutation moved embedding away = position important)
4. Average over 100 sequences -> saturation map (19, 4).
5. Compare with PWM information content per position (Pearson r).
"""
import os
os.environ["HF_HOME"] = "/home/ekorshunov/hf-cache"

import json
import argparse
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from transformers import AutoTokenizer, AutoModel, AutoConfig


# JASPAR MA0139.1 (CTCF, Homo sapiens) — counts of A/C/G/T at each of 19 positions.
CTCF_COUNTS = np.array([
    [87,  167, 281, 56],
    [167, 46,  18,  166],
    [97,  18,  281, 13],
    [0,   1,   388, 6],
    [47,  4,   228, 116],
    [1,   389, 4,   1],
    [0,   100, 4,   291],
    [27,  295, 31,  41],
    [1,   294, 7,   91],
    [28,  6,   357, 2],
    [125, 30,  4,   245],
    [136, 100, 75,  35],
    [89,  133, 154, 102],
    [24,  336, 16,  19],
    [33,  3,   368, 0],
    [0,   361, 26,  18],
    [72,  21,  297, 25],
    [29,  10,  329, 51],
    [85,  72,  264, 0],
], dtype=np.float64)
CTCF_PWM = (CTCF_COUNTS + 1) / (CTCF_COUNTS.sum(axis=1, keepdims=True) + 4)
NUCLEOTIDES = ["A", "C", "G", "T"]
MOTIF_LEN = CTCF_PWM.shape[0]


def sample_from_pwm(pwm, rng):
    return "".join(NUCLEOTIDES[rng.choice(4, p=pwm[i])] for i in range(len(pwm)))


def random_dna(length, rng):
    return "".join(NUCLEOTIDES[i] for i in rng.choice(4, size=length))


def info_content(pwm, bg=0.25):
    return np.array([np.sum(p * np.log2(p / bg + 1e-12)) for p in pwm])


def batch_motif_embeddings(seqs, motif_start, motif_end, tokenizer, model, device,
                           batch_size=32, max_len=512):
    """For each seq, return mean-pooled last hidden state over tokens covering [motif_start, motif_end)."""
    out_emb = []
    for batch_start in range(0, len(seqs), batch_size):
        batch = seqs[batch_start:batch_start + batch_size]
        enc = tokenizer(batch, return_tensors="pt", padding=True, truncation=True,
                        max_length=max_len, return_offsets_mapping=True)
        input_ids = enc["input_ids"].to(device)
        attn = enc["attention_mask"].to(device)
        offsets = enc["offset_mapping"]  # (B, T, 2) tensor

        with torch.no_grad():
            outputs = model(input_ids=input_ids, attention_mask=attn)

        # last hidden state
        if hasattr(outputs, "hidden_states") and outputs.hidden_states is not None:
            h = outputs.hidden_states[-1]
        elif hasattr(outputs, "last_hidden_state"):
            h = outputs.last_hidden_state
        else:
            h = outputs[0]
        h = h.cpu().float()  # (B, T, H)
        attn_cpu = attn.cpu()

        for i in range(len(batch)):
            offsets_i = offsets[i].tolist()
            # Find tokens overlapping motif region (skip pad/CLS/SEP which have (0,0))
            mask = torch.zeros(h.shape[1])
            for j, (s, e) in enumerate(offsets_i):
                if attn_cpu[i, j].item() == 0:
                    continue
                if e == 0 and s == 0:  # special tokens
                    continue
                if not (e <= motif_start or s >= motif_end):
                    mask[j] = 1.0
            if mask.sum() == 0:
                out_emb.append(torch.zeros(h.shape[-1]).numpy())
                continue
            emb = (h[i] * mask.unsqueeze(-1)).sum(0) / mask.sum()
            out_emb.append(emb.numpy())

    return np.stack(out_emb, axis=0)  # (N, H)


def main(n_seqs=100, flank_len=200, model_name="AIRI-Institute/gena-lm-bert-base-t2t",
         results_dir="/home/ekorshunov/gena-lm-airi/results", seed=42, batch_size=32):
    print(f"=== E5 CTCF saturation mutagenesis (embedding-distance) ===", flush=True)
    print(f"n_seqs={n_seqs}, flank={flank_len}, batch={batch_size}", flush=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}", flush=True)

    print("\nLoading model (AutoModel + custom architecture via trust_remote_code)...", flush=True)
    tok = AutoTokenizer.from_pretrained(model_name)
    config = AutoConfig.from_pretrained(model_name, trust_remote_code=True)
    config.output_hidden_states = True   # CRUCIAL: without this, custom model returns MLM logits as outputs[0] (32000-dim, not 768)
    model = AutoModel.from_pretrained(model_name, trust_remote_code=True, config=config).to(device).eval()
    print(f"  {type(model).__name__}, {sum(p.numel() for p in model.parameters())/1e6:.1f}M params", flush=True)
    print(f"  hidden_size: {config.hidden_size}, num_layers: {config.num_hidden_layers}", flush=True)

    ic = info_content(CTCF_PWM)
    print(f"\nCTCF PWM IC per position: {ic.round(2)}", flush=True)
    print(f"Most informative positions (IC>1.0): {[i+1 for i, v in enumerate(ic) if v > 1.0]}", flush=True)

    rng = np.random.default_rng(seed)
    motif_start = flank_len
    motif_end = flank_len + MOTIF_LEN

    # Generate sequences
    print(f"\nGenerating {n_seqs} synthetic sequences with motif at [{motif_start}, {motif_end})...", flush=True)
    seqs, motifs = [], []
    for _ in range(n_seqs):
        seq = random_dna(flank_len, rng) + sample_from_pwm(CTCF_PWM, rng) + random_dna(flank_len, rng)
        seqs.append(seq)
        motifs.append(seq[motif_start:motif_end])
    print(f"  example: ...{seqs[0][motif_start-3:motif_end+3]}...", flush=True)

    # Build all mutant sequences
    print(f"\nBuilding mutant sequences...", flush=True)
    all_seqs = []
    indexing = []
    for seq_idx, seq in enumerate(seqs):
        all_seqs.append(seq)
        indexing.append((seq_idx, -1, -1))
        for p in range(MOTIF_LEN):
            pos_in_seq = motif_start + p
            orig_nuc = seq[pos_in_seq]
            for alt_idx, alt in enumerate(NUCLEOTIDES):
                if alt == orig_nuc:
                    continue
                mutant = seq[:pos_in_seq] + alt + seq[pos_in_seq + 1:]
                all_seqs.append(mutant)
                indexing.append((seq_idx, p, alt_idx))
    print(f"  total: {len(all_seqs)} = {n_seqs} * (1 + {MOTIF_LEN}*3)", flush=True)

    # Compute motif-region embeddings
    print(f"\nComputing motif-region embeddings...", flush=True)
    embeddings = batch_motif_embeddings(all_seqs, motif_start, motif_end, tok, model, device, batch_size=batch_size)
    print(f"  embedding tensor: {embeddings.shape}", flush=True)

    # Cosine-distance saturation: 1 - cos(emb_orig, emb_mut)
    # First: index baselines by seq_idx
    baseline_emb = {}
    for idx, (seq_idx, p, alt_idx) in enumerate(indexing):
        if p == -1:
            baseline_emb[seq_idx] = embeddings[idx]

    # Compute effect
    saturation = np.zeros((n_seqs, MOTIF_LEN, 4), dtype=np.float64)
    for idx, (seq_idx, p, alt_idx) in enumerate(indexing):
        if p == -1:
            continue
        e_orig = baseline_emb[seq_idx]
        e_mut = embeddings[idx]
        cos_sim = float(np.dot(e_orig, e_mut) / (np.linalg.norm(e_orig) * np.linalg.norm(e_mut) + 1e-12))
        saturation[seq_idx, p, alt_idx] = 1.0 - cos_sim  # in [0, 2], higher = more changed

    mean_effect = saturation.mean(axis=0)  # (19, 4)
    pos_mean_effect = mean_effect.sum(axis=1) / 3.0  # average over 3 substitutions

    pearson_pos = float(np.corrcoef(pos_mean_effect, ic)[0, 1])
    print(f"\nPearson r (position-wise: mean substitution effect vs PWM IC): {pearson_pos:.3f}", flush=True)

    # Per-element: substitution effect should be larger where (1-PWM) is larger (off-consensus subs)
    inv_pwm = 1.0 - CTCF_PWM
    mask = mean_effect != 0
    pearson_full = float(np.corrcoef(mean_effect[mask], inv_pwm[mask])[0, 1])
    print(f"Pearson r (full off-consensus map vs (1-PWM_prob)): {pearson_full:.3f}", flush=True)

    print(f"\nPer-position mean effect:", flush=True)
    for p in range(MOTIF_LEN):
        print(f"  pos {p+1:2d}: IC={ic[p]:.2f}  effect={pos_mean_effect[p]:.4f}", flush=True)

    # Plot
    results_dir = Path(results_dir)
    (results_dir / "figures").mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(3, 1, figsize=(14, 9), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1, 1]})

    im0 = axes[0].imshow(CTCF_PWM.T, aspect="auto", cmap="Greys", vmin=0, vmax=1)
    axes[0].set_yticks(range(4))
    axes[0].set_yticklabels(NUCLEOTIDES)
    axes[0].set_title("Reference: JASPAR MA0139.1 CTCF PWM (probabilities)")
    plt.colorbar(im0, ax=axes[0], fraction=0.02)

    vmax = float(np.abs(mean_effect).max())
    im1 = axes[1].imshow(mean_effect.T, aspect="auto", cmap="Reds", vmin=0, vmax=vmax)
    axes[1].set_yticks(range(4))
    axes[1].set_yticklabels(NUCLEOTIDES)
    axes[1].set_title(f"GENA-LM saturation mutagenesis: mean 1 - cos(emb_orig, emb_mut)\n"
                      f"(higher = mutation moved motif embedding more = position important)")
    plt.colorbar(im1, ax=axes[1], fraction=0.02)

    axes[2].bar(range(1, MOTIF_LEN + 1), ic, alpha=0.55, color="grey", label="PWM IC (bits)")
    ax2b = axes[2].twinx()
    ax2b.plot(range(1, MOTIF_LEN + 1), pos_mean_effect, marker="o", color="red", label="Mean mutation effect")
    axes[2].set_xlabel("Position in motif (1-19)")
    axes[2].set_ylabel("Information content (bits)")
    ax2b.set_ylabel("Mean cosine distance")
    axes[2].set_title(f"Position-wise: PWM information content vs GENA-LM mutation sensitivity\n"
                      f"Pearson r = {pearson_pos:.3f}")
    axes[2].legend(loc="upper left")
    ax2b.legend(loc="upper right")
    axes[2].set_xticks(range(1, MOTIF_LEN + 1))

    plt.tight_layout()
    fig_path = results_dir / "figures" / "e5_ctcf_saturation.png"
    plt.savefig(fig_path, dpi=200, bbox_inches="tight")
    print(f"\nSaved: {fig_path}", flush=True)

    metrics = {
        "experiment": "E5 CTCF saturation mutagenesis (embedding-distance)",
        "model": model_name,
        "n_sequences": n_seqs,
        "motif_len": int(MOTIF_LEN),
        "flank_len": flank_len,
        "metric": "1 - cosine_similarity(emb_orig, emb_mut) at motif region",
        "pearson_r_pos_IC_vs_effect": pearson_pos,
        "pearson_r_full_map_off_consensus": pearson_full,
        "info_content_per_position": ic.tolist(),
        "mean_effect_per_position": pos_mean_effect.tolist(),
        "saturation_map_19x4_ACGT": mean_effect.tolist(),
        "PWM_19x4_ACGT": CTCF_PWM.tolist(),
    }
    with open(results_dir / "e5_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"Saved: {results_dir / 'e5_metrics.json'}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--flank", type=int, default=200)
    ap.add_argument("--bs", type=int, default=32)
    ap.add_argument("--model", type=str, default="AIRI-Institute/gena-lm-bert-base-t2t")
    args = ap.parse_args()
    main(n_seqs=args.n, flank_len=args.flank, model_name=args.model, batch_size=args.bs)
