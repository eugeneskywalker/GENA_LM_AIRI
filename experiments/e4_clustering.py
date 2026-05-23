"""E4: Embedding clustering — UMAP of GENA-LM hidden states for functional DNA classes.

Uses Genomic Benchmarks datasets (ready-to-use, no manual data prep):
- human_nontata_promoters (positive class = promoter, negative = random-ish)
- human_enhancers_cohn (enhancers)
- human_ocr_ensembl (open chromatin = regulatory)
"""
import os
os.environ["HF_HOME"] = "/home/ekorshunov/hf-cache"

import json
import argparse
from pathlib import Path

import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from datasets import load_dataset
from transformers import AutoTokenizer, AutoModel, AutoConfig
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, normalized_mutual_info_score, adjusted_rand_score
import umap


# (dataset_name, split, label_value_to_keep)
DATASETS = {
    "promoter":   ("katarinagresova/Genomic_Benchmarks_human_nontata_promoters", "test", 1),
    "enhancer":   ("katarinagresova/Genomic_Benchmarks_human_enhancers_cohn",     "test", 1),
    "ocr":        ("katarinagresova/Genomic_Benchmarks_human_ocr_ensembl",        "test", 1),
    "negative":   ("katarinagresova/Genomic_Benchmarks_human_nontata_promoters", "test", 0),
}


def load_model(model_name, device):
    """Load tokenizer + model with hidden states enabled."""
    tok = AutoTokenizer.from_pretrained(model_name)
    config = AutoConfig.from_pretrained(model_name, trust_remote_code=True)
    config.output_hidden_states = True
    model = AutoModel.from_pretrained(model_name, trust_remote_code=True, config=config).to(device).eval()
    print(f"  class: {type(model).__name__}, layers: {config.num_hidden_layers}, hidden: {config.hidden_size}", flush=True)
    return tok, model


def get_embeddings(seqs, tok, model, device, batch_size=16, max_len=512, layer_idx=-1):
    """Mean-pooled embeddings from given layer of hidden_states (-1 = last layer)."""
    embs = []
    n = len(seqs)
    for i in range(0, n, batch_size):
        batch = seqs[i:i+batch_size]
        inp = tok(batch, return_tensors="pt", padding=True, truncation=True, max_length=max_len).to(device)
        with torch.no_grad():
            out = model(**inp)
        hidden = out.hidden_states[layer_idx]  # (B, T, H)
        mask = inp["attention_mask"].unsqueeze(-1).float()
        emb = (hidden * mask).sum(1) / mask.sum(1).clamp(min=1)
        embs.append(emb.cpu().float().numpy())
        if (i // batch_size) % 5 == 0:
            print(f"    batch {i//batch_size+1}/{(n+batch_size-1)//batch_size}", flush=True)
    return np.concatenate(embs, axis=0)


def main(n_per_class=500, model_name="AIRI-Institute/gena-lm-bert-base-t2t",
         results_dir="/home/ekorshunov/gena-lm-airi/results", seed=42, batch_size=16):
    print(f"=== E4 clustering | n_per_class={n_per_class} | model={model_name} ===", flush=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}, gpu: {torch.cuda.get_device_name(0) if device=='cuda' else '-'}", flush=True)

    print("Loading model...", flush=True)
    tok, model = load_model(model_name, device)

    rng = np.random.default_rng(seed)
    X_list, y_list, classes = [], [], []
    for cls_name, (ds_name, split, target_label) in DATASETS.items():
        print(f"\n--- {cls_name}: {ds_name} (split={split}, label={target_label}) ---", flush=True)
        ds = load_dataset(ds_name, split=split)
        labels = np.array(ds["label"])
        seqs_all = np.array(ds["seq"])[labels == target_label]
        if len(seqs_all) < n_per_class:
            print(f"  WARN: only {len(seqs_all)} sequences available", flush=True)
        idx = rng.choice(len(seqs_all), size=min(n_per_class, len(seqs_all)), replace=False)
        seqs = list(seqs_all[idx])
        avg_len = int(np.mean([len(s) for s in seqs]))
        print(f"  sampled {len(seqs)}, avg len={avg_len} bp", flush=True)

        emb = get_embeddings(seqs, tok, model, device, batch_size=batch_size)
        X_list.append(emb)
        y_list.extend([cls_name] * len(seqs))
        classes.append(cls_name)

    X = np.concatenate(X_list, axis=0)
    y = np.array(y_list)
    y_int = np.array([classes.index(lab) for lab in y])
    print(f"\nEmbedding matrix: {X.shape}", flush=True)

    print("UMAP...", flush=True)
    reducer = umap.UMAP(n_components=2, random_state=seed, n_neighbors=15, min_dist=0.1, metric="cosine")
    X_2d = reducer.fit_transform(X)

    print("Metrics...", flush=True)
    kmeans = KMeans(n_clusters=len(classes), random_state=seed, n_init=10).fit(X)
    sil = float(silhouette_score(X, y_int, metric="cosine"))
    nmi = float(normalized_mutual_info_score(y_int, kmeans.labels_))
    ari = float(adjusted_rand_score(y_int, kmeans.labels_))
    print(f"  silhouette (cosine, true labels): {sil:.3f}", flush=True)
    print(f"  NMI (kmeans vs true):             {nmi:.3f}", flush=True)
    print(f"  ARI (kmeans vs true):             {ari:.3f}", flush=True)

    results_dir = Path(results_dir)
    (results_dir / "figures").mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 7))
    palette = sns.color_palette("Set1", n_colors=len(classes))
    for i, lab in enumerate(classes):
        m = y == lab
        ax.scatter(X_2d[m, 0], X_2d[m, 1], c=[palette[i]], label=f"{lab} (n={m.sum()})",
                   s=20, alpha=0.7, edgecolors="none")
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    ax.set_title(f"GENA-LM-base-t2t embeddings (mean-pooled, last layer)\n"
                 f"silhouette={sil:.3f} | NMI={nmi:.3f} | ARI={ari:.3f}")
    ax.legend(loc="best")
    plt.tight_layout()
    fig_path = results_dir / "figures" / "e4_umap.png"
    plt.savefig(fig_path, dpi=200, bbox_inches="tight")
    print(f"\nSaved figure: {fig_path}", flush=True)

    np.savez(results_dir / "e4_embeddings.npz", X=X, y=y, X_2d=X_2d, classes=np.array(classes))

    metrics = {
        "experiment": "E4 embedding clustering",
        "model": model_name,
        "n_per_class": n_per_class,
        "classes": classes,
        "metrics": {"silhouette_cosine": sil, "nmi": nmi, "ari": ari},
        "embedding_dim": int(X.shape[1]),
        "datasets": {k: v[0] for k, v in DATASETS.items()},
    }
    with open(results_dir / "e4_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"Saved metrics + embeddings.", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--bs", type=int, default=16)
    ap.add_argument("--model", type=str, default="AIRI-Institute/gena-lm-bert-base-t2t")
    args = ap.parse_args()
    main(n_per_class=args.n, batch_size=args.bs, model_name=args.model)
