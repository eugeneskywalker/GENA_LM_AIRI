"""E4 v2: Embedding clustering across multiple layers.

Compare GENA-LM embeddings from 3 layers (embedding, mid, last) for unsupervised
clustering of 4 functional DNA classes. Tests the hypothesis (suggested by E1
results) that middle layers contain more semantically informative representations
than the last layer.
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


DATASETS = {
    "promoter":   ("katarinagresova/Genomic_Benchmarks_human_nontata_promoters", "test", 1),
    "enhancer":   ("katarinagresova/Genomic_Benchmarks_human_enhancers_cohn",     "test", 1),
    "ocr":        ("katarinagresova/Genomic_Benchmarks_human_ocr_ensembl",        "test", 1),
    "negative":   ("katarinagresova/Genomic_Benchmarks_human_nontata_promoters", "test", 0),
}

LAYERS = [0, 4, 12]  # embedding, mid (peak from E1), last


def load_model(model_name, device):
    tok = AutoTokenizer.from_pretrained(model_name)
    config = AutoConfig.from_pretrained(model_name, trust_remote_code=True)
    config.output_hidden_states = True
    model = AutoModel.from_pretrained(model_name, trust_remote_code=True, config=config).to(device).eval()
    return tok, model, config


def extract_layer_embeddings(seqs, tok, model, device, layers, hidden_size, batch_size=16, max_len=512):
    """Return dict: layer_idx -> (N, H) mean-pooled embeddings."""
    n = len(seqs)
    out = {l: np.zeros((n, hidden_size), dtype=np.float32) for l in layers}
    for i in range(0, n, batch_size):
        batch = seqs[i:i+batch_size]
        inp = tok(batch, return_tensors="pt", padding=True, truncation=True, max_length=max_len).to(device)
        with torch.no_grad():
            res = model(**inp)
        mask = inp["attention_mask"].unsqueeze(-1).float()
        denom = mask.sum(1).clamp(min=1)
        for l in layers:
            h = res.hidden_states[l]
            emb = (h * mask).sum(1) / denom
            out[l][i:i+len(batch)] = emb.cpu().float().numpy()
        if (i // batch_size) % 5 == 0:
            print(f"    batch {i//batch_size+1}/{(n+batch_size-1)//batch_size}", flush=True)
    return out


def cluster_metrics(X, y_int, n_clusters, seed=42):
    kmeans = KMeans(n_clusters=n_clusters, random_state=seed, n_init=10).fit(X)
    return {
        "silhouette": float(silhouette_score(X, y_int, metric="cosine")),
        "nmi": float(normalized_mutual_info_score(y_int, kmeans.labels_)),
        "ari": float(adjusted_rand_score(y_int, kmeans.labels_)),
    }


def main(n_per_class=500, model_name="AIRI-Institute/gena-lm-bert-base-t2t",
         results_dir="/home/ekorshunov/gena-lm-airi/results", seed=42, batch_size=16):
    print(f"=== E4 v2 multilayer clustering | n_per_class={n_per_class} ===", flush=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}", flush=True)

    tok, model, config = load_model(model_name, device)
    print(f"  layers: {config.num_hidden_layers}, hidden: {config.hidden_size}", flush=True)

    rng = np.random.default_rng(seed)
    all_X = {l: [] for l in LAYERS}
    y_list, classes = [], []
    for cls_name, (ds_name, split, target_label) in DATASETS.items():
        print(f"\n--- {cls_name}: {ds_name} (label={target_label}) ---", flush=True)
        ds = load_dataset(ds_name, split=split)
        labels = np.array(ds["label"])
        seqs_all = np.array(ds["seq"])[labels == target_label]
        idx = rng.choice(len(seqs_all), size=min(n_per_class, len(seqs_all)), replace=False)
        seqs = list(seqs_all[idx])
        print(f"  sampled {len(seqs)}, avg len={int(np.mean([len(s) for s in seqs]))} bp", flush=True)

        emb_dict = extract_layer_embeddings(seqs, tok, model, device, LAYERS, config.hidden_size, batch_size)
        for l in LAYERS:
            all_X[l].append(emb_dict[l])
        y_list.extend([cls_name] * len(seqs))
        classes.append(cls_name)

    y = np.array(y_list)
    y_int = np.array([classes.index(lab) for lab in y])

    # Build figure
    results_dir = Path(results_dir)
    (results_dir / "figures").mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    palette = sns.color_palette("Set1", n_colors=len(classes))

    all_metrics = {}
    for ax, layer_idx in zip(axes, LAYERS):
        X = np.concatenate(all_X[layer_idx], axis=0)
        m = cluster_metrics(X, y_int, n_clusters=len(classes), seed=seed)
        all_metrics[f"layer_{layer_idx}"] = m
        print(f"\nlayer {layer_idx}: silhouette={m['silhouette']:.3f}  NMI={m['nmi']:.3f}  ARI={m['ari']:.3f}", flush=True)

        reducer = umap.UMAP(n_components=2, random_state=seed, n_neighbors=15, min_dist=0.1, metric="cosine")
        X_2d = reducer.fit_transform(X)
        for i, lab in enumerate(classes):
            mask = y == lab
            ax.scatter(X_2d[mask, 0], X_2d[mask, 1], c=[palette[i]], label=f"{lab} (n={mask.sum()})",
                       s=18, alpha=0.7, edgecolors="none")
        layer_label = {0: "embedding (L0)", 4: "mid (L4, E1 peak)", 12: "last (L12)"}[layer_idx]
        ax.set_xlabel("UMAP-1")
        ax.set_ylabel("UMAP-2")
        ax.set_title(f"Layer {layer_idx} - {layer_label}\nsil={m['silhouette']:.3f} | NMI={m['nmi']:.3f} | ARI={m['ari']:.3f}")
        if layer_idx == LAYERS[0]:
            ax.legend(loc="best", fontsize=9)

    plt.suptitle(f"GENA-LM-bert-base-t2t: embedding clustering across layers ({n_per_class}/class)")
    plt.tight_layout()
    fig_path = results_dir / "figures" / "e4_multilayer_umap.png"
    plt.savefig(fig_path, dpi=200, bbox_inches="tight")
    print(f"\nSaved figure: {fig_path}", flush=True)

    out = {
        "experiment": "E4 v2 multilayer clustering",
        "model": model_name,
        "n_per_class": n_per_class,
        "classes": classes,
        "layers_compared": LAYERS,
        "metrics_per_layer": all_metrics,
        "datasets": {k: v[0] for k, v in DATASETS.items()},
    }
    with open(results_dir / "e4_multilayer_metrics.json", "w") as f:
        json.dump(out, f, indent=2)
    print("Saved metrics.", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--bs", type=int, default=16)
    ap.add_argument("--model", type=str, default="AIRI-Institute/gena-lm-bert-base-t2t")
    args = ap.parse_args()
    main(n_per_class=args.n, batch_size=args.bs, model_name=args.model)
