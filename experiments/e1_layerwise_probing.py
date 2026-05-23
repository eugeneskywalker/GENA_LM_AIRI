"""E1: Layer-wise linear probing for GENA-LM.

For each of 4 binary regulatory tasks, extract mean-pooled embeddings from
EVERY hidden layer (embedding + 12 transformer layers = 13 layers total) and
train logistic regression on each layer separately. The resulting F1/MCC/AUC
vs layer curves reveal at which depth functional info is encoded.

Tasks (all from Genomic Benchmarks on HF Hub):
1. promoter:  human_nontata_promoters
2. enhancer:  human_enhancers_cohn
3. ocr:       human_ocr_ensembl
4. coding:    demo_coding_vs_intergenomic_seqs
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
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import f1_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import train_test_split


TASKS = {
    "promoter": "katarinagresova/Genomic_Benchmarks_human_nontata_promoters",
    "enhancer": "katarinagresova/Genomic_Benchmarks_human_enhancers_cohn",
    "ocr":      "katarinagresova/Genomic_Benchmarks_human_ocr_ensembl",
    "coding":   "katarinagresova/Genomic_Benchmarks_demo_coding_vs_intergenomic_seqs",
}


def load_model(model_name, device):
    tok = AutoTokenizer.from_pretrained(model_name)
    config = AutoConfig.from_pretrained(model_name, trust_remote_code=True)
    config.output_hidden_states = True
    model = AutoModel.from_pretrained(model_name, trust_remote_code=True, config=config).to(device).eval()
    return tok, model, config


def extract_all_layer_embeddings(seqs, tok, model, device, n_layers, hidden_size,
                                 batch_size=16, max_len=512):
    n = len(seqs)
    out_arr = np.zeros((n, n_layers + 1, hidden_size), dtype=np.float32)
    for i in range(0, n, batch_size):
        batch = seqs[i:i+batch_size]
        inp = tok(batch, return_tensors="pt", padding=True, truncation=True, max_length=max_len).to(device)
        with torch.no_grad():
            out = model(**inp)
        mask = inp["attention_mask"].unsqueeze(-1).float()
        denom = mask.sum(1).clamp(min=1)
        for layer_idx, h in enumerate(out.hidden_states):
            emb = (h * mask).sum(1) / denom
            out_arr[i:i+len(batch), layer_idx] = emb.cpu().float().numpy()
        if (i // batch_size) % 10 == 0:
            print(f"      batch {i//batch_size+1}/{(n+batch_size-1)//batch_size}", flush=True)
    return out_arr


def sample_balanced(ds, n_per_class, rng):
    labels = np.array(ds["label"])
    seqs_all = np.array(ds["seq"])
    seqs, lbs = [], []
    for lbl in (0, 1):
        idx_pool = np.where(labels == lbl)[0]
        n = min(n_per_class, len(idx_pool))
        idx = rng.choice(idx_pool, size=n, replace=False)
        seqs.extend(seqs_all[idx].tolist())
        lbs.extend([lbl] * n)
    return seqs, np.array(lbs)


def probe_layer(X_layer, y, n_seeds=3):
    f1s, mccs, aucs = [], [], []
    for seed in range(n_seeds):
        Xtr, Xte, ytr, yte = train_test_split(X_layer, y, test_size=0.3, random_state=seed, stratify=y)
        clf = make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=2000, C=1.0, random_state=seed),
        )
        clf.fit(Xtr, ytr)
        pred = clf.predict(Xte)
        prob = clf.predict_proba(Xte)[:, 1]
        f1s.append(f1_score(yte, pred))
        mccs.append(matthews_corrcoef(yte, pred))
        aucs.append(roc_auc_score(yte, prob))
    return {
        "f1_mean": float(np.mean(f1s)), "f1_std": float(np.std(f1s)),
        "mcc_mean": float(np.mean(mccs)), "mcc_std": float(np.std(mccs)),
        "auc_mean": float(np.mean(aucs)), "auc_std": float(np.std(aucs)),
    }


def main(n_per_class=500, model_name="AIRI-Institute/gena-lm-bert-base-t2t",
         results_dir="/home/ekorshunov/gena-lm-airi/results", seed=42, batch_size=16):
    print(f"=== E1 layer-wise probing | n_per_class={n_per_class} | model={model_name} ===", flush=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}, gpu: {torch.cuda.get_device_name(0) if device=='cuda' else '-'}", flush=True)

    print("Loading model...", flush=True)
    tok, model, config = load_model(model_name, device)
    n_layers = config.num_hidden_layers
    hidden_size = config.hidden_size
    print(f"  layers: {n_layers}, hidden: {hidden_size}", flush=True)

    rng = np.random.default_rng(seed)
    all_results = {}

    for task_name, ds_name in TASKS.items():
        print(f"\n=== TASK: {task_name} | {ds_name} ===", flush=True)
        ds = load_dataset(ds_name, split="test")
        seqs, y = sample_balanced(ds, n_per_class, rng)
        print(f"  sampled {len(seqs)} balanced, avg len={int(np.mean([len(s) for s in seqs]))} bp", flush=True)

        print(f"  Extracting embeddings from all {n_layers+1} layers...", flush=True)
        X_all = extract_all_layer_embeddings(seqs, tok, model, device, n_layers, hidden_size, batch_size=batch_size)
        print(f"  Embedding tensor: {X_all.shape}", flush=True)

        print(f"  Probing each layer (LogReg, 3 seeds)...", flush=True)
        task_results = []
        for layer_idx in range(n_layers + 1):
            r = probe_layer(X_all[:, layer_idx], y, n_seeds=3)
            r["layer"] = layer_idx
            task_results.append(r)
            print(f"    layer {layer_idx:2d}: F1={r['f1_mean']:.3f}+/-{r['f1_std']:.3f}  MCC={r['mcc_mean']:.3f}  AUC={r['auc_mean']:.3f}", flush=True)
        all_results[task_name] = task_results

    results_dir = Path(results_dir)
    (results_dir / "figures").mkdir(parents=True, exist_ok=True)
    (results_dir / "tables").mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True)
    metric_keys = [("f1_mean", "f1_std", "F1"), ("mcc_mean", "mcc_std", "MCC"), ("auc_mean", "auc_std", "ROC-AUC")]
    palette = sns.color_palette("tab10", n_colors=len(TASKS))

    for ax, (m_mean, m_std, m_label) in zip(axes, metric_keys):
        for i, (task, results) in enumerate(all_results.items()):
            layers = [r["layer"] for r in results]
            means = [r[m_mean] for r in results]
            stds = [r[m_std] for r in results]
            ax.errorbar(layers, means, yerr=stds, label=task, marker="o", color=palette[i], capsize=3)
        ax.set_xlabel("Layer index (0 = embedding, 12 = last)")
        ax.set_ylabel(m_label)
        ax.set_title(f"{m_label} per layer (mean +/- std, 3 seeds)")
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.suptitle("GENA-LM-bert-base-t2t - layer-wise linear probing on regulatory DNA tasks")
    plt.tight_layout()
    fig_path = results_dir / "figures" / "e1_layerwise_probing.png"
    plt.savefig(fig_path, dpi=200, bbox_inches="tight")
    print(f"\nSaved figure: {fig_path}", flush=True)

    metrics_out = {
        "experiment": "E1 layer-wise probing",
        "model": model_name,
        "n_per_class": n_per_class,
        "n_layers": n_layers + 1,
        "hidden_size": hidden_size,
        "tasks": TASKS,
        "results": all_results,
    }
    with open(results_dir / "e1_metrics.json", "w") as f:
        json.dump(metrics_out, f, indent=2)

    import csv
    with open(results_dir / "tables" / "e1_results.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["task", "layer", "f1_mean", "f1_std", "mcc_mean", "auc_mean"])
        for task, results in all_results.items():
            for r in results:
                w.writerow([task, r["layer"], f"{r['f1_mean']:.4f}", f"{r['f1_std']:.4f}", f"{r['mcc_mean']:.4f}", f"{r['auc_mean']:.4f}"])

    print("Saved metrics + CSV.", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--bs", type=int, default=16)
    ap.add_argument("--model", type=str, default="AIRI-Institute/gena-lm-bert-base-t2t")
    args = ap.parse_args()
    main(n_per_class=args.n, batch_size=args.bs, model_name=args.model)
