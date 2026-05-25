"""E1-large: Layer-wise linear probing for gena-lm-bert-LARGE-t2t (24 layers, 336M).

Tests whether the «mid-layer peak» pattern observed for bert-base (12 layers,
peak at L4–L5) holds, shifts deeper, or breaks in the larger 24-layer variant.

Same protocol as E1: 4 binary regulatory tasks, frozen mean-pool, LogReg,
3 seeds, mean ± std per layer.
"""
import os
os.environ["HF_HOME"] = "/home/ekorshunov/hf-cache"

import json
import argparse
import time
from pathlib import Path

import numpy as np
import torch
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
    tok = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    config = AutoConfig.from_pretrained(model_name, trust_remote_code=True)
    config.output_hidden_states = True
    model = AutoModel.from_pretrained(model_name, trust_remote_code=True, config=config).to(device).eval()
    return tok, model, config


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


def get_all_layer_embeddings(seqs, tok, model, device, n_layers, batch_size=4, max_len=512):
    """Return shape (n_layers+1, N, hidden_dim) — all 25 layers (0=embedding + 24 transformer)."""
    embs = [[] for _ in range(n_layers + 1)]
    for i in range(0, len(seqs), batch_size):
        batch = seqs[i:i+batch_size]
        inp = tok(batch, return_tensors="pt", padding=True, truncation=True, max_length=max_len).to(device)
        with torch.no_grad():
            out = model(**inp)
        # out.hidden_states is a tuple of (n_layers+1) tensors
        mask = inp["attention_mask"].unsqueeze(-1).float()
        for L in range(n_layers + 1):
            h = out.hidden_states[L]
            emb = (h * mask).sum(1) / mask.sum(1).clamp(min=1)
            embs[L].append(emb.cpu().float().numpy())
    return [np.concatenate(layer_embs, axis=0) for layer_embs in embs]


def probe(X, y, n_seeds=3):
    f1s, mccs, aucs = [], [], []
    for seed in range(n_seeds):
        Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.3, random_state=seed, stratify=y)
        clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0, random_state=seed))
        clf.fit(Xtr, ytr)
        pred = clf.predict(Xte)
        prob = clf.predict_proba(Xte)[:, 1]
        f1s.append(f1_score(yte, pred))
        mccs.append(matthews_corrcoef(yte, pred))
        aucs.append(roc_auc_score(yte, prob))
    return {"f1_mean": float(np.mean(f1s)), "f1_std": float(np.std(f1s)),
            "mcc_mean": float(np.mean(mccs)), "auc_mean": float(np.mean(aucs))}


def main(n_per_class=300, results_dir="/home/ekorshunov/gena-lm-airi/results",
         model_name="AIRI-Institute/gena-lm-bert-large-t2t", seed=42, batch_size=4):
    print(f"=== E1-large: layer-wise probing on {model_name} ===", flush=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}", flush=True)

    tok, model, config = load_model(model_name, device)
    n_layers = config.num_hidden_layers
    n_params = sum(p.numel() for p in model.parameters())
    print(f"  {type(model).__name__}, {n_params/1e6:.1f}M params, {n_layers} transformer layers", flush=True)

    rng = np.random.default_rng(seed)
    results = {"model": model_name, "n_params_M": round(n_params/1e6, 1),
               "n_layers": int(n_layers), "n_per_class": n_per_class, "tasks": {}}
    for task_name, ds_name in TASKS.items():
        print(f"\n--- task: {task_name} ---", flush=True)
        ds = load_dataset(ds_name, split="test")
        seqs, y = sample_balanced(ds, n_per_class, rng)
        print(f"  {len(seqs)} sequences, computing all-layer embeddings...", flush=True)
        t0 = time.time()
        all_embs = get_all_layer_embeddings(seqs, tok, model, device, n_layers=n_layers, batch_size=batch_size)
        print(f"  extracted {len(all_embs)} layers in {time.time()-t0:.1f}s", flush=True)

        per_layer = []
        for L, X in enumerate(all_embs):
            r = probe(X, y)
            r["layer"] = int(L)
            per_layer.append(r)
        f1s = [r["f1_mean"] for r in per_layer]
        best_L = int(np.argmax(f1s))
        results["tasks"][task_name] = {"per_layer": per_layer, "best_layer": best_L,
                                       "best_f1": f1s[best_L]}
        print(f"  best layer: L{best_L}, F1={f1s[best_L]:.3f}", flush=True)
        # Print full per-layer curve
        for L, r in enumerate(per_layer):
            marker = " ★" if L == best_L else ""
            print(f"    L{L:2d}: F1={r['f1_mean']:.3f}±{r['f1_std']:.3f}{marker}", flush=True)

    # Save
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    out_file = results_dir / "e1_large_metrics.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_file}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=300, help="seqs per class (smaller than e1 base because large model is heavier)")
    ap.add_argument("--bs", type=int, default=4)
    ap.add_argument("--model", type=str, default="AIRI-Institute/gena-lm-bert-large-t2t")
    args = ap.parse_args()
    main(n_per_class=args.n, model_name=args.model, batch_size=args.bs)
