"""E3 v2: GENA-LM at best-per-task layer (from E1) vs HyenaDNA-tiny.

Closes the «cherry-picked worst layer» risk: original E3 used last-layer (L12)
embeddings for GENA-LM, but E1 showed best layer = L4–L5 for promoter/coding.
Here we extract GENA-LM embeddings at the E1-best layer for each task and
re-run the head-to-head with HyenaDNA-tiny (last hidden state, as before).

Best layers per task from E1 (results/e1_metrics.json):
- coding:   L5
- promoter: L4
- enhancer: L1   (near-chance; included for completeness)
- ocr:      L2   (near-chance; included for completeness)

Same probing protocol as E3 — frozen embeddings, LogReg, 3 seeds, mean ± std.
"""
import os
os.environ["HF_HOME"] = "/home/ekorshunov/hf-cache"

import json
import time
import argparse
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

# Best layer per task from E1 (results/e1_metrics.json)
GENA_BEST_LAYER = {
    "promoter": 4,
    "enhancer": 1,
    "ocr":      2,
    "coding":   5,
}

GENA_MODEL = "AIRI-Institute/gena-lm-bert-base-t2t"
HYENA_MODEL = "LongSafari/hyenadna-tiny-1k-seqlen-hf"


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


def get_embeddings_at_layer(seqs, tok, model, device, layer_idx, batch_size=8, max_len=512):
    """Mean-pool hidden states at a specific layer (0=embedding, 1..12=transformer layers)."""
    embs = []
    for i in range(0, len(seqs), batch_size):
        batch = seqs[i:i+batch_size]
        try:
            inp = tok(batch, return_tensors="pt", padding=True, truncation=True, max_length=max_len).to(device)
        except Exception:
            tokenized = [tok(s, return_tensors="pt", truncation=True, max_length=max_len)["input_ids"][0] for s in batch]
            max_t = max(t.shape[0] for t in tokenized)
            pad_id = getattr(tok, "pad_token_id", 0) or 0
            ids = torch.full((len(tokenized), max_t), pad_id, dtype=torch.long)
            am = torch.zeros((len(tokenized), max_t), dtype=torch.long)
            for j, t in enumerate(tokenized):
                ids[j, :t.shape[0]] = t
                am[j, :t.shape[0]] = 1
            inp = {"input_ids": ids.to(device), "attention_mask": am.to(device)}

        with torch.no_grad():
            try:
                out = model(**inp)
            except TypeError:
                out = model(input_ids=inp["input_ids"])

        if hasattr(out, "hidden_states") and out.hidden_states is not None:
            h = out.hidden_states[layer_idx]  # SPECIFIC layer
        else:
            raise RuntimeError("Model output does not expose hidden_states — set output_hidden_states=True.")

        if "attention_mask" in inp:
            mask = inp["attention_mask"].unsqueeze(-1).float()
            emb = (h * mask).sum(1) / mask.sum(1).clamp(min=1)
        else:
            emb = h.mean(1)
        embs.append(emb.cpu().float().numpy())
    return np.concatenate(embs, axis=0)


def get_embeddings_last(seqs, tok, model, device, batch_size=8, max_len=512):
    """Generic last-hidden-state mean-pool (for HyenaDNA)."""
    embs = []
    for i in range(0, len(seqs), batch_size):
        batch = seqs[i:i+batch_size]
        try:
            inp = tok(batch, return_tensors="pt", padding=True, truncation=True, max_length=max_len).to(device)
        except Exception:
            tokenized = [tok(s, return_tensors="pt", truncation=True, max_length=max_len)["input_ids"][0] for s in batch]
            max_t = max(t.shape[0] for t in tokenized)
            pad_id = getattr(tok, "pad_token_id", 0) or 0
            ids = torch.full((len(tokenized), max_t), pad_id, dtype=torch.long)
            am = torch.zeros((len(tokenized), max_t), dtype=torch.long)
            for j, t in enumerate(tokenized):
                ids[j, :t.shape[0]] = t
                am[j, :t.shape[0]] = 1
            inp = {"input_ids": ids.to(device), "attention_mask": am.to(device)}

        with torch.no_grad():
            try:
                out = model(**inp)
            except TypeError:
                out = model(input_ids=inp["input_ids"])

        if hasattr(out, "last_hidden_state"):
            h = out.last_hidden_state
        elif isinstance(out, (tuple, list)):
            h = out[0]
        else:
            h = out

        if "attention_mask" in inp:
            mask = inp["attention_mask"].unsqueeze(-1).float()
            emb = (h * mask).sum(1) / mask.sum(1).clamp(min=1)
        else:
            emb = h.mean(1)
        embs.append(emb.cpu().float().numpy())
    return np.concatenate(embs, axis=0)


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
    return {
        "f1_mean": float(np.mean(f1s)), "f1_std": float(np.std(f1s)),
        "mcc_mean": float(np.mean(mccs)),
        "auc_mean": float(np.mean(aucs)),
    }


def main(n_per_class=500, results_dir="/home/ekorshunov/gena-lm-airi/results", seed=42, batch_size=8):
    print(f"=== E3-v2: GENA-LM @ best-per-task layer vs HyenaDNA-tiny ===", flush=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}", flush=True)

    rng = np.random.default_rng(seed)
    data = {}
    for task_name, ds_name in TASKS.items():
        ds = load_dataset(ds_name, split="test")
        seqs, y = sample_balanced(ds, n_per_class, rng)
        data[task_name] = {"seqs": seqs, "y": y}
        print(f"  {task_name}: {len(seqs)} seqs", flush=True)

    results = {}

    # ---- GENA-LM at best layer per task ----
    print(f"\n========== GENA-LM @ best-per-task layer ==========", flush=True)
    tok = AutoTokenizer.from_pretrained(GENA_MODEL, trust_remote_code=True)
    config = AutoConfig.from_pretrained(GENA_MODEL, trust_remote_code=True)
    config.output_hidden_states = True
    model = AutoModel.from_pretrained(GENA_MODEL, trust_remote_code=True, config=config).to(device).eval()
    n_params = sum(p.numel() for p in model.parameters())
    print(f"  loaded: {type(model).__name__}, {n_params/1e6:.1f}M params", flush=True)

    gena_res = {"params_M": round(n_params/1e6, 1), "best_layer_per_task": GENA_BEST_LAYER, "tasks": {}}
    for task_name, layer_idx in GENA_BEST_LAYER.items():
        print(f"\n  -- task: {task_name} | layer L{layer_idx} --", flush=True)
        t1 = time.time()
        X = get_embeddings_at_layer(data[task_name]["seqs"], tok, model, device, layer_idx=layer_idx, batch_size=batch_size)
        y = data[task_name]["y"]
        r = probe(X, y)
        r["layer"] = int(layer_idx)
        r["extract_time_s"] = round(time.time() - t1, 1)
        gena_res["tasks"][task_name] = r
        print(f"    F1={r['f1_mean']:.3f}+/-{r['f1_std']:.3f}  MCC={r['mcc_mean']:.3f}  AUC={r['auc_mean']:.3f}", flush=True)
    results["GENA-LM-base @ best-per-task layer"] = gena_res
    del model
    torch.cuda.empty_cache()

    # ---- HyenaDNA last hidden state ----
    print(f"\n========== HyenaDNA-tiny (last hidden state) ==========", flush=True)
    tok = AutoTokenizer.from_pretrained(HYENA_MODEL, trust_remote_code=True)
    config = AutoConfig.from_pretrained(HYENA_MODEL, trust_remote_code=True)
    model = AutoModel.from_pretrained(HYENA_MODEL, trust_remote_code=True, config=config).to(device).eval()
    n_params = sum(p.numel() for p in model.parameters())
    print(f"  loaded: {type(model).__name__}, {n_params/1e6:.1f}M params", flush=True)

    hyena_res = {"params_M": round(n_params/1e6, 2), "tasks": {}}
    for task_name in TASKS:
        print(f"\n  -- task: {task_name} (HyenaDNA last) --", flush=True)
        t1 = time.time()
        X = get_embeddings_last(data[task_name]["seqs"], tok, model, device, batch_size=batch_size)
        y = data[task_name]["y"]
        r = probe(X, y)
        r["layer"] = "last"
        r["extract_time_s"] = round(time.time() - t1, 1)
        hyena_res["tasks"][task_name] = r
        print(f"    F1={r['f1_mean']:.3f}+/-{r['f1_std']:.3f}  MCC={r['mcc_mean']:.3f}  AUC={r['auc_mean']:.3f}", flush=True)
    results["HyenaDNA-tiny"] = hyena_res

    # ---- Save ----
    results_dir = Path(results_dir)
    (results_dir / "tables").mkdir(parents=True, exist_ok=True)
    with open(results_dir / "e3_v2_metrics.json", "w") as f:
        json.dump({"experiment": "E3-v2: GENA-LM @ best-per-task layer vs HyenaDNA-tiny",
                   "n_per_class": n_per_class, "results": results}, f, indent=2)

    # Print final table
    print(f"\n\n=== FINAL COMPARISON (GENA-LM @ best layer) ===", flush=True)
    print(f"{'Model / config':<40}{'Task':<12}{'F1 (mean±std)':<18}{'MCC':<8}{'AUC':<8}")
    for model_name, mres in results.items():
        for task_name, r in mres["tasks"].items():
            layer_tag = f"L{r['layer']}" if isinstance(r['layer'], int) else str(r['layer'])
            f1str = f"{r['f1_mean']:.3f}±{r['f1_std']:.3f}"
            print(f"{model_name+' ('+layer_tag+')':<40}{task_name:<12}{f1str:<18}{r['mcc_mean']:<8.3f}{r['auc_mean']:<8.3f}")

    print(f"\nSaved: {results_dir / 'e3_v2_metrics.json'}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--bs", type=int, default=8)
    args = ap.parse_args()
    main(n_per_class=args.n, batch_size=args.bs)
