"""E6: DART-Eval-style 3-way comparison.

Adds a supervised CNN baseline trained from scratch on raw DNA, alongside the
frozen-embedding probing of GENA-LM and HyenaDNA from E3. This is the core
question of DART-Eval (Patel et al., NeurIPS 2024): do DNA foundation models
beat lightweight supervised baselines on regulatory DNA classification?

Task: cCRE-style discrimination on the Genomic Benchmarks `human_ocr_ensembl`
(Open Chromatin Regions, 1:1 vs random control).

Models compared:
1. GENA-LM-base-t2t  (110M, frozen, last hidden state mean-pooled, LogReg)
2. HyenaDNA-tiny     (0.4M, frozen, last hidden state mean-pooled, LogReg)
3. Tiny 1D-CNN       (~50k params, trained from scratch on one-hot DNA)

Plus standard sklearn baselines:
4. 3-mer counts + LogReg (no DL at all)

The cleanest DART-Eval-style comparison we can run without Synapse / hg38
reference. Same protocol as Patel et al. Task 1 / 2: AUROC, AUPRC, params.
"""
import os
os.environ["HF_HOME"] = "/home/ekorshunov/hf-cache"

import json
import time
import argparse
from pathlib import Path
from collections import Counter
from itertools import product

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
from datasets import load_dataset
from transformers import AutoTokenizer, AutoModel, AutoConfig
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import f1_score, matthews_corrcoef, roc_auc_score, average_precision_score
from sklearn.model_selection import train_test_split


# DART-Eval-style Task: Open Chromatin Regions vs random (regulatory DNA classification)
DATASET = "katarinagresova/Genomic_Benchmarks_human_ocr_ensembl"
NUC2IDX = {"A": 0, "C": 1, "G": 2, "T": 3}

# Random seeds used for all probes (LogReg + TinyCNN) — reported as mean ± std.
SEEDS = (42, 7, 123)


# ============================================================
# Data
# ============================================================
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


def one_hot_encode(seqs, max_len=500):
    """Convert list of strings to (N, 4, L) one-hot tensor."""
    n = len(seqs)
    out = np.zeros((n, 4, max_len), dtype=np.float32)
    for i, s in enumerate(seqs):
        s = s[:max_len].upper()
        for j, nuc in enumerate(s):
            if nuc in NUC2IDX:
                out[i, NUC2IDX[nuc], j] = 1.0
    return out


def kmer_counts(seqs, k=3):
    """k-mer count feature vectors, normalised."""
    vocab = ["".join(t) for t in product("ACGT", repeat=k)]
    vocab_idx = {kmer: i for i, kmer in enumerate(vocab)}
    out = np.zeros((len(seqs), len(vocab)), dtype=np.float32)
    for i, s in enumerate(seqs):
        s = s.upper()
        for j in range(len(s) - k + 1):
            kmer = s[j:j+k]
            if kmer in vocab_idx:
                out[i, vocab_idx[kmer]] += 1
        if out[i].sum() > 0:
            out[i] /= out[i].sum()
    return out


# ============================================================
# Foundation-model embeddings (same as E3)
# ============================================================
def load_foundation_model(model_name, device):
    tok = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    config = AutoConfig.from_pretrained(model_name, trust_remote_code=True)
    if hasattr(config, "output_hidden_states"):
        config.output_hidden_states = True
    model = AutoModel.from_pretrained(model_name, trust_remote_code=True, config=config).to(device).eval()
    return tok, model, config


def get_embeddings_generic(seqs, tok, model, device, batch_size=8, max_len=512):
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

        h = None
        if hasattr(out, "hidden_states") and out.hidden_states is not None:
            h = out.hidden_states[-1]
        elif hasattr(out, "last_hidden_state"):
            h = out.last_hidden_state
        elif isinstance(out, (tuple, list)):
            h = out[0]
        else:
            h = out
        if h.dim() == 3 and h.shape[-1] > 10000:
            raise RuntimeError(f"Got {h.shape} — looks like MLM logits, not hidden states.")

        if "attention_mask" in inp:
            mask = inp["attention_mask"].unsqueeze(-1).float()
            emb = (h * mask).sum(1) / mask.sum(1).clamp(min=1)
        else:
            emb = h.mean(1)
        embs.append(emb.cpu().float().numpy())
    return np.concatenate(embs, axis=0)


# ============================================================
# Supervised CNN (from scratch, no foundation model)
# ============================================================
class TinyCNN(nn.Module):
    """One-hot DNA -> 2x Conv1D -> GlobalMaxPool -> Dense (~50k params)."""
    def __init__(self, max_len=500, n_filters_1=64, n_filters_2=128, kernel_1=21, kernel_2=11):
        super().__init__()
        self.conv1 = nn.Conv1d(4, n_filters_1, kernel_size=kernel_1, padding=kernel_1//2)
        self.pool1 = nn.MaxPool1d(4)
        self.conv2 = nn.Conv1d(n_filters_1, n_filters_2, kernel_size=kernel_2, padding=kernel_2//2)
        self.gap = nn.AdaptiveMaxPool1d(1)
        self.dropout = nn.Dropout(0.3)
        self.fc1 = nn.Linear(n_filters_2, 64)
        self.fc2 = nn.Linear(64, 1)

    def forward(self, x):
        x = F.relu(self.conv1(x))
        x = self.pool1(x)
        x = F.relu(self.conv2(x))
        x = self.gap(x).squeeze(-1)
        x = self.dropout(x)
        x = F.relu(self.fc1(x))
        return self.fc2(x).squeeze(-1)


def train_cnn(X_train, y_train, X_val, y_val, device, max_len=500,
              epochs=20, batch_size=32, lr=1e-3, patience=5, seed=42):
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = TinyCNN(max_len=max_len).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    loss_fn = nn.BCEWithLogitsLoss()

    tr_ds = TensorDataset(torch.from_numpy(X_train), torch.from_numpy(y_train).float())
    tr_dl = DataLoader(tr_ds, batch_size=batch_size, shuffle=True)
    val_x = torch.from_numpy(X_val).to(device)
    val_y = torch.from_numpy(y_val).float().to(device)

    best_val_auc = 0.0
    best_state = None
    patience_count = 0
    history = []
    for epoch in range(epochs):
        model.train()
        for xb, yb in tr_dl:
            xb, yb = xb.to(device), yb.to(device)
            opt.zero_grad()
            logits = model(xb)
            loss = loss_fn(logits, yb)
            loss.backward()
            opt.step()
        # val
        model.eval()
        with torch.no_grad():
            val_logits = model(val_x).cpu().numpy()
            val_y_cpu = val_y.cpu().numpy()
            val_auc = roc_auc_score(val_y_cpu, val_logits)
        history.append({"epoch": epoch+1, "val_auc": float(val_auc)})
        if val_auc > best_val_auc:
            best_val_auc = val_auc
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            patience_count = 0
        else:
            patience_count += 1
            if patience_count >= patience:
                break
    model.load_state_dict(best_state)
    return model, history


# ============================================================
# Eval pipeline
# ============================================================
def eval_predictions(y_true, y_score):
    y_pred = (y_score > 0.5).astype(int)
    return {
        "f1": float(f1_score(y_true, y_pred)),
        "mcc": float(matthews_corrcoef(y_true, y_pred)),
        "auroc": float(roc_auc_score(y_true, y_score)),
        "auprc": float(average_precision_score(y_true, y_score)),
    }


def probe_with_logreg(X_train, X_val, X_test, y_train, y_val, y_test, seed=42):
    """Train LogReg on train, validate on val, report on test."""
    clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0, random_state=seed))
    clf.fit(X_train, y_train)
    val_proba = clf.predict_proba(X_val)[:, 1]
    test_proba = clf.predict_proba(X_test)[:, 1]
    return {
        "val": eval_predictions(y_val, val_proba),
        "test": eval_predictions(y_test, test_proba),
    }


def probe_with_logreg_multiseed(X_train, X_val, X_test, y_train, y_val, y_test, seeds=(42, 7, 123)):
    """Train LogReg across multiple seeds; report mean ± std of test/val metrics."""
    runs = []
    for s in seeds:
        runs.append(probe_with_logreg(X_train, X_val, X_test, y_train, y_val, y_test, seed=s))
    # Aggregate test/val metrics across seeds
    keys = ("auroc", "auprc", "f1", "mcc")
    agg = {"val": {}, "test": {}, "seeds": list(seeds), "per_seed": runs}
    for split in ("val", "test"):
        for k in keys:
            vals = [r[split][k] for r in runs]
            agg[split][k] = float(np.mean(vals))
            agg[split][k + "_std"] = float(np.std(vals))
    return agg


def train_cnn_multiseed(X_train, y_train, X_val, y_val, X_test, y_test, device, max_len=500,
                        epochs=20, batch_size=32, lr=1e-3, patience=5, seeds=(42, 7, 123)):
    """Train TinyCNN across multiple seeds; report mean ± std of val/test metrics."""
    runs = []
    last_model = None
    for s in seeds:
        model, history = train_cnn(X_train, y_train, X_val, y_val, device,
                                   max_len=max_len, epochs=epochs, batch_size=batch_size,
                                   lr=lr, patience=patience, seed=s)
        model.eval()
        with torch.no_grad():
            val_score = torch.sigmoid(model(torch.from_numpy(X_val).to(device))).cpu().numpy()
            test_score = torch.sigmoid(model(torch.from_numpy(X_test).to(device))).cpu().numpy()
        runs.append({
            "val": eval_predictions(y_val, val_score),
            "test": eval_predictions(y_test, test_score),
            "epochs_trained": len(history),
        })
        last_model = model
    keys = ("auroc", "auprc", "f1", "mcc")
    agg = {"val": {}, "test": {}, "seeds": list(seeds), "per_seed": runs}
    for split in ("val", "test"):
        for k in keys:
            vals = [r[split][k] for r in runs]
            agg[split][k] = float(np.mean(vals))
            agg[split][k + "_std"] = float(np.std(vals))
    return agg, last_model


# ============================================================
# Main
# ============================================================
def main(n_per_class=2500, max_len=500, results_dir="/home/ekorshunov/gena-lm-airi/results",
         seed=42, batch_size=8, cnn_epochs=20):
    print(f"=== E6 DART-Eval-style 3-way comparison ===", flush=True)
    print(f"task: cCRE-like discrimination on Open Chromatin Regions", flush=True)
    print(f"n_per_class={n_per_class}, max_len={max_len}", flush=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}", flush=True)

    # 1) Load and split data
    print(f"\nLoading {DATASET}...", flush=True)
    rng = np.random.default_rng(seed)
    ds = load_dataset(DATASET, split="test")
    seqs, y = sample_balanced(ds, n_per_class, rng)
    print(f"  total: {len(seqs)} balanced, avg len={int(np.mean([len(s) for s in seqs]))} bp", flush=True)

    # train/val/test 60/20/20 — stratified
    idx_tr, idx_test, ytr, ytest = train_test_split(np.arange(len(seqs)), y, test_size=0.2, random_state=seed, stratify=y)
    idx_tr, idx_val, ytr2, yval = train_test_split(idx_tr, ytr, test_size=0.25, random_state=seed, stratify=ytr)
    seqs_train = [seqs[i] for i in idx_tr]
    seqs_val = [seqs[i] for i in idx_val]
    seqs_test = [seqs[i] for i in idx_test]
    y_train = ytr2
    y_val = yval
    y_test = ytest
    print(f"  splits: train={len(seqs_train)} val={len(seqs_val)} test={len(seqs_test)}", flush=True)

    results = {}

    # 2) k-mer + LogReg baseline
    print(f"\n--- BASELINE: 3-mer + LogReg ---", flush=True)
    t0 = time.time()
    X_train = kmer_counts(seqs_train, k=3)
    X_val = kmer_counts(seqs_val, k=3)
    X_test = kmer_counts(seqs_test, k=3)
    feat_time = time.time() - t0
    print(f"  feature dim: {X_train.shape[1]}, feature_time={feat_time:.1f}s", flush=True)
    t0 = time.time()
    rmets = probe_with_logreg_multiseed(X_train, X_val, X_test, y_train, y_val, y_test, seeds=SEEDS)
    rmets["train_time_s"] = round(time.time() - t0, 1)
    rmets["params"] = X_train.shape[1] + 1
    rmets["embedding_dim"] = X_train.shape[1]
    results["3mer+LogReg"] = rmets
    print(f"  test: AUROC={rmets['test']['auroc']:.3f}±{rmets['test']['auroc_std']:.3f}  AUPRC={rmets['test']['auprc']:.3f}±{rmets['test']['auprc_std']:.3f}", flush=True)

    # 3) Supervised CNN from scratch
    print(f"\n--- SUPERVISED: TinyCNN (from scratch) ---", flush=True)
    t0 = time.time()
    X_train_cnn = one_hot_encode(seqs_train, max_len=max_len)
    X_val_cnn = one_hot_encode(seqs_val, max_len=max_len)
    X_test_cnn = one_hot_encode(seqs_test, max_len=max_len)
    encode_time = time.time() - t0
    print(f"  one-hot shape: {X_train_cnn.shape}, encode_time={encode_time:.1f}s", flush=True)

    t0 = time.time()
    cnn_mets, last_cnn = train_cnn_multiseed(
        X_train_cnn, y_train, X_val_cnn, y_val, X_test_cnn, y_test,
        device, max_len=max_len, epochs=cnn_epochs, seeds=SEEDS,
    )
    train_time = time.time() - t0
    cnn_params = sum(p.numel() for p in last_cnn.parameters())
    cnn_mets["params"] = cnn_params
    cnn_mets["train_time_s"] = round(train_time, 1)
    results["TinyCNN-supervised"] = cnn_mets
    print(f"  params: {cnn_params/1000:.1f}k, total train_time={train_time:.1f}s (across {len(SEEDS)} seeds)", flush=True)
    print(f"  test: AUROC={cnn_mets['test']['auroc']:.3f}±{cnn_mets['test']['auroc_std']:.3f}  AUPRC={cnn_mets['test']['auprc']:.3f}±{cnn_mets['test']['auprc_std']:.3f}", flush=True)

    # 4) Foundation models — frozen embeddings + LogReg
    foundation_models = {
        "GENA-LM-base":  "AIRI-Institute/gena-lm-bert-base-t2t",
        "HyenaDNA-tiny": "LongSafari/hyenadna-tiny-1k-seqlen-hf",
    }
    for name, model_id in foundation_models.items():
        print(f"\n--- FROZEN: {name} ({model_id}) ---", flush=True)
        t0 = time.time()
        try:
            tok, model, config = load_foundation_model(model_id, device)
        except Exception as e:
            print(f"  ERROR loading: {e}", flush=True)
            results[name] = {"status": "load_failed", "error": str(e)}
            continue
        n_params = sum(p.numel() for p in model.parameters())
        print(f"  {type(model).__name__}, {n_params/1e6:.1f}M params, load_time={time.time()-t0:.1f}s", flush=True)

        t0 = time.time()
        Xtr = get_embeddings_generic(seqs_train, tok, model, device, batch_size=batch_size)
        Xva = get_embeddings_generic(seqs_val, tok, model, device, batch_size=batch_size)
        Xte = get_embeddings_generic(seqs_test, tok, model, device, batch_size=batch_size)
        ext = time.time() - t0
        print(f"  embeddings: train={Xtr.shape} extract_time={ext:.1f}s", flush=True)

        t0 = time.time()
        mets = probe_with_logreg_multiseed(Xtr, Xva, Xte, y_train, y_val, y_test, seeds=SEEDS)
        probe_time = time.time() - t0
        mets["params"] = n_params
        mets["embedding_dim"] = int(Xtr.shape[1])
        mets["extract_time_s"] = round(ext, 1)
        mets["probe_time_s"] = round(probe_time, 1)
        results[name] = mets
        print(f"  test: AUROC={mets['test']['auroc']:.3f}±{mets['test']['auroc_std']:.3f}  AUPRC={mets['test']['auprc']:.3f}±{mets['test']['auprc_std']:.3f}", flush=True)

        del model
        torch.cuda.empty_cache()

    # 5) Print and save final table
    print(f"\n\n=== FINAL COMPARISON (test set, n={len(SEEDS)} seeds, mean ± std) ===", flush=True)
    header = f"{'Model':<25}{'Params':<12}{'AUROC':<18}{'AUPRC':<18}{'F1':<18}{'MCC':<18}"
    print(header, flush=True)
    print("-" * len(header), flush=True)
    for name, r in results.items():
        if "test" not in r:
            print(f"{name:<25}FAILED")
            continue
        p_str = f"{r['params']/1e6:.1f}M" if r["params"] > 1e6 else f"{r['params']/1e3:.1f}k"
        t = r["test"]
        def cell(k):
            return f"{t[k]:.3f}±{t.get(k+'_std', 0.0):.3f}"
        print(f"{name:<25}{p_str:<12}{cell('auroc'):<18}{cell('auprc'):<18}{cell('f1'):<18}{cell('mcc'):<18}", flush=True)

    # Save
    results_dir = Path(results_dir)
    (results_dir / "tables").mkdir(parents=True, exist_ok=True)
    out = {
        "experiment": "E6 DART-Eval-style comparison",
        "dataset": DATASET,
        "n_per_class": n_per_class,
        "splits": {"train": len(seqs_train), "val": len(seqs_val), "test": len(seqs_test)},
        "results": results,
    }
    with open(results_dir / "e6_metrics.json", "w") as f:
        json.dump(out, f, indent=2)

    import csv
    with open(results_dir / "tables" / "e6_results.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "params", "AUROC_mean", "AUROC_std", "AUPRC_mean", "AUPRC_std",
                    "F1_mean", "F1_std", "MCC_mean", "MCC_std", "n_seeds"])
        for name, r in results.items():
            if "test" not in r:
                continue
            t = r["test"]
            w.writerow([name, r["params"],
                        t["auroc"], t.get("auroc_std", 0.0),
                        t["auprc"], t.get("auprc_std", 0.0),
                        t["f1"],    t.get("f1_std", 0.0),
                        t["mcc"],   t.get("mcc_std", 0.0),
                        len(SEEDS)])

    print(f"\nSaved: {results_dir / 'e6_metrics.json'}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2500, help="sequences per class")
    ap.add_argument("--max_len", type=int, default=500)
    ap.add_argument("--bs", type=int, default=8)
    ap.add_argument("--cnn_epochs", type=int, default=20)
    args = ap.parse_args()
    main(n_per_class=args.n, max_len=args.max_len, batch_size=args.bs, cnn_epochs=args.cnn_epochs)
