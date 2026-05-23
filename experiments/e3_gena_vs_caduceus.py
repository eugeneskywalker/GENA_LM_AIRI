"""E3: GENA-LM vs Caduceus head-to-head on 4 binary regulatory DNA tasks.

Frozen-embedding linear probing. For each model:
1. Extract mean-pooled embeddings (last layer) on test split of each Genomic Benchmark task.
2. Train LogReg (3 seeds, 70/30 stratified split).
3. Report F1 / MCC / AUC.

Result: comparison table for the Research Proposal.

Models compared:
- GENA-LM:  AIRI-Institute/gena-lm-bert-base-t2t  (110M, BPE transformer)
- Caduceus: kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16
            (~7M, BiMamba + RC-equivariance, single-nucleotide)

Both are frozen — only the linear probe is trained.
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

MODELS = {
    "GENA-LM-base":     "AIRI-Institute/gena-lm-bert-base-t2t",     # 110M, BPE transformer
    "DNABERT-2":        "zhihan1996/DNABERT-2-117M",                # 117M, BPE transformer (direct competitor)
    "HyenaDNA-tiny":    "LongSafari/hyenadna-tiny-1k-seqlen-hf",    # 1.7M, Hyena, single-nucleotide (alt arch)
    # "Caduceus-ph-7M": "kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16",  # needs mamba_ssm
}


def load_model_safe(model_name, device):
    """Try to load model; return (tok, model, config) or None on failure."""
    try:
        tok = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
        config = AutoConfig.from_pretrained(model_name, trust_remote_code=True)
        if hasattr(config, "output_hidden_states"):
            config.output_hidden_states = True
        model = AutoModel.from_pretrained(model_name, trust_remote_code=True, config=config).to(device).eval()
        return tok, model, config
    except Exception as e:
        print(f"  ERROR loading {model_name}: {type(e).__name__}: {e}", flush=True)
        return None


def get_embeddings_generic(seqs, tok, model, device, batch_size=8, max_len=512):
    """Mean-pool last hidden state. Handles BERT-style and Mamba-style models."""
    embs = []
    for i in range(0, len(seqs), batch_size):
        batch = seqs[i:i+batch_size]
        try:
            inp = tok(batch, return_tensors="pt", padding=True, truncation=True, max_length=max_len).to(device)
        except Exception:
            # Caduceus tokenizer may not accept padding directly — pad per-sample
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
                # Some models don't take attention_mask
                out = model(input_ids=inp["input_ids"])

        # Extract last hidden state. PRIORITY:
        # 1. out.hidden_states[-1]  (BertForMaskedLM with output_hidden_states=True)
        # 2. out.last_hidden_state  (BertModel, encoder-only)
        # 3. out[0] tuple                (Caduceus / raw output)
        # AVOID out.logits — it's the MLM head, vocab_size dim (32000), not 768.
        h = None
        if hasattr(out, "hidden_states") and out.hidden_states is not None:
            h = out.hidden_states[-1]
        elif hasattr(out, "last_hidden_state"):
            h = out.last_hidden_state
        elif isinstance(out, (tuple, list)):
            h = out[0]
        else:
            # raw tensor (e.g. some Mamba-based models return tensor directly)
            h = out
        # Sanity: h should be (B, T, H) with H << vocab_size
        if h.dim() == 3 and h.shape[-1] > 10000:
            raise RuntimeError(f"Got {h.shape} — likely MLM logits, not hidden states. Check output_hidden_states config.")

        # Mean pool
        if "attention_mask" in inp:
            mask = inp["attention_mask"].unsqueeze(-1).float()
            emb = (h * mask).sum(1) / mask.sum(1).clamp(min=1)
        else:
            emb = h.mean(1)
        embs.append(emb.cpu().float().numpy())
        if (i // batch_size) % 10 == 0:
            print(f"    batch {i//batch_size+1}/{(len(seqs)+batch_size-1)//batch_size}", flush=True)
    return np.concatenate(embs, axis=0)


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
    print(f"=== E3 GENA-LM vs Caduceus | n_per_class={n_per_class} ===", flush=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}", flush=True)

    rng = np.random.default_rng(seed)
    # Pre-sample sequences (same for both models — fair comparison)
    print("\nLoading datasets...", flush=True)
    data = {}
    for task_name, ds_name in TASKS.items():
        ds = load_dataset(ds_name, split="test")
        seqs, y = sample_balanced(ds, n_per_class, rng)
        data[task_name] = {"seqs": seqs, "y": y, "avg_len": int(np.mean([len(s) for s in seqs]))}
        print(f"  {task_name}: {len(seqs)} seqs, avg len={data[task_name]['avg_len']} bp", flush=True)

    results = {}
    for model_name, model_id in MODELS.items():
        print(f"\n========== MODEL: {model_name} ({model_id}) ==========", flush=True)
        t0 = time.time()
        loaded = load_model_safe(model_id, device)
        if loaded is None:
            print(f"  SKIP {model_name} — failed to load", flush=True)
            results[model_name] = {"status": "load_failed"}
            continue
        tok, model, config = loaded
        n_params = sum(p.numel() for p in model.parameters())
        print(f"  loaded {type(model).__name__}, {n_params/1e6:.1f}M params, in {time.time()-t0:.1f}s", flush=True)

        model_res = {"params_M": round(n_params/1e6, 1), "tasks": {}}
        for task_name in TASKS:
            print(f"\n  -- task: {task_name} --", flush=True)
            t1 = time.time()
            X = get_embeddings_generic(data[task_name]["seqs"], tok, model, device, batch_size=batch_size)
            y = data[task_name]["y"]
            extract_time = time.time() - t1
            print(f"    embeddings: {X.shape}, extract_time={extract_time:.1f}s", flush=True)

            r = probe(X, y)
            r["extract_time_s"] = round(extract_time, 1)
            r["embedding_dim"] = int(X.shape[1])
            model_res["tasks"][task_name] = r
            print(f"    F1={r['f1_mean']:.3f}+/-{r['f1_std']:.3f}  MCC={r['mcc_mean']:.3f}  AUC={r['auc_mean']:.3f}", flush=True)

        results[model_name] = model_res

        # free GPU mem
        del model
        torch.cuda.empty_cache()

    # Save
    results_dir = Path(results_dir)
    (results_dir / "tables").mkdir(parents=True, exist_ok=True)
    with open(results_dir / "e3_metrics.json", "w") as f:
        json.dump({"experiment": "E3 GENA-LM vs Caduceus", "n_per_class": n_per_class, "results": results}, f, indent=2)

    # CSV
    import csv
    with open(results_dir / "tables" / "e3_results.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "params_M", "task", "f1_mean", "f1_std", "mcc_mean", "auc_mean", "extract_time_s"])
        for model_name, mres in results.items():
            if "tasks" not in mres:
                continue
            for task_name, r in mres["tasks"].items():
                w.writerow([model_name, mres["params_M"], task_name,
                            f"{r['f1_mean']:.4f}", f"{r['f1_std']:.4f}",
                            f"{r['mcc_mean']:.4f}", f"{r['auc_mean']:.4f}",
                            r["extract_time_s"]])

    # Print final table
    print("\n\n=== FINAL COMPARISON ===", flush=True)
    print(f"{'Model':<25}{'Params':<10}{'Task':<12}{'F1':<8}{'MCC':<8}{'AUC':<8}")
    for model_name, mres in results.items():
        if "tasks" not in mres:
            print(f"{model_name:<25}LOAD FAILED")
            continue
        for task_name, r in mres["tasks"].items():
            print(f"{model_name:<25}{mres['params_M']:<10}{task_name:<12}{r['f1_mean']:<8.3f}{r['mcc_mean']:<8.3f}{r['auc_mean']:<8.3f}")

    print(f"\nSaved: {results_dir / 'e3_metrics.json'}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--bs", type=int, default=8)
    args = ap.parse_args()
    main(n_per_class=args.n, batch_size=args.bs)
