"""E7: Fine-tune GENA-LM-base on promoter task — one calibration data point.

Addresses the «frozen vs fine-tuned» critique of E3/E6 head-to-head.
Trains a classification head on top of `gena-lm-bert-base-t2t` with the
backbone unfrozen, on `human_nontata_promoters`. Reports final test F1.

Expected: fine-tuned F1 substantially higher than frozen-probing F1 (≈ 0.78)
— illustrating that our frozen-regime conclusions do NOT extrapolate to
fine-tuned usage of GENA-LM.

Protocol:
- AdamW, lr 2e-5, batch 16, 3 epochs (early stop on val F1)
- BertForSequenceClassification head added manually on top of AutoModel
- Same 70/30 stratified split as E1/E3, 3 seeds for stability
"""
import os
os.environ["HF_HOME"] = "/home/ekorshunov/hf-cache"

import json
import argparse
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
from datasets import load_dataset
from transformers import AutoTokenizer, AutoModel, AutoConfig, get_linear_schedule_with_warmup
from sklearn.metrics import f1_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import train_test_split


MODEL = "AIRI-Institute/gena-lm-bert-base-t2t"
TASK_DS = "katarinagresova/Genomic_Benchmarks_human_nontata_promoters"


class SeqDataset(Dataset):
    def __init__(self, seqs, labels, tok, max_len=512):
        self.seqs = seqs
        self.labels = labels
        self.tok = tok
        self.max_len = max_len

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, idx):
        enc = self.tok(self.seqs[idx], return_tensors="pt", padding="max_length",
                       truncation=True, max_length=self.max_len)
        return {
            "input_ids": enc["input_ids"].squeeze(0),
            "attention_mask": enc["attention_mask"].squeeze(0),
            "label": torch.tensor(int(self.labels[idx]), dtype=torch.long),
        }


class GenaClassifier(nn.Module):
    """GENA-LM backbone + mean-pool + linear classification head (2-way)."""
    def __init__(self, backbone, hidden_size, n_classes=2, dropout=0.1):
        super().__init__()
        self.backbone = backbone
        self.dropout = nn.Dropout(dropout)
        self.classifier = nn.Linear(hidden_size, n_classes)

    def forward(self, input_ids, attention_mask):
        out = self.backbone(input_ids=input_ids, attention_mask=attention_mask)
        # Use hidden_states[-1] for the custom GENA-LM modeling (avoid out.logits which is MLM head)
        if hasattr(out, "hidden_states") and out.hidden_states is not None:
            h = out.hidden_states[-1]
        elif hasattr(out, "last_hidden_state"):
            h = out.last_hidden_state
        else:
            h = out[0]
        mask = attention_mask.unsqueeze(-1).float()
        pooled = (h * mask).sum(1) / mask.sum(1).clamp(min=1)
        return self.classifier(self.dropout(pooled))


def evaluate(model, loader, device):
    model.eval()
    probs, preds, ys = [], [], []
    with torch.no_grad():
        for batch in loader:
            ids = batch["input_ids"].to(device)
            am = batch["attention_mask"].to(device)
            logits = model(ids, am)
            prob = torch.softmax(logits, dim=-1)[:, 1].cpu().numpy()
            pred = logits.argmax(dim=-1).cpu().numpy()
            probs.extend(prob.tolist())
            preds.extend(pred.tolist())
            ys.extend(batch["label"].numpy().tolist())
    return {
        "f1": float(f1_score(ys, preds)),
        "mcc": float(matthews_corrcoef(ys, preds)),
        "auc": float(roc_auc_score(ys, probs)),
    }


def train_one_seed(seqs_train, y_train, seqs_val, y_val, seqs_test, y_test,
                   tok, device, seed, batch_size=16, epochs=3, lr=2e-5, max_len=512):
    torch.manual_seed(seed)
    np.random.seed(seed)

    config = AutoConfig.from_pretrained(MODEL, trust_remote_code=True)
    config.output_hidden_states = True
    backbone = AutoModel.from_pretrained(MODEL, trust_remote_code=True, config=config)
    model = GenaClassifier(backbone, hidden_size=config.hidden_size, n_classes=2).to(device)

    train_ds = SeqDataset(seqs_train, y_train, tok, max_len=max_len)
    val_ds   = SeqDataset(seqs_val,   y_val,   tok, max_len=max_len)
    test_ds  = SeqDataset(seqs_test,  y_test,  tok, max_len=max_len)
    train_dl = DataLoader(train_ds, batch_size=batch_size, shuffle=True, num_workers=2, pin_memory=True)
    val_dl   = DataLoader(val_ds,   batch_size=batch_size, shuffle=False, num_workers=2, pin_memory=True)
    test_dl  = DataLoader(test_ds,  batch_size=batch_size, shuffle=False, num_workers=2, pin_memory=True)

    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=0.01)
    total_steps = len(train_dl) * epochs
    sched = get_linear_schedule_with_warmup(opt, num_warmup_steps=int(0.1*total_steps), num_training_steps=total_steps)
    loss_fn = nn.CrossEntropyLoss()

    best_val_f1 = 0.0
    best_test = None
    history = []
    for ep in range(epochs):
        model.train()
        t0 = time.time()
        for batch in train_dl:
            ids = batch["input_ids"].to(device)
            am = batch["attention_mask"].to(device)
            labels = batch["label"].to(device)
            opt.zero_grad()
            logits = model(ids, am)
            loss = loss_fn(logits, labels)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            sched.step()
        val_m = evaluate(model, val_dl, device)
        test_m = evaluate(model, test_dl, device)
        history.append({"epoch": ep+1, "val": val_m, "test": test_m, "time_s": round(time.time()-t0, 1)})
        print(f"    epoch {ep+1}/{epochs}: val F1={val_m['f1']:.3f} test F1={test_m['f1']:.3f} ({time.time()-t0:.1f}s)", flush=True)
        if val_m["f1"] > best_val_f1:
            best_val_f1 = val_m["f1"]
            best_test = test_m
    return best_test, history


def main(n_per_class=2000, batch_size=16, epochs=3, lr=2e-5, max_len=512,
         results_dir="/home/ekorshunov/gena-lm-airi/results", seeds=(42, 7, 123)):
    print(f"=== E7 Fine-tune GENA-LM-base on promoter ===", flush=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}, n_per_class={n_per_class}, epochs={epochs}, batch={batch_size}", flush=True)

    tok = AutoTokenizer.from_pretrained(MODEL, trust_remote_code=True)

    print(f"\nLoading {TASK_DS}...", flush=True)
    rng = np.random.default_rng(42)
    ds = load_dataset(TASK_DS, split="test")
    labels = np.array(ds["label"])
    seqs_all = np.array(ds["seq"])
    seqs, y = [], []
    for lbl in (0, 1):
        idx_pool = np.where(labels == lbl)[0]
        n = min(n_per_class, len(idx_pool))
        idx = rng.choice(idx_pool, size=n, replace=False)
        seqs.extend(seqs_all[idx].tolist())
        y.extend([lbl] * n)
    y = np.array(y)
    print(f"  total: {len(seqs)} balanced (avg len ≈ {int(np.mean([len(s) for s in seqs]))} bp)", flush=True)

    # Hold out 10% test once (same across seeds), then 70/20 train/val varies per seed
    idx_trainval, idx_test, ytv, ytest = train_test_split(np.arange(len(seqs)), y, test_size=0.1, random_state=42, stratify=y)

    per_seed = []
    for s in seeds:
        print(f"\n--- seed {s} ---", flush=True)
        idx_train, idx_val, ytr, yv = train_test_split(idx_trainval, ytv, test_size=0.2, random_state=s, stratify=ytv)
        seqs_train = [seqs[i] for i in idx_train]
        seqs_val   = [seqs[i] for i in idx_val]
        seqs_test  = [seqs[i] for i in idx_test]
        print(f"  splits: train={len(seqs_train)} val={len(seqs_val)} test={len(seqs_test)}", flush=True)
        best, hist = train_one_seed(seqs_train, ytr, seqs_val, yv, seqs_test, ytest,
                                    tok, device, seed=s, batch_size=batch_size,
                                    epochs=epochs, lr=lr, max_len=max_len)
        per_seed.append({"seed": s, "best_test": best, "history": hist})
        torch.cuda.empty_cache()

    # Aggregate
    test_f1s  = [r["best_test"]["f1"]  for r in per_seed]
    test_mccs = [r["best_test"]["mcc"] for r in per_seed]
    test_aucs = [r["best_test"]["auc"] for r in per_seed]
    out = {
        "experiment": "E7 Fine-tune GENA-LM-base on promoter (nontata)",
        "model": MODEL,
        "task": TASK_DS,
        "n_per_class": n_per_class,
        "seeds": list(seeds),
        "epochs": epochs,
        "lr": lr,
        "batch_size": batch_size,
        "max_len": max_len,
        "test_f1_mean":  float(np.mean(test_f1s)),
        "test_f1_std":   float(np.std(test_f1s)),
        "test_mcc_mean": float(np.mean(test_mccs)),
        "test_auc_mean": float(np.mean(test_aucs)),
        "per_seed": per_seed,
    }

    Path(results_dir).mkdir(parents=True, exist_ok=True)
    out_path = Path(results_dir) / "e7_fine_tune_metrics.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    print(f"\n=== FINAL (n={len(seeds)} seeds, mean ± std) ===", flush=True)
    print(f"  test F1  = {out['test_f1_mean']:.3f} ± {out['test_f1_std']:.3f}", flush=True)
    print(f"  test MCC = {out['test_mcc_mean']:.3f}", flush=True)
    print(f"  test AUC = {out['test_auc_mean']:.3f}", flush=True)
    print(f"Saved: {out_path}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--bs", type=int, default=16)
    ap.add_argument("--epochs", type=int, default=3)
    ap.add_argument("--lr", type=float, default=2e-5)
    ap.add_argument("--max_len", type=int, default=512)
    args = ap.parse_args()
    main(n_per_class=args.n, batch_size=args.bs, epochs=args.epochs, lr=args.lr, max_len=args.max_len)
