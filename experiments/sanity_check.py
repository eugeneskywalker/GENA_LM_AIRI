"""GENA-LM sanity check: load model, forward pass on dummy DNA, get hidden states."""
import os
os.environ["HF_HOME"] = "/home/ekorshunov/hf-cache"

import torch
from transformers import AutoTokenizer, AutoModel, AutoConfig

print("=" * 60)
print(f"PyTorch: {torch.__version__}")
print(f"CUDA: {torch.version.cuda}")
print(f"GPU: {torch.cuda.get_device_name(0)}")
print(f"GPU mem total: {torch.cuda.get_device_properties(0).total_memory/1e9:.1f} GB")
print("=" * 60)

MODEL = "AIRI-Institute/gena-lm-bert-base-t2t"

print("\n[1/3] Loading tokenizer + model (with output_hidden_states=True)...")
tok = AutoTokenizer.from_pretrained(MODEL)
config = AutoConfig.from_pretrained(MODEL, trust_remote_code=True)
config.output_hidden_states = True
model = AutoModel.from_pretrained(MODEL, trust_remote_code=True, config=config).to("cuda").eval()
print(f"  vocab size: {tok.vocab_size}")
print(f"  model class: {type(model).__name__}")
print(f"  num hidden layers: {config.num_hidden_layers}")
print(f"  hidden size: {config.hidden_size}")
print(f"  parameters: {sum(p.numel() for p in model.parameters())/1e6:.1f}M")

print("\n[2/3] Forward pass on 1000 bp dummy DNA...")
seq = "ATCG" * 250  # 1000 nucleotides
inputs = tok(seq, return_tensors="pt").to("cuda")
print(f"  input_ids shape: {inputs['input_ids'].shape}")

with torch.no_grad():
    out = model(**inputs)

print(f"\n[3/3] Inspecting output...")
print(f"  output type: {type(out).__name__}")
# hidden_states is a tuple of (n_layers+1,) tensors — embedding + each layer
print(f"  hidden_states: tuple of {len(out.hidden_states)} tensors")
print(f"  hidden_states[0].shape (embedding): {out.hidden_states[0].shape}")
print(f"  hidden_states[-1].shape (last layer): {out.hidden_states[-1].shape}")
print(f"  hidden_states[-1].dtype: {out.hidden_states[-1].dtype}")
print(f"  GPU mem peak: {torch.cuda.max_memory_allocated()/1e9:.2f} GB")

# Compute mean-pooled embedding (used for E4 clustering)
mask = inputs["attention_mask"].unsqueeze(-1).float()
mean_emb = (out.hidden_states[-1] * mask).sum(1) / mask.sum(1)
print(f"\n  mean-pooled embedding shape: {mean_emb.shape}")

print("\nSanity check PASSED. GENA-LM ready for experiments E1, E4.")
