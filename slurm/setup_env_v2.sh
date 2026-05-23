#!/bin/bash
set -e
source /opt/conda/etc/profile.d/conda.sh
conda activate gena

# Use /scratch for tmp (3.5TB available, /tmp only 2GB)
mkdir -p /scratch/ekorshunov/tmp
export TMPDIR=/scratch/ekorshunov/tmp
export PIP_CACHE_DIR=/scratch/ekorshunov/pip-cache

which python
python --version

echo "=== pip install torch ==="
pip install --no-cache-dir torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121 2>&1 | tail -3

echo "=== pip install transformers + ML stack ==="
pip install --no-cache-dir transformers==4.36.2 datasets accelerate scikit-learn umap-learn matplotlib seaborn pandas numpy biopython loguru tqdm 2>&1 | tail -3

echo "=== DONE ==="
python -c "import torch, transformers, sklearn, umap; print('torch:', torch.__version__); print('transformers:', transformers.__version__); print('sklearn:', sklearn.__version__); print('umap:', umap.__version__)"
