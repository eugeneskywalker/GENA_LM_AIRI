#!/bin/bash
set -e
source /opt/conda/etc/profile.d/conda.sh
conda activate solo_prac

echo "=== [1/4] mamba create gena ==="
mamba create -n gena python=3.11 pip -c conda-forge -y

echo "=== [2/4] activate gena ==="
conda deactivate
conda activate gena
which python
python --version

echo "=== [3/4] pip install torch (with CUDA 12.1) ==="
pip install --no-cache-dir torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121

echo "=== [4/4] pip install transformers + ML stack ==="
pip install --no-cache-dir transformers==4.36.2 datasets accelerate scikit-learn umap-learn matplotlib seaborn pandas numpy biopython pybedtools pysam loguru tqdm

echo "=== DONE ==="
python -c "import torch, transformers, sklearn, umap; print('torch:', torch.__version__); print('transformers:', transformers.__version__); print('sklearn:', sklearn.__version__); print('umap:', umap.__version__)"
