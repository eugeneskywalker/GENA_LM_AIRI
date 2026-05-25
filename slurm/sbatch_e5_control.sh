#!/bin/bash
#SBATCH --job-name=gena-e5-ctrl
#SBATCH --gres=gpu:1
#SBATCH --constraint=gpu
#SBATCH --partition=short
#SBATCH --time=01:30:00
#SBATCH --output=/home/ekorshunov/gena-lm-airi/logs/e5_control_%j.out
#SBATCH --error=/home/ekorshunov/gena-lm-airi/logs/e5_control_%j.err

source /opt/conda/etc/profile.d/conda.sh
conda activate gena
export HF_HOME=/home/ekorshunov/hf-cache
cd /home/ekorshunov/gena-lm-airi
echo "node: $(hostname)"
echo "started: $(date)"

# Negative control: mutate MOTIF_LEN random positions OUTSIDE the CTCF motif.
# Expected: overall_mean_effect for control ≪ motif-mode value (≈ in e5_metrics.json)
# if GENA-LM has internalised CTCF; comparable if effect is BPE-tokenisation artefact.
python -u experiments/e5_saturation_mutagenesis.py --n 100 --bs 32 --control

echo "finished: $(date)"
