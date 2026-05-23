#!/bin/bash
#SBATCH --job-name=gena-e3
#SBATCH --gres=gpu:1
#SBATCH --constraint=gpu
#SBATCH --partition=short
#SBATCH --time=01:30:00
#SBATCH --output=/home/ekorshunov/gena-lm-airi/logs/e3_%j.out
#SBATCH --error=/home/ekorshunov/gena-lm-airi/logs/e3_%j.err

source /opt/conda/etc/profile.d/conda.sh
conda activate gena
export HF_HOME=/home/ekorshunov/hf-cache
cd /home/ekorshunov/gena-lm-airi
echo "node: $(hostname)"
echo "started: $(date)"
python -u experiments/e3_gena_vs_caduceus.py --n 500 --bs 8
echo "finished: $(date)"
