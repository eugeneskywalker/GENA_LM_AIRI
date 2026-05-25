#!/bin/bash
#SBATCH --job-name=gena-e6
#SBATCH --gres=gpu:1
#SBATCH --constraint=gpu
#SBATCH --partition=short
#SBATCH --time=01:30:00
#SBATCH --output=/home/ekorshunov/gena-lm-airi/logs/e6_%j.out
#SBATCH --error=/home/ekorshunov/gena-lm-airi/logs/e6_%j.err

source /opt/conda/etc/profile.d/conda.sh
conda activate gena
export HF_HOME=/home/ekorshunov/hf-cache
cd /home/ekorshunov/gena-lm-airi
echo "node: $(hostname)"
echo "started: $(date)"
python -u experiments/e6_dart_eval_style.py --n 2500 --bs 8
echo "finished: $(date)"
