#!/bin/bash
#SBATCH --job-name=gena-e7
#SBATCH --gres=gpu:1
#SBATCH --constraint=gpu
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --output=/home/ekorshunov/gena-lm-airi/logs/e7_%j.out
#SBATCH --error=/home/ekorshunov/gena-lm-airi/logs/e7_%j.err

source /opt/conda/etc/profile.d/conda.sh
conda activate gena
export HF_HOME=/home/ekorshunov/hf-cache
cd /home/ekorshunov/gena-lm-airi
echo "node: $(hostname)"
echo "started: $(date)"
python -u experiments/e7_fine_tune_promoter.py --n 2000 --bs 16 --epochs 3 --lr 2e-5
echo "finished: $(date)"
