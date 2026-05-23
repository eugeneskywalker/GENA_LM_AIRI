#!/bin/bash
#SBATCH --job-name=gena-e4v2
#SBATCH --gres=gpu:1
#SBATCH --constraint=gpu
#SBATCH --partition=short
#SBATCH --time=01:30:00
#SBATCH --output=/home/ekorshunov/gena-lm-airi/logs/e4v2_%j.out
#SBATCH --error=/home/ekorshunov/gena-lm-airi/logs/e4v2_%j.err

source /opt/conda/etc/profile.d/conda.sh
conda activate gena
export HF_HOME=/home/ekorshunov/hf-cache
cd /home/ekorshunov/gena-lm-airi
echo "node: $(hostname)"
echo "started: $(date)"
python -u experiments/e4_clustering_multilayer.py --n 500 --bs 16
echo "finished: $(date)"
