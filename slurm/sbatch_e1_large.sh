#!/bin/bash
#SBATCH --job-name=gena-e1L
#SBATCH --gres=gpu:1
#SBATCH --constraint=gpu
#SBATCH --partition=short
#SBATCH --time=01:30:00
#SBATCH --output=/home/ekorshunov/gena-lm-airi/logs/e1_large_%j.out
#SBATCH --error=/home/ekorshunov/gena-lm-airi/logs/e1_large_%j.err

source /opt/conda/etc/profile.d/conda.sh
conda activate gena
export HF_HOME=/home/ekorshunov/hf-cache
cd /home/ekorshunov/gena-lm-airi
echo "node: $(hostname)"
echo "started: $(date)"
python -u experiments/e1_bert_large.py --n 300 --bs 4
echo "finished: $(date)"
