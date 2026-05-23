#!/bin/bash
#SBATCH --job-name=gena-e5
#SBATCH --gres=gpu:1
#SBATCH --constraint=gpu
#SBATCH --partition=short
#SBATCH --time=01:30:00
#SBATCH --output=/home/ekorshunov/gena-lm-airi/logs/e5_%j.out
#SBATCH --error=/home/ekorshunov/gena-lm-airi/logs/e5_%j.err

source /opt/conda/etc/profile.d/conda.sh
conda activate gena
export HF_HOME=/home/ekorshunov/hf-cache
cd /home/ekorshunov/gena-lm-airi
echo "node: $(hostname)"
echo "started: $(date)"
python -u experiments/e5_saturation_mutagenesis.py --n 100 --bs 32
echo "finished: $(date)"
