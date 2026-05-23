#!/bin/bash
#SBATCH --job-name=gena-sanity
#SBATCH --gres=gpu:1
#SBATCH --constraint=gpu
#SBATCH --partition=short
#SBATCH --time=00:20:00
#SBATCH --output=/home/ekorshunov/gena-lm-airi/logs/sanity_%j.out
#SBATCH --error=/home/ekorshunov/gena-lm-airi/logs/sanity_%j.err

source /opt/conda/etc/profile.d/conda.sh
conda activate gena
cd /home/ekorshunov/gena-lm-airi
echo "node: $(hostname)"
echo "started: $(date)"
python -u sanity_check.py
echo "finished: $(date)"
