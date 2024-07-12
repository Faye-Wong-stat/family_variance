#!/bin/bash
#SBATCH -D /home/wang9418
#SBATCH --partition=med2
#SBATCH --job-name=select_best_parents_gametes
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/home/wang9418/family_variance/codes/logs/run_select_best_parents_gametes.out.txt
#SBATCH --error=/home/wang9418/family_variance/codes/logs/run_select_best_parents_gametes.err.txt
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=faywang@ucdavis.edu

source ~/.bashrc
conda activate strawberry
Rscript family_variance/codes/select_best_parents_gametes.R 