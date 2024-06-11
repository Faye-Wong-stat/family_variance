#!/bin/bash
#SBATCH -D /home/wang9418/family_variance
#SBATCH --partition=med
#SBATCH --job-name=simulate_crosses_best_parents_gametes
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --output=/home/wang9418/family_variance/codes/logs/run_simulate_crosses_best_parents_gametes-%A_%a.out.txt
#SBATCH --error=/home/wang9418/family_variance/codes/logs/run_simulate_crosses_best_parents_gametes-%A_%a.err.txt
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=faywang@ucdavis.edu
#SBATCH --array=1-300



echo $SLURM_ARRAY_TASK_ID

name=($(head -n $SLURM_ARRAY_TASK_ID select_best_parents_gametes/best_parent_mean.txt | tail -n 1 | tr -d '"'))
echo ${name[@]}

module load R/4.2.2
echo finished loading R
Rscript codes/simulate_crosses_best_parents_gametes.R ${name[@]}

# run_simulate_crosses_best_parents_gametes.sh