#!/bin/bash
#SBATCH -D /home/wang9418/family_variance
#SBATCH --partition=med
#SBATCH --job-name=simulate_phenotypes_gametes
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --output=/home/wang9418/family_variance/codes/logs/run_simulate_phenotypes_gametes-%A_%a.out.txt
#SBATCH --error=/home/wang9418/family_variance/codes/logs/run_simulate_phenotypes_gametes-%A_%a.err.txt
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=faywang@ucdavis.edu
#SBATCH --array=1-10000



echo $SLURM_ARRAY_TASK_ID

name=$(head -n $SLURM_ARRAY_TASK_ID simulate_gametes/file_names.txt | tail -n 1)
echo $name

source ~/.bashrc
conda activate strawberry
Rscript codes/simulate_phenotypes_gametes.R $name


