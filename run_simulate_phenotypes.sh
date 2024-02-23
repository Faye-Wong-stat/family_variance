#!/bin/bash
#SBATCH -D /home/wang9418
#SBATCH --partition=med
#SBATCH --job-name=simulate_phenotypes
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --output=/home/wang9418/family_variance/codes/logs/simulate_phenotypes.out
#SBATCH --error=/home/wang9418/family_variance/codes/logs/simulate_phenotypes.err
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=faywang@ucdavis.edu

module load R/4.2.2
Rscript family_variance/codes/simulate_phenotypes.R
