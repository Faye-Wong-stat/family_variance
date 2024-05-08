#!/bin/bash
#SBATCH -D /home/wang9418/family_variance
#SBATCH --partition=med2
#SBATCH --job-name=run_simulate_crosses
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/home/wang9418/family_variance/codes/logs/run_simulate_crosses.out
#SBATCH --error=/home/wang9418/family_variance/codes/logs/run_simulate_crosses.err
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=faywang@ucdavis.edu

module load R/4.2.2
Rscript codes/simulate_crosses.R 

