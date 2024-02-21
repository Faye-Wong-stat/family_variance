#!/bin/bash
#SBATCH -D /home/wang9418
#SBATCH --partition=med
#SBATCH --job-name=get_phase_error
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --output=/home/wang9418/family_variance/codes/logs/run_get_phase_error.out
#SBATCH --error=/home/wang9418/family_variance/codes/logs/run_get_phase_error.err
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=faywang@ucdavis.edu

module load R/4.2.2
Rscript family_variance/codes/get_phase_error.R