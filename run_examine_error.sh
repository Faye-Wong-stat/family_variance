#!/bin/bash
#SBATCH -D /home/wang9418
#SBATCH --partition=med2
#SBATCH --job-name=examine_error
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/home/wang9418/family_variance/codes/logs/run_examine_error.out.txt
#SBATCH --error=/home/wang9418/family_variance/codes/logs/run_examine_error.err.txt
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=faywang@ucdavis.edu

module load R/4.2.2
Rscript family_variance/codes/examine_error.R 