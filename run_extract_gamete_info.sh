#!/bin/bash
#SBATCH -D /home/wang9418
#SBATCH --partition=med2
#SBATCH --job-name=extract_gamete_info
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/home/wang9418/family_variance/codes/logs/run_extract_gamete_info.out.txt
#SBATCH --error=/home/wang9418/family_variance/codes/logs/run_extract_gamete_info.err.txt
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=faywang@ucdavis.edu

module load R/4.2.2
Rscript family_variance/codes/extract_gamete_info.R 