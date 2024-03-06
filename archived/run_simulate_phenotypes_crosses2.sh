#!/bin/bash
#SBATCH -D /home/wang9418/family_variance
#SBATCH --partition=med
#SBATCH --job-name=simulate_phenotypes_crosses2
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --output=/home/wang9418/family_variance/codes/logs/run_simulate_phenotypes_crosses2-%A_%a.out
#SBATCH --error=/home/wang9418/family_variance/codes/logs/run_simulate_phenotypes_crosses2-%A_%a.err
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=faywang@ucdavis.edu
#SBATCH --array=1-520

#while read name
#do 
# echo job_started
# echo $SGE_TASK_ID
echo $SLURM_ARRAY_TASK_ID
name=$(head -n $SLURM_ARRAY_TASK_ID simulate_crosses/file_names.txt | tail -n 1)
#name=$(sed -n -e "$SGE_TASK_ID p" simulate_crosses/file_names.txt)
echo $name
# echo job_done
module load R/4.2.2
Rscript codes/simulate_phenotypes_crosses2.R $name
#done < simulated_data/file_names.txt
