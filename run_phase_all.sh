#!/bin/bash
#SBATCH -D /home/wang9418
#SBATCH --partition=med2
#SBATCH --job-name=run_phase_all
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=10:00:00
#SBATCH --output=/home/wang9418/family_variance/codes/logs/run_phase_all.out
#SBATCH --error=/home/wang9418/family_variance/codes/logs/run_phase_all.err
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=faywang@ucdavis.edu

module load jdk/
java -Xmx1g -jar beagle.05May22.33a.jar gt=~/family_variance/generate_vcffiles/genotype.vcf.gz out=~/family_variance/phased_data/phased seed=0 nthreads=2

