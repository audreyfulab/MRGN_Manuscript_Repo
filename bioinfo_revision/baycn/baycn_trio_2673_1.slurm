#!/bin/bash
#SBATCH --job-name baycn_2673_1
#SBATCH -q primary
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ht2699@wayne.edu
#SBATCH -o output_%j.out
#SBATCH -e errors_%j.err
#SBATCH -t 04:00:00

module purge
module load gnu9/9.1.0
module load r/4.5.0

Rscript GTEx_wholeblood_baycn_save_data.R 1e6 400 2673 1

