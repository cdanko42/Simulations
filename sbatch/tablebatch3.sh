#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --exclusive
#SBATCH --ntasks=40
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=20
#SBATCH --mem=1024
#SBATCH --output=regress_stdout.txt
#SBATCH --error=regress_stderr.txt
#SBATCH --time=5:00
#SBATCH --job-name=table3
#SBATCH --mail-user=christopher.p.danko-1@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/danko/Simulations/R
#
#################################################
module load R
Rscript tablegen3.R > table3.txt