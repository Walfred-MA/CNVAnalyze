#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=200:00:00
#SBATCH --account=mchaisso_100
#SBATCH --partition=qcb


aldy profile /project/mchaisso_100/datasets/1kg_phase3_related/$1 > profiles/"$1".profile
