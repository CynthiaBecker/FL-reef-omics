#! /bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=quast
#SBATCH --mail-type=ALL			# When you get mail
#SBATCH --mail-user=cbecker@whoi.edu	# where the mail is sent
#SBATCH --ntasks=1			# Max 85 compute nodes. Number of CPUs or how many "tasks"; nodes??
#SBATCH --cpus-per-task=36		# Max 36, number of CPU cores per ntask for MULTITHREADING
#SBATCH --mem=100G			# 192gb limit - Job memory request, with "gb" for gigabyte
#SBATCH --time=10:00:00			# hr:min:sec
#SBATCH --output=logs/quast_%j.log	# standard output saved to what file
#export OMP_NUM_THREADS=1

## usage from FLK2019NextSeq folder: sbatch scripts/quast.sh 

## The quast conda environment must be active for this script to run
## If you have not created or activated it, please do the following:
## Create the environment:
## $ conda env create -f envs/quast.yml
## Activate the environment:
## $ conda activate quast

python quast.py -o output/quast_megahit -t 36 FLK2019_assembly2/final.contigs.fa

## -o specifies the output folder
## The metagenome assembly (final.contigs.fa_) is the thing quast is evaluating the metrics of.

