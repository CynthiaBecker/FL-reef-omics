#! /bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=fastqc
#SBATCH --mail-type=ALL			# When you get mail
#SBATCH --mail-user=cbecker@whoi.edu	# where the mail is sent
#SBATCH --ntasks=1			# Max 85 compute nodes. Number of CPUs or how many "tasks"; nodes??
#SBATCH --cpus-per-task=6		# Max 36, number of CPU cores per ntask for MULTITHREADING
#SBATCH --mem=100G			# 192gb limit - Job memory request, with "gb" for gigabyte
#SBATCH --time=4:00:00			# hr:min:sec
#SBATCH --output=logs/fastqc_%j.log	# standard output saved to what file
#export OMP_NUM_THREADS=1

## usage from FLK2019NextSeq folder: sbatch scripts/fastqc.sh 

## The qc conda environment must be active for this script to run
## If you have not created or activated it, please do the following:
## Create the environment:
## $ conda env create -f envs/fastqc.yml
## Activate the environment:
## $ conda activate qc

fastqc RawFastq/*fastq.gz -o output/fastqc/

## -o specifies the output folder. Must be present prior to running command
