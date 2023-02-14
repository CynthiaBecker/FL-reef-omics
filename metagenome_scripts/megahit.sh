#! /bin/bash

#SBATCH --partition=bigmem
#SBATCH --job-name=megahit
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbecker@whoi.edu
#SBATCH --ntasks=1				# number of tasks (in this case, 1)
#SBATCH --cpus-per-task=60			# number of CPU cores (max 36) per ntask for multithreading
#SBATCH --mem=2000gb				# max 3 tb i think 
#SBATCH --time=24:00:00				# max 24 hr, hr:min:sec
#SBATCH --output=logs/megahit_%j.log
#export OMP_NUM_THREADS=60

## usage from FLK2019NextSeq folder: sbatch scripts/megahit.sh

## NOTES: conda environment "megahit" must be active before running the command
## install the environment with the megahit.yml file
## $ conda env create -f envs/megahit.yml
## Activate with `conda activate megahit` #megahit v1.2.9

## make variables that are comma-separated lists of F and R reads
R1=$(cat data/trimmedFastq/r1)
R2=$(cat data/trimmedFastq/r2)

echo $R1 ##make sure look correct
echo $R2 

megahit -1 $R1 -2 $R2 --min-contig-len 500 --continue -o FLK2019_assembly2 -t 60

## --min-contig-len is the minimum length of the contigs. The meren lab tutorial used a length of 1000, but I will do half that, at 500. I believe the default is 200
## -o is the output folder that megahit creates
## -t is the number of threads to use, so number of cpus available on the computer system 
## --continue Was added because the large coassembly timed out. I am hoping this will allow it to restart in the same spot. 

