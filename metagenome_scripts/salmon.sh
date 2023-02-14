#! /bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=salmon
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbecker@whoi.edu
#SBATCH --ntasks=1				# number of tasks (in this case, 1)
#SBATCH --cpus-per-task=36			# number of CPU cores (max 36 for compute; 80 for bigmem) per ntask for multithreading
#SBATCH --mem=150gb				# max 3 tb for bigmem and 192 gb for compute 
#SBATCH --qos=unlim				# max usually 24 hr. To use my exception, do --qos=unlim
#SBATCH --time=12:00:00			# must have time. Here days-hours:min:sec
#SBATCH --output=logs/salmon_%j.log
#export OMP_NUM_THREADS=36

## usage from FLK2019NextSeq folder: sbatch scripts/salmon.sh

## NOTES: conda environment "salmon" must be active using `conda activate salmon`
## If environment isn't active yet, use `conda env create -f envs/salmon.yml` from the FLK2019NextSeq folder

cd /vortexfs1/home/cbecker/FLK2019NextSeq/output/salmonquant/

## salmon index -t BacteriaMG.ffn -i MG_index -k 31

for file in *_1.fastq.gz
do
tail1=_1.fastq.gz
tail2=_2.fastq.gz
BASE=${file/$tail1/}
salmon quant --meta -i MG_index --libType A \
        -1 $BASE$tail1 -2 $BASE$tail2 -o $BASE.quant \
        -p 36
done

## INFO ON FLAGS
## -k is length of kmers
## -i is index name
## --perfectHash: the --prefectHash cuts the construction memory for the hash function by ~40%. Good if you have a large metagenomic reference set
## Salmon quant info:
## --meta changes the slmon quant command to be more specific for metagenomic data rather than transcriptomic data
## --libType IU specifies an unstranded library that is "inward", meaning the reads face each other. 
##  --libType A tells the program to automatically decide the library type 
