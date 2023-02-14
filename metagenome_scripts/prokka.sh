#! /bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=prokka
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbecker@whoi.edu
#SBATCH --ntasks=1				# number of tasks (in this case, 1)
#SBATCH --cpus-per-task=36			# number of CPU cores (max 36 for compute; 80 for bigmem) per ntask for multithreading
#SBATCH --mem=150gb				# max 3 tb for bigmem and 192 gb for compute 
#SBATCH --qos=unlim				# max usually 24 hr. To use my exception, do --qos=unlim
#SBATCH --time=5-12:00:00			# must have time. Here days-hours:min:sec
#SBATCH --output=logs/prokka_%j.log
#export OMP_NUM_THREADS=36

## usage from FLK2019NextSeq folder: sbatch scripts/prokka.sh

## NOTES: conda environment "prokka" must be active using `conda activate prokka`
## If environment isn't active yet, use `conda env create -f envs/prokka.yml` from the FLK2019NextSeq folder

prokka FLK2019_assembly2/final.contigs.1000plus.fa --outdir output/prokka3 --prefix BacteriaMG --norrna --notrna --metagenome --cpus 36

## prokka FLK2019_assembly2/final.contigs.1000plus.fa --kingdom Archaea --outdir output/prokka --prefix ArchaeaMG --norrna --notrna --cpus 36

## --kingdom [X]      Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
## --outdir specify the directory to put output files. It will create the filename specified
## --prefix [X]       Filename output prefix [auto] (default '') 
## --metagenome       Improve gene predictions for highly fragmented genomes (default OFF)
## Be sure to include the contigs file that has been subset (mine was subset to only include contigs of 1000+ bp
## --norrna and --notrna Tells Prokka to not run  the rRNA or tRNA searches

