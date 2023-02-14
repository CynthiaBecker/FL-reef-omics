#! /bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=trimming
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbecker@whoi.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=150gb
#SBATCH --time=10:00:00
#SBATCH --output=trimmed_%j.log
export OMP_NUM_THREADS=4

# usage: sbatch trim.sh
# submit from FLK2019NextSeq folder
# conda environment with trimmomatic must be active. If not yet created and activated, you can set it up with the following code using the yaml file in the envs folder
# conda env create -f envs/trim.yml
# conda activat trim

cd RawFastq

for infile in *_R1_001.fastq.gz
do
   base=$(basename ${infile} _R1_001.fastq.gz)
   trimmomatic PE -version -threads 4 ${infile} ${base}_R2_001.fastq.gz \
                ${base}_1.trim.fastq.gz ${base}_1.untrim.fastq.gz \
                ${base}_2.trim.fastq.gz ${base}_2.untrim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:50
done

## Sliding window trimming means over an average of 4 bases, if quality drops below 25, the read will be trimmed there
## Minlen means if a read is below 100 bp it will be dropped. 
## -version prints the version
## -threads will parallelize the trimming over 4 tasks. I think the multithreading is enabled with the cpus-per-task and the OMP_NUM_THREADS part

