#!/bin/bash -e

#SBATCH --account       project-code                 						# project-code
#SBATCH --job-name      fastqc_raw           						        # job-name
#SBATCH --time          00:45:00                    					        # hh:mm:ss time allocated 
#SBATCH --mem           2100MB                         						# memory allocated
#SBATCH --cpus-per-task 6                           						# 1 thread per file.
#SBATCH --output        ./fastqc_raw.%j.output
#SBATCH --error         ./fastqc_raw.%j.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/taxon/

module purge
module load FastQC/0.11.9

time fastqc -f fastq -o ./results/fastqc/raw_reads -t ${SLURM_CPUS_PER_TASK} \
./path/to/reads/raw/taxon_1_F.fq.gz \
./path/to/reads/raw/taxon_1_R.fq.gz \
./path/to/reads/raw/taxon_2_F.fq.gz \
./path/to/reads/raw/taxon_2_R.fq.gz \
./path/to/reads/raw/taxon_3_F.fq.gz \
./path/to/reads/raw/taxon_3_R.fq.gz


