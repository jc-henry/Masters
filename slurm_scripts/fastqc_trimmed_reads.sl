#!/bin/bash -e

#SBATCH --account       project-code                 						        # project-code
#SBATCH --job-name      fastqc_trimmed           					                # job-name
#SBATCH --time          00:45:00                    							# hh:mm:ss time allocated 
#SBATCH --mem           2100MB                         							# memory allocated
#SBATCH --cpus-per-task 6                           						        # 1 thread per file.
#SBATCH --output        ./fastqc_trimmed.%j.output
#SBATCH --error         ./fastqc_trimmed.%j.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/taxon/

module purge
module load FastQC/0.11.9

time fastqc -f fastq -o ./results/fastqc/trimmed_reads -t ${SLURM_CPUS_PER_TASK} \
./path/to/reads/trimmed/taxon_1_R1_val_1.fq.gz \
./path/to/reads/trimmed/taxon_1_R2_val_2.fq.gz \
./path/to/reads/trimmed/taxon_2_R1_val_1.fq.gz \
./path/to/reads/trimmed/taxon_2_R2_val_2.fq.gz \
./path/to/reads/trimmed/taxon_3_R1_val_1.fq.gz \
./path/to/reads/trimmed/taxon_3_R2_val_2.fq.gz


