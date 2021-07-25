#!/bin/bash -e

#SBATCH --account       project-code                 								# project-code
#SBATCH --job-name      taxon_sceleratus_bowtie2_local          		                                # job-name
#SBATCH --time          01:00:00                    								# hh:mm:ss time allocated
#SBATCH --mem           1G                           							        # memory allocated 
#SBATCH --cpus-per-task 12                           								# CPUs allocated
#SBATCH --array			1-3										# run as an array
#SBATCH --output        ./taxon_sceleratus_bowtie2_local.%A_%a.out
#SBATCH --error         ./taxon_sceleratus_bowtie2_local.%A_%a.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz


# bowtie2 index already created using:
#time bowtie2-build --threads ${SLURM_CPUS_PER_TASK} ./ranunculus_sceleratus_chloroplast.fasta bowtie2_sceleratus

module purge
module load Bowtie2/2.3.5-GCC-7.4.0
module load SAMtools/1.10-GCC-9.2.0

DIR=/nesi/nobackup/project-code/taxon/data/reads/short/trimmed

time bowtie2 -q --very-sensitive-local --no-unal --threads ${SLURM_CPUS_PER_TASK} \
-x ./bowtie2_sceleratus \
-1 $DIR/taxon_${SLURM_ARRAY_TASK_ID}_R1_val_1.fq.gz \
-2 $DIR/taxon_${SLURM_ARRAY_TASK_ID}_R2_val_2.fq.gz \
| samtools view -b - | samtools sort -o - > taxon_${SLURM_ARRAY_TASK_ID}_sceleratus_sorted.bam


