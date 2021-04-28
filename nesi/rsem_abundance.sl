#!/bin/bash -e

#SBATCH --account       project-code                 								# project-code
#SBATCH --job-name      rsem_abundance				                                # job-name
#SBATCH --time          18:00:00                    								# hh:mm:ss time allocated
#SBATCH --mem           40G                           							    # memory allocated
#SBATCH --cpus-per-task 16                           								# threads allocated
#SBATCH --array			1-3                                                         # run as an array
#SBATCH --output        /nesi/nobackup/project-code/path/to/rsem/taxon_rsem_abundance_%A_%a.output
#SBATCH --error         /nesi/nobackup/project-code/path/to/rsem/taxon_rsem_abundance_%A_%a.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/path/to/rsem/

module purge
module load Trinity/2.11.0-gimkl-2020a
module load Bowtie2/2.4.1-GCC-9.2.0
module load RSEM/1.3.3-gimkl-2020a
module load SAMtools/1.10-GCC-9.2.0

READ_DIR=/nesi/nobackup/project-code/taxon/path/to/reads/trimmed
ASS_DIR=/nesi/nobackup/project-code/taxon/path/to/trinity_assembly

time align_and_estimate_abundance.pl \
--transcripts $ASS_DIR/Trinity_taxon.fasta \
--seqType fq \
--left $READ_DIR/taxon_${SLURM_ARRAY_TASK_ID}_R1_val_1.fq.gz \
--right $READ_DIR/taxon_${SLURM_ARRAY_TASK_ID}_R2_val_2.fq.gz \
--est_method RSEM \
--output_dir ./taxon_${SLURM_ARRAY_TASK_ID}.out \
--aln_method bowtie2 \
--SS_lib_type RF \
--thread_count ${SLURM_CPUS_PER_TASK} \
--trinity_mode

