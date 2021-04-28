#!/bin/bash -e

#SBATCH --account       project-code                 								# project-code
#SBATCH --job-name      rsem_contigs				                                # job-name
#SBATCH --time          00:03:00                    								# hh:mm:ss time allocated
#SBATCH --mem           1G                           							    # memory allocated
#SBATCH --cpus-per-task 1                           								# single threaded
#SBATCH --output        /nesi/nobackup/project-code/taxon/path/to/trinity_assembly/rsem_contigs.%j.out
#SBATCH --error         /nesi/nobackup/project-code/taxon/path/to/trinity_assembly/rsem_contigs.%j.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/taxon/path/to/trinity_assembly/

module purge
module load Trinity/2.11.0-gimkl-2020a

ASS_DIR=/nesi/nobackup/project-code/taxon/path/to/trinity_assembly
MAT_DIR=/nesi/nobackup/project-code/taxon/path/to/rsem

time filter_low_expr_transcripts.pl \
--matrix $MAT_DIR/taxon_rsem_contigs_TPM.txt \
--transcripts $ASS_DIR/Trinity_taxon.fasta \
--highest_iso_only \
--trinity_mode > high_expression_contigs_taxon.fasta


