#!/bin/bash -e

#SBATCH --account       project-code                 								# project-code
#SBATCH --job-name      taxon_transdecoder         						            # job-name
#SBATCH --time          00:20:00                    								# hh:mm:ss time allocated
#SBATCH --mem           1500MB                               						# memory allocated
#SBATCH --cpus-per-task 1                           								# single-threaded process
#SBATCH --output        ./taxon_transdecoder.%j.output
#SBATCH --error         ./taxon_transdecoder.%j.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/taxon/path/to/transdecoder/

module purge
module load TransDecoder/5.5.0-GCC-9.2.0-Perl-5.30.1

ASSEMBLY_DIR=/nesi/nobackup/project-code/taxon/path/to/trinity_assembly

time TransDecoder.LongOrfs -t $ASSEMBLY_DIR/high_expression_contigs_taxon.fasta -S -m 100 \
--output_dir ./

time TransDecoder.Predict -t $ASSEMBLY_DIR/high_expression_contigs_taxon.fasta \
--single_best_only --output_dir ./

# Note strand specific flag (-S)

