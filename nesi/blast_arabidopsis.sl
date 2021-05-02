#!/bin/bash -e

#SBATCH --account       project-code                 								# project-code
#SBATCH --job-name      taxon_blast_arabidopsis            				            # job-name
#SBATCH --time          01:00:00                    								# hh:mm:ss time allocated
#SBATCH --mem           400MB                               						# memory allocated
#SBATCH --cpus-per-task 30                           								# threads allocated
#SBATCH --output        /nesi/nobackup/project-code/taxon/path/to/blast_arabidopsis/taxon_blast_arabidopsis.%j.output
#SBATCH --error         /nesi/nobackup/project-code/taxon/path/to/blast_arabidopsis/taxon_blast_arabidopsis.%j.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/taxon/

module purge
module load BLAST/2.10.0-GCC-9.2.0

time blastp -query ./path/to/transdecoder/high_expression_contigs_taxon.fasta.transdecoder.pep \
-db /nesi/project/project-code/resources/blast_databases/Araport11_genes.201606.pep \
-out ./path/to/blast_arabidopsis/blastp_taxon_arabidopsis.outfmt6 \
-evalue 1e-5 -num_threads ${SLURM_CPUS_PER_TASK} -max_target_seqs 5 -max_hsps 1 -outfmt 6




