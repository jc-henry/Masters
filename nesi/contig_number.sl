#!/bin/bash -e

#SBATCH --account       project-code                 								# project-code
#SBATCH --job-name      blast_test_contig_number_taxon           		                                # job-name
#SBATCH --time          01:00:00                    								# hh:mm:ss time allocated
#SBATCH --mem           1500MB                               						        # memory allocated
#SBATCH --cpus-per-task 20                           								# run as an array
#SBATCH --array			1-3
#SBATCH --output         ./blast_test_contig_number_taxon.%A_%a.out
#SBATCH --error          ./blast_test_contig_number_taxon.%A_%a.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz

module purge
module load BLAST/2.10.0-GCC-9.2.0

blastp -query /nesi/nobackup/project-code/taxon/results/taxon_${SLURM_ARRAY_TASK_ID}_trinity/transdecoder/taxon_${SLURM_ARRAY_TASK_ID}_trinity_rep_contigs.fasta.transdecoder.pep \
-db /nesi/project/project-code/resources/blast_databases/Araport11_genes.201606.pep \
-out ./blastp_taxon_${SLURM_ARRAY_TASK_ID}_arabidopsis.outfmt6 \
-evalue 1e-5 -num_threads ${SLURM_CPUS_PER_TASK} -max_target_seqs 1 -max_hsps 1 -outfmt 6

blastp -query /nesi/nobackup/project-code/taxon/results/taxon_${SLURM_ARRAY_TASK_ID}_trinity/transdecoder/taxon_${SLURM_ARRAY_TASK_ID}_trinity_rep_contigs.fasta.transdecoder.pep \
-db /nesi/project/project-code/resources/blast_databases/Araport11_genes.201606.pep \
-out ./blastp_taxon_${SLURM_ARRAY_TASK_ID}_arabidopsis.outfmt6 \
-evalue 1e-5 -num_threads ${SLURM_CPUS_PER_TASK} -max_target_seqs 1 -max_hsps 1 -outfmt 6

blastp -query /nesi/nobackup/project-code/taxon/results/taxon_${SLURM_ARRAY_TASK_ID}_trinity/transdecoder/taxon_${SLURM_ARRAY_TASK_ID}_trinity_rep_contigs.fasta.transdecoder.pep \
-db /nesi/project/project-code/resources/blast_databases/Araport11_genes.201606.pep \
-out ./blastp_taxon_${SLURM_ARRAY_TASK_ID}_arabidopsis.outfmt6 \
-evalue 1e-5 -num_threads ${SLURM_CPUS_PER_TASK} -max_target_seqs 1 -max_hsps 1 -outfmt 6






