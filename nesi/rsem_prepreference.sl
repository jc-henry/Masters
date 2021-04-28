#!/bin/bash -e

#SBATCH --account       project-code                 								# project-code
#SBATCH --job-name      rsem_prepreference				                            # job-name
#SBATCH --time          00:45:00                    								# hh:mm:ss time allocated
#SBATCH --mem           3G                           							    # memory allocated
#SBATCH --cpus-per-task 1                           								# appears single threaded
#SBATCH --output        /nesi/nobackup/project-code/path/to/rsem/taxon_rsem_prepreference.%j.output
#SBATCH --error         /nesi/nobackup/project-code/path/to/rsem/taxon_rsem_prepreference.%j.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/path/to/rsem/

module purge
module load Trinity/2.11.0-gimkl-2020a
module load Bowtie2/2.4.1-GCC-9.2.0
module load RSEM/1.3.3-gimkl-2020a
  
ASS_DIR=/nesi/nobackup/project-code/taxon/path/to/trinity_assembly

time align_and_estimate_abundance.pl \
--transcripts $ASS_DIR/Trinity_taxon.fasta \
--est_method RSEM \
--aln_method bowtie2 \
--thread_count ${SLURM_CPUS_PER_TASK} \
--trinity_mode \
--prep_reference

# the index and gene to transcript map will all be made in the "trinity_assembly" dir (not the rsem dir)
