#!/bin/bash -e

#SBATCH --account       project-code                 								# project-code
#SBATCH --job-name      taxon_busco           						                # job-name
#SBATCH --time          24:00:00                    								# hh:mm:ss time allocated
#SBATCH --mem           2G                          							    # memory allocated
#SBATCH --cpus-per-task 1                           								# appears single threaded
#SBATCH --output        ./taxon_busco.%j.output
#SBATCH --error         ./taxon_busco.%j.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/taxon/path/to/busco/

module purge
module load BUSCO/4.1.4-gimkl-2020a

time busco --in /nesi/nobackup/project-code/taxon/path/to/trinity_assembly/Trinity_taxon.fasta \
--cpu ${SLURM_CPUS_PER_TASK} \
--out taxon_trinity --out_path ./ \
--mode transcriptome  --lineage_dataset embryophyta_odb10 --force --quiet



