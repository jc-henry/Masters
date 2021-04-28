#!/bin/bash -e

#SBATCH --account       project-code                  		# project-code
#SBATCH --job-name      taxon_x_salmon                      # job-name
#SBATCH --time          00:45:00                     		# hh:mm:ss time allocated
#SBATCH --mem           3G                         	        # memory allocated
#SBATCH --cpus-per-task 6                            		# would potentially be faster with more threads
#SBATCH --output        ./taxon_x_salmon.%j.output
#SBATCH --error         ./taxon_x_salmon.%j.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/taxon/path/to/taxon_x_trinity/salmon/

module purge
module load Salmon/1.3.0-gimkl-2020a

time salmon index --index ./salmon_index_taxon_x \
--transcripts /nesi/nobackup/project-code/taxon/path/to/taxon_x_trinity/trinity_assembly/Trinity_taxon_x.fasta -k 31

READ_DIR=/nesi/nobackup/project-code/taxon/path/to/reads/trimmed/

time salmon quant --index ./salmon_index_taxon_x --libType ISR \
-1 $READ_DIR/taxon_x_R1_val_1.fq.gz \
-2 $READ_DIR/taxon_x_R2_val_2.fq.gz \
-p ${SLURM_CPUS_PER_TASK} --gcBias --validateMappings --dumpEq taxon_x \
--output ./taxon_x.out

