#!/bin/bash -e

#SBATCH --account       project-code                  		# project-code
#SBATCH --job-name      trimgalore   	        	        # job-name
#SBATCH --time          07:00:00                     		# hh:mm:ss time allocated
#SBATCH --mem           80MB                                    # memory allocated
#SBATCH --cpus-per-task 1                            		# reads are already uncompressed so a single thread is sufficient
#SBATCH --array=1-3                                             # run as an array
#SBATCH --output        ./path/to/reads/trimmed/taxon_trimgalore_%A_%a.out
#SBATCH --error         ./path/to/reads/trimmed/taxon_trimgalore_%A_%a.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/taxon/

module purge
module load TrimGalore/0.6.4-gimkl-2018b

time trim_galore --length 50 --paired --basename taxon_${SLURM_ARRAY_TASK_ID} \
-o ./path/to/reads/trimmed --retain_unpaired --gzip \
./path/to/reads/rcorrected/unfixrm_taxon_${SLURM_ARRAY_TASK_ID}_F.cor.fq \
./path/to/reads/rcorrected/unfixrm_taxon_${SLURM_ARRAY_TASK_ID}_R.cor.fq 

# trims at default phred of 20
