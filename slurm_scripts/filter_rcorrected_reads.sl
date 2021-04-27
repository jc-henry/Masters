#!/bin/bash -e

#SBATCH --account       project-code                  		# project-code
#SBATCH --job-name      filter_rcorrected   	 	        # job-name
#SBATCH --time          02:00:00                     		# hh:mm:ss time allocated
#SBATCH --mem           50MB                            	# memory allocated
#SBATCH --cpus-per-task 1                            		# single threaded
#SBATCH --array		1-3                                     # run as an array
#SBATCH --output        filter_rcorrected_%A_%a.out
#SBATCH --error         filter_rcorrected_%A_%a.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/massey02848/project-code/path/to/reads/rcorrected/

module purge
module load Python/2.7.18-gimkl-2020a

time python /path/to/docs/FilterUncorrectabledPEfastq.py \
-1 ./taxon_${SLURM_ARRAY_TASK_ID}_F.cor.fq.gz \
-2 ./taxon_${SLURM_ARRAY_TASK_ID}_R.cor.fq.gz \
-s taxon_${SLURM_ARRAY_TASK_ID}_filtered


