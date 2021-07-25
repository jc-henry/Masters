#!/bin/bash -e

#SBATCH --account       project-code                 								# project-code
#SBATCH --job-name      rcorrector           						        	        # job-name
#SBATCH --time          24:00:00                    								# hh:mm:ss time allocated
#SBATCH --mem           21G                           								# memory allocated
#SBATCH --cpus-per-task 8                           								# threads allocated
#SBATCH --output        ./path/to/reads/rcorrected/rcorrector.%j.output
#SBATCH --error         ./path/to/reads/rcorrected/rcorrector.%j.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/taxon/

module purge
module load Rcorrector/1.0.4-gimkl-2020a

time run_rcorrector.pl \
-1 ./path/to/reads/raw/taxon_1_F.fq.gz,./path/to/reads/raw/taxon_2_F.fq.gz,./path/to/reads/raw/taxon_3_F.fq.gz \
-2 ./path/to/reads/raw/taxon_1_R.fq.gz,./path/to/reads/raw/taxon_2_R.fq.gz,./path/to/reads/raw/taxon_3_R.fq.gz \
-k 23 -t ${SLURM_CPUS_PER_TASK} -od ./path/to/reads/rcorrected


