#!/bin/bash -e

#SBATCH --account       project-code                  		# project-code
#SBATCH --job-name      hutt_monroi_1_rnaspades_corset      # job-name
#SBATCH --time          00:20:00                     		# hh:mm:ss time allocated
#SBATCH --mem           600MB                            	# memory allocated
#SBATCH --cpus-per-task 1                            		# single-threaded
#SBATCH --output        ./hutt_monroi_1_corset.%j.output
#SBATCH --error         ./hutt_monroi_1_corset.%j.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/path/to/hutt_monroi_1_rnaspades/corset/

module purge
module load Corset/1.09-GCC-9.2.0

time corset -m 0 -n hutt_monroi_1 -p hutt_monroi_1 -i salmon_eq_classes ../salmon_corset/*.out/aux_info/eq_classes.txt

# -m 0 means Corset won't remove transcripts with fewer than 10 reads mapping

