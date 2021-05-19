#!/bin/bash -e

#SBATCH --job-name			hutt_monroi_1_rnaspades
#SBATCH --account			project-code   						# project-code
#SBATCH --time				6:00:00       						# hh:mm:ss time allocated
#SBATCH --cpus-per-task		12   		   	   				    # number of threads allocated
#SBATCH --mem				20G            						# memory allocated
#SBATCH --output        	./rnaspades_assembly/rnaspades.%j.output
#SBATCH --error         	./rnaspades_assembly/rnaspades.%j.err
#SBATCH --mail-type     	END
#SBATCH --mail-user     	j.henry2@massey.ac.nz
#SBATCH --chdir         	/nesi/nobackup/project-code/path/to/hutt_monroi_1_rnaspades/

module purge
module load SPAdes/3.14.0-gimkl-2020a

ASS_DIR=/nesi/nobackup/project-code/path/to/reads/trimmed

rnaspades.py \
--pe-1 1 $ASS_DIR/hutt_monroi_1_R1_val_1.fq --pe-2 1 $ASS_DIR/hutt_monroi_1_R2_val_2.fq \
-o ./rnaspades_assembly/ \
--ss rf \
--threads ${SLURM_CPUS_PER_TASK} \
--memory ${SLURM_MEM_PER_NODE}



  

