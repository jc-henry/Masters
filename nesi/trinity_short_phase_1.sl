#!/bin/bash -e

#SBATCH --job-name			trinity-phase1
#SBATCH --account			project-code   						# project-code
#SBATCH --time				12:00:00       						# hh:mm:ss time allocated
#SBATCH --ntasks			1            						# always 1
#SBATCH --cpus-per-task		18   		   	   				    # number of threads to use for Trinity 
#SBATCH --mem				170G            					# memory allocated
#SBATCH --hint				nomultithread  						# disable hyper-threading
#SBATCH --output        	./path/to/trinity_assembly/trinity-phase1.%j.output
#SBATCH --error         	./path/to/trinity_assembly/trinity-phase1.%j.err
#SBATCH --mail-type     	END
#SBATCH --mail-user     	j.henry2@massey.ac.nz
#SBATCH --chdir         	/nesi/nobackup/project-code/taxon/

module purge
module load Trinity/2.11.0-gimkl-2020a

# run trinity, stop before phase 2

srun Trinity --no_distributed_trinity_exec \
--CPU ${SLURM_CPUS_PER_TASK} \
--max_memory 160G \
--seqType fq \
--left ./path/to/reads/trimmed/taxon_1_R1_val_1.fq.gz,./path/to/reads/trimmed/taxon_2_R1_val_1.fq.gz,./path/to/reads/trimmed/taxon_3_R1_val_1.fq.gz \
--right .path/to/reads/trimmed/taxon_1_R2_val_2.fq.gz,./path/to/reads/trimmed/taxon_2_R2_val_2.fq.gz,./path/to/reads/trimmed/taxon_3_R2_val_2.fq.gz \
--output ./path/to/trinity_assembly/ \
--SS_lib_type RF \
--min_contig_length 300 \
--include_supertranscripts 

  
# max_memory should be same (or maybe slightly lower, so you have a small buffer) than the value specified with the sbatch option --mem

