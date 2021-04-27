#!/bin/bash -e
#SBATCH --job-name					trinity-phase2_grid
#SBATCH --account					project-code  				# project-code
#SBATCH --time						05:00:00      				# hh:mm:ss time allocated
#SBATCH --ntasks					1           				# always 1 - this is the master process
#SBATCH --cpus-per-task				1    					    # always 1
#SBATCH --mem						6G            				# memory requirements for master process
#SBATCH --hint						nomultithread				# disable hyper-threading
#SBATCH --output        	        ./path/to/trinity_assembly/trinity-phase1.%j.output
#SBATCH --error         	        ./path/to/trinity_assembly/trinity-phase1.%j.err
#SBATCH --mail-type     			END
#SBATCH --mail-user     			j.henry2@massey.ac.nz
#SBATCH --chdir         			/nesi/nobackup/project-code/taxon/

module purge

# load Trinity and HPC GridRunner
module load Trinity/2.11.0-gimkl-2020a
module load HpcGridRunner/20181005

# run Trinity - this will be the master HPC GridRunner process that handles
#   submitting sub-jobs (batches of commands) to the Slurm queue
srun Trinity --CPU ${SLURM_CPUS_PER_TASK} --max_memory 6G \
 --grid_exec "hpc_cmds_GridRunner.pl --grid_conf ${SLURM_SUBMIT_DIR}/SLURM.conf -c" \
--seqType fq \
--left ./path/to/reads/trimmed/taxon_1_R1_val_1.fq.gz,./path/to/reads/trimmed/taxon_2_R1_val_1.fq.gz,./path/to/reads/trimmed/taxon_3_R1_val_1.fq.gz \
--right .path/to/reads/trimmed/taxon_1_R2_val_2.fq.gz,./path/to/reads/trimmed/taxon_2_R2_val_2.fq.gz,./path/to/reads/trimmed/taxon_3_R2_val_2.fq.gz \
--output ./path/to/trinity_assembly/ \
--SS_lib_type RF \
--min_contig_length 300 \
--include_supertranscripts 

# Note --CPU and --max_memory in Trinity commands don't do anything when using HPC GridRunner