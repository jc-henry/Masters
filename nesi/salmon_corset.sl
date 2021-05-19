#!/bin/bash -e

#SBATCH --account       project-code                  		# project-code
#SBATCH --job-name      hutt_monroi_1_salmon_corset         # job-name
#SBATCH --time          00:45:00                     		# hh:mm:ss time allocated 
#SBATCH --mem           2500MB                         	    # memory allocated
#SBATCH --cpus-per-task 6                            		# threads allocated
#SBATCH --output        ./hutt_monroi_1_salmon_corset.%j.output
#SBATCH --error         ./hutt_monroi_1_salmon_corset.%j.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/path/to/hutt_monroi_1_rnaspades/salmon_corset/

module purge
module load Salmon/1.3.0-gimkl-2020a

time salmon index --index ./salmon_index_hutt_monroi_1 \
--transcripts /nesi/nobackup/project-code/path/to/hutt_monroi_1_rnaspades/rnaspades_assembly/transcripts.fasta -k 31


READ_DIR=/nesi/nobackup/project-code/path/to/reads/trimmed/

time salmon quant --index ./salmon_index_hutt_monroi_1 --libType ISR \
-1 $READ_DIR/hutt_monroi_1_R1_val_1.fq.gz \
-2 $READ_DIR/hutt_monroi_1_R2_val_2.fq.gz \
-p ${SLURM_CPUS_PER_TASK} --gcBias --hardFilter --dumpEq hutt_monroi_1_ \
--output ./hutt_monroi_1.out

# Need --hardFilter flag if using for Corset clustering

# Corset can't work with zipped files
gunzip ./hutt_monroi_1.out/aux_info/eq_classes.txt.gz