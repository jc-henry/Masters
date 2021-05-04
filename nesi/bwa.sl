#!/bin/bash -e

#SBATCH --account       project-code                 								# project-code
#SBATCH --job-name      taxon_x_trinity_bwa           				                # job-name
#SBATCH --time          03:00:00                    								# hh:mm:ss time allocated
#SBATCH --mem           7G                           							    # memory allocated
#SBATCH --cpus-per-task 10                           								# threads allocated
#SBATCH --output        /nesi/nobackup/project-code/taxon/path/to/taxon_x_trinity/bwa/taxon_x_bwa.%j.output
#SBATCH --error         /nesi/nobackup/project-code/taxon/path/to/taxon_x_trinity/bwa/taxon_x_bwa.%j.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/taxon/path/to/taxon_x_trinity/bwa/

module purge
module load BWA/0.7.17-gimkl-2017a
module load SAMtools/1.10-GCC-9.2.0

READ_DIR=/nesi/nobackup/project-code/taxon/path/to/reads/trimmed
ASS_DIR=/nesi/nobackup/project-code/taxon/path/to/taxon_x_trinity/trinity_assembly

time bwa index \
-a is \
-p bwa_highexpression_taxon_x_trinity \
$ASS_DIR/taxon_x_trinity_rep_contigs.fasta

time bwa mem -t ${SLURM_CPUS_PER_TASK} -v 3 \
bwa_highexpression_taxon_x_trinity \
$READ_DIR/taxon_x_R1_val_1.fq.gz \
$READ_DIR/taxon_x_R2_val_2.fq.gz \
| samtools view -b - | samtools sort -o - > bwa_highexpression_taxon_x_sorted.bam


