#!/bin/bash -e

#SBATCH --account       project-code                 								# project-code
#SBATCH --job-name      bowtie2          		                                    # job-name
#SBATCH --time          15:00:00                    								# hh:mm:ss time allocated
#SBATCH --mem           3G                           							    # memory allocated 
#SBATCH --cpus-per-task 36                           								# threads allocated 
#SBATCH --output        ./taxon_bowtie2.%j.out
#SBATCH --error         ./taxon_bowtie2.%j.err
#SBATCH --mail-type     END
#SBATCH --mail-user     j.henry2@massey.ac.nz
#SBATCH --chdir         /nesi/nobackup/project-code/taxon/path/to/bowtie2/

module purge
module load Bowtie2/2.3.5-GCC-7.4.0
module load SAMtools/1.10-GCC-9.2.0

time bowtie2-build \
--threads ${SLURM_CPUS_PER_TASK} \
/path/to/trinity_assembly/Trinity_taxon.fasta \
bowtie2_taxon

time bowtie2 -q --sensitive --no-unal --threads ${SLURM_CPUS_PER_TASK} --fr --nofw \
-x ./bowtie2_taxon \
-1 /path/to/reads/trimmed/taxon_1_R1_val_1.fq.gz,/path/to/reads/trimmed/taxon_2_R1_val_1.fq.gz,/path/to/reads/trimmed/taxon_3_R1_val_1.fq.gz \
-2 /path/to/reads/trimmed/taxon_1_R2_val_2.fq.gz,/path/to/reads/trimmed/taxon_2_R2_val_2.fq.gz,/path/to/reads/trimmed/taxon_3_R2_val_2.fq.gz \
| samtools view -b - | samtools sort -o - > bowtie2_taxon_sorted.bam

# note the strand-specific options "--fr --nofw"
# end to end alignment is default