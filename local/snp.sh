#/bin/bash

# extract read pairs mapped on same contig and in proper orientation (samtools -f 3 flag)
samtools view -h -b -f 3 bwa_highexpression_taxon_x_sorted.bam > bwa_highexpression_taxon_x_sorted_concordant.bam

# index the BAM file
samtools index bwa_highexpression_taxon_x_sorted_concordant.bam

# get pileup summary only for contigs part of SCO clusters and output bcf file 
bcftools mpileup --regions-file ../contigs_of_interest_nogaps70/clean_taxon_x.txt \
-Ou -f ../trinity_best_supported_contig_assemblies/taxon_x_trinity_rep_contigs.fasta \
./bwa_highexpression_taxon_x_sorted_concordant.bam > taxon_x.bcf

# perform variant calling on the bcf file
# -m for multiallelic calling (output all sites - not just variant sites -v)
bcftools call -m -Ov taxon_x.bcf  > taxon_x_unfiltered.vcf

# perform filtering with vcfutils.pl
# minimum of depth of 10, maximum depth of 100 reads for variant calling
vcfutils.pl varFilter -Q 30 -d 10 -D 100 taxon_x_unfiltered.vcf > taxon_x_vcfutils_filtered.vcf

# further filtering with vcftools
# allow only 2 alleles called
# minimum minor frequency 25 %
vcftools --vcf taxon_x_vcfutils_filtered.vcf --maf 0.25 --max-alleles 2 --min-alleles 2 --minQ 30 --out taxon_x_vcftools_filtered \
--recode --stdout | gzip -c > taxon_x_vcfutils_vcftools_filtered.vcf.gz

# unzip the final output
gunzip taxon_x_vcfutils_vcftools_filtered.vcf.gz

# report the statistics
rtg vcfstats taxon_x_vcfutils_vcftools_filtered.vcf > stats_taxon_x_vcfutils_vcftools_filtered.vcf
