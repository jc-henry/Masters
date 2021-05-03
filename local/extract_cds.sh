#!/bin/bash
# Would need to alter file paths
# make sure seqtk installed
# need to create the concatenated .cds (all_taxa_nuc) file
for filename in /path/to/OrthoFinder/Results_date/Single_Copy_Orthologue_Sequences/*.fa 
do
    base=$(basename ${filename} .fa)
    grep ">" ${filename} | sed 's/^>\(.*\)/\1/' > ${base}.txt # sed removes the ">" from the ID. File created for each cluster with IDs
    seqtk subseq /path/to/transdecoder_cds/all_taxa_nuc ${base}.txt >> /path/to/nucleotide_orthologues/${base}.fa # search the concatenated file for the IDs and export to the final .fa cluster files
    rm ${base}.txt # Remove the temporary ID file
done

