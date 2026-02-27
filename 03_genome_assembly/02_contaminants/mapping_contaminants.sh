#!/bin/bash

# Mapping short reads on polished assembly for contaminants detection. Mutate SAM into BAM

# Short reads
minimap2 -ax sr --MD -t 8 Anoste_pol.fasta SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq | samtools view -Sb - > Anoste_pol_sr.bam
samtools sort -@8 -o Anoste_pol_sr_sorted.bam Anoste_pol_sr.bam
samtools index Anoste_pol_sr_sorted.bam
rm Anoste_pol_sr.bam


