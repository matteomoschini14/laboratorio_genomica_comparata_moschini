# Decontamination

## Preliminary steps

The first step requires re-mapping the short reads to the polished genome.

```bash
#[assembly]
minimap2 -ax sr --MD -t 8 Anoste_pol.fasta SRR11672503_1_paired_fastq SRR11672503_2_paired_fastq | samtools view -Sb - > Anoste-pol_sr.bam
samtools sort -@8 -o Anoste_pol_sr_sorted.bam Anoste_pol_sr.bam
samtools index Anoste_pol_sr_sorted.bam
rm Anoste_pol_sr.bam
```
In this case, the analysis involving long reads is no longer required.

Second step is the taxonomic annotation of the contigs using `blastn`

```bash
blastn -query <ASSEMBLY> -db <PATH/TO/nt/> -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -max_target_seqs 25 -max_hsps 1 -num_threads 25 -evalue 1e-25 -out <OUTFILE>
```
> Due to time constraints, this step was not performed; instead, a link to a pre-existing assembly has been provided.

----

## Contaminants detection

To identify potential contaminants, we utilized [*BlobTools*](https://github.com/DRL/blobtools). This software enables the analysis of genome assembly composition by integrating taxonomic classification with GC content and coverage metrics. Such visualization is instrumental in detecting exogenous sequences or assembly artifacts.

```bash
#[sequence]
blobtools create -i Anoste_pol.fasta -b Anoste_pol_sorted.bam -t Anoste_blast.tsv -o Anoste_blob
blobtools view -i Anoste_blob.blobDB.json -o Anoste
blobtools plot -i Anoste_blob.blobDB.json -o Anoste
```
> The `.json file` is generated as the output of the initial `blobtools create` command.

At this stage, the primary output of interest is the `_table.txt` file; from this file, we aim to retain only the records that specify the target taxonomy—Arthropoda—in the third column.

```bash
grep 'Arthropoda' Anoste.Anoste_blob.blobDB_table.txt | column -t -s $'\t' | less`
```

Subsequently, the following script was used to isolate all previously identified sequences of interest within a single file; notably, this script is compatible with multi-line FASTA formats.

```bash
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' Anoste_pol.fasta | grep -w -Ff <(cut -f1 contig_arthropoda.tsv) - | tr "\t" "\n" > Anoste_decontaminated.fasta
#The first part ensures the FASTA is in oneline form and reformat the file obtaing ">header\tsequence". In this way it is easier to process the file with grep. Then filters out everything is contained inside the patterns file, outputting a fasta file oneline with only desired contigs.
```

The assembly obtained at this stage has been successfully decontaminated.
