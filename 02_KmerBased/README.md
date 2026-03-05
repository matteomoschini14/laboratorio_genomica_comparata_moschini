# *Kmer* - Based Genome Survery

## Reads quality check

We used [*FastQC*](https://github.com/s-andrews/FastQC) to assess read quality. It is important to provide files in FASTA format.
```bash
c#[assembly]
fastqc SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz
```

The *FastQC* results are provided as two .html files located in this [folder](./FastQC_results).

-----

## Trimming

The trimming phase is performed to be able to discard sequence portions that do not meet the specified quality criteria. The trimming process was performed on raw FASTQ files before analysis.

```bash
#[assembly]
trimmomatic PE -threads 8 -phred33 SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz SRR11672503_1_paired.fastq SRR11672503_1_unpaired.fastq SRR11672503_2_paired.fastq SRR11672503_2_unpaired.fastq ILLUMINACLIP:/opt/miniforge3/envs/assembly/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> stats_trimmomatic.log
```

[*Trimmomatic*](https://github.com/usadellab/Trimmomatic) yields two output files for each sample: one containing paired reads and another for unpaired reads.

-----

## Compute *k-mer* frequency

Focusing exclusively on the paired reads obtained after trimming, we proceeded with [*KAT*](https://github.com/EarlhamInst/KAT) to further investigate the raw assembly characteristics through k-mer analysis. Specifically, *KAT* utilizes k-mer distribution to provide key insights into the genome assembly.
We used a k-mer size of 27.

```bash
#[kat]
kat hist -t 8 -m 27 -o Anoste SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq
```

The acronym `Anoste` was used to refer to *Anopheles stephensi*, the species whose genome was assembled and annotated in the study.
The primary results of interest are the `Anoste.hist` and `Anoste.png` files, both of which are located in this [folder](./KAT_results).
`Anoste.hist` was uploaded to [GenomeScope](http://genomescope.org/genomescope2.0/) for evaluation: the obtained graph represents the coverage of kmers on the genome.es, both of which are located in this folder. Anoste.hist was uploaded to GenomeScope for evaluation: the obtained graph represents the coverage of kmers on the genome.
