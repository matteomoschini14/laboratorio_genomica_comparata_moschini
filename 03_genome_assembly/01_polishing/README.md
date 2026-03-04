 # Polishing

### Mapping short and long reads

The first step in the read polishing process is mapping, a procedure in which reads are realigned to a specific assembly. This allows us to determine the positioning of short reads on the assembly and evaluate their alignment accuracy.
The mapping process is performed using a script based on the `minimap2` command, and it is applied to both short and long reads.

#### Mapping short reads

```bash
#Short reads
#[assembly]
minimap2 -ax sr --MD -t 6 Anoste_raw.fasta SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq > Anoste_raw_sr.sam
samtools view -Sb Anoste_raw_sr.sam > Anoste_raw.bam
rm Anoste_raw_sr.sam
samtools sort -@6 -o Anoste_raw_sr_sorted.bam Anoste_raw_sr.bam
samtools index Anoste_raw_sr.bam
rm Anoste_raw_sr.bam
```

#### Mapping long reads

```bash
#Long reads
#[assembly]
minimap2 -ax map-pb --MD -t 6 Anoste_raw.fasta SRR11672506.fastq.gz > Anoste_raw_lr.sam
samtools view -Sb Anoste_raw_lr.sam > Anoste_raw_lr.bam
rm Anoste_raw_lr.sam
samtools sort -@6 -o Anoste_raw_lr_sorted.bam Anoste_raw_lr.bam
samtools index Anoste_raw_lr_sorted.bam
rm Anoste_raw_lr.bam
```

### Hypo

The actual polishing is performed using [*HyPo*](https://github.com/kensung-lab/hypo); however, providing the short-read coverage information is a prerequisite for running the tool.

#### Mean coverage of short reads

[*Mosdepth*](https://github.com/brentp/mosdepth) was used to obtain the missing coverage information required to proceed.

```bash
mosdepth -n --fast-mode --by 500 Anoste_raw_sr Anoste_raw_sr_sorted.bam
```

#### Hypo
```bash
#[assembly]
echo -e "$R1\n$R2" > Sr.Path
hypo -d Anoste_raw.fasta -r @Sr.path -s 227m -c 136 -b Anoste_raw_sr_sorted.bam -B Anoste_raw_lr_sorted.bam -t 6
```
> `-s` is the approximate genome size, data taken from GenomeScope.

> `@Sr.path` is a file that contains the paths of the reads (pasted with `echo` command).

Upon completion, the process generates the file hypo_anoste_raw.fasta, which was subsequently renamed to anoste_pol.fasta.

### Quality assessment of the polished assembly

Finally, a quality control assessment of the polished assembly was conducted using the same methods described in section [00_Assembly_raw](../00_Assembly_raw/README.md): N50, BUSCO, and Spectra-cn analysis.

```bash
#[assembly]
assembly-stats Anoste_pol.fasta > Anoste_pol.stats

#[sequence]
export NUMEXPR_MAX_THREADS=80
busco -m geno -l $Busco/culicidae_odb12 -c 8 -o Anoste_pol_busco -i Anoste_pol.fasta

#[kat]
kat comp -t 8 -o Anoste_pol 'SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq' Anoste_pol.fasta
```
