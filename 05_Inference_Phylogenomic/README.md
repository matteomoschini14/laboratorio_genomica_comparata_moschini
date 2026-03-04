# Ortholgy inference and Phylogenomic

This represents the comparative analysis proper, where the functional portion of the genome and selective pressures are examined to link specific evolutionary patterns to corresponding phenotypic patterns.
The analysis aims to identify genes that exhibit similar selective pressures across species sharing the same phenotype.

## Dataset creation

To perform the comparison, it is essential to have a comprehensive genomic dataset, in this case provided in the [mosquitoes.tsv](../00_practice/00_dataset/mosquitoes.tsv) file.

To retrieve the genomes selected for the comparative analysis, a dedicated script was developed, as shown below:

```bash
#!/bin/bash

#This script downloads genomes and relative gff from NCBI datasets creating ready to use folders

AN2name=$1

mkdir 00_genome
mkdir 01_gff

while IFS=$'\t' read -r AN sname ID; do
	echo $AN
	#download specifying the name of the folder
	datasets download genome accession "$AN" --filename "$ID".zip --include genome,gff3
	#unzip specifying the name of the folder
	unzip "$ID".zip -d "$ID"
	#rename the two file of interest
	mv "$ID"/ncbi_dataset/data/"$AN"/*.fna 00_genome/"$ID".fna
	mv "$ID"/ncbi_dataset/data/"$AN"/*.gff 01_gff/"$ID".gff
	#delete the folder
	rm -r "$ID"/
done < "$AN2name"
```

### Longest isoform and translation

Once the genomes have been downloaded, the raw data consists of nucleotide sequences. However, in comparative genomics, analyses are typically conducted at the amino acid level rather than using raw nucleotides. To perform this conversion, the AGAT toolkit is utilized. Specifically, once the script results are obtained, a for-cicle can be executed to process the files through this command.

```bash
#[sequence]
for gff in *.gff; do agat_sp_keep_longest_isoform.pl --gff "$gff" -o ${gff/.gff/_longest.gff}; done
```

This initial script is designed to select only the longest isoform for each gene. Including all alternative isoforms would artificially inflate the proteome size; by retaining only the longest sequence, we ensure that the maximum number of exons is represented for each locus.

```bash
#[sequence]
for gff in *_longest.gff; do agat_sp_extract_sequences.pl -g "$gff" -f ../00_genome/${gff/_longest.gff/.fna} -t cds -p --cfs --output ../02_raw_proteoms/${gff/_longest.gff/.faa}; done
```

This second script is designed to extract CDS (Coding Sequence) features and translate them into proteins using the standard eukaryotic genetic code.

### Removal of pseudogenes

Despite being reference-grade sequences, these genomes may still contain pseudogenes (defined here as any sequence containing internal premature stop codons). Consequently, these must be removed. The following script was employed for this purpose:

```bash
#!/bin/bash/

# list each plausible pseudogene present. 

mkdir raw_proteomes
mv *.faa raw_proteomes/

#make proteome one line
cd raw_proteomes
for proteome in *.faa; do
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' "$proteome" > ../${proteome/.faa}".faa"
done
cd ..

#Extract all pseudogene names
mkdir 00_pseudogene_name

for proteome in *.faa; do
	species=$(basename -s .faa "$proteome")
	grep -B1 '*' "$proteome" | grep ">" >> 00_pseudogene_name/"$species"_pseudogenes_name.txt
done

#removes sequences identified as pseudogenes

for pseudo_file in 00_pseudogene_name/*_pseudogenes_name.txt; do
	species=$(basename -s _pseudogenes_name.txt "$pseudo_file")
	while IFS=$'\t' read -r header; do
		sed -E -i "/${header}/{N;d;}" "$species".faa # N option loads the next line found after the pattern and put it into pattern space too; d delete the pattern space
	done < "$pseudo_file" 
done

mv 00_pseudogene_name ../00_genome
```
### Header modification

At this stage, the FASTA headers must be modified to ensure they are clear and functional for downstream analyses. Specifically, the headers were truncated to retain only the species identifier and the unique sequence ID. This step is indispensable for the subsequent stages of the analysis.

```bash
for prot in *.faa; do ID=$(basename -s .faa "$prot"); sed -i.old -E "s/>(rna-XM_[0-9]+\.[0-9]) (gene=gene-.[^ ]+) name=(.[^ ]+) .+$/ >${ID}\|\3/" "$prot"; done
```

----

## Orthofinder

[*OrthoFinder*](https://github.com/davidemms/OrthoFinder) was employed to identify orthologous groups. The software identifies orthologs through a pipeline that searches for maximum sequence similarity to infer orthology. It utilizes a tree-based method, defining orthology through speciation events; consequently, it distinguishes between duplication events, which create paralogs, and speciation events, which result in orthologs, based on the reconstructed gene trees.

```bash
#[orthofinder]
orthofinder -t 8 -a 8 -f .
```

The software generates a comprehensive output organized into several directories, each containing specific types of evolutionary and genomic information.

----

## Paralog filtering

To refine the dataset, it was necessary to filter the orthogroups to remove paralogs. For this purpose, we utilized [*DISCO*](https://github.com/JSdoubleL/DISCO), which decomposes complex orthogroups into distinct, orthology-consistent sub-clusters. DISCO script requires a specific header syntax to properly work. Luckily, the supported syntax is the one we already implemented and used so far (SPECIES|SEQUENCE_ID).

```bash
#[tree]
while IFS=' ' read -r OG tree; do python3 ~/GG_Laboratorio/99_scripts/disco.py -i <(echo "$tree") -o ../../../01_disco/${OG/:/}.nwk -d "|" -m 4 --remove_in_paralogs --keep-labels --verbose >> ../../../01_disco/disco.log; done < <(sed -E 's/[A-Z][a-z]{5}_//g; s/\)n[0-9]*+/\)/g' Resolved_Gene_Trees.txt)
```

> The `disco.py` script is available at the following [link](https://github.com/jsdoublel/DISCO/blob/master/disco.py)

To optimize computational resources and storage, empty DISCO output files (zero-size files) were identified and removed from the dataset.

```bash
find . -size 0 -print > empty_disco.txt
find . -size 0 -delete
```

Following the filtering process, the amino acid sequences associated with each gene in the identified orthogroups were reconstructed. We then proceeded to partition the phylogenetic trees (in Newick format, .nwk) of each orthogroup into individual amino acid sequence files (.faa). This step was performed using a dedicated script (split_disco_output.sh), which extracts the phylogenetic trees from the DISCO output and maps them to their corresponding protein sequences retrieved from the OrthoFinder results.


```bash
#!/bin/bash

echo "Splitting DISCO outputs ..."

mkdir split

# Split the multiline output in multiple outputs
for disco in *.nwk; do
	OGname=$(basename -s .nwk "$disco")
	split -d -a 2 -l 1 --additional-suffix=.nwk "$disco" split/"$OGname"_
done

mkdir ../disco_OG

# Define original sequence folder with $1
OG_folder=$(realpath $1)

# recreate orthogroups using original ones. It is MANDATORY to substitute all '-' with '_' because it is te only way to grap exactly
# the meant sequence.

echo "Recreating orthogroups ..."

cd split

for tree in *_*.nwk; do
    name="${tree%.nwk}"
    OG="${name%%_*}"  # Extract the part before the number
    output="../../disco_OG/$name.faa"

    # Extract unique sequence names from tree
    grep -o -E "[A-Z][a-z]{5}\.[^:]+" "$tree" | sort -u > tmp_seqs.txt

    # Extract sequences in one go using awk
    awk -v seqs=tmp_seqs.txt '
        BEGIN {
            while ((getline < seqs) > 0) wanted[$1] = 1
        }
        /^>/ {
            seq = substr($0, 2)
            keep = (seq in wanted)
        }
        keep { print }
    ' "$OG_folder/$OG.faa" >> "$output"
done

# Clean up temporary file
rm -f tmp_seqs.txt
```

----

## Alignments and trimming

Following *DISCO*, the pipeline requires multiple sequence alignment (MSA) and trimming. These steps are crucial because evolutionary analysis depends not just on orthology, but on establishing positional information to identify conserved and divergent sites.

### Alignment

The alignment process requires single-copy orthologs, which are located in the 'Single_Copy_Orthologue_Sequence' directory generated by OrthoFinder. From this dataset, 200 sequences were randomly selected and exported into a .txt file for downstream processing.

```bash
ls *.fa | shuf -n 200 > species_tree_OG.txt
```

At this stage, the [*MAFFT*](https://github.com/GSLBiotech/mafft) aligner can be executed to perform the multiple sequence alignment.

```bash
#[sequence]
for OG in $(cat species_tree_OG.txt); do mafft --auto --anysymbol "$OG" > ../../../03_aligned/${OG/.faa/_aligned.faa}; done
```

To align every sequence obtained after DISCO, we used the following command line:

```bash
#[sequence]
for OG in *.faa; do mafft --auto --anysymbol "$OG" > ../../../03_aligned/${OG/.faa/_aligned.faa}; done
```

### Trimming

The pipeline continues with trimming using *BMGE*, which enhances alignment quality by identifying and removing poorly aligned regions. This software employs an entropy-based approach to filter out non-informative sites. The objective of this phase is to clean the dataset and eliminate sequences that appear erroneous or misaligned following the alignment process.
 
```bash
conda deactivate 
for OG in *; do bmge -i "$OG" -t AA -m BLOSUM62 -e 0.5 -g 0.4 -of ../04_trimmed/${OG/_aligned.faa/_trimmed.faa}; done
```

> It is important to note that BMGE must be executed outside of any Conda environment to avoid dependency conflicts and ensure proper functionality.

----

## Concatenation

At this stage, all the individual sequences—already aligned and trimmed—from each orthogroup are merged into a single large-scale concatenation. This process was carried out using the [AMAS.py](https://github.com/marekborowiec/AMAS/blob/master/amas/AMAS.py) script, resulting in a supermatrix that represents the combined evolutionary signal of the dataset.
Before running `AMAS.py`, it is crucial to modify the FASTA headers so they contain only the species name. This is because the script uses the header as a key to match sequences across different orthogroups, concatenating all sequences with the same identifier into a single taxonomic entry.

```bash
sed -E 's/\|.+$//' * #Preliminary step to standardize headers for the subsequent concatenation process
~/GG_Laboratorio/99_scripts/AMAS.py concat -y nexus -i *.faa -f fasta -d aa -t conc_species_tree
```

----

## Species tree reconstruction

All necessary inputs are now available to initiate the phylogenetic analysis for the species tree reconstruction. In this study, we employed *IQ-TREE*, utilizing the concatenated sequences as the primary dataset for the inference.

```bash
#[tree]
iqtree -m TESTNEW -b 100 -s conc_species_tree --prefix species_tree -nt 9
```
