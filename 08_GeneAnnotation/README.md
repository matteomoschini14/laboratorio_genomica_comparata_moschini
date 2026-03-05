# Gene Annotation and Functional Enrichment

## Annotation

The annotation process involves identifying genomic features and assigning biological functions to these sequences. The most widely used tools for identifying sequence similarities are BLAST and DIAMOND:

[Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi): The industry standard for comparing primary biological sequence information. It calculates the statistical significance of matches to identify local regions of similarity between sequences.

[Diamond](https://github.com/bbuchfink/diamond): A high-throughput alternative to BLASTX and BLASTP, designed for the analysis of vast datasets. It offers significant speed improvements (up to 20,000x) over BLAST while maintaining a high degree of sensitivity

While HMMER is the annotation tool that utilizes probabilistic models:

[HMMER](http://hmmer.org/): A software suite used for searching sequence databases for homologs and making sequence alignments. It implements methods using profile hidden Markov models (profile HMMs). Unlike BLAST, which compares sequences directly, HMMER constructs a statistical model of a multiple sequence alignment, making it significantly more sensitive in detecting remote homologs and defining protein domains.

DIAMOND was selected for this analysis due to its ability to compare sequences against massive databases while maintaining a high level of sensitivity.

----

## Input File

The first thing to functionally explore data is creating a file that can be annotate. The best one is a list of protein choosen by selecting the longest sequence in each trimmed orthogroup without taking into consideration gaps. This secure the use of the longest protein that is for sure clustered in an orthogroup and implemented in the analyses. The following script does that:

```bash
#!/bin/bash

# this script looks in every OGs looking for the longest protein in order to create a fasta file ready to be used to infer GO and KEGG terms

> longest_protein_OGs.txt

for orthogroup in *_trimmed.fa; do
	orthoname=$(basename -s _trimmed.fa "$orthogroup")
	maxlen=0
	maxname=""
	#create an orthogroup file without any gap and while read it load two different variables that correspond to consecutive lines
	while IFS= read -r name && IFS= read -r sequence; do
		#get the length of the sequence
		lenprotein=${#sequence}
		#if the measured length is greater to the one already stored, substitute the max length and sequence name variable
		if (( lenprotein > maxlen )); then
			maxname=${name#>}
			maxlen=${lenprotein}
		fi
	done < <(sed 's/-//g' "$orthogroup")
	#retrive che complete sequence from where original orthgroups are stored
	sequence=$(grep -A1 "$maxname" <original_orthogroups_folder>${orthoname}.fa | tail -n1)
	printf ">%s@%s\n%s\n" "$orthoname" "$maxname" "$sequence" >> longest_protein_OGs.txt
done
```
----

## Databases

* Nr: Non-redundant protein sequences from GenPept, Swissprot, PIR, PDF, PDB, and NCBI RefSeq
* Nt: nucleotide sequences
* Swiss-Prot: manually annotated and reviewd proteins ([UniProt](https://www.uniprot.org/))
* [Pfam](http://pfam.xfam.org/): is a large collection of protein families, each represented by multiple sequence alignments and hidden Markov models
* See also [InterPro consortium](http://www.ebi.ac.uk/interpro/)

----

## Diamond 

[Diamond](https://github.com/bbuchfink/diamond) is optimized for large input files of >1 million proteins. DIAMOND was utilized to perform the all-vs-all sequence alignment required for orthogroup inference. Its speed allows for the efficient processing of large-scale proteomic datasets to identify homologous sequences prior to clustering.

The program may use quite a lot of memory and also temporary disk space. Should the program fail due to running out of either one, you need to set a lower value for the block size parameter -b.


```bash
diamond makedb --in /home/STUDENTI/giovanni.galassi3/GG_Laboratorio/05_OG.Inference_Phylogenomic/04_trimmed/longest_protein_OGs.txt --db ./nr_diamond
```

> Due to server issues we coudn't be able to run DIAMOND.

----

## Gene Onthology Annotation

To link our protein sequences to biological functions, we use Gene Ontology terms (GO terms).

GOterms annotation can be performes with a graet variety of programs, such as [PANNZER](http://ekhidna2.biocenter.helsinki.fi/sanspanz/) and [eggNOG mapper](http://eggnog-mapper.embl.de/). They are often reduntant, so the annotation of the proteome is performed using the command line program [InterProScan](https://interproscan-docs.readthedocs.io/en/latest/index.html).

InterProScan is a software package that allows sequences to be scanned against the InterPro database's member consortium. It integrates predictive models (HMMs, regular expressions, and profiles) from multiple source databases into a single analysis pipeline.

```bash
/home/PERSONALE/dbs/interproscan-5.65-97.0/interproscan.sh -i longest_protein_OGs.txt -goterms -pa -b longest_giovanni.tsv -cpu <N_CPUS>
```

> Unfortunately, due to server issue, we couldn't be able to run InterProScan.

----

## GO Enrichment Analysis

To see if our genes of interest show an enrichment in some functional annotation. GO Enrichment Analysis is a statistical method used to interpret the biological significance of the gene lists identified by CAFE. This analysis determines whether specific Gene Ontology terms—representing biological processes, molecular functions, or cellular components—are over-represented in the subset of interest compared to the genomic background.

The core objective is to identify functional "signals" that rise above the background "noise" of the entire genome. We compare a background dataset, composed by every orthogroup that has been assigned at least one GO term, with a foreground one, the subset of "interesting" genes, identified by CAFE as significantly expanded, which are suspected to be under positive selection or involved in a specific phenotypic trait. 

[topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html) was used for this analysis and requires:

1. **Gene universe** (the complete list of genes)
2. **Genes of interest**
3. **GO annotation**


### Background preparation

To proceed with the functional enrichment analysis, a background file ([go_back.tsv](./go_back.tsv)) must first be generated.

The background serves as the reference universe, representing the complete list of OrthoFinder-identified Orthogroups that have been assigned at least one functional annotation (GO term).

To ensure compatibility with the R script, the raw input file (longest_giovanni.tsv) requires preprocessing. This involves removing extraneous metadata, pipe separators (|), and textual descriptions.

```bash
awk -F'\t' '{
  gsub(/@.*/,"",$1); gsub(/\([^)]*\)/,"",$2); gsub(/\|/,",",$2);
  split($2, a, ",");
  for(i in a) if(a[i]!="") seen[$1,a[i]]=1
}
END {
  for(k in seen){
    split(k, b, SUBSEP)
    groups[b[1]] = (groups[b[1]] ? groups[b[1]] "," b[2] : b[2])
  }
  for(g in groups) print g "\t" groups[g]
}' <(cut -f1,14 longest_giovanni.tsv) | grep -v "-" > go_back.tsv
```

A filtering step was applied to the Gamma_asr.tre (derived from the best-fitting model identified via AIC/BIC analysis: `00_2L/2K`) output to retain only significant gene families. Trees were extracted based on specific asterisk (*) combinations indicating significant expansions or contractions.

This operation resulted in three separate files, created by filtering the main tree file for different evolutionary scenarios using grep.

```bash
grep '<1>\*' Gamma_asr.tre | grep '<2>\*' | grep '<3>_' | grep '<4>_' | grep '<5>_' | grep '<6>_' > OG_malaria.txt

grep '<1>_' Gamma_asr.tre | grep '<2>_' | grep '<3>\*' | grep '<4>\*' | grep '<5>_' | grep '<6>_' > OG_dengue.txt

grep '<1>_' Gamma_asr.tre | grep '<2>_' | grep '<3>_' | grep '<4>_' | grep '<5>\*' | grep '<6>\*' > OG_westnile.txt
```

The resulting files are:

+ [OG_malaria.txt](../07_GeneFamilies_Evolution/CAFE_interesting_OG/OG_malaria.txt): Significance in Malaria Vectors;
+ [OG_dengue.txt](../07_GeneFamilies_Evolution/CAFE_interesting_OG/OG_dengue.txt): Significance in Dengue Vectors;
+ [OG_westnile.txt](../07_GeneFamilies_Evolution/CAFE_interesting_OG/OG_westnile.txt): Significance in West Nile Vectors.

### Enrichment Analysis

#### topGO

The functional enrichment analysis is performed in RStudio, utilizing the input files generated in the previous steps. The execution script is provided below:

Prior to executing the R script, the functional annotation data underwent further formatting to ensure compatibility with the required input structure. Specifically, the go_back file was reformatted by collapsing all GOterms assigned to the same orthogroup into a single entry per line. 

```bash
awk -F'\t' '{split($1,a,"_"); n=split($2,g,","); for(i=1;i<=n;i++) if(!seen[a[1],g[i]]++) out[a[1]]=(out[a[1]]?out[a[1]]","g[i]:g[i])} END{for(o in out) print o"\t"out[o]}' go_back.ts
```

Then, the script was launched on RStudio

```bash
#!/usr/bin/env Rscript

library(tidyverse)
library(topGO)

gene_universe <- readMappings(file =
                                "go_back_collapsed")
geneUniverse <- names(gene_universe)

genesOfInterest <- read.table("interesting.txt",header=FALSE)
list_interest1 <- list( "name_interest" = genesOfInterest)

#upload of gene of interest
GOenrichment <- function(trait, trait_name) {
  if (!dir.exists("01_enrichment")) {
    dir.create("01_enrichment")
  }
  
  genesOfInterest <- as.character(trait[[1]]) #as vector not character 
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  
  print(trait_name)
  
  ontology_values = c("BP", "MF", "CC")
  
  GOdata_list <- lapply(ontology_values, function(ontology_value) {
    GOdata_name <- paste("GOdata_", ontology_value, sep = "")
    # annot = annFUN.gene2GO this imparts the program which annotation it should use. In this case, it is specified that it will be in gene2GO format and provided by the user.
    # gene2GO = gene_universe is the argument used to tell where is the annotation
    assign(GOdata_name, new("topGOdata", ontology=ontology_value, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = gene_universe))
  })
  
  elim_list <- lapply(seq_along(ontology_values), function(i) {
    elim_name <- paste("elim_", ontology_values[i], sep="")
    assign(elim_name, runTest(GOdata_list[[i]], algorithm="elim", statistic="fisher"))
  })
  
  results_elim <- function(GO_data, elim_data) {
    num_nodes <- min(1000, length(elim_data@score))
    resulte <- GenTable(GO_data, Classic_Fisher = elim_data,
                        orderBy = "Classic_Fisher", topNodes=num_nodes, numChar=1000)
    resulte$Classic_Fisher <- as.numeric(resulte$Classic_Fisher)
    resulte <- subset(resulte, Classic_Fisher < 0.05)
    return(resulte)
  }
  
  results_elim_list <- lapply(seq_along(ontology_values), function(i) {
    resulte_name <- paste("resulte_", ontology_values[i], sep="")
    assign(resulte_name, envir = .GlobalEnv, results_elim(GOdata_list[[i]], elim_list[[i]]))
  })
  
  write_elim_results <- function(result, ontology_value, trait_name) {
    table_name <- paste("01_enrichment/topGOe_", trait_name, "_", ontology_value, ".txt", sep="")
    write.table(result, file=table_name, quote=F, sep = "\t", row.names = F)
  }
  
  lapply(seq_along(ontology_values), function(i) {
    write_elim_results(results_elim_list[[i]], ontology_values[i], trait_name)
  })
}

#Complete function to perform enrichment
##this particular syntax has been necessary since it was impossible to give the function the trait name it was computing.
GO_enrichment <- function(list) {
  lapply(seq_along(list), function(i) {
  GOenrichment(list[[i]], names(list)[i])
  })
}

#Run the complete function
GO_enrichment(list_interest1)

#search for specific gene
genesOfInterest <- as.character(notext_biggest$V1) #as vector not character 
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
topGO <- new("topGOdata", ontology="BP", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = gene_universe)
print(intersect(genesInTerm(topGO, "GO:0051289")$'GO:0051289', notext_biggest$V1)) #to investigate which genes are annotated with a particular GO term
```

The enrichment results are organized into three topGO folders, categorized by specific disease.

The various files were formatted to ensure compatibility with REVIGO for the subsequent analysis step.

```bash
for i in *; do cut -f1,6 $i >> final_revigo.tsv; done
```
#### REVIGO

Finally, we utilized REVIGO to visualize and summarize the functional enrichment results obtained from topGO. This tool allowed us to reduce the redundancy of the enriched Gene Ontology terms by clustering them based on semantic similarity, thereby facilitating the identification of the most representative biological themes.

The REVIGO results are available in the [Revigo_results](./Revigo_results) folder.
