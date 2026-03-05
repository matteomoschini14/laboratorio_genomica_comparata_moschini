# Genome Annotation 

The annotation phase of a newly assembled genome allows for the identification of all its functional characteristics. This process relies on the .gff format, which is essential for mapping the coordinates of various features across the genome.
Functional annotation involves assigning specific biological functions to nucleotide sequences, effectively mapping functional roles to precise genomic coordinates.

## Annotation of repetitive elements

Repetitive elements can significantly interfere with genome analysis, especially when focusing on the coding regions. The objective of this phase is to mask these elements, ensuring they are excluded from consideration by downstream annotation software.
To achieve a high-quality annotation, it is recommended to use a transcriptome (or proteome) as a reference. This provides additional biological evidence, as it contains crucial information regarding isoform diversity.

For a quick and dirty transposon annotation, we employ [*RepeatModeler*](https://github.com/Dfam-consortium/RepeatModeler) and [*RepeatMasker*](https://github.com/Dfam-consortium/RepeatMasker). Briefly, RepeatModeler identifies a set of consensus sequences representative of the repeats present in the genome. This library is then used by RepeatMasker to annotate the assembly.

> Due to computational constraints, RepeatModeler has been pre-run, and the consensus library can be found at the following path: /home/PERSONALE/mirko.martini3/01_2024/00_Data/02_annotation/Anoste_RepeatModeler_library.fa

## First Annotation round (MAKER)

The software used for the annotation is *MAKER*. To function correctly, it requires several control and configuration files, which are generated using the command below.

```bash
maker -CTL
```

This command generates the files required for MAKER to operate. The program relies on configuration files rather than standard command-line arguments; it is designed to parse inputs and settings from these pre-filled files. All required parameters have been correctly configured.

```bash
    #-----Genome (these are always required)
    genome= Anoste_pol.decontaminated #genome sequence (fasta file or fasta embeded in GFF3 file)
    organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic
    
    #-----Re-annotation Using MAKER Derived GFF3
    maker_gff= #MAKER derived GFF3 file
    est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
    altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
    protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
    rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
    model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
    pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
    other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no
    
    #-----EST Evidence (for best results provide a file for at least one)
    est= #set of ESTs or assembled mRNA-seq in fasta format
    altest= #EST/cDNA sequence file in fasta format from an alternate organism
    est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
    altest_gff= #aligned ESTs from a closly relate species in GFF3 format
    
    #-----Protein Homology Evidence (for best results provide a file for at least one)
    protein= <PROTEOME> #protein sequence file in fasta format (i.e. from mutiple oransisms)
    protein_gff=  #aligned protein homology evidence from an external GFF3 file
    
    #-----Repeat Masking (leave values blank to skip repeat masking)
    model_org=<EMPTY> #select a model organism for RepBase masking in RepeatMasker
    rmlib= <RepeatModeler_library> #provide an organism specific repeat library in fasta format for RepeatMasker
    repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
    rm_gff= #pre-identified repeat elements from an external GFF3 file
    prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
    softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)
    
    #-----Gene Prediction
    snaphmm= #SNAP HMM file
    gmhmm= #GeneMark HMM file
    augustus_species= #Augustus gene prediction species model
    fgenesh_par_file= #FGENESH parameter file
    pred_gff= #ab-initio predictions from an external GFF3 file
    model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
    est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
    protein2genome=<0 OR 1> #infer predictions from protein homology, 1 = yes, 0 = no
    trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
    snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
    unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
    
    #-----Other Annotation Feature Types (features MAKER doesn't recognize)
    other_gff= #extra features to pass-through to final MAKER generated GFF3 file
    
    #-----External Application Behavior Options
    alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
    cpus=<CPUS> #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)
    
    #-----MAKER Behavior Options
    max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
    min_contig=1 #skip genome contigs below this length (under 10kb are often useless)
    
    pred_flank=200 #flank for extending evidence clusters sent to gene predictors
    pred_stats=<0 or 1> #report AED and QI statistics for all predictions as well as models
    AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
    min_protein=50 #require at least this many amino acids in predicted proteins
    alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
    always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
    map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
    keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)
    
    split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
    single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
    single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
    correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes
    
    tries=2 #number of times to try a contig if there is a failure for some reason
    clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
    clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
    TMP= #specify a directory other than the system default temporary directory for temporary files
```

MAKER was then launched

```bash
maker -base <OUTPUT PREFIX>
```

Once the process was complete, it was necessary to merge the resulting .fasta and .gff files.

```bash
fasta_merge -d <DATASTORE INDEX FILE>
gff3_merge -d <DATASTORE INDEX FILE>
```

The summary statistics can then be evaluated; the most commonly analyzed metric is *AGAT*, that concerns the repetitive element content.

```bash
agat_sp_statistics.pl --gff file.gff -o <output_file>    # Gene models summary statistics.
agat_sq_repeats_analyzer.pl -i <input_file> -o <output_file>   # Repeats' summary statistics. 
```

This concludes the first round of MAKER annotation.

-----

## Second Annotation round (SNAP and Augustus)

SNAP is a gene prediction tool designed to identify protein-coding genes within DNA sequences, particularly in de novo assembled genomes. It utilizes a  HMM (Hidden Markov Model) probabilistic approach and can be trained using a set of known genes to significantly enhance its prediction accuracy.
SNAP was run using a custom script with the following sequence of commands:

```bash
maker2zff -c 0 -e 0  -l 80 -x 0.1 -d <datastore_index> #To extract gene models based on mutiple filter criterion. It transforms maker output into .zff files, the format for SNAP
fathom <.ANN FILE> <.DNA FILE> -gene-stats #Print some summary statistics of the selected gene models
fathom <.ANN FILE> <.DNA FILE> -validate #Validate gene models and print summary statistics. output is in the stdout and cannot be redirected
fathom <.ANN FILE> <.DNA FILE> -categorize 1000 #Extract gene modeles together with 1000 bp at both ends for training
fathom <UNI.ANN FILE> <UNI.DNA FILE> -export 1000 -plus #Export and convert unique genes to strand plus. Tipically, only the unique genes are sued for training and testing. The value limits the intergenic sequence at the ends.
forge <export.ann> <export.dna> #call the parameter estimation program, better in another directory
hmm-assembler.pl <NAME> <FORGE DIRECTORY> > <OUTPUT HMM FILE>
```
> The proteome-to-genome results must be provided to SNAP in order to identify and define the protein prediction models.

Augustus is a widely used gene prediction tool that identifies protein-coding genes in eukaryotic genomes. It employs an ab initio approach based on a generalized Hidden Markov Model (HMM), which allows it to predict gene structures, even in the absence of external transcriptomic evidence.
Running Augustus as a standalone tool is highly computationally intensive. However, BUSCO allows for a more time-efficient execution of the software.

```bash
busco -i ../../03_GenomeAssembly/03_scaffolding/Anoste_chr.fasta -c 30 -l $BUSCO/culicidae_odb12 --augustus --long -m genome --out Anoste_cu --augustus_parameters='--progress=true'
```

At this stage, the second round of MAKER can be initiated, ensuring that the parameters are correctly updated to incorporate the newly trained models.

> The annotation quality improves significantly when transitioning from the first to the second round.
