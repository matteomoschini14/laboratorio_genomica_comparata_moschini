# Gene Families Evolution

## CAFE 

*CAFE* is one of the most widely used software tools based on gene and species trees. Its primary purpose is to analyze changes in gene family sizes, accounting for phylogenetic history to provide a statistical foundation for evolutionary inferences.

The program employs a birth-and-death process to model gene gain and loss across a user-specified phylogenetic tree. The distribution of family sizes generated under this model serves as a baseline to assess the statistical significance of observed size differences among taxa. In brief, CAFE estimates the ancestral states of gene families by comparing the tip states (extant family sizes) across the phylogeny.

----

### Data Preparation

*CAFE* requires two main inputs: the ultrametric tree (generated in the previous section) and a gene count table. The table must be formatted as follows:
NONE Orthogroup ID ID ID ID
NONE OG0000000  93 0  5  9 

*OrthoFinder* generates a similar table, located in the results folder as `Orthogroups.GeneCount.`. However, this file requires preprocessing to match CAFE's input format:

-The first column required by CAFE is missing.

-There is an extra column that is unnecessary and must be removed.

This file was processed to meet CAFE's input requirements, ensuring the correct column structure and formatting.

```bash
sed $'s/^/NONE\t/g' Orthogroups.GeneCount.tsv | rev Orthogroups.GeneCount.tsv  | cut -f 2- | rev > ../../../../07_GeneFamilies_Evolution/GeneCount_CAFE.tsv 
```
Finally, please verify that your phylogenetic tree is formatted as Newick rather than Nexus before proceeding with the analysis.

----

### Generating the Error Model

While the program can immediately infer one or more birth-death parameters (lambda), it is standard practice to first estimate a general error model. 
This statistical model describes the probability that the observed gene count in a family differs from the true count due to technical artifacts (such as genome assembly or annotation errors). The statistics derived from this model will be incorporated into subsequent analysis steps to ensure more accurate evolutionary inferences.

```bash
#[tree]
cafe5 -i GenCount_CAFE.tsv -t timetree.nwk -o Error_model -e
```

----

### Single-Lambda Analysis
Then we can run the program specifying the error model just inferred. An initial *CAFE* analysis was performed using a single $\lambda$ (lambda) parameter for the entire tree. This model assumes a uniform turnover rate, implying that the tendency for gene families to expand or contract remains constant across all species.

```bash
#[tree]
for k in {1..5}; do for n in {1..10}; do mkdir -p 00_1L/${k}K/${n}N; cafe5 -i GenCount_ -t timetree.nwk -o 00_1L/${k}K/${n}N -e./Error_model/Base_error_model.txt -k ${k}; done; done
```
----

### Two-Lambda Analysis

To allow for rate heterogeneity across the phylogeny, we performed a two-lambda analysis. This model assigns distinct turnover rates to specific lineages. To compute multiple $\lambda$ rates, you must define the number of distinct rate categories and assign them to specific species or clades.These assignments are specified in a separate tree file (distinct from the main ultrametric tree) passed via the -y flag. In this file, integers are used to map specific branches to their corresponding $\lambda$ category.

```bash
(Anofun:1,(Anoste:1,((Culqui:2,Culpip:2):1,(Aedalb:1,Aedaeg:1):1):1):1);
```

The execution script for the two-lambda analysis is provided below:

```bash
#[tree]
for k in {1..5}; do for n in {1..10}; do mkdir -p 00_2L/${k}K/${n}N; cafe5 -i GeneCount_ -t timetree.nwk -o 00_2L/${k}K/${n}N -y timetree_2Lambda.nwk -e./Error_model/Base_error_model.txt -k ${k}; done; done
```

----

### Key Parameters

- lambda: maximum-likelihood estimation of a global or local gene family evolutionary rates. Since the tree used is a time tree, branch lengths are in million of years, the value describe the turnover per gene per million years. With a reasonable evolution hypothesis, different species can be grouped into different lambdas.
- gamma: even if all species could be described by the same λ, this does not mean that all gene families evolve at the same rate. Γ divides gene families into rate categories that share the same mean λ (based on how this parameter has been distributed). Reasonable values for gamma are between two and five.
- alpha (present ony if gamma >= 2): defines the shape of the gamma distribution build around λ. Small α means that some gene families evolve very fast and others very slow, while big α the opposite.
- epsilon (present ony if gamma >= 2): it represents the number of static families that do not show any turnover.

----

## Model Selection

This project employs *AIC* and *BIC* to determine the optimal model.

#### AIC (Akaike Information Criterion)
*AIC* estimates the relative amount of information lost by a given model. It is generally preferred when the primary goal is prediction. A lower AIC value indicates a better model. Rewards models that fit the data well while penalizing unnecessary complexity.

$$AIC = 2k - 2\ln(\hat{L})$$

#### BIC (Bayesian Information Criterion)
*BIC* penalizes model complexity (the number of parameters) more strictly than AIC, especially for large datasets. It is preferred when the goal is to identify the true underlying model or when a more parsimonious solution is desired.

$$BIC = k\ln(n) - 2\ln(\hat{L})$$ 

> $\hat{L}$ stands for Likelihood


### Model estimation

Likelihood values from all 10 technical replicates of the single-lambda CAFE run were extracted and saved to `sum_results.txt`. This procedure was repeated for each directory representing the 5 different Gamma values (ranging from 1K to 5K).

 ```bash
# The following command applies exclusively to the 1-Gamma analysis (Base_results.txt).
for folder in */; do lnL=$(grep "lnL" ${folder}/Base_results.txt | grep -oE "[0-9]*([\.,][0-9]*)?"); L=$(grep "Lambda" ${folder}/Base_results.txt | grep -oE "[0-9]*\.[0-9]*"); E=$(grep "Epsilon" ${folder}/Base_results.txt | grep -oE "[0-9]*\.[0-9]*"); echo -e "$lnL\t$L\t$E" >> sum_results.tsv; done

# The following command covers the remaining analyses for 2 to 5-Gamma categories (Gamma_results.txt)
for i in */; do cd $i; for folder in */; do lnL=$(grep "lnL" ${folder}/Gamma_results.txt | grep -oE "[0-9]*([\.,][0-9]*)?"); L=$(grep "Lambda" ${folder}/Gamma_results.txt | grep -oE "[0-9]*\.[0-9]*"); E=$(grep "Epsilon" ${folder}/Gamma_results.txt | grep -oE "[0-9]*\.[0-9]*"); A=$(grep "Alpha" ${folder}/Gamma_results.txt | grep -oE "[0-9]*\.[0-9]*"); echo -e "$lnL\t$L\t$E\t$A" >> sum_results.tsv; done; cd ..; done
```

Upon completion of this step, the best likelihood values were extracted from each generated `sum_results.txt` file and consolidated into a new file named `all_L.txt`.
```bash
for f in */; do cut -f1 "$f"/sum_results.tsv | sort -n | head -n1; done > all_L.txt
```

The same procedure was also applied to the 2-lambda analysis.

In all_L.txt files, the first likelihood value belongs to single-Gamma value run, the second to the run with two values, and so on. To perform the test, it is necessary to manually define the number of 'free parameters' for each run, which are determined by summing possible values assumed by Lambda and Gamma. 

```bash
00_1L/all_L.txt
100007  2
97485.3 3
98129.7 4
98290.9 5 
98428.3 6

00_2L/all_L.txt
98345.6 3
96613.3 4
97170.6 5
96840.9 6
97544.1 7
```

It is also necessary to define natural logarithm ($\ln$) of trees' number reported in `Gamma_asr.tre` (or `Base_asr.tre`).

```bash
grep "OG" Gamma_asr.tre | wc -l 
> 10526

ln(10526) = 9.26
```

Finally, we perform the calculation of the two indicators, which will be carried out separately for both the `00_1L` and `00_2L` runs, to then compare them and choose the best-fitting model. 

```bash
paste --delimiters=$"\t" all_L.txt <(while IFS=$'\t' read -r L k; do echo "2*$k + 2*$L" | bc; done < all_L.txt) <(while IFS=$'\t' read -r L k; do echo "$k*9.21 + 2*$L" | bc; done < all_L.txt) | sort -k4,4n > AIC_BIC.tsv
```

This generates the [AIC_BIC.tsv](./Model_selection/AIC_BIC.tsv) file, which contains the best values from the 1-lambda and 2-lambda runs. Within the file, entries are sorted so that the first one represents the best-fitting model. In this case, it corresponds to `00_2L/2K`.
 
----
