# CPAP
### Polyadenylation site prediction in *Chlamydomonas reinhardtii*

[![Language Perl](https://img.shields.io/badge/Language-Perl-3876BB)](Scripts/CPAP.pl)
[![Project Bachelor thesis 2015](https://img.shields.io/badge/Project-Bachelor_thesis_2015-8062AB)](#cpap)
[![License MIT](https://img.shields.io/badge/License-MIT-137D87)](LICENSE)

Chlamydomonas Polyadenylation Prediction Tool (CPAP) : A machine learning tool to accurately predict polyadenylation (polyA) sites across all genomic regions (3' UTR, 5' UTR, CDS, introns, and intergenic regions) in *Chlamydomonas reinhardtii*.

I developed CPAP as my Bachelor of Science thesis project in Fall 2015. This repository contains the historical Perl prediction script, supporting model files, and a summary of the thesis methods and results.

## Overview of Polyadenylation

Polyadenylation is a vital post-transcriptional modification in eukaryotes that occurs in the nucleus to protect messenger RNA (mRNA) from degradation during cytoplasmic transport.

### Key Mechanism
The process follows a distinct **two-step pathway**:
1. **Cleavage:** The mRNA strand is cleaved at a specific site.
2. **Tail Addition:** A poly(A) tail is appended directly to the newly formed 3'-end.

## Dataset

Positive examples were assembled from published poly(A) evidence. The documented workflow maps source sequences to the *C. reinhardtii* reference genome, clusters nearby sites, and extracts sequence windows. Negative controls were extracted from genomic background regions and screened against positive-site evidence.

**Source studies listed in the original repository**

- https://doi.org/10.1534/g3.114.010249
- https://doi.org/10.1534/genetics.108.088971

[Positive Dataset](assets/Positive%20Dataset.png) | [Negative Dataset](assets/Negative%20Dataset.png)

## Model Development 

Following were the 4 components of model development:
1. Feature Generation
2. Feature Selection and Model Development 
3. Model Performance Evaluation 
4. Performance Evaluation and Comparison with PASS 2.0 

### Feature Generation

Total of 65 features were generated based on the following data points:

1. Finding novel motifs in the NUE, FUE, DSE, and CE regions around poly (A) sites in addition to the previously reported motifs
   | Region | Full name | Relative base offsets | Length |
   | :---: | :---: | :---: | :---: |
   | FUE | Far upstream element | −140 to −41 | 100 nt |
   | NUE | Near upstream element | −40 to −1 | 40 nt |
   | CE | Cleavage element | −20 to +19 | 40 nt |
   | DSE | Downstream element | 0 to +49 | 50 nt |
3. Using multiple physiochemical and conformational properties of the poly (A) site regions
4. Cleavage site composition
5. Combinations of various K-mer patterns in NUE, FUE, CE, and DSE regions
   
[Feature Generation Pipeline](assets/Feature%20Generation.png)

### Feature Selection and Model Development

Feature selection resulted in 56 features to be used in model development and scoring. 

[Feature Selection Pipeline](assets/Feature%20Selection.png)

Final model was developed with the remaining 56 features.

### Model Performance Evaluation 

Model was evaluated using sensitivity and specificity. Following formulas were used for calculation:

$$ \text{Sensitivity} = \frac{\text{True Positives (TP)}}{\text{True Positives (TP)} + \text{False Negatives (FN)}} $$

$$ \text{Specificity} = \frac{\text{True Negatives (TN)}}{\text{True Negatives (TN)} + \text{False Positives (FP)}} $$

The error margin was ±20 i.e., if the predicted poly (A) site is within ±20 nts of a specified location then it was considered as positive prediction. CPAP had a sensitivity and specificity of **79%** and **41%** on testing dataset. 

### Performance Evaluation and Comparison with PASS 2.0 

CPAP was compared with PASS2.0 by scoring on the same test dataset. 

CPAP had better sensitivity then PASS 2.0 at both overall level and individual region level.

Lower specificity showed higher number of false positives. To counter that, prediction error margin was reduced and both models were re-evaluated.

[Sensitivity Comparison](assets/Sensitivity%20Comparison.png) | [Reduced Error Margin Evaluation](assets/Performance%20Evaluation.png)

It was observed that individual performance of both models dropped, however, if by pooling the predictions of both models together 70% sensitivity was still recovered. 

## Conclusion

* **New Prediction Tool**: Developed CPAP, a new program for poly (A) sites prediction in *Chlamydomonas reinhardtii* using logistic regression for feature selection.
* **Performance Boost**: CPAP outperformed PASS 2.0 across key metrics:
  - Sensitivity (Sn): 6% higher on the combined positive dataset.
  - Specificity (Sp): 5% better on coding sequences (CDS) and 21% better on intronic control sequences.
* **Complementary Approach**: Combining correct predictions from both CPAP and PASS 2.0 yielded the best results. The distinct features used by each program led to a large number of valuable, unique predictions.
* **Proposed Solution**: Because both tools individually showed lower specificity, a complementary pipeline was proposed to reduce the error margin and secure high-quality poly (A) sites with better Sp.

