# CPAP
Chlamydomonas Polyadenylation Prediction Tool (CPAP) : A machine learning tool to accurately predict polyadenylation (polyA) sites across all genomic regions (3' UTR, 5' UTR, CDS, introns, and intergenic regions) in *Chlamydomonas reinhardtii*.

This was my Bachelor of Science Thesis Project. 

## Overview of Polyadenylation

Polyadenylation is a vital post-transcriptional modification in eukaryotes that occurs in the nucleus to protect messenger RNA (mRNA) from degradation during cytoplasmic transport.

### Key Mechanism
The process follows a distinct **two-step pathway**:
1. **Cleavage:** The mRNA strand is cleaved at a specific site.
2. **Tail Addition:** A poly(A) tail is appended directly to the newly formed 3'-end.

## Dataset

Both Positive and Negative datasets were collected separately.

### Positive Dataset

Positive Dataset was extracted as per the following workflow:

![alt text](https://github.com/shahrukh27/CPAP/blob/main/Positive%20Dataset.png)

#### Data Sources: 

https://doi.org/10.1534/g3.114.010249
https://doi.org/10.1534/genetics.108.088971

### Negative Dataset

Negative Dataset was extracted as per the following workflow:

![alt text](https://github.com/shahrukh27/CPAP/blob/main/Negative%20Dataset.png)

## Model Development 

Following were the 4 components of model development:
1. Feature Generation
2. Feature Selection and Model Development 
3. Model Performance Evaluation 
4. Performance Evaluation and Comparison with PASS 2.0 

### Feature Generation

Total of 65 features were generated based on the following data points:

1. Finding novel motifs in the NUE, FUE, DSE, and CE regions around poly (A) sites in addition to the previously reported motifs
  • Near Upstream Element (NUE) region [101:140]
  • Far Upstream Element (FUE) region [1:100]
  • Cleavage Element (CE) region [121:160]
  • Downstream Element (DSE) region [141:190]
3. Using multiple physiochemical and conformational properties of the poly (A) site regions
4. Cleavage site composition
5. Combinations of various K-mer patterns in NUE, FUE, CE, and DSE regions
   
![alt_text](https://github.com/shahrukh27/CPAP/blob/main/Feature%20Generation.png)

### Feature Selection and Model Development

Feature selection resulted in 56 features to be used in model development. 

![alt_text](https://github.com/shahrukh27/CPAP/blob/main/Feature%20Selection.png)

Final model was developed with the remaining 56 features.

### Model Performance Evaluation 

Model was evaluated using sensitivity and specificity. Following formulas were used for calculation:

$$ \text{Sensitivity} = \frac{\text{True Positives (TP)}}{\text{True Positives (TP)} + \text{False Negatives (FN)}} $$

$$ \text{Specificity} = \frac{\text{True Negatives (TN)}}{\text{True Negatives (TN)} + \text{False Positives (FP)}} $$

The error margin was ±20 i.e., if the predicted poly (A) site is within ±20 nts of a specified location then it was considered as positive prediction. CPAP had a sensitivity and specificity of **79%** and **41%** on testing dataset. 

### Performance Evaluation and Comparison with PASS 2.0 

PASS 2.0 was also scored on the same testing dataset. Comparison was done on model evaluations of both CPAP and PASS 2.0. 

![alt_text](https://github.com/shahrukh27/CPAP/blob/main/Sensitivity%20Comparison.png)

CPAP performed 6% better than PASS 2.0 on overall testing dataset. It performed better in individual regions as well. 

Lower specificity showed higher number of false positives. To counter that, we reduced the prediction error margin ad evaluated both models again as follows:

![alt_text](https://github.com/shahrukh27/CPAP/blob/main/Performance%20Evaluation.png)

We observed that individually performance of both models dropped, however, if we can pool the predictions of both models together we can still achieve 70% sensitivty.

## Conclusion

• **New Prediction Tool**: Developed CPAP, a new program for poly (A) sites prediction in *Chlamydomonas reinhardtii* using logistic regression for feature selection.
• **Performance Boost**: CPAP outperformed PASS 2.0 across key metrics:
	• Sensitivity (Sn): 6% higher on the combined positive dataset.
	• Specificity (Sp): 5% better on coding sequences (CDS) and 21% better on intronic control sequences.
• **Complementary Approach**: Combining correct predictions from both CPAP and PASS 2.0 yielded the best results. The distinct features used by each program led to a large number of valuable, unique predictions.
• **Proposed Solution**: Because both tools individually showed lower specificity, a complementary pipeline was proposed to reduce the error margin and secure high-quality poly (A) sites with better Sp.

