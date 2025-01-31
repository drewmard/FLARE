# FLARE
FLARE is a context-specific functional genomic model of constraint that helps prioritize impactful rare non-coding variants.

<img src="img/FLARE_schematic1.png" alt="FLARE Schematic1" width="500">

FLARE is a lasso regression model that integrates deep learning predictions with evolutionary conservation. Evolutionary conservation is a key predictor of disease risk, influenced by diverse variant mechanisms. FLARE aims to predict PhyloP conservation scores using TSS distance, nearest gene constraint, peak overlap, ChromBPNet scores, and ChromBPNet scores conditional on the variant residing within a peak. Since PhyloP scores are, by definition, not context-specific, we expect FLARE to model the relationship between genomic context, predicted regulatory effects, and evolutionary conservation specifically in cell contexts where regulation is highly relevant to conservation. 

![FLARE Schematic2](img/FLARE_schematic_2.png)

Thus, FLARE:

- (i) disentangles the contributions of context-specific accessibility and regulatory effects to conservation, 
- (ii) provides an intuitive framework for integrating multiple functional genomic features into a unified model, and 
- (iii) captures variants’ regulatory potential across multiple cell types through its predictions. 

A powerful advantage of FLARE is that the model can be trained on any context of interest, such as  different tissues and developmental contexts, by using ChromBPNet predictions across 8,757,029 ultra-rare variants in 1KG.

In this repository, we provide code for FLARE model training and predictions. 

## Step 1. Create input feature matrix

```

model="baseline + fb peaks + cbp"
./FLARE_preprocess.R -i $input -o $output -m $model
```

## Step 2. Train FLARE models

FLARE fits a lasso regression model on the input feature matrix to predict PhyloP conservation scores. Regularization hyperparameters are tuned with 4-fold cross-validation. For each chromosome, FLARE scores are predicted from a model trained on the remaining chromosomes.



## Step 3. Compute FLARE scores for new variants

## Example application to de novo mutations in autism families.

