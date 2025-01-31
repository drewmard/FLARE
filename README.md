<div align="center">
<img src="img/FLARE_logo_text.png" alt="FLARE_logo_text" width="150">
</div>

FLARE integrates ChromBPNet deep learning predictions with evolutionary conservation to help prioritize impactful rare non-coding variants. Evolutionary conservation is a key predictor of disease risk, influenced by diverse variant mechanisms. FLARE aims to predict PhyloP conservation scores using TSS distance, nearest gene constraint, peak overlap, ChromBPNet scores, and ChromBPNet scores conditional on the variant residing within a peak. Since PhyloP scores are, by definition, not context-specific, we expect FLARE to model the relationship between genomic context, predicted regulatory effects, and evolutionary conservation specifically in cell contexts where regulation is highly relevant to conservation. 

<div align="center">
<img src="img/FLARE_schematic_2.png" alt="FLARE Schematic2" width="375">
</div>

Overall, FLARE:

- (i) disentangles the contributions of context-specific accessibility and regulatory effects to conservation, 
- (ii) provides an intuitive framework for integrating multiple functional genomic features into a unified model, and 
- (iii) captures variants’ regulatory potential across multiple cell types through its predictions. 

A powerful advantage of FLARE is that the model can be trained on any context of interest, such as  different tissues and developmental contexts, by using ChromBPNet predictions across 8,757,029 ultra-rare variants in 1KG.

<div align="center">
<img src="img/FLARE_schematic_1.png" alt="FLARE Schematic1" width="500">
</div>

In this repository, we provide code for FLARE model training and FLARE predictions. Below, we provide an example workflow using ASD de novo mutations for training and prediction since this is a much smaller set of variants to test whether the code is working. However, the results described in the manuscript from Marderstein, Kundu et al. are done by training on the rare variants matrix, and making predictions in the ASD set.

## Step 0.

Make sure that the following packages can be loaded in `R`:

```
library(data.table)
library(optparse)
library(glmnet)
```

Conda environments can be used, e.g. `conda activate r`.

Next, download the input data from the Synapse directory: https://www.synapse.org/Synapse:syn64693551/files/ 

Finally, you need a lot of compute to read the rare variant matrices into R! So, request it. The uncompressed rare variant file is over 40G! This should be plenty of memory, requested interactively:

```
srun --account=smontgom --partition=batch --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
```

## Step 1. Create input feature matrix

From the FLARE root directory:

```
cd scripts
input="/oak/stanford/groups/smontgom/amarder/synapse/predictions/asd.all_dataset.K562_bias.annot2.txt.gz"
output="/oak/stanford/groups/smontgom/amarder/FLARE/data/ASD.FLARE-fb.txt"
model="fetal_brain"
./FLARE_Preprocess.R -i $input -o $output -m $model
```

7 different types of models can be specified by default:

```
"baseline" (baseline model)
"fetal_brain_peaksonly", (FLARE-fetal brain excluding ChromBPNet predictions)
"fetal_brain", (FLARE-fetal brain)
"adult_brain" (FLARE-adult brain)
"brain", (FLARE-brain)
"heart", (FLARE-heart)
"all" (FLARE-all)
```

## Step 2. Train FLARE models

FLARE fits a lasso regression model on the input feature matrix to predict PhyloP conservation scores. Regularization hyperparameters are tuned with 4-fold cross-validation. For each chromosome, FLARE scores are predicted from a model trained on the remaining chromosomes.

```
input="/oak/stanford/groups/smontgom/amarder/FLARE/data/ASD.FLARE-fb.txt"
outdir="/oak/stanford/groups/smontgom/amarder/FLARE/models/ASD.FLARE-fb"
model="fetal_brain"
./FLARE_Training.R -i $input -o $outdir
```

## Step 3. Compute FLARE scores for new variants

```
input="/oak/stanford/groups/smontgom/amarder/FLARE/data/ASD.FLARE-fb.txt"
output="/oak/stanford/groups/smontgom/amarder/FLARE/predictions/ASD.FLARE-fb.ASD.txt"
modelpath="/oak/stanford/groups/smontgom/amarder/FLARE/models/ASD.FLARE-fb"
./FLARE_Predict.R -i $input -o $output -m $modelpath
```

## Example application to de novo mutations in autism families.

