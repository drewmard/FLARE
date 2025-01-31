#!/usr/bin/env Rscript

# srun --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# conda activate r

library(glmnet)
library(data.table)
library(optparse)

FLARE_Training = function(f.input,outdir) {
  # input matrix must include "chr" and "phylop" column
  # outdir must be a valid path, where the final folder will be a newly created folder
  
  # Throw error if "snp_id", "chr", or "phylop" not in feature matrix
  df = fread(f.input,data.table = F,stringsAsFactors = F)
  if (!("snp_id" %in% colnames(df))) {stop('"snp_id" column is missing!')}
  if (!("chr" %in% colnames(df))) {stop('"chr" column is missing!')}
  if (!("phylop" %in% colnames(df))) {stop('"phylop" column is missing!')}
  
  for (chrNum in 1:22) {
    
    cat("Training FLARE using leave-one-chromosome-out schema (chr ",chrNum,")...\n",sep = '')
    
    # Curate data:
    ind_chr_exclude = df$chr!=paste0("chr",chrNum)
    cols_exclude = !(colnames(df) %in% c("chr","phylop","snp_id"))
    x = as.matrix(df[ind_chr_exclude,cols_exclude])
    y = df$phylop[ind_chr_exclude]
    
    ############################################################################
    # # Lasso
    set.seed(123)  # Set the seed for reproducibility:
    
    # Set peak, chrombpnet, and gene constraint info to have non-negative constraints.
    constrained_cols = unique(c(
      grep("abs_logfc.mean",colnames(x),value = TRUE),
      grep("cbp",colnames(x),value = TRUE),
      grep("peak",colnames(x),value = TRUE),
      grep("int_",colnames(x),value=TRUE),
      grep("s_het_1",colnames(x),value = TRUE)
    ))
    constrained_limits = rep(-Inf,ncol(x))
    constrained_limits[colnames(x) %in% constrained_cols] = 0
    
    # L1 loss with non-negative constraints - identify lambda values
    alphaUse = 1
    cv_mod <- cv.glmnet(x, y, alpha = alphaUse, lower.limits = constrained_limits, nfolds = 4)  # 4-fold cross-validation
    
    # Fit model with optimal lambda
    best_lambda <- cv_mod$lambda.min
    final_mod <- glmnet(x, y, alpha = alphaUse, lower.limits = constrained_limits, lambda = best_lambda)
    
    # Save model
    # Example output path directory: "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/FLARE/models/flare_fb"
    dir.create(outdir)
    f.out = paste0(outdir,"/flare.chr",chrNum,".rds")
    saveRDS(final_mod,file=f.out)
  }
}

# Define the command-line arguments
option_list = list(
  make_option(c("-i", "--input"), type = "character",
              help = "File path for input matrix."),
  make_option(c("-o", "--outdir"), type = "character",
              help = "File path for models.")
)

# Parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

FLARE_Training(opt$input,opt$outdir)


