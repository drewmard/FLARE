#!/usr/bin/env Rscript

# srun --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# conda activate r

library(glmnet)
library(data.table)
library(optparse)

FLARE_Predict = function(f.input,f.output,f.modelpath) {
  
  # Throw error if "snp_id", "chr", or "phylop" not in feature matrix
  df = fread(f.input,data.table = F,stringsAsFactors = F)
  if (!("snp_id" %in% colnames(df))) {stop('"snp_id" column is missing!')}
  if (!("chr" %in% colnames(df))) {stop('"chr" column is missing!')}
  if (!("phylop" %in% colnames(df))) {stop('"phylop" column is missing!')}
  
  predictions_df.all = list(); cor_result = list()
  
  for (chrNum in 1:22) {
    cat("Making FLARE predictions using LOCO schema (chr ",chrNum,")...\n",sep = '')
    
    # Subset to chr of interest
    ind_chr_include = df$chr==paste0("chr",chrNum)
    
    # Build feature matrix
    cols_exclude = !(colnames(df) %in% c("chr","phylop","snp_id"))
    x = as.matrix(df[ind_chr_include,cols_exclude])

    # Load model
    f = paste0(f.modelpath,"/flare.chr",chrNum,".rds")
    final_mod = readRDS(f)
    
    # Making FLARE predictions
    predictions_lasso <- as.numeric(predict(final_mod, newx = x)[,1])
    
    # Store predictions
    predictions_df = data.frame(snp_id = df[ind_chr_include,"snp_id"],FLARE=predictions_lasso)
    
    # Store results
    predictions_df.all[[chrNum]] = predictions_df
  }
  
  # Save predictions
  predictions_df = as.data.frame(do.call(rbind,predictions_df.all))
  fwrite(predictions_df,f.output,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}

################################################################################

# Define the command-line arguments
option_list = list(
  make_option(c("-i", "--input"), type = "character",
              help = "File path for input matrix."),
  make_option(c("-o", "--output"), type = "character",
              help = "File path for predictions."),
  make_option(c("-m","--modelpath"),type = "character",
              help = "File path for saved models.")
)

# Parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Run:
FLARE_Predict(opt$input,opt$output,opt$modelpath)

