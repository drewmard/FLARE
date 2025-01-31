#!/usr/bin/env Rscript

# srun --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# conda activate r

# specify model as one of:
# "baseline"
# "baseline + fb peaks",
# "baseline + fb peaks + cbp",
# "baseline + ab peaks + cbp"
# "baseline + brain peaks + cbp",
# "baseline + heart peaks + cbp",
# "complete"

# load
library(data.table)
library(optparse)

################################################################################

# Define the command-line arguments
option_list = list(
  make_option(c("-i", "--input_file"), type = "character",
              help = "Input file."),
  make_option(c("-o", "--output_file"), type = "character",
              help = "Output file."),
  make_option(c("-m", "--model"), type = "character", default = "baseline + fb peaks + cbp",
              help = "Model type.")
)

# Parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

################################################################################

# Process input data
df = initial_data_load(f.input_file)
df = make_FLARE_input_data(df,model)

# Save
fwrite(df,opt$output_file,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)



