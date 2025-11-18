#!/usr/bin/env Rscript

# srun --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# conda activate r

# specify model as one of:
# "baseline"
# "fetal_brain_peaksonly",
# "fetal_brain",
# "adult_brain"
# "brain",
# "heart",
# "all"

# load
library(data.table)
library(optparse)

################################################################################

# Functions:
snp_identifier = "snp_id"
snp_identifier = "variant_id"

initial_data_load = function(f.input_file,
                             path_to_shet = "/data1/offitk/mardera1/data/neuro_variants_notsynapse/s_het_hgnc_data.txt"
                             ) {
  # read dataframe
  cat("Reading data...\n")
  df = fread(f.input_file,data.table = F,stringsAsFactors = F)
  
  # Use log10 distance to nearest TSS in the model:
  df$gene_distance_1.log10 = log10(df$gene_distance_1+1)
  
  # THESE ARE STEPS REQUIRED FOR TRAINING, OTHERWISE IGNORED:
  if ("phylop" %in% colnames(df)) {
    # Filter SNPs with missing PhyloP values
    cat("Filter SNPs with missing PhyloP values... \n")
    ind = !is.na(df$phylop)
    df = df[ind,]
  }
  
  if ("Consequence" %in% colnames(df)) {
    # Create dummy variables for the categorical variable
    ind = !grepl("splice",df$Consequence)
    cat("Removing ",nrow(df) - sum(ind)," splice-related SNPs...\n")
    df = df[ind,]
  }
  
  if (!("s_het_1") %in% colnames(df)) {
    shet = fread(path_to_shet,data.table = F,stringsAsFactors = F)
    df = merge(df,shet,by.x="closest_gene_1",by.y="hgnc",all.x=TRUE)
    df$s_het_1[is.na(df$s_het_1)] = mean(shet$s_het_1,na.rm=TRUE)
  }
  
  # return data
  return(df)
}


make_FLARE_input_data = function(df,model) {
  
  # select model
  print("Model")
  print("...")
  
  baseline_cols = c("s_het_1","gene_distance_1.log10")
  
  if (model=="baseline") {
    cols = baseline_cols
  } else if (model=="fetal_brain_peaksonly") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    peak_cols = grep("trevino_2021|domcke_2020",peak_cols,value = TRUE)
    peak_cols <- grep("heart", peak_cols, invert = TRUE, value = TRUE)
    summary_cols = c("num_peaks_fb")
    # cols = c(baseline_cols,peak_cols,summary_cols)
    cols = c(baseline_cols,peak_cols)
  } else if (model=="fetal_brain") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    peak_cols = grep("trevino_2021|domcke_2020",peak_cols,value = TRUE)
    peak_cols <- grep("heart", peak_cols, invert = TRUE, value = TRUE)
    cbp_cols <- grep("abs_logfc.mean", colnames(df), value = TRUE)
    cbp_cols <- grep("pval", cbp_cols, invert = TRUE, value = TRUE)
    cbp_cols = grep("trevino_2021|domcke_2020",cbp_cols,value = TRUE)
    cbp_cols <- grep("heart", cbp_cols, invert = TRUE, value = TRUE)
    cols = c(baseline_cols,peak_cols,cbp_cols)
  } else if (model=="adult_brain") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    peak_cols = grep("corces_2020",peak_cols,value = TRUE)
    cbp_cols <- grep("abs_logfc.mean", colnames(df), value = TRUE)
    cbp_cols <- grep("pval", cbp_cols, invert = TRUE, value = TRUE)
    cbp_cols = grep("corces_2020",cbp_cols,value = TRUE)
    summary_cols = c("num_cbp_ab","num_peaks_ab","cbp_max_score_ab")
    # cols = c(baseline_cols,peak_cols,cbp_cols,summary_cols)
    cols = c(baseline_cols,peak_cols,cbp_cols)
  } else if (model=="brain") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    peak_cols = grep("corces_2020|trevino_2021|domcke_2020",peak_cols,value = TRUE)
    peak_cols <- grep("heart", peak_cols, invert = TRUE, value = TRUE)
    cbp_cols <- grep("abs_logfc.mean", colnames(df), value = TRUE)
    cbp_cols <- grep("pval", cbp_cols, invert = TRUE, value = TRUE)
    cbp_cols = grep("corces_2020|trevino_2021|domcke_2020",cbp_cols,value = TRUE)
    cbp_cols <- grep("heart", cbp_cols, invert = TRUE, value = TRUE)
    summary_cols = grep("num_cbp_|num_peaks_|cbp_max_score_",colnames(df),value = TRUE)
    # cols = c(baseline_cols,peak_cols,cbp_cols,summary_cols)
    cols = c(baseline_cols,peak_cols,cbp_cols)
  } else if (model=="all") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    cbp_cols <- grep("abs_logfc.mean", colnames(df), value = TRUE)
    cbp_cols <- grep("pval", cbp_cols, invert = TRUE, value = TRUE)
    summary_cols = grep("num_cbp|num_peaks|cbp_max_score",colnames(df),value = TRUE)
    # cols = c(baseline_cols,peak_cols,cbp_cols,summary_cols)
    cols = c(baseline_cols,peak_cols,cbp_cols)
  } else if (model=="heart") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    peak_cols = grep("domcke_2020|ameen_2022|encode_2024",peak_cols,value = TRUE)
    peak_cols <- grep("brain", peak_cols, invert = TRUE, value = TRUE)
    cbp_cols <- grep("abs_logfc.mean", colnames(df), value = TRUE)
    cbp_cols <- grep("pval", cbp_cols, invert = TRUE, value = TRUE)
    cbp_cols = grep("domcke_2020|ameen_2022|encode_2024",cbp_cols,value = TRUE)
    cbp_cols <- grep("brain", cbp_cols, invert = TRUE, value = TRUE)
    summary_cols = grep("num_cbp_|num_peaks_|cbp_max_score_",colnames(df),value = TRUE)
    # cols = c(baseline_cols,peak_cols,cbp_cols,summary_cols)
    cols = c(baseline_cols,peak_cols,cbp_cols)
  }
  
  ##############################################################################
  # Interaction terms:
  cat("Creating interaction terms...\n")
  
  # This line excludes all columns in df whose names start with int_. 
  # The code resets df in each iteration to a version of df without those columns.
  df = df[,!(colnames(df) %in% grep("^int_",colnames(df),value=TRUE))] # reset in a for loop
  if (!(model %in% c("baseline","fetal_brain_peaksonly"))) {
    for (k in 1:length(peak_cols)) {
      # cat(k,"/",length(peak_cols),"\n",sep = '')
      df[,paste0("int_",k)] = df[,cbp_cols[k]] * df[,peak_cols[k]]
    }
    cols = c(cols,paste0("int_",1:length(peak_cols)))
  }
  
  if ("phylop" %in% colnames(df)) {
    return(df[,c(snp_identifier,"chr","phylop",cols)])
  }
  else {
    return(df[,c(snp_identifier,"chr",cols)])
  }
}


################################################################################

# Define the command-line arguments
option_list = list(
  make_option(c("-i", "--input_file"), type = "character",
              help = "Input file."),
  make_option(c("-o", "--output_file"), type = "character",
              help = "Output file."),
  make_option(c("-m", "--model"), type = "character", default = "fetal_brain",
              help = "Model type.")
)

# Parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

################################################################################

# Process input data
df = initial_data_load(opt$input_file)
df = make_FLARE_input_data(df,opt$model)

# Save
fwrite(df,opt$output_file,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)



