#!/usr/bin/env Rscript

library(data.table)
library(glmnet)
library(optparse)

split_arg = function(x) {
  if (is.null(x) || is.na(x) || x == "") {
    return(character(0))
  }
  trimws(unlist(strsplit(x, ",")))
}

read_named_arg = function(x) {
  values = split_arg(x)
  if (length(values) == 0) {
    return(setNames(character(0), character(0)))
  }

  split_pos = regexpr("=", values, fixed = TRUE)
  if (all(split_pos < 1)) {
    names_out = basename(normalizePath(values, mustWork = FALSE))
    names(values) = names_out
    return(values)
  }
  if (any(split_pos < 1)) {
    stop("Do not mix plain paths and name=path values in --model-dirs.")
  }

  names_out = substr(values, 1, split_pos - 1)
  paths_out = substr(values, split_pos + 1, nchar(values))
  if (any(names_out == "") || any(paths_out == "")) {
    stop("Model names and paths cannot be empty.")
  }

  names(paths_out) = names_out
  paths_out
}

feature_type = function(feature) {
  out = rep("other", length(feature))
  out[feature %in% c("s_het_1", "gene_distance_1.log10")] = "base"
  out[grepl("^peak_overlap\\.", feature)] = "peak"
  out[grepl("logfc\\.mean|abs_logfc\\.mean", feature)] = "cbp"
  out[grepl("^int_", feature)] = "int"
  out
}

read_model_weights = function(model_dir, chr_values) {
  beta_lst = list()

  for (chr_num in chr_values) {
    f = file.path(model_dir, paste0("flare.chr", chr_num, ".rds"))
    if (!file.exists(f)) {
      stop(paste0("Missing model file: ", f))
    }
    final_mod = readRDS(f)
    if (!is.null(final_mod$beta)) {
      beta = as.matrix(final_mod$beta)
    } else {
      beta = as.matrix(coef(final_mod))[-1, , drop = FALSE]
    }
    beta_lst[[as.character(chr_num)]] = beta
    cat("Read ", basename(f), ": ", nrow(beta), " features\n", sep = "")
  }

  all_features = unique(unlist(lapply(beta_lst, rownames)))
  if (length(all_features) == 0) {
    stop(paste0("No coefficient rows found in model directory: ", model_dir))
  }

  beta_mat = matrix(0, nrow = length(all_features), ncol = length(beta_lst),
                    dimnames = list(all_features, names(beta_lst)))

  for (chr_name in names(beta_lst)) {
    beta = beta_lst[[chr_name]]
    beta_mat[rownames(beta), chr_name] = beta[, 1]
  }

  data.frame(
    feature = rownames(beta_mat),
    mean_weight = rowMeans(beta_mat),
    nonzero_chromosomes = rowSums(beta_mat != 0),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

option_list = list(
  make_option(c("-m", "--model-dirs"), type = "character",
              dest = "model_dirs",
              help = "Model dir path, comma-separated model dir paths, or comma-separated name=path model dirs."),
  make_option(c("-o", "--output"), type = "character",
              help = "Output lasso weights TSV."),
  make_option(c("--chromosomes"), type = "character", default = paste(1:22, collapse = ","),
              help = "Comma-separated chromosome model suffixes to read. Default: 1,2,...,22.")
)

opt = parse_args(OptionParser(option_list = option_list))

model_dirs = read_named_arg(opt$model_dirs)
chr_values = split_arg(opt$chromosomes)

if (length(model_dirs) == 0) {
  stop("No model directories were provided. Use --model-dirs name=/path/to/models.")
}
if (length(chr_values) == 0) {
  stop("No chromosomes were provided.")
}

res_lst = list()
for (model_name in names(model_dirs)) {
  cat("Reading model directory '", model_name, "': ", model_dirs[[model_name]], "\n", sep = "")
  weights = read_model_weights(model_dirs[[model_name]], chr_values)
  weights$model = model_name
  weights$type = feature_type(weights$feature)
  res_lst[[model_name]] = weights[, c("model", "feature", "type", "mean_weight", "nonzero_chromosomes")]
}

res = rbindlist(res_lst, fill = TRUE)
cat("Writing ", nrow(res), " feature rows to ", opt$output, "\n", sep = "")
fwrite(res, opt$output, quote = FALSE, na = "NA", sep = "\t", row.names = FALSE, col.names = TRUE)
