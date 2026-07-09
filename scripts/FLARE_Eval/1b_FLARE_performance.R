#!/usr/bin/env Rscript

library(data.table)
library(optparse)

split_arg = function(x) {
  if (is.null(x) || is.na(x) || x == "") {
    return(character(0))
  }
  trimws(unlist(strsplit(x, ",")))
}

infer_file_label = function(path) {
  label = basename(path)
  label = sub("\\.gz$", "", label)
  label = sub("\\.bgz$", "", label)
  label = sub("\\.txt$", "", label)
  label = sub("\\.tsv$", "", label)
  label = sub("\\.FLARE\\.pred$", "", label)
  label = sub("\\.FLARE$", "", label)
  label = sub("\\.pred$", "", label)

  parts = unlist(strsplit(label, "\\."))
  if (length(parts) > 1) {
    label = parts[length(parts)]
  }
  label
}

read_path_arg = function(x) {
  values = split_arg(x)
  if (length(values) == 0) {
    return(setNames(character(0), character(0)))
  }

  split_pos = regexpr("=", values, fixed = TRUE)
  if (all(split_pos < 1)) {
    names(values) = vapply(values, infer_file_label, character(1))
    return(values)
  }
  if (any(split_pos < 1)) {
    stop("Do not mix plain paths and name=path values.")
  }

  names_out = substr(values, 1, split_pos - 1)
  paths_out = substr(values, split_pos + 1, nchar(values))
  if (any(names_out == "") || any(paths_out == "")) {
    stop("Names and paths cannot be empty.")
  }
  names(paths_out) = names_out
  paths_out
}

safe_cor_test = function(df, truth_col, predictor_col) {
  keep = complete.cases(df[, c(truth_col, predictor_col)])
  if (sum(keep) < 3) {
    return(data.frame(r = NA_real_, l = NA_real_, h = NA_real_, n = sum(keep)))
  }

  cor_result = tryCatch(
    cor.test(df[[truth_col]][keep], df[[predictor_col]][keep]),
    error = function(e) NULL
  )
  if (is.null(cor_result)) {
    return(data.frame(r = NA_real_, l = NA_real_, h = NA_real_, n = sum(keep)))
  }
  conf_int = cor_result$conf.int
  if (length(conf_int) < 2) {
    conf_int = c(NA_real_, NA_real_)
  }
  data.frame(
    r = unname(cor_result$estimate),
    l = conf_int[1],
    h = conf_int[2],
    n = sum(keep)
  )
}

performance_table = function(prediction_files,
                              output_file,
                              truth_files = NULL,
                              predictors = NULL,
                              truth_col = "phylop",
                              id_col = "snp_id",
                              mean_output_file = NULL) {
  res_lst = list()
  mean_lst = list()
  k = 0

  for (variant_set in names(prediction_files)) {
    df = fread(prediction_files[[variant_set]], data.table = FALSE, stringsAsFactors = FALSE)

    if (!(truth_col %in% colnames(df))) {
      if (is.null(truth_files) || !(variant_set %in% names(truth_files))) {
        stop(paste0("'", truth_col, "' is missing for ", variant_set,
                    ". Provide it in the prediction file or pass --truth-files ", variant_set, "=path."))
      }
      truth_df = fread(truth_files[[variant_set]], data.table = FALSE, stringsAsFactors = FALSE)
      if (!(id_col %in% colnames(df)) || !(id_col %in% colnames(truth_df))) {
        stop(paste0("'", id_col, "' is required to merge truth values for ", variant_set, "."))
      }
      if (!(truth_col %in% colnames(truth_df))) {
        stop(paste0("'", truth_col, "' is missing from truth file for ", variant_set, "."))
      }
      df = merge(df, truth_df[, c(id_col, truth_col)], by = id_col, all.x = TRUE)
    }

    predictor_cols = predictors
    if (is.null(predictor_cols) || length(predictor_cols) == 0) {
      excluded = c(id_col, truth_col, "chr", "variant_id")
      predictor_cols = names(df)[vapply(df, is.numeric, logical(1))]
      predictor_cols = setdiff(predictor_cols, excluded)
    }
    missing_predictors = setdiff(predictor_cols, colnames(df))
    if (length(missing_predictors) > 0) {
      stop(paste0("Missing predictor columns for ", variant_set, ": ",
                  paste(missing_predictors, collapse = ", ")))
    }

    for (predictor in predictor_cols) {
      k = k + 1
      stats = safe_cor_test(df, truth_col, predictor)
      res_lst[[k]] = data.frame(
        variantSet = variant_set,
        model = predictor,
        stats,
        row.names = NULL
      )
    }

    if (!is.null(mean_output_file)) {
      mean_cols = c(truth_col, predictor_cols)
      y = vapply(df[, mean_cols, drop = FALSE], mean, numeric(1), na.rm = TRUE)
      y2 = vapply(df[, mean_cols, drop = FALSE], sd, numeric(1), na.rm = TRUE)
      n = vapply(df[, mean_cols, drop = FALSE], function(x) sum(!is.na(x)), integer(1))
      mean_lst[[variant_set]] = data.frame(
        variantSet = variant_set,
        metric = names(y),
        mu = y,
        se = y2 / sqrt(n),
        n = n,
        row.names = NULL
      )
    }
  }

  res = as.data.frame(do.call(rbind, res_lst))
  fwrite(res, output_file, quote = FALSE, na = "NA", sep = "\t", row.names = FALSE, col.names = TRUE)

  if (!is.null(mean_output_file)) {
    means = as.data.frame(do.call(rbind, mean_lst))
    fwrite(means, mean_output_file, quote = FALSE, na = "NA", sep = "\t", row.names = FALSE, col.names = TRUE)
  }
}

option_list = list(
  make_option(c("-p", "--prediction-files"), type = "character",
              dest = "prediction_files",
              help = "Prediction file path, comma-separated prediction file paths, or comma-separated name=path prediction files."),
  make_option(c("-o", "--output"), type = "character",
              help = "Output performance TSV."),
  make_option(c("-t", "--truth-files"), type = "character", default = NULL,
              dest = "truth_files",
              help = "Optional file path, comma-separated file paths, or comma-separated name=path files containing snp_id and phylop."),
  make_option(c("--predictors"), type = "character", default = NULL,
              help = "Optional comma-separated predictor columns. Defaults to numeric columns other than IDs/truth."),
  make_option(c("--truth-col"), type = "character", default = "phylop",
              dest = "truth_col",
              help = "Observed truth column."),
  make_option(c("--id-col"), type = "character", default = "snp_id",
              dest = "id_col",
              help = "Variant ID column used for truth merges."),
  make_option(c("--mean-output"), type = "character", default = NULL,
              dest = "mean_output",
              help = "Optional output TSV for means and standard errors.")
)

opt = parse_args(OptionParser(option_list = option_list))

prediction_files = read_path_arg(opt$prediction_files)
truth_files = NULL
if (!is.null(opt$truth_files)) {
  truth_files = read_path_arg(opt$truth_files)
}

performance_table(
  prediction_files = prediction_files,
  output_file = opt$output,
  truth_files = truth_files,
  predictors = split_arg(opt$predictors),
  truth_col = opt$truth_col,
  id_col = opt$id_col,
  mean_output_file = opt$mean_output
)
