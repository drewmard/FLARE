#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(optparse)

split_arg = function(x) {
  if (is.null(x) || is.na(x) || x == "") {
    return(character(0))
  }
  trimws(unlist(strsplit(x, ",")))
}

pretty_variant_set = function(x) {
  out = gsub("_", " ", x)
  out = tools::toTitleCase(tolower(out))
  out[out == "Asd"] = "ASD"
  out
}

option_list = list(
  make_option(c("-i", "--input"), type = "character",
              help = "Performance TSV from 1b_FLARE_performance.R."),
  make_option(c("-o", "--output"), type = "character",
              help = "Output PDF path."),
  make_option(c("--model-order"), type = "character", default = NULL,
              dest = "model_order",
              help = "Optional comma-separated model order."),
  make_option(c("--variant-order"), type = "character", default = NULL,
              dest = "variant_order",
              help = "Optional comma-separated variantSet order."),
  make_option(c("--width"), type = "double", default = 5.6,
              help = "PDF width."),
  make_option(c("--height"), type = "double", default = 3.2,
              help = "PDF height.")
)

opt = parse_args(OptionParser(option_list = option_list))

df = fread(opt$input, data.table = FALSE, stringsAsFactors = FALSE)
required_cols = c("variantSet", "model", "r", "l", "h")
missing_cols = setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop(paste0("Missing required columns: ", paste(missing_cols, collapse = ", ")))
}

model_order = split_arg(opt$model_order)
if (length(model_order) > 0) {
  df$model = factor(df$model, levels = model_order)
}

df$variantSetLabel = pretty_variant_set(df$variantSet)
variant_order = split_arg(opt$variant_order)
if (length(variant_order) > 0) {
  df$variantSetLabel = factor(df$variantSetLabel, levels = pretty_variant_set(variant_order))
}

p = ggplot(df, aes(x = model, y = r, fill = model)) +
  geom_col(color = "black", width = 0.75) +
  geom_errorbar(aes(ymin = l, ymax = h), width = 0.2) +
  facet_wrap(. ~ variantSetLabel, ncol = min(3, length(unique(df$variantSetLabel)))) +
  labs(x = "Model", y = expression("PhyloP Corr. (Pred vs Obs)"), title = "Model Comparison", fill = "Model") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

pdf(opt$output, width = opt$width, height = opt$height)
print(p)
dev.off()
