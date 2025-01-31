make_FLARE_input_data = function(df,model) {
  
  # select model
  print("Model")
  print("...")
  
  baseline_cols = c("s_het_1","gene_distance_1.log10")
  
  if (model=="baseline") {
    cols = baseline_cols
  } else if (model=="baseline + fb peaks") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    peak_cols = grep("trevino_2021|domcke_2020",peak_cols,value = TRUE)
    peak_cols <- grep("heart", peak_cols, invert = TRUE, value = TRUE)
    summary_cols = c("num_peaks_fb")
    # cols = c(baseline_cols,peak_cols,summary_cols)
    cols = c(baseline_cols,peak_cols)
  } else if (model=="baseline + fb peaks + cbp") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    peak_cols = grep("trevino_2021|domcke_2020",peak_cols,value = TRUE)
    peak_cols <- grep("heart", peak_cols, invert = TRUE, value = TRUE)
    cbp_cols <- grep("abs_logfc.mean", colnames(df), value = TRUE)
    cbp_cols <- grep("pval", cbp_cols, invert = TRUE, value = TRUE)
    cbp_cols = grep("trevino_2021|domcke_2020",cbp_cols,value = TRUE)
    cbp_cols <- grep("heart", cbp_cols, invert = TRUE, value = TRUE)
    cols = c(baseline_cols,peak_cols,cbp_cols)
  } else if (model=="baseline + ab peaks + cbp") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    peak_cols = grep("corces_2020",peak_cols,value = TRUE)
    cbp_cols <- grep("abs_logfc.mean", colnames(df), value = TRUE)
    cbp_cols <- grep("pval", cbp_cols, invert = TRUE, value = TRUE)
    cbp_cols = grep("corces_2020",cbp_cols,value = TRUE)
    summary_cols = c("num_cbp_ab","num_peaks_ab","cbp_max_score_ab")
    # cols = c(baseline_cols,peak_cols,cbp_cols,summary_cols)
    cols = c(baseline_cols,peak_cols,cbp_cols)
  } else if (model=="baseline + brain peaks + cbp") {
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
  } else if (model=="complete") {
    peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
    cbp_cols <- grep("abs_logfc.mean", colnames(df), value = TRUE)
    cbp_cols <- grep("pval", cbp_cols, invert = TRUE, value = TRUE)
    summary_cols = grep("num_cbp|num_peaks|cbp_max_score",colnames(df),value = TRUE)
    # cols = c(baseline_cols,peak_cols,cbp_cols,summary_cols)
    cols = c(baseline_cols,peak_cols,cbp_cols)
  } else if (model=="baseline + heart peaks + cbp") {
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
  if (!(model %in% c("baseline","baseline + fb peaks"))) {
    for (k in 1:length(peak_cols)) {
      # cat(k,"/",length(peak_cols),"\n",sep = '')
      df[,paste0("int_",k)] = df[,cbp_cols[k]] * df[,peak_cols[k]]
    }
    cols = c(cols,paste0("int_",1:length(peak_cols)))
  }
  
  return(df[,c("snp_id","chr","phylop",cols)])
}
