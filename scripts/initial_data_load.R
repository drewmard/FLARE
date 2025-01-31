initial_data_load = function(f.input_file) {
  # read dataframe
  cat("Reading data...\n")
  df = fread(f.input_file,data.table = F,stringsAsFactors = F)
  
  # Filter SNPs with missing PhyloP values
  cat("Filter SNPs with missing PhyloP values... \n")
  ind = !is.na(df$phylop)
  df = df[ind,]
  
  cat("Performing misc. filtering... \n")
  
  # Use log10 distance to nearest TSS in the model:
  df$gene_distance_1.log10 = log10(df$gene_distance_1+1)
  
  # Create dummy variables for the categorical variable
  ind = !grepl("splice",df$Consequence)
  cat("Removing ",nrow(df) - sum(ind)," splice-related SNPs...\n")
  df = df[ind,]
  
  # return data
  return(df)
}
