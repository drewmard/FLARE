# srun --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# conda activate r

# load
library(data.table)

variantSet="asd"
# variantSet="rare"
bias="K562_bias"
i=3  

initial_data_load = function(variantSet) {
  # read dataframe
  cat("Reading data...\n")
  f=paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/",variantSet,".","all_dataset",".",bias,".annot2.txt")
  df = fread(f,data.table = F,stringsAsFactors = F)
  return(df)
}

for (variantSet in c("asd","common","rare")) {
  cat(variantSet)
  df = initial_data_load(variantSet)
  peak_cols = grep("peak_overlap.", colnames(df), value = TRUE)
  peak_cols = grep("trevino_2021|domcke_2020",peak_cols,value = TRUE)
  peak_cols <- grep("heart", peak_cols, invert = TRUE, value = TRUE)
  cat("num_peaks_fb...")
  df$num_peaks_fb = apply(df[,peak_cols],1,sum)
  cat("saving...")
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/etc/num_peaks_fb.",variantSet,".txt")
  fwrite(df[,c("snp_id","num_peaks_fb")],f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}





library(data.table)

for (variantSet in c("asd","rare","common")) { 
  print(variantSet)
  df1 = fread(paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/",variantSet,".FLARE.txt"),data.table = F,stringsAsFactors = F)
  # df2 = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/rare.FLARE.txt",data.table = F,stringsAsFactors = F)
  # df3 = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/common.FLARE.txt",data.table = F,stringsAsFactors = F)
  
  f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/etc/num_peaks_fb.",variantSet,".txt")
  df1.meta = fread(f,data.table = F,stringsAsFactors = F)
  # f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/etc/num_peaks_fb.rare.txt")
  # df2.meta = fread(f,data.table = F,stringsAsFactors = F)
  
  # cor_result = cor.test(df1$phylop,df1[,"FLARE_fb"])
  # cor_result
  tmp = list()
  for (num_cells in 0:(ceiling(0.75*max(df1.meta$num_peaks_fb)))) {
    print(num_cells)
    df1.sub = subset(df1,snp_id %in% subset(df1.meta,num_peaks_fb >= num_cells)$snp_id)
    cor_result1 = cor(df1.sub$phylop,df1.sub[,"FLARE_baseline"])
    cor_result1
    cor_result2 = cor(df1.sub$phylop,df1.sub[,"FLARE_fb_peaks"])
    cor_result2
    cor_result3 = cor(df1.sub$phylop,df1.sub[,"FLARE_fb"])
    cor_result3
    tmp[[num_cells+1]] = data.frame(num_cells,FLARE_baseline=cor_result1,FLARE_fb_peaks=cor_result2,FLARE_fb=cor_result3)
  }
  out = as.data.frame(do.call(rbind,tmp))
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/flare/num_peaks_fb.",variantSet,".txt")
  fwrite(out,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  
}


