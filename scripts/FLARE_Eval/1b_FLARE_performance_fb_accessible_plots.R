library(data.table)
library(ggplot2)
df = fread("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/flare/num_peaks_fb.asd.txt",data.table = F,stringsAsFactors = F)
df = df[,-2]
df.melt = melt(df,id.vars="num_cells")
f.out = ("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/flare/flare_predictive_accuracy.r2.pdf")
pdf(f.out,width = 5,height = 5)
ggplot(df.melt,aes(x=100*num_cells/30,y=value^2,col=variable)) +
  # geom_line() +
  geom_point() +
  geom_smooth(fill='lightgrey',alpha=0.25) +
  scale_color_manual(values=c("FLARE_baseline" = "#66BBBB",
                              "FLARE_fb_peaks" = "black",
                              "FLARE_fb" =  "#E0CA70"),
                     labels=c("FLARE_baseline" = "Baseline",
                              "FLARE_fb_peaks" = "Fetal brain (peaks)",
                              "FLARE_fb" =  "Fetal brain")
  ) +
  ggpubr::theme_pubr() +
  scale_x_continuous(breaks=c(0,25,50,75)) +
  labs(x="Minimum % fb accessible contexts",y=expression("Estimated R"^2),title = "SNP conservation prediction") +
  ylim(0,0.05) +
  theme(plot.title = element_text(hjust=0.5),  legend.position = c(0.05, 0.95),
        legend.justification = c("left", "top"),legend.title = element_blank(),legend.background = element_blank())
dev.off()

f.out = ("/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/flare/flare_predictive_accuracy.r.pdf")
pdf(f.out,width = 5,height = 5)
ggplot(df.melt,aes(x=100*num_cells/30,y=value,col=variable)) +
  # geom_line() +
  geom_point() +
  geom_smooth(fill='lightgrey',alpha=0.25) +
  scale_color_manual(values=c("FLARE_baseline" = "#66BBBB",
                              "FLARE_fb_peaks" = "black",
                              "FLARE_fb" =  "#E0CA70"),
                     labels=c("FLARE_baseline" = "Baseline",
                              "FLARE_fb_peaks" = "Fetal brain (peaks)",
                              "FLARE_fb" =  "Fetal brain")
  ) +
  ggpubr::theme_pubr() +
  scale_x_continuous(breaks=c(0,25,50,75)) +
  labs(x="Minimum % fb accessible contexts",y=expression("PhyloP Corr. (Pred vs Obs)"),title = "SNP conservation prediction") +
  # ylim(0,0.05) +
  theme(plot.title = element_text(hjust=0.5),  legend.position = c(0.05, 0.95),
        legend.justification = c("left", "top"),legend.title = element_blank(),legend.background = element_blank())
dev.off()


df$R2_ratio = (df[,3]^2) / (df[,2]^2)
