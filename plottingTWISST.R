library(cowplot)

dev.off()
for (i in c(1:17, "Z")) {
  name.windows <- paste0("phyml__allspp_25kb_SUPER_",i,".data.tsv")
  name.weights <- paste0("by_chrom/SUPER_",i,"_NJ_wurs_25kb.weights.csv")
  windows <- read_delim(name.windows, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  weights <- read_table(name.weights, skip = 10395)

  weights_clean <- data.frame(
    species_tree = weights$topo8203,
    ILS = weights$topo8202,
    mitochondrial = weights$topo3589,
    Vseo_intro = weights$topo3772 + weights$topo9358 + weights$topo3043 + weights$topo9486 + weights$topo3007 + weights$topo9353 + weights$topo10198 + weights$topo3409)
  df <- weights_clean / 256 ##divided by 256 to have weigthings up to 1. If the weigthing is different from 1, I won't take it into account for the plot.
  df[df < 1] <- NA
    # remove NA rows
  df <- replace(df, is.na(df), 0)
  df$mid <- windows$mid 
  df$mid <- df$mid / 1000000

  df$combi <- ifelse(df$species_tree == 1, "species_tree",
                     ifelse(df$mitochondrial == 1, "mitochondrial",
                            ifelse(df$Vseo_intro == 1, "Vseo_intro",
                                   ifelse(df$ILS == 1, "ILS", "others"))))
  # Plot
  chromosome <- ggplot() +
    labs(x = "Chromosome position (Mb)") +
    geom_bar(data = df, aes(x = mid, y = Vseo_intro, fill = factor(Vseo_intro)), stat = "identity", width = 0.0001, alpha=0.7) +
    scale_fill_manual(values = c("black", "#E95635", "#E95635", "#E95635", "#E95635", "#E95635", "#E95635", "#E95635", "#E95635", "#E95635", "#E95635"), labels = c("0", "1")) +
    theme_classic() + guides(fill = FALSE) + 
    scale_y_continuous(breaks = c(0, 1), name = paste0("Chr. ", i)) +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5, size = 10), 
          axis.title.x = if (i == "Z") element_text(size = 11,margin = margin(t = 4.5)) else element_blank(), 
          axis.text.x = if (i == "Z") element_text(size = 8) else element_blank(), 
          axis.ticks.x = if (i == "Z") element_line() else element_blank(),  
          axis.line.x = if (i == "Z") element_line() else element_blank(), 
          axis.text.y = element_text(size = 6.5)) +
    xlim(0, 370.562366) 
  chromosome2 <- chromosome + geom_rect(data = df, aes(xmin = 0, xmax = max(mid), ymin = 0, ymax = 1), 
                                            fill = NA, color = "black", size = 0.6)  
    
  assign(paste0("chr_vseo_", i), chromosome2)
}

plots = list(chr_vseo_1, chr_vseo_2, chr_vseo_3, chr_vseo_4, chr_vseo_5, chr_vseo_6,chr_vseo_7, chr_vseo_8, chr_vseo_9, chr_vseo_10, chr_vseo_11, chr_vseo_12, chr_vseo_13, chr_vseo_14,chr_vseo_15, chr_vseo_16, chr_vseo_17, chr_vseo_Z)
altura_por_grafico <- rep(1/19, length(plots))
altura_por_grafico[length(plots)] <- 1.82/19  


plot_grid(plotlist = plots, ncol = 1, align = "v", rel_heights = altura_por_grafico, axis = "lr")
