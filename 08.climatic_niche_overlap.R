#Script by Fernando Martínez-Feiría 2024
library(devtools)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(pca3d)


#upload / creation species/ groups' PCs

clim_VLA<-read.csv("VLA_main_groups.csv",h=T,sep=";")

clim_VSE<-read.csv("VSE_main_groups.csv",h=T,sep=";")

clim_VLAg_int<-read.csv("VLA_gad_x_seoaneiW.csv",h=T,sep=";")

clim_VSEW_int<-read.csv("VSE_seoaneiW_x_gad.csv",h=T,sep=";")

clim_VLAnoad <-read.csv("VLA_no-admix.csv",h=T,sep=";")

clim_VSEnoad <-read.csv("VSE_no-admix.csv",h=T,sep=";")


subsp_VLA <- clim_VLA$group
subsp_VSE <- clim_VSE$group


clim_lat <- subset(clim_VLA, subsp_VLA == "latastei")
clim_gad  <- subset(clim_VLA, subsp_VLA == "gaditana")

clim_seo_W <- subset(clim_VSE, subsp_VSE == "seoanei_W")
clim_seo_E <- subset(clim_VSE, subsp_VSE == "seoanei_E")

clim_VLA_VSE <- rbind(clim_VLA,clim_VSE)
clim_VLA_VSEnoad <- rbind(clim_VLAnoad,clim_VSEnoad)

clim_gad_seo_W <- rbind(clim_gad,clim_seo_W)
clim_lat_seo_E <- rbind(clim_lat,clim_seo_E)

VLAgad_int <- subset(clim_VLAg_int, clim_VLAg_int$group == "gaditana_seoW")
VLAgad_NO_int <- subset(clim_VLAg_int, clim_VLAg_int$group == "gaditana")

VSEW_int <- subset(clim_VSEW_int, clim_VSEW_int$group == "seoanei_W_gad")
VSEW_NO_int <- subset(clim_VSEW_int, clim_VSEW_int$group == "seoanei_W")


VLAg_VSEW_NO_int <- rbind(VLAgad_NO_int,VSEW_NO_int)



## 2D and 3D PCA plots for the comparisons between groups/species

#latastei - seoanei
VLA_VSE_pca1_2 <- ggplot(clim_VLA_VSE, aes(x = pcac1, y = pcac2, color = species)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
VLA_VSE_pca1_3 <- ggplot(clim_VLA_VSE, aes(x = pcac1, y = pcac3, color = species)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
VLA_VSE_pca2_3 <- ggplot(clim_VLA_VSE, aes(x = pcac2, y = pcac3, color = species)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")

legend <- cowplot::get_legend(ggplot(clim_VLA_VSE, aes(x = pcac1, y = pcac2, color = species)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw())


figure <- ggarrange(VLA_VSE_pca1_2, VLA_VSE_pca1_3, VLA_VSE_pca2_3, legend , align = "v", ncol = 2, nrow = 2)
annotate_figure(figure,
                top = text_grob("V. latastei - V. seoanei", color = "black", face = "bold", size = 14))


group<-clim_VLA_VSE$species
pc1<-clim_VLA_VSE$pcac1
pc2<-clim_VLA_VSE$pcac2
pc3<-clim_VLA_VSE$pcac3

data1<-cbind(pc1, pc2, pc3)

pca3d(data1, group = group, palette = c("red", "green"), axes.color = "black", show.plane = FALSE, bg = "white", radius = 1.2, show.ellipses = TRUE, ellipse.ci = 0.95, shape=1, show.shapes=TRUE)

snapshotPCA3d(file="3D_latastei-seoanei.png")


#latastei - seoanei NO ADMIXT

VLA_VSEnoad_pca1_2 <- ggplot(clim_VLA_VSEnoad, aes(x = pcac1, y = pcac2, color = species)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
VLA_VSEnoad_pca1_3 <- ggplot(clim_VLA_VSEnoad, aes(x = pcac1, y = pcac3, color = species)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
VLA_VSEnoad_pca2_3 <- ggplot(clim_VLA_VSEnoad, aes(x = pcac2, y = pcac3, color = species)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")

legend <- cowplot::get_legend(ggplot(clim_VLA_VSEnoad, aes(x = pcac1, y = pcac2, color = species)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw())


figure <- ggarrange(VLA_VSEnoad_pca1_2, VLA_VSEnoad_pca1_3, VLA_VSEnoad_pca2_3, legend , align = "v", ncol = 2, nrow = 2)
annotate_figure(figure,
                top = text_grob("V. latastei - V. seoanei, no admixed", color = "black", face = "bold", size = 14))


group<-clim_VLA_VSEnoad$species
pc1<-clim_VLA_VSEnoad$pcac1
pc2<-clim_VLA_VSEnoad$pcac2
pc3<-clim_VLA_VSEnoad$pcac3

data1<-cbind(pc1, pc2, pc3)

pca3d(data1, group = group, palette = c("red", "green"), axes.color = "black", show.plane = FALSE, bg = "white", radius = 1.2, show.ellipses = TRUE, ellipse.ci = 0.95, shape=1, show.shapes=TRUE)

snapshotPCA3d(file="3D_latastei-seoanei_no_admixed.png")


#latastei-gaditana

lat_gad_pca1_2 <- ggplot(clim_VLA, aes(x = pcac1, y = pcac2, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
lat_gad_pca1_3 <- ggplot(clim_VLA, aes(x = pcac1, y = pcac3, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
lat_gad_pca2_3 <- ggplot(clim_VLA, aes(x = pcac2, y = pcac3, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")

legend <- cowplot::get_legend(ggplot(clim_VLA, aes(x = pcac1, y = pcac2, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw())

figure <- ggarrange(lat_gad_pca1_2, lat_gad_pca1_3, lat_gad_pca2_3, legend , align = "v", ncol = 2, nrow = 2)
annotate_figure(figure,
                top = text_grob("V. l. gaditana - V. l. latastei", color = "black", face = "bold", size = 14))

group<-clim_VLA$group
pc1<-clim_VLA$pcac1
pc2<-clim_VLA$pcac2
pc3<-clim_VLA$pcac3

data2<-cbind(pc1, pc2, pc3)

pca3d(data2, group = group, palette = c("red", "orange"), axes.color = "black", show.plane = FALSE, bg = "white", radius = 1.2, show.ellipses = TRUE, ellipse.ci = 0.95, shape=1, show.shapes=TRUE)

snapshotPCA3d(file="3D_latastei-gaditana.png")


#seoanei_W - E

seo_W_E_pca1_2 <- ggplot(clim_VSE, aes(x = pcac1, y = pcac2, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
seo_W_E_pca1_3 <- ggplot(clim_VSE, aes(x = pcac1, y = pcac3, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
seo_W_E_pca2_3 <- ggplot(clim_VSE, aes(x = pcac2, y = pcac3, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")

legend <- cowplot::get_legend(ggplot(clim_VSE, aes(x = pcac2, y = pcac3, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw())

figure <- ggarrange(seo_W_E_pca1_2, seo_W_E_pca1_3, seo_W_E_pca2_3, legend , align = "v", ncol = 2, nrow = 2)
annotate_figure(figure,
                top = text_grob("seoanei E - seoanei W", color = "black", face = "bold", size = 14))


group<-clim_VSE$group
pc1<-clim_VSE$pcac1
pc2<-clim_VSE$pcac2
pc3<-clim_VSE$pcac3

data3<-cbind(pc1, pc2, pc3)

pca3d(data3, group = group, palette = c("green", "darkgreen"), axes.color = "black", show.plane = FALSE, bg = "white", radius = 1.2, show.ellipses = TRUE, ellipse.ci = 0.95, shape=1, show.shapes=TRUE)

snapshotPCA3d(file="3D_seoanei_W-E.png")


#gaditana - seoanei W
gad_seoW_pca1_2 <- ggplot(clim_gad_seo_W, aes(x = pcac1, y = pcac2, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
gad_seoW_pca1_3 <- ggplot(clim_gad_seo_W, aes(x = pcac1, y = pcac3, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
gad_seoW_pca2_3 <- ggplot(clim_gad_seo_W, aes(x = pcac2, y = pcac3, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")

legend <- cowplot::get_legend(ggplot(clim_gad_seo_W, aes(x = pcac1, y = pcac3, color = species)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw())

figure <- ggarrange(gad_seoW_pca1_2, gad_seoW_pca1_3, gad_seoW_pca2_3, legend , align = "v", ncol = 2, nrow = 2)
annotate_figure(figure,
                top = text_grob("gaditana - seoanei W", color = "black", face = "bold", size = 14))

group<-clim_gad_seo_W$group
pc1<-clim_gad_seo_W$pcac1
pc2<-clim_gad_seo_W$pcac2
pc3<-clim_gad_seo_W$pcac3


data4<-cbind(pc1, pc2, pc3)

pca3d(data4, group = group, palette = c("red", "darkgreen"), axes.color = "black", show.plane = FALSE, bg = "white", radius = 1.2, show.ellipses = TRUE, ellipse.ci = 0.95, shape=1, show.shapes=TRUE)

snapshotPCA3d(file="3D_seoanei_W-gaditana.png")



#latastei - seoanei E
lat_seoE_pca1_2 <- ggplot(clim_lat_seo_E, aes(x = pcac1, y = pcac2, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
lat_seoE_pca1_3 <- ggplot(clim_lat_seo_E, aes(x = pcac1, y = pcac3, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
lat_seoE_pca2_3 <- ggplot(clim_lat_seo_E, aes(x = pcac2, y = pcac3, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")

legend <- cowplot::get_legend(ggplot(clim_lat_seo_E, aes(x = pcac1, y = pcac3, color = species)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw())


figure <- ggarrange(lat_seoE_pca1_2, lat_seoE_pca1_3, lat_seoE_pca2_3, legend , align = "v", ncol = 2, nrow = 2)
annotate_figure(figure,
                top = text_grob("latastei - seoanei E", color = "black", face = "bold", size = 14))


group<-clim_lat_seo_E$group
pc1<-clim_lat_seo_E$pcac1
pc2<-clim_lat_seo_E$pcac2
pc3<-clim_lat_seo_E$pcac3

data5<-cbind(pc1, pc2, pc3)

pca3d(data5, group = group, palette = c("orange", "green"), axes.color = "black", show.plane = FALSE, bg = "white", radius = 1.2, show.ellipses = TRUE, ellipse.ci = 0.95, shape=1, show.shapes=TRUE)

snapshotPCA3d(file="3D_seoanei_E-latastei.png")



#gaditana - seoanei W - NO ADMIX
p1 <- ggplot(VLAg_VSEW_NO_int, aes(x = pcac1, y = pcac2, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
gad_seoW_NINT_pca1_2 <- p1 + geom_point(data = VLAgad_int, aes(x = pcac1, y = pcac2, color = group), size=3) + geom_point(data = VSEW_int, aes(x = pcac1, y = pcac2, color = group), size=3) + theme(legend.position = "none")

p2 <- ggplot(VLAg_VSEW_NO_int, aes(x = pcac1, y = pcac3, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
gad_seoW_NINT_pca1_3 <- p2 + geom_point(data = VLAgad_int, aes(x = pcac1, y = pcac2, color = group), size=3) + geom_point(data = VSEW_int, aes(x = pcac1, y = pcac2, color = group), size=3) + theme(legend.position = "none")

p3 <- ggplot(VLAg_VSEW_NO_int, aes(x = pcac2, y = pcac3, color = group)) + geom_point(size=2) + stat_ellipse(level = 0.95) + theme_bw() + theme(legend.position = "none")
gad_seoW_NINT_pca2_3 <- p3 + geom_point(data = VLAgad_int, aes(x = pcac1, y = pcac2, color = group), size=3) + geom_point(data = VSEW_int, aes(x = pcac1, y = pcac2, color = group), size=3) + theme(legend.position = "none")


legend <- cowplot::get_legend(ggplot(VLAg_VSEW_NO_int, aes(x = pcac1, y = pcac2, color = group)) + geom_point(size=2) + geom_point(data = VLAgad_int, aes(x = pcac1, y = pcac2, color = group), size=3) + geom_point(data = VSEW_int, aes(x = pcac1, y = pcac2, color = group), size=3))

figure <- ggarrange(gad_seoW_NINT_pca1_2, gad_seoW_NINT_pca1_3, gad_seoW_NINT_pca2_3, legend , align = "v", ncol = 2, nrow = 2)
annotate_figure(figure,
                top = text_grob("gaditana - seoanei W, admixt", color = "black", face = "bold", size = 14))


group<-VLAg_VSEW_NO_int$group
pc1<-VLAg_VSEW_NO_int$pcac1
pc2<-VLAg_VSEW_NO_int$pcac2
pc3<-VLAg_VSEW_NO_int$pcac3


data6<-cbind(pc1, pc2, pc3)

pca3d(data6, group = group, palette = c("red", "darkgreen"), axes.color = "black", show.plane = FALSE, bg = "white", radius = 1.2, show.ellipses = TRUE, ellipse.ci = 0.95, shape=1, show.shapes=TRUE)

snapshotPCA3d(file="3D_seoanei_W-gaditana_no-admix.png")




#NICHE OVERLAP TESTS - HYPERVOLUMES
library(hypervolume)
library(alphahull)

### comparisons: 
#1 - clim_VLA - clim_VSE
#2 - clim_VLAnoad - clim_VSEnoad
#3 - clim_lat - clim_gad
#4 - clim_seo_W - clim_seo_E
#5 - clim_lat - clim_seo_E
#6 - clim_gad - clim_seo_W
#7 - VLAgad_NO_int - VSEW_NO_int

#for each comparison change the followng parts:

## species 1 -PCs
pc1 <- clim_VLAnoad$pcac1
pc2 <- clim_VLAnoad$pcac2
pc3 <- clim_VLAnoad$pcac3

sp1 <- cbind(pc1,pc2,pc3)


## species 2 -PCs
pc1<-clim_VSEnoad$pcac1
pc2<-clim_VSEnoad$pcac2
pc3<-clim_VSEnoad$pcac3

sp2 <- cbind(pc1,pc2,pc3)


sorenson.coef <- function(X, Y, I) {
  ## Calculates the Soerenson coefficient
  ## X - number of species in sample X
  ## Y - number of species in sample Y
  ## I - number of species common to both samples
  ## NOTE: the formula is equivalent to (2a)/(2a+b+c) but avoids to
  ## include the intersection at the denominator and it is simpler
  ## to calculate the maximum Sorenson.
  (2*I)/(X+Y)
}


test.hv <- function(sp1, sp2, reps = 99, quantile = 0.95,
                    repsperpoint = 1000, verbose = TRUE, ...)  {
  ### Test hv siginificance with Soerenson and overlap index
  
  
  ## Build the observed hypervolumes for each species
  hv.sp1 <- hypervolume(sp1, quantile.requested=quantile,
                        samples.per.point = repsperpoint,
                        kde.bandwidth = estimate_bandwidth(sp1),
                        verbose = verbose > 1, ...)
  
  hv.sp2<- hypervolume(sp2, quantile.requested=quantile,
                       samples.per.point = repsperpoint,
                       kde.bandwidth = estimate_bandwidth(sp2),
                       verbose = verbose > 1, ...)
  
  hv.set <- hypervolume_set(hv.sp1, hv.sp2, check.memory=FALSE,
                            verbose = verbose > 1)
  
  volumes <- get_volume(hv.set)
  
  ## merge both species for a common niche/null distribution
  sp <- rbind(sp1, sp2)
  
  sp1.n <- nrow(sp1)
  sp2.n <- nrow(sp2)
  sp.n <- nrow(sp)
  
  null.dist <- matrix(NA, reps+1,3)
  colnames(null.dist) <- c("Soerensen", "MaxSorensen", "OverlapIndex")
  for (rep in 1:reps) {
    if (verbose) cat(sprintf("%.1f", rep/reps*100), "%\r")
    
    ## Sample the common niche 
    samples <- sample(1:sp.n, sp1.n)
    sp1.rep <- sp[samples,]
    sp2.rep <- sp[-samples,]
    
    sp1.rep.bw <- estimate_bandwidth(sp1.rep)
    sp2.rep.bw <- estimate_bandwidth(sp2.rep) 
    
    sp1.hv.rep <- hypervolume(sp1.rep, quantile.requested=quantile,
                              samples.per.point = repsperpoint,
                              kde.bandwidth = sp1.rep.bw,
                              verbose = verbose > 1, ...)
    
    sp2.hv.rep <- hypervolume(sp2.rep, quantile.requested=quantile,
                              samples.per.point = repsperpoint,
                              kde.bandwidth = sp2.rep.bw,
                              verbose = verbose > 1, ...)
    
    
    msg <- capture.output(
      hv.rep.set <- hypervolume_set(sp1.hv.rep, sp2.hv.rep,
                                    verbose = FALSE,
                                    check.memory=FALSE))
    
    v.rep <- get_volume(hv.rep.set)
    
    null.dist[rep, 1] <- sorenson.coef(v.rep[1], v.rep[2], v.rep[3])
    null.dist[rep, 2] <- sorenson.coef(v.rep[1], v.rep[2],
                                       min(v.rep[1:2]))
  }
  
  
  S.obs <- sorenson.coef(volumes[1], volumes[2], volumes[3])
  null.dist[reps+1, 1] <- S.obs
  
  Smax.obs <- sorenson.coef(volumes[1], volumes[2],
                            min(volumes[1:2]))
  null.dist[reps+1, 2] <- Smax.obs
  
  null.dist[,3] <- null.dist[,1] / null.dist[,2]
  
  OI.obs <- null.dist[reps+1, 3]
  
  p.value.sor <- sum(null.dist[,1] <= S.obs) / (reps+1)
  p.value.oi <- sum(null.dist[,3] <= OI.obs) / (reps+1)
  
  results <- data.frame(sorensen = c(S.obs, p.value.sor),
                        maxsorensen = c(Smax.obs, NA),
                        overlapIndex = c(OI.obs, p.value.oi))
  rownames(results) <- c("Observed", "p.value")
  
  list(results=results, null.dist = null.dist)
}



### Plots hypervolume wihtout test
hv.sp1 <- hypervolume(sp1, name="sp1", quantile.requested=0.95,
                      samples.per.point = 100)

hv.sp2<- hypervolume(sp2, name="sp2", quantile.requested=0.95,
                     samples.per.point = 100)

hv.set <- hypervolume_set(hv.sp1, hv.sp2, check.memory=FALSE)
volumes <- get_volume(hv.set)
x11()
plot(hv.sp1)
x11()
plot(hv.sp2)
x11()
plot(hv.set)

volumes



#Soerenson test with significance

spstat <- test.hv(sp1, sp2, reps=99, quantile=0.95, repsperpoint = 100) #It might take some time to run 99 reps

### Results
spstat$results


### plot Sorensen and null distribution
x11()
hist(spstat$null.dist[,1], col='gray', breaks=50,
     xlab="Sorensen index", main = "")
arrows(spstat$results[1,1], 100, spstat$results[1,1], 0, col='red')
### plot Overlap Index and null distribution
x11()
hist(spstat$null.dist[,3], col='gray', breaks=50,
     xlab="Overlap index", main="")
arrows(spstat$results[1,3], 100, spstat$results[1,3], 0, col='red')

save.image("hypervol.RData")
