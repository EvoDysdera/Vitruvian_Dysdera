rm(list = ls())
dev.off()

library(geomorph)
library(factoextra)
library(plyr)

load("data_p4.bin")

spec.sel <- all.data$specs
spec.sel$sp_island <- revalue(spec.sel$sp_island,
                              c("maa_t"="mac_t","mae_t"="mac_t","maw_t"="mac_t"))

comb.lm.Q <- combine.subsets(Q1 = all.data$gpa.resampled$Q1, 
                           Q2 = all.data$gpa.resampled$Q2,
                           Q3 = all.data$gpa.resampled$Q3,
                           #C1 = all.data$gpa.resampled$C1,
                           #C3 = all.data$gpa.resampled$C3,
                           #L11 = all.data$gpa.resampled$L11,
                           #L41 = all.data$gpa.resampled$L41,
                           #L14 = all.data$gpa.resampled$L14,
                           #L44 = all.data$gpa.resampled$L44, 
                           gpa = TRUE)

comb.lm.L <- combine.subsets(#Q1 = all.data$gpa.resampled$Q1, 
                             #Q2 = all.data$gpa.resampled$Q2,
                             #Q3 = all.data$gpa.resampled$Q3,
                             #C1 = all.data$gpa.resampled$C1,
                             #C3 = all.data$gpa.resampled$C3,
                             L11 = all.data$gpa.resampled$L11,
                             L41 = all.data$gpa.resampled$L41,
                             L14 = all.data$gpa.resampled$L14,
                             L44 = all.data$gpa.resampled$L44, 
                             gpa = TRUE)

comb.lm.C <- combine.subsets(#Q1 = all.data$gpa.resampled$Q1, 
                               #Q2 = all.data$gpa.resampled$Q2,
                               #Q3 = all.data$gpa.resampled$Q3,
                               C1 = all.data$gpa.resampled$C1,
                               C3 = all.data$gpa.resampled$C3,
                               #L11 = all.data$gpa.resampled$L11,
                               #L41 = all.data$gpa.resampled$L41,
                               #L14 = all.data$gpa.resampled$L14,
                               #L44 = all.data$gpa.resampled$L44, 
                               gpa = TRUE)

comb.lm.all <- combine.subsets(Q1 = all.data$gpa.resampled$Q1, 
                             Q2 = all.data$gpa.resampled$Q2,
                             Q3 = all.data$gpa.resampled$Q3,
                             C1 = all.data$gpa.resampled$C1,
                             C3 = all.data$gpa.resampled$C3,
                             L11 = all.data$gpa.resampled$L11,
                             L41 = all.data$gpa.resampled$L41,
                             L14 = all.data$gpa.resampled$L14,
                             L44 = all.data$gpa.resampled$L44, 
                             gpa = TRUE)

### PROBLEMATIC INDIVS FOR REVIEW (C1 y C3)
#ind.outliers <- plotOutliers(comb.lm.C$coords, inspect.outliers = F)
#spec.sel[ind.outliers[1:6],]


# Individual PCAs ###
pca.Q <- gm.prcomp(comb.lm.Q$coords)
pca.L <- gm.prcomp(comb.lm.L$coords)
pca.C <- gm.prcomp(comb.lm.C$coords)
pca.all <- gm.prcomp(comb.lm.all$coords)


spec.sel$colours <- as.character(revalue(spec.sel$sp_island,
        c("amb_t"="gray48", "bre_t"="darkkhaki", "cur_t"="aquamarine3", 
          "ins_t"="gold", "mac_t"="darkolivegreen3", "ram_g"="mediumorchid3", 
          "rat_p"="white", 
          "sil_g"="darkorange2", "sil_h"= "darkorange2",  
          "sil_p"="darkorange2", "ung_t"="black", "ver_t"="darkred")))
gp.col <- as.character(revalue(unique(spec.sel$sp_island),
                               c("amb_t"="gray48", "bre_t"="darkkhaki", "cur_t"="aquamarine3", 
                                 "ins_t"="gold", "mac_t"="darkolivegreen3", "ram_g"="mediumorchid3", 
                                 "rat_p"="white", 
                                 "sil_g"="darkorange2", "sil_h"= "darkorange2",  
                                 "sil_p"="darkorange2", "ung_t"="black", "ver_t"="darkred")))
gp.col <- data.frame(gp.col)       
row.names(gp.col) <- unique(spec.sel$sp_island)
gp.col$sp <- c("D. ratonensis", "D. unguimmanis", "D. ramblae", "D. silvatica H",
               "D. insulana", "D. brevisetae", "D. macra", "D. verneui",
               "D. curviseta", "D. ambulotenta", "D. silvatica G", "D. silvatica P")
gp.col <- gp.col[c(10, 1, 2, 11, 12, 4, 8, 5, 3, 6, 7, 9),] 

pdf("PCA_indivs.pdf", family = "Times")        
layout(matrix(c(1, 2, 5, 3, 4, 5), nrow = 2, byrow = T), widths = c(3, 3, 1))
par(mar = c(3, 3, 0.75, 0.5), mgp = c(1.75, 0.5, 0))
plot(pca.all, 
     cex = 1.25, cex.lab = 1.25, font.lab = 2,
     pch = ifelse(spec.sel$sp_island=="sil_g", 22, 
                  ifelse(spec.sel$sp_island=="sil_h", 23,
                         ifelse(spec.sel$sp_island=="sil_p", 24, 21))), 
     bg = spec.sel$colours)
legend("topleft", legend = "A: All subsets", text.font = 2, cex = 1.25, bty = "n", 
       adj = c(0.25, 0))

plot(pca.Q, 
     cex = 1.25, cex.lab = 1.25, font.lab = 2,
     pch = ifelse(spec.sel$sp_island=="sil_g", 22, 
                  ifelse(spec.sel$sp_island=="sil_h", 23,
                         ifelse(spec.sel$sp_island=="sil_p", 24, 21))), 
     bg = spec.sel$colours)
legend("topright", legend = "B: Chelicera", text.font = 2, cex = 1.25, bty = "n", 
       adj = c(0, 0))

plot(pca.C, 
     cex = 1.25, cex.lab = 1.25, font.lab = 2,
     pch = ifelse(spec.sel$sp_island=="sil_g", 22, 
                  ifelse(spec.sel$sp_island=="sil_h", 23,
                         ifelse(spec.sel$sp_island=="sil_p", 24, 21))), 
     bg = spec.sel$colours)
legend("topleft", legend = "C: Carapace", text.font = 2, cex = 1.25, bty = "n", 
       adj = c(0.25, 0))

plot(pca.L, 
     cex = 1.25, cex.lab = 1.25, font.lab = 2,
     pch = ifelse(spec.sel$sp_island=="sil_g", 22, 
                  ifelse(spec.sel$sp_island=="sil_h", 23,
                         ifelse(spec.sel$sp_island=="sil_p", 24, 21))), 
     bg = spec.sel$colours)
legend("topright", legend = "D: Legs", text.font = 2, cex = 1.25, bty = "n", 
       adj = c(0, 0))



par(xpd = T, mar = rep(0, 4))
plot.new()
legend("left", pch = c(rep(21, 3), 22, 24, 23, rep(21, 6)), 
       pt.bg = gp.col$gp.col, bty = "n", legend = gp.col$sp,
       text.font = 4, pt.cex = 1.5)
dev.off()

###############################################################
# IndivPCAs - grids ####
load("linksQ.bin")
load("linksC.bin")
load("linksL.bin")
load("linksALL.bin")

subset.pts_Q <- c(0, 14, 38, 51)
names(subset.pts_Q)<-c("", "Q1", "Q2", "Q3")
sh.pc1.Q <- list()
for(i in 2:length(subset.pts_Q)){
        sh.pc1.Q[[names(subset.pts_Q)[i]]] <- lapply(pca.Q$shapes$shapes.comp1, function(x){
                x[(subset.pts_Q[i-1]+1):subset.pts_Q[i],]
                })
}

subset.pts_C <- c(0, 32, 54)
names(subset.pts_C)<-c("", "C1", "C3")
sh.pc1.C <- list()
for(i in 2:length(subset.pts_C)){
        sh.pc1.C[[names(subset.pts_C)[i]]] <- lapply(pca.C$shapes$shapes.comp1, function(x){
                x[(subset.pts_C[i-1]+1):subset.pts_C[i],]
                })
}

subset.pts_L <- c(0, 14, 28, 44, 60)
names(subset.pts_L)<-c("", "L11", "L41", "L14", "L44")
sh.pc1.L <- list()
for(i in 2:length(subset.pts_L)){
        sh.pc1.L[[names(subset.pts_L)[i]]] <- lapply(pca.L$shapes$shapes.comp1, function(x){
                x[(subset.pts_L[i-1]+1):subset.pts_L[i],]
                })
}

subset.pts_all <- c(0, 14, 38, 51, 83, 105, 119, 133, 149, 165)
names(subset.pts_all)<-c("", "Q1", "Q2", "Q3", "C1", "C3", "L11", "L41", "L14", "L44")
sh.pc1.all <- list()
for(i in 2:length(subset.pts_all)){
        sh.pc1.all[[names(subset.pts_all)[i]]] <- lapply(pca.all$shapes$shapes.comp1, function(x){
                x[(subset.pts_all[i-1]+1):subset.pts_all[i],]
                })
}

# Plot with all analyses ####
pdf("PCshapes_indiv.pdf", family = "Times")
layout(matrix(1:40, nrow = 10, byrow = F), heights = c(2, rep(5, 9)))
par(mar = rep(0,4))
plot.new()
text(0.5, 0.5, "All subsets", font = 2, cex = 1.5)
lapply(1:length(sh.pc1.all), function(x){
        plotRefToTarget(sh.pc1.all[[x]]$min, sh.pc1.all[[x]]$max, links = linksALL[[x]],
                        gridPars = gridPar(grid.col = "grey30"))
})

plot.new()
text(0.5, 0.5, "Chelicera", font = 2, cex = 1.5)
lapply(1:length(sh.pc1.Q), function(x){
        plotRefToTarget(sh.pc1.Q[[x]]$max, sh.pc1.Q[[x]]$min, links = linksq[[x]],
                        gridPars = gridPar(grid.col = "grey30"))
})
for(i in 1:6) {plot.new()}

plot.new()
text(0.5, 0.5, "Carapace", font = 2, cex = 1.5)
for(i in 1:3) {plot.new()}

lapply(1:length(sh.pc1.C), function(x){
        plotRefToTarget(sh.pc1.C[[x]]$min, sh.pc1.C[[x]]$max, links = linksc[[x]],
                        gridPars = gridPar(grid.col = "grey30"))
})
for(i in 1:4) {plot.new()}

plot.new()
text(0.5, 0.5, "Legs", font = 2, cex = 1.5)
for(i in 1:5) {plot.new()}

lapply(1:length(sh.pc1.L), function(x){
        plotRefToTarget(sh.pc1.L[[x]]$max, sh.pc1.L[[x]]$min, links = linksl[[x]],
                        gridPars = gridPar(grid.col = "grey30"))
})
dev.off()


###############################################################
# Species PCAs ####
all.comb <- list(all = comb.lm.all, Q = comb.lm.Q, L = comb.lm.L, C = comb.lm.C)
sp.pca <- lapply(all.comb, function(x) {
        gr.data <- coords.subset(x$coords, spec.sel$sp_island)
        gr.means <- lapply(gr.data, mshape)
        gr.means <- simplify2array(gr.means)
        pca.gr <- gm.prcomp(gr.means)
}
)

pdf("PCA_sps.pdf", family = "Times")        
layout(matrix(c(1, 2, 5, 3, 4, 5), nrow = 2, byrow = T), widths = c(3, 3, 1))
par(mar = c(3, 3, 0.75, 0.5), mgp = c(1.75, 0.5, 0))
plot(sp.pca$all, 
     cex = 2, cex.lab = 1.25, font.lab = 2,
     pch = ifelse(rownames(sp.pca$all$x)=="sil_g", 22, 
                  ifelse(rownames(sp.pca$all$x)=="sil_h", 23,
                         ifelse(rownames(sp.pca$all$x)=="sil_p", 24, 21))), 
     bg = gp.col[rownames(sp.pca$all$x),"gp.col"])
legend("topleft", legend = "A: All subsets", text.font = 2, cex = 1.25, bty = "n", 
       adj = c(0.25, 0))

plot(sp.pca$Q, 
     cex = 2, cex.lab = 1.25, font.lab = 2,
     pch = ifelse(rownames(sp.pca$all$x)=="sil_g", 22, 
                  ifelse(rownames(sp.pca$all$x)=="sil_h", 23,
                         ifelse(rownames(sp.pca$all$x)=="sil_p", 24, 21))), 
     bg = gp.col[rownames(sp.pca$all$x),"gp.col"])
legend("topright", legend = "B: Chelicera", text.font = 2, cex = 1.25, bty = "n", 
       adj = c(0, 0))
plot(sp.pca$C, 
     cex = 2, cex.lab = 1.25, font.lab = 2,
     pch = ifelse(rownames(sp.pca$all$x)=="sil_g", 22, 
                  ifelse(rownames(sp.pca$all$x)=="sil_h", 23,
                         ifelse(rownames(sp.pca$all$x)=="sil_p", 24, 21))), 
     bg = gp.col[rownames(sp.pca$all$x),"gp.col"])
legend("topleft", legend = "C: Carapace", text.font = 2, cex = 1.25, bty = "n", 
       adj = c(0.25, 0))

plot(sp.pca$L, 
     cex = 2, cex.lab = 1.25, font.lab = 2,
     pch = ifelse(rownames(sp.pca$all$x)=="sil_g", 22, 
                  ifelse(rownames(sp.pca$all$x)=="sil_h", 23,
                         ifelse(rownames(sp.pca$all$x)=="sil_p", 24, 21))), 
     bg = gp.col[rownames(sp.pca$all$x),"gp.col"])
legend("topright", legend = "D: Legs", text.font = 2, cex = 1.25, bty = "n", 
       adj = c(0, 0))



par(xpd = T, mar = rep(0, 4))
plot.new()
legend("left", pch = c(rep(21, 3), 22, 24, 23, rep(21, 6)), 
       pt.bg = gp.col$gp.col, bty = "n", legend = gp.col$sp,
       text.font = 4, pt.cex = 1.5)
dev.off()

###############################################################
# speciesPCAs - grids ####
sp.pc1.Q <- list()
for(i in 2:length(subset.pts_Q)){
        sp.pc1.Q[[names(subset.pts_Q)[i]]] <- lapply(sp.pca$Q$shapes$shapes.comp1, function(x){
                x[(subset.pts_Q[i-1]+1):subset.pts_Q[i],]
        })
}

sp.pc1.C <- list()
for(i in 2:length(subset.pts_C)){
        sp.pc1.C[[names(subset.pts_C)[i]]] <- lapply(sp.pca$C$shapes$shapes.comp1, function(x){
                x[(subset.pts_C[i-1]+1):subset.pts_C[i],]
        })
}

sp.pc1.L <- list()
for(i in 2:length(subset.pts_L)){
        sp.pc1.L[[names(subset.pts_L)[i]]] <- lapply(sp.pca$L$shapes$shapes.comp1, function(x){
                x[(subset.pts_L[i-1]+1):subset.pts_L[i],]
        })
}

sp.pc1.all <- list()
for(i in 2:length(subset.pts_all)){
        sp.pc1.all[[names(subset.pts_all)[i]]] <- lapply(sp.pca$all$shapes$shapes.comp1, function(x){
                x[(subset.pts_all[i-1]+1):subset.pts_all[i],]
        })
}

# Plot with all analyses ####
pdf("PCshapes_SP.pdf", family = "Times")
layout(matrix(1:40, nrow = 10, byrow = F), heights = c(2, rep(5, 9)))
par(mar = rep(0,4))
plot.new()
text(0.5, 0.5, "All subsets", font = 2, cex = 1.5)
lapply(1:length(sp.pc1.all), function(x){
        plotRefToTarget(sp.pc1.all[[x]]$max, sp.pc1.all[[x]]$min, links = linksALL[[x]],
                        gridPars = gridPar(grid.col = "grey30"))
})

plot.new()
text(0.5, 0.5, "Quelicera", font = 2, cex = 1.5)
lapply(1:length(sp.pc1.Q), function(x){
        plotRefToTarget(sp.pc1.Q[[x]]$min, sp.pc1.Q[[x]]$max, links = linksq[[x]],
                        gridPars = gridPar(grid.col = "grey30"))
})
for(i in 1:6) {plot.new()}

plot.new()
text(0.5, 0.5, "Carapace", font = 2, cex = 1.5)
for(i in 1:3) {plot.new()}
lapply(1:length(sp.pc1.C), function(x){
        plotRefToTarget(sp.pc1.C[[x]]$max, sp.pc1.C[[x]]$min, links = linksc[[x]],
                        gridPars = gridPar(grid.col = "grey30"))
})
for(i in 1:4) {plot.new()}

plot.new()
text(0.5, 0.5, "Limbs", font = 2, cex = 1.5)
for(i in 1:5) {plot.new()}
lapply(1:length(sp.pc1.L), function(x){
        plotRefToTarget(sp.pc1.L[[x]]$min, sp.pc1.L[[x]]$max, links = linksl[[x]],
                        gridPars = gridPar(grid.col = "grey30"))
})
dev.off()


