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



# Individual PCAs ####
pca.Q <- gm.prcomp(comb.lm.Q$coords)
pca.L <- gm.prcomp(comb.lm.L$coords)
pca.C <- gm.prcomp(comb.lm.C$coords)
pca.all <- gm.prcomp(comb.lm.all$coords)

# Species PCAs ####
all.comb <- list(all = comb.lm.all, Q = comb.lm.Q, L = comb.lm.L, C = comb.lm.C)
sp.pca <- lapply(all.comb, function(x) {
        gr.data <- coords.subset(x$coords, spec.sel$sp_island)
        gr.means <- lapply(gr.data, mshape)
        gr.means <- simplify2array(gr.means)
        pca.gr <- gm.prcomp(gr.means)
}
)

# PCA.all vs body parts ####
# Ind PCAs
ed.Q <- as.matrix(dist(pca.Q$x))
ed.L <- as.matrix(dist(pca.L$x))
ed.C <- as.matrix(dist(pca.C$x))
ed.ALL <- as.matrix(dist(pca.all$x))

library(vegan)
mantel(ed.Q, ed.ALL)
# Mantel statistic r: 0.2572 
# Significance: 0.001 

mantel(ed.C, ed.ALL)
# Mantel statistic r: 0.3943 
# Significance: 0.001 

mantel(ed.L, ed.ALL)
# Mantel statistic r: 0.4469 
# Significance: 0.001 

# Species PCAs
ed.sp.Q <- as.matrix(dist(sp.pca$Q$x))
ed.sp.L <- as.matrix(dist(sp.pca$L$x))
ed.sp.C <- as.matrix(dist(sp.pca$C$x))
ed.sp.ALL <- as.matrix(dist(sp.pca$all$x))

mantel(ed.sp.Q, ed.sp.ALL)
# Mantel statistic r: 0.4964 
# Significance: 0.004 

mantel(ed.sp.C, ed.sp.ALL)
# Mantel statistic r: 0.4765 
# Significance: 0.028

mantel(ed.sp.L, ed.sp.ALL)
# Mantel statistic r: 0.6766 
# Significance: 0.012

############################################################################
# PC direction comparisons ####
# Compare the direction of PC1 between the individual and species levels
source("PCcompare.R")

pc.comp.Q <- pc.compare(x = two.d.array(comb.lm.Q$coords),
                           gr = spec.sel$sp_island, pc = 1)
pc.comp.Q[c(1, 3)]
# theta.obs = 15.86286
# p.val = 0.001

pc.comp.C <- pc.compare(x = two.d.array(comb.lm.C$coords),
                        gr = spec.sel$sp_island, pc = 1)
pc.comp.C[c(1, 3)]
# theta.obs = 53.21916
# p.val = 0.029

pc.comp.L <- pc.compare(x = two.d.array(comb.lm.L$coords),
                        gr = spec.sel$sp_island, pc = 1)
pc.comp.L[c(1, 3)]
# theta.obs = 8.84914
# p.val = 1

pc.comp.ALL <- pc.compare(x = two.d.array(comb.lm.all$coords),
                        gr = spec.sel$sp_island, pc = 1)
pc.comp.ALL[c(1, 3)]
# theta.obs = 10.13414
# p.val = 0.001
