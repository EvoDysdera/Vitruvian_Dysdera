rm(list = ls())
dev.off()

library(geomorph)
library(plyr)

load("data_p4.bin")

spec.sel <- all.data$specs
spec.sel$sp_island <- revalue(spec.sel$sp_island,
                              c("maa_t"="mac_t","mae_t"="mac_t","maw_t"="mac_t"))
spec.sel$colours <- as.character(revalue(spec.sel$sp_island,
                                         c("amb_t"="gray32", "bre_t"="brown", "cur_t"="aquamarine", 
                                           "ins_t"="darkgoldenrod1", "mac_t"="brown1", "ram_g"="chartreuse", 
                                           "rat_p"="cornflowerblue", 
                                           "sil_g"="darkorange2", "sil_h"= "darkorange3",  
                                           "sil_p"="darkorange1", "ung_t"="blue4", "ver_t"="darkorchid")))

# Pairwise integration of structures - ORI ####
dt.ori <- lapply(all.data$gpa.resampled, function(x){
  dt <- x$coords
  dimnames(dt)[[3]] <- spec.sel$code
  return(dt)
})

pls.ORI <- list()
for(i in 1:(length(dt.ori)-1)){
  for(j in (i+1):length(dt.ori)){
    temp <- integration.test(dt.ori[[i]], dt.ori[[j]],
                             iter = 999, print.progress = F)
    comp.name <- paste(names(dt.ori)[i], names(dt.ori)[j], sep = "_")
    pls.ORI[[comp.name]] <- temp
  }
}

ori.rPLS <- unlist(lapply(1:length(pls.ORI), function(x) {pls.ORI[[x]][1]}))
ori.pvalues <- unlist(lapply(1:length(pls.ORI), function(x) {pls.ORI[[x]][3]}))
names(ori.rPLS) <- names(ori.pvalues) <- names(pls.ORI)
ori.rPLS
ori.pvalues

# Compare levels of integration (pairwise - total chaos)
pls.comp.ORI <- compare.pls(pls.ORI)

# Review zvalues across structures 
zval <- pls.comp.ORI$sample.z[c(1, 2, 9, 22, 31:36, 3:8, 10:21, 23:30)]
plot(hclust(dist(zval)))

# ANOVA on integration zvalues. 
comps <- as.factor(c(rep("within", 10), rep("between", 26)))
comps.lm <- lm.rrpp(zval ~ comps)
anova(comps.lm)
boxplot(zval ~ comps)
barplot(zval, las = 2, col = c(rep("grey", 10), rep("white", 26)))
