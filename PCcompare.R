## AK, 4 November 2015

### Functions for calculating PC axis for individuals, and group means,
# calculating their angle, and then running a randomization test to compare them.

pc.angle <- function(x, gr, pc){
  group <- list(gr)
  # Observed means, pcs, and angles
  ind.pc <- prcomp(x)$rotation[,pc] # Indiv PC1
  
  lsm <- as.matrix(aggregate(x, mean, by=group)[,-1])
  gr.pc <- prcomp(lsm)$rotation[,pc]
  
  theta <- (180/pi)*acos(t(ind.pc)%*%gr.pc)
  if(theta>90) theta <- 180-theta
  
  return(theta)
}

# Randomization function
pc.compare <- function(x, gr, pc, nperm=999, show.plot=F){
  x <- as.matrix(x)
  all.theta <- vector("numeric", nperm+1)
  all.theta[1] <- pc.angle(x, gr, pc)
  for (p in 1:nperm){
    yhat <- predict(lm(x~gr))
    rand.res <- sample(resid(lm(x~gr)))
    rand.x <- yhat + rand.res
    all.theta[p+1] <- pc.angle(rand.x, gr, pc)
  }
  res <- list()
  res[["theta.obs"]] <- all.theta[1]
  res[["rand.thetas"]] <- all.theta[-1]
  res[["p.val"]] <- ((nperm+1) - (which(sort(all.theta)==all.theta[1])-1))/(nperm+1)
  if(show.plot==T) {
    hist(res[["rand.thetas"]], xlim=range(all.theta))
    abline(v=res[["theta.obs"]], col="red", lwd=2)
  }
  return(res)
}
