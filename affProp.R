############################################################
# Affinity Propagation
# Clustering by Passing Messages Between Data Points
# B.J. Frey and D. Dueck
# Science 315, 972
# 2007
############################################################
affProp <- function(S, pref=NA, lambda = 0.5, iter.max = 100)
{
  # init
  dimnames(S) <- NULL
  n <- nrow(S)
  A <- matrix(0,nr=n,nc=n)
  R <- matrix(0,nr=n,nc=n)
  # set preferences
  if(is.na(pref[1]))
  {
    diag(S) <- median(c(S[lower.tri(S)],S[upper.tri(S)]))
  }
  else
  {
    diag(S) <- pref
  }
  # remove degeneracies in S
  S <- S + 1e-12 * matrix(rnorm(n),nc=n,nr=n) * (max(S) - min(S))
  # compute
  for(i in 1:iter.max)
  {
    # responsibilities
    R.old <- R
    AS <- A + S
    maxAS <- apply(AS,1,max)
    maxAS.idx <- apply(AS,1,which.max)
    for(i in 1:n)
    {
      AS[i,maxAS.idx[i]] <- -.Machine$double.xmax
    }
    maxAS2 <- apply(AS,1,max)
    maxAS2.idx <- apply(AS,1,which.max)
    R <- S - matrix(maxAS,nr=n,nc=n)
    for(i in 1:n)
    {
      R[i,maxAS.idx[i]] <- S[i,maxAS.idx[i]] - maxAS2[i]
    }
    # dampen responsibilities
    R <- (1 - lambda) * R + lambda * R.old
    # availabilities
    A.old <- A
    Rp <- apply(R,c(1,2),function(u) max(u,0))
    for(k in 1:n)
    {
      Rp[k,k] <- R[k,k]
    }
    A <- matrix(colSums(Rp),nc=n,nr=n,byrow=T) - Rp
    dA <- diag(A)
    A <- apply(A,c(1,2),function(u) min(u,0))
    for(k in 1:n)
    {
      A[k,k] <- dA[k]
    }
    # dampen availabilities
    A <- (1 - lambda) * A + lambda * A.old
  }
  # pseudomarginals
  E <- R + A
  # indices of exemplars
  I <- which(diag(E)>0)
  # number of exemplars
  K <- length(I)
  # assignments
  tmp <- apply(S[,I],1,max)
  c <- apply(S[,I],1,which.max)
  c[I] <- 1:K
  idx <- I[c]
  ce <- unique(idx)
  return(list(k=K,cluster=idx,centers=ce,comb.evidence=E))
}

test.affProp <- function()
{
  # kmeans
  x <- rbind(matrix(rnorm(20, sd = 0.05), ncol = 2), matrix(rnorm(100, mean = 1, sd = 0.2), ncol = 2))
  colnames(x) <- c("x", "y")
  k=2
  cl <- kmeans(x, k)
  # affprop
  S <- as.matrix(dist(x)^2)
  S <- -S
  # pref <- median(c(S[lower.tri(S)],S[upper.tri(S)]))
  pref <- -2
  ap <- affProp(S,pref)
  # plot
  par(mfrow=c(1,2))
  plot(x, col = cl$cluster,main="result from kmeans")
  points(cl$centers, col = 1:k, pch = 8, cex=2)
  plot(x, col = ap$cluster,main="result from affprop")
  points(x[ap$centers,],col=ap$centers,pch=8,cex=2)
  # return
  return(list(data=x,sim=S,result.cl=cl,result.ap=ap))
}

d2test.affProp <- function(x)
{
  S <- -as.matrix(dist(x)^2)
  K <- c()
  I <- c()
  for(i in seq(-10,-0.1,0.1))
  {
    pref <- i
    k <- affProp(S,pref)$k
    if(!is.element(k,K))
    {
      K <- c(K,k)
      I <- c(I,i)
    }
  }
  return(list(K=K,I=I))
}

# parallel <- rbind(matrix(c(runif(100,1,5),runif(100,4,5)),100),matrix(c(runif(100,1,5),runif(100,1,2)),100))
# boxes <- rbind(matrix(c(runif(100,2,5),runif(100,2,5)),nr=100),matrix(c(runif(20,7,8),runif(20,3,4)),nr=20),matrix(c(runif(10,5,7),runif(10,3.4,3.6)),nr=10))
