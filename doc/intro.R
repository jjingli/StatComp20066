## ----eval=FALSE---------------------------------------------------------------
#  Ra <- function(x, sigma) {
#    if (any(x < 0)) return (0)
#    stopifnot(sigma > 0)
#    return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  MetropolisRa<-function(m,sigma){
#    x <- numeric(m)
#    x[1] <- rchisq(1, df=1)
#    k <- 0
#    u <- runif(m)
#    for (i in 2:m) {
#      xt <- x[i-1]
#      y <- rchisq(1, df = xt)
#      num <- Ra(y, sigma) * dchisq(xt, df = y)
#      den <- Ra(xt, sigma) * dchisq(y, df = xt)
#      if (u[i] <= num/den) x[i] <- y else {
#        x[i] <- xt
#        k <- k+1
#      }
#    }
#    return(list(x=x,k=k))
#  }

## ----eval=TRUE----------------------------------------------------------------
indep_norm <- function(m,a,b,p,n,mu,sigma){
  xt <- numeric(m)
  # generate the observed sample
  i <- sample(1:2, size=n, replace=TRUE, prob=c(p, 1-p))
  x <- rnorm(n, mu[i], sigma[i])
  # generate the independence sampler chain
  u <- runif(m)
  y <- rbeta(m, a, b) #proposal distribution
  xt[1] <- .5
  for (i in 2:m) {
    fy <- y[i] * dnorm(x, mu[1], sigma[1]) +
      (1-y[i]) * dnorm(x, mu[2], sigma[2])
    fx <- xt[i-1] * dnorm(x, mu[1], sigma[1]) +
      (1-xt[i-1]) * dnorm(x, mu[2], sigma[2])
    r <- prod(fy / fx) *
      (xt[i-1]^(a-1) * (1-xt[i-1])^(b-1)) /
      (y[i]^(a-1) * (1-y[i])^(b-1))
    if (u[i] <= r) xt[i] <- y[i] else
      xt[i] <- xt[i-1]
  }
  return(list(xt=xt,p_hat=mean(xt)))
}


