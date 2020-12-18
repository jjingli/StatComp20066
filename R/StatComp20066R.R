#' @title The Rayleigh density function
#' @description The Rayleigh density function
#' @param x the variable
#' @param sigma the parameter of Rayleigh density function
#' @return The Rayleigh density function
#' @examples
#' \dontrun{
#' ra <- Ra(x,4)
#' }
#' @export
Ra <- function(x, sigma) {
  if (any(x < 0)) return (0)
  stopifnot(sigma > 0)
  return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}

#' @title A Metropolis-Hastings sampler using R
#' @description A Metropolis-Hastings sampler using R
#' @importFrom stats rchisq runif dchisq
#' @param m the number of samples
#' @param sigma the parameter of Rayleigh distribution
#' @return a random sample of size and the number of reject 
#' @examples
#' \dontrun{
#' MH.Ra <- MetropolisRa(10000,4)
#' }
#' @export
MetropolisRa<-function(m,sigma){
  x <- numeric(m)
  x[1] <- rchisq(1, df=1)
  k <- 0
  u <- runif(m)
  for (i in 2:m) {
    xt <- x[i-1]
    y <- rchisq(1, df = xt)
    num <- Ra(y, sigma) * dchisq(xt, df = y)
    den <- Ra(xt, sigma) * dchisq(y, df = xt)
    if (u[i] <= num/den) x[i] <- y else {
      x[i] <- xt
      k <- k+1 
    }
  }
  return(list(x=x,k=k))
}

#' @title A Independence sampler using R
#' @description A Independence sampler using R
#' @importFrom stats rnorm runif rbeta dnorm
#' @param m the length of chain
#' @param a parameter of Beta(a,b) proposal dist.
#' @param b parameter of Beta(a,b) proposal dist.
#' @param p mixing parameter
#' @param n sample size
#' @param mu parameters of the normal densities include two numbers
#' @param sigma parameters of the normal densities include two numbers
#' @return the estimation of mix parameter p
#' @examples
#' \dontrun{
#' I<-indep_norm(10000,1,1,0.2,50,c(0, 5),c(1, 1))
#' print(I$xt)
#' print(I$p_hat)
#' plot(xt, type="l", ylab="p")
#' hist(xt[101:m], main="", xlab="p", prob=TRUE)
#' }
#' @export
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
