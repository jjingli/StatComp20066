## ----eval=FALSE---------------------------------------------------------------
#  install.packages("knitr")
#  install.packages("ggplot2")

## -----------------------------------------------------------------------------
library(ggplot2)
data(diamonds)

## ----fig.cap = "Scatter map with cut as color scale"--------------------------
qplot(carat,price,data = diamonds,colour = cut)

## ----fig.cap = "Scatter map with cut as color scale"--------------------------
qplot(carat,price,data = diamonds,colour = clarity)

## ----results='asis'-----------------------------------------------------------
library(knitr)
kable(diamonds[1:10,],caption = "A knitr kable",align = "c")

## -----------------------------------------------------------------------------
set.seed(16273)
n <- 100
u <- runif(n)   #产生100个服从均匀分布的随机数
x <- 2*exp(-((log(1-u)))/2) # F(x) = 1-(2/x)^2,x>=2,生成服从Pareto(2,2)的随机数x

## -----------------------------------------------------------------------------
hist(x, prob = TRUE,breaks = 10)  #绘制x的频率分布直方图
y <- seq(2, 20, .01)    #绘制Pareto(2,2)的概率密度函数
lines(y, 8/y^3)

## -----------------------------------------------------------------------------
set.seed(26667)
n <-1000
u1 <- runif(n,-1,1)
u2 <- runif(n,-1,1)
u3 <- runif(n,-1,1)
for (i in 1:n) {
    if ((abs(u3[i]) >= abs(u2[i]) & abs(u3[i]) >= abs(u1[i]))){
      x[i] <- u2[i]
    }else{
      x[i] <- u3[i]
    }
}                 #根据题中所给算法从目标分布中产生随机数

hist(x, prob = TRUE,breaks = 10)      #构建频率分布直方图


## -----------------------------------------------------------------------------
set.seed(14256)
n <- 1e5
x <- runif(n, min=0, max=pi/3)   #生成模拟的随机数
I <- mean(sin(x)) * pi/3   #根据蒙特卡洛方法计算积分值
print(c(I,1/2,I*2))   #打印估计值，理论值，以及他们的比值（越接近1模拟地越准确）

## -----------------------------------------------------------------------------
MC.Phi <- function(x, R, antithetic = TRUE) {
  u <- runif(R/2,0,x)
  if (!antithetic) v <- runif(R/2) else v <- 1 - u
  u <- c(u, v)   #生成随机数
  g <- x * exp(u )
  cdf<- mean(g) 
}      #计算e^x在[0,x]上的积分值
set.seed(1234)
m <- 1000
MC1 <- MC2 <- numeric(m)
x <- 1   #计算e^x在[0,1]上的积分值
for (i in 1:m) {
  MC1[i] <- MC.Phi(x, 1000, anti = FALSE)  #the simple Monte Carlo method
  MC2[i] <- MC.Phi(x, 1000)   #the antithetic variate approach
}
format(c(var(MC1),var(MC2),1-var(MC2)/var(MC1)),scientific = FALSE)   #输出模拟的简单蒙特卡洛方差，对偶方法的方差，以及使用对偶方法方差减少的百分比


## -----------------------------------------------------------------------------
x <- seq(1,1000,0.1)
y <- x^2/sqrt(2*pi)*exp(-x^2/2)
plot(x,y,type="l", lwd=1,col="DarkTurquoise",lty=1)
lines(x^(-2),col="DeepPink",lty=2)
lines(x*exp(-(x^2-1)/2),col="RosyBrown",lty=3)
legend(700,0.25,c("g(x)","f1","f2"),col=c("DarkTurquoise","DeepPink","RosyBrown"),text.col=c("DarkTurquoise","DeepPink","RosyBrown"),lty=c(1,2,3))


## -----------------------------------------------------------------------------
n <- 1000
m <- 100
theta_1 <- theta_2 <- numeric(m)
for (i in 1:100) {
  x0<- runif(n)
  x1 <- 1/(1-x0)  #用逆变换法求服从密度函数为f1的随机数x1
  x2 <- sqrt(1-2*log(1-x0))  #用逆变换法求服从密度函数为f2的随机数x2
  theta_1[i] <- mean((x1)^4/sqrt(2*pi)*exp((-(x1)^2/2)))
  theta_2[i] <- mean((x2)*exp(-1/2)/sqrt(2*pi))  #分别用f1,f2估计积分值theta
}
print(format(c(var(theta_1),var(theta_2),var(theta_1)/var(theta_2))))
mean(theta_1)
mean(theta_2)

## -----------------------------------------------------------------------------
M <- 10000
k <- 5 
r <- M/k #replicates per stratum
N <- 50 #number of times to repeat the estimation
T <- numeric(k)
est <- numeric(N)
g<-function(x)exp(-x)/(1+x^2)
f3<-function(x)exp(-x)/(1-exp(-1))
for (i in 1:N) {
  for(j in 1:k){
    x0 <- runif(M/k,0,1)
    x1 <- -log(1-x0*(1-exp(-1)))
    T[j]<-mean(g(x1)/f3(x1)*(x1>(j-1)/k)*(x1<j/k))
    }
  est[i] <- sum(T)   
}
print(c(mean(est),sd(est)))

## -----------------------------------------------------------------------------
m<-10000
n<-20
alpha<-0.05
count<-0
lcl<-ucl<-numeric(n)
for (i in 1:m) {
  for (j in 1:n) {
    x <- rlnorm(n, meanlog = 0, sdlog = 2)
    y <- log(x)
  }
  lcl[i] <- mean(y)-qt(1-alpha/2,n-1)*sd(y)/sqrt(n)
  ucl[i] <- mean(y)+qt(1-alpha/2,n-1)*sd(y)/sqrt(n)
  if (lcl[i]<0 & ucl[i]>0)
     {count<-count+1}
  else
     {count<-count}
}
conf_level<-count/m
print(c(conf_level))

## -----------------------------------------------------------------------------
m<-10000
n<-20
alpha<-0.05
count<-0
df<-2
lcl<-ucl<-numeric(n)
for (i in 1:m) {
  for (j in 1:n) {
    x <- rchisq(n, df)
  }
  lcl[i] <- mean(x)-qt(1-alpha/2,n-1)*sd(y)/sqrt(n)
  ucl[i] <- mean(x)+qt(1-alpha/2,n-1)*sd(y)/sqrt(n)
  if (lcl[i]<df & ucl[i]>df)
     {count<-count+1}
  else
     {count<-count}
}
conf_level<-count/m
print(conf_level)

## -----------------------------------------------------------------------------
sk <- function(x) {
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}
alpha <- .1
n <- 30
m <- 2500
beta_alpha <- seq(1,100,5)
N <- length(beta_alpha)
pwr <- numeric(N)
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { 
  b <- beta_alpha[j]
  sktests <- numeric(m)
  for (i in 1:m) { 
    x <- rbeta(n, b, b)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
plot(beta_alpha, pwr, type = "b", ylim = c(0,0.2), col = "red")
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) 
lines(beta_alpha, pwr+se, lty = 3, col = "blue")
lines(beta_alpha, pwr-se, lty = 3, col = "blue")

## -----------------------------------------------------------------------------
sk <- function(x) {
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}
alpha <- 0.1
n <- 30
m <- 2500
v <- seq(1,20,2)
N <- length(v)
pwr <- numeric(N)
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) {
  b <- v[j]
  sktests <- numeric(m)
  for (i in 1:m) {
    x <- rt(n, v)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
plot(v, pwr, type = "b",xlab = bquote(beta_alpha), ylim = c(0,1), col = "red")
abline(h = 0.1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) 
lines(v, pwr+se, lty = 3, col = "blue")
lines(v, pwr-se, lty = 3, col = "blue")


## -----------------------------------------------------------------------------
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}
ftest <- function(x,y,alpha,n){
  fvalue <- var(x)/var(y)
  u <- as.integer(fvalue < qf(1-alpha/2,n,n) | fvalue> qf(1-alpha/2,n,n))
  return(u)
}
m <-10000
# generate samples under H1 to estimate power
sigma1 <- 1
sigma2 <- 1.5
alpha <- 0.055
n <- c(5,10,20,50,100,200,500,1000)
results <- matrix(0,length(n),3)
results[,1] <- n
power1 <- power2 <- numeric(length(n))
for (i in 1:length(n)) {
  power1[i] <- mean(replicate(m, expr={
    x <- rnorm(n[i], 0, sigma1)
    y <- rnorm(n[i], 0, sigma2)
    count5test(x, y)
  }))
  power2[i] <- mean(replicate(m, expr={
    x <- rnorm(n[i], 0, sigma1)
    y <- rnorm(n[i], 0, sigma2)
    ftest(x,y,alpha,n[i])
  }))
  results[i,2] <- power1[i]
  results[i,3] <- power2[i]
}
print(results)
plot(results[,1], results[,2], col= "blue", ylim = c(0, 1), type = "l",
     xlab = bquote(n), ylab = "power")
lines(results[,1], results[,3], lty = 2, col = "green")
abline(h = alpha, lty = 3, col = "red")
legend("right", 1, c("count5test", "ftest"),
       lty = c(1,2),col = c("blue","green"), inset = .02)


## -----------------------------------------------------------------------------
library(MASS)
n <- 10 #num. of sample
d <- 5 #dim. of sample
m <- 10000 #num. of experiment
alpha <- 0.05
u <- numeric(m)
for (k in 1:m) {
  mean <- c(numeric(d))
  sigma <- diag(d)
  x <- mvrnorm(n,mean,sigma) #generate sample
  hat_sigma <- cov(x)
  inf_hat_sigma <- solve(hat_sigma)
  x_bar <- apply(x,2,mean)
  cv <- qchisq(1-alpha,d*(d+1)*(d+2)/6)
  for (i in 1:n) {
    for (j in 1:n) {
      b <- (t(x[i,]-x_bar) %*% inf_hat_sigma %*% (x[j,]-x_bar))^3/n^2
    }
  }
  u[k] <- as.integer(b < cv | b > cv)
}
print(mean(u))  

## ----warning = FALSE----------------------------------------------------------
library(bootstrap)
x <- law
n <- nrow(x)
b.cor <- function(x,i) cor(x[i,1],x[i,2])
theta.hat <- b.cor(x,1:n)
theta.jack <- numeric(n)
for(i in 1:n){
  theta.jack[i] <- b.cor(x,(1:n)[-i])
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(original=theta.hat,bias.jack=bias.jack,
        se.jack=se.jack),3)

## ----warning = FALSE----------------------------------------------------------
library(boot);set.seed(12345)
dat<-cbind(aircondit$hours)
boot.mean <- function(x,ind) mean(x[ind,1])
de <- boot(data=dat,statistic=boot.mean, R = 1000)
ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
print(ci)

## ----warning = FALSE----------------------------------------------------------
library(bootstrap)
x<-scor
n <- nrow(x)
b.cov <- function(x,i,n){
  M<-((n-1)/n)*cov(x[i,])
  ev <- eigen(M)
  s<-sum(ev$values)
  theta_hat<-max(ev$values)/s
  return(theta_hat)
}
theta.hat <- b.cov(x,(1:n),n)
theta.jack <- numeric(n)
for(i in 1:n){
  theta.jack[i] <- b.cov(x,(1:n)[-i],n)
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(original=theta.hat,bias.jack=bias.jack,
        se.jack=se.jack),6)

## ----warning = FALSE----------------------------------------------------------
library(plyr)
CVgroup <- function(k,datasize,seed){
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]    #将数据分成K份，并生成的完成数据集n
  temp <- sample(n,datasize)   #把n打乱
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])  #dataseq中随机生成k个随机有序数据列
  return(cvlist)
}
library(DAAG); 
data<-ironslag
k <- 26
datasize <- length(data$magnetic)-1   #样本有53个，去掉1个样本，用剩下的52个样本做26折交叉验证，每个subset中包含2个样本
cvlist <- CVgroup(k = k,datasize = datasize,seed = 1206)
cvlist   #生成样本分组情况

## ----warning = FALSE----------------------------------------------------------
e1 <- e2 <- e3 <- e4 <- matrix(0,k,2)
# for k-fold cross validation
# fit models on leave-two-out samples
for (k in 1:26) {
  y <- data$magnetic[-cvlist[[k]]]
  x <- data$chemical[-cvlist[[k]]]
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * data$chemical[cvlist[[k]]]
  e1[k,] <- data$magnetic[cvlist[[k]]] - yhat1
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * data$chemical[cvlist[[k]]] + J2$coef[3] * data$chemical[cvlist[[k]]]^2
  e2[k,] <- data$magnetic[cvlist[[k]]] - yhat2
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * data$chemical[cvlist[[k]]]
  yhat3 <- exp(logyhat3)
  e3[k,] <- data$magnetic[cvlist[[k]]] - yhat3
  J4 <- lm(log(y) ~ log(x))
  logyhat4 <- J4$coef[1] + J4$coef[2] * log(data$chemical[cvlist[[k]]])
  yhat4 <- exp(logyhat4)
  e4[k,] <- data$magnetic[cvlist[[k]]] - yhat4
}

c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))   #计算均方误差

## -----------------------------------------------------------------------------
n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m=1000
K=1:50
R=999
reps=numeric(R)
set.seed(123)
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  z=c(X,Y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  t0=max(c(outx, outy))
  z=c(X,Y)
  # return 1 (reject) or 0 (do not reject H0)
  for(i in 1:R){
    k=sample(K,size=25,replace = FALSE)
    x1=z[k]
    y1=z[-k]
    X <- x1 - mean(x1)
    Y <- y1 - mean(y1)
    outx <- sum(X > max(Y)) + sum(X < min(Y))
    outy <- sum(Y > max(X)) + sum(Y < min(X))
    reps[i]=max(c(outx, outy))
  }
  return(mean(abs(c(t0,reps))>5))
}
alpha= mean(replicate(m, expr={
  x <- rnorm(n1, mu1, sigma1)
  y <- rnorm(n2, mu2, sigma2)
  x <- x - mean(x)
  y <- y - mean(y)
  count5test(x, y)
}))
alpha


## -----------------------------------------------------------------------------
library(RANN) 
library(boot)
library(Ball)
library(energy)
set.seed(123)
mu1=0
mu2=0
sigma1=1
sigma2=1.5
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}
m <- 1e3; k<-3; p<-2; mu <- 0.3; set.seed(12345)
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- matrix(rnorm(n1*p,mu1,sigma1),ncol=p);
y <- matrix(rnorm(n1*p,mu2,sigma2),ncol=p);
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.1;
pow <- colMeans(p.values<alpha)
pow


## -----------------------------------------------------------------------------
library(RANN) 
library(boot)
library(Ball)
library(energy)
set.seed(123)
mu1=0
mu2=0.5
sigma1=1
sigma2=1.5
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}
m <- 1e3; k<-3; p<-2; mu <- 0.3; set.seed(12345)
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- matrix(rnorm(n1*p,mu1,sigma1),ncol=p);
y <- matrix(rnorm(n1*p,mu2,sigma2),ncol=p);
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.1;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
library(RANN) 
library(boot)
library(Ball)
library(energy)
set.seed(123)
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}
m <- 1e3; k<-3; p<-2; mu <- 0.3; set.seed(12345)
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- matrix(rt(n1*p,1),ncol=p)
y <- matrix(rt(n1*p,1),ncol=p)
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.1;
pow1 <- colMeans(p.values<alpha)
pow1
for(i in 1:m){
sigma= sample(c(1, 10), size = n1*p,
replace = TRUE, prob = c(0.4,0.6))
x <- matrix(rnorm(n1*p, 0, sigma),ncol=p);
sigma= sample(c(1, 10), size = n1*p,
replace = TRUE, prob = c(0.4,0.6))
y <- matrix(rnorm(n1*p, 0, sigma),ncol=p);
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.1;
pow2 <- colMeans(p.values<alpha)
pow2

## -----------------------------------------------------------------------------
library(RANN) 
library(boot)
library(Ball)
library(energy)
set.seed(123)
mu1=0
mu2=0
sigma1=1
sigma2=1
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}
m <- 1e3; k<-3; p<-2; mu <- 0.3; set.seed(12345)
n1 <-10; n2 <- 100; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- matrix(rnorm(n1*p,mu1,sigma1),ncol=p);
y <- matrix(rnorm(n2*p,mu2,sigma2),ncol=p);
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.1;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
laplace<-function(x){1/2*exp(-abs(x))}
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (laplace(y) / laplace(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}
set.seed(1234)
N <- 15000
sigma <- c(1, 2, 3, 4)
x0 <- 25
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)
print(c(rw1$k, rw2$k, rw3$k, rw4$k))
print(c(rw1$k/N,rw2$k/N,rw3$k/N,rw4$k/N)) # compute the acceptance rates of each chain

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}
laplace<-function(x){1/2*exp(-abs(x))}
#normal.chain <- function(sigma, x0, N) {
  #generates a Metropolis chain for Normal(0,1)
  #with Normal(X[t], sigma) proposal distribution
  #and starting value X
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (laplace(y) / laplace(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
      }
  }
  return(x)
}
set.seed(1234)
sigma <- c(1, 2, 3, 4) #parameter of proposal distribution
k <- 4 #number of chains to generate
n <- 15000 #length of chains
b <- 100 #burn-in length
#choose overdispersed initial values
x0 <- 25
#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
  X[i, ] <- rw.Metropolis(sigma[i], x0,n)
#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))
#plot psi for the four chains
#par(mfrow=c(2,2))
#for (i in 1:k)
  #plot(psi[i, (b+1):n], type="l",
       #xlab=i, ylab=bquote(psi))
#par(mfrow=c(1,1)) #restore default
#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
#plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
#abline(h=1.1, lty=2)

## -----------------------------------------------------------------------------
kopt <- c(4:25,100,500,1000)
Ak <- matrix(0,length(kopt),1)
set.seed(13425)
for(i in 1:length(kopt)){
  k <- kopt[i]
  f <- function(a) pt(sqrt(a^2*(k-1)/(k-a^2)),df=k-1)-pt(sqrt(a^2*k/(k+1-a^2)),df=k)
  res <- uniroot(f,c(0,sqrt(k))-1e-10,maxiter = 1e6)
  res <- as.numeric(unlist(res)[1])
  if(res>0 & res<sqrt(k)) Ak[i] <- res else Ak[i] <- NA
}
print(Ak)

## -----------------------------------------------------------------------------
na <- 444
nb <- 132
noo <- 361
nab <- 63
n <- na + nb + noo +nab
k <- 10    #times of iters
p <- q <- r <- E <- numeric(k)
p[1] <- q[1] <- r[1] <- 1/3
for (i in 2:k){
  p_hat <- p[i-1]
  q_hat <- q[i-1]
  r_hat <- r[i-1]
  Eaa <- na*p_hat/(p_hat+2*r_hat)
  Ebb <- nb*q_hat/(q_hat+2*r_hat)
  p[i] <- (Eaa+na+nab)/2/n
  q[i] <- (Ebb+nb+nab)/2/n
  r[i] <- 1-p[i]-q[i]
  E[i] <- na*p[i]/(p[i]+2*r[i])*log(p[i]/2/r[i])+na*log(2*p[i]*r[i])+
          nb*q[i]/(q[i]+2*r[i])*log(q[i]/2/r[i])+nb*log(2*q[i]*r[i])+
          2*noo*log(r[i])+nab*log(2*p[i]*q[i])
}
print(cbind(p,q,r,E))

## -----------------------------------------------------------------------------
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)

## -----------------------------------------------------------------------------
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
#use lapply() to fit linear models
res_1 <- lapply(formulas, function(x) lm(x,mtcars))
#use loops to fit linear models
res_2 <- list(NULL)
length(res_2) <- 4
for (i in 1:4) {
  res_2[[i]] <- (lm(formulas[[i]],mtcars))
}
#output the results
print(res_1)
print(res_2)

## -----------------------------------------------------------------------------
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify  =  FALSE
)

## -----------------------------------------------------------------------------
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
#use sapply() and an anonymous function to extract the p-value
res_1 <- sapply(trials,function(x) x$p.value)
#get rid of the anonymous function by using [[ directly
res_2 <- sapply(trials,"[[",3)
# print the results
print(res_1)
print(res_2)

## -----------------------------------------------------------------------------
meapply <- function(x,f,f.value){
  vapply(x,function(x) Map(f,x),f.value)
}

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction('List rwcc(double sigma, double x0, int N) {
NumericVector x(N) ;
x[0] = x0;
NumericVector u ;
u= runif(N);
int k = 0;
for (int i=1;i<N; i+=1) {
NumericVector y = rnorm(1, x[i-1], sigma);
if (u[i] <= (exp(-abs(y[0])) / exp(-abs(x[i-1]))))
{x[i] = y[0] ;
}else  {
x[i] = x[i-1];
k = k + 1;
} }
return List::create(Named("x") =x, Named("k") = k);
}')

## -----------------------------------------------------------------------------
set.seed(1321)
library(microbenchmark)
library(Rcpp)
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (exp(-abs(y)) / exp(-abs(x[i-1]))))
    x[i] <- y else {
      x[i] <- x[i-1]
      k <- k + 1
      } 
    }
  return(list(x=x, k=k))
}
N <- 2000
sigma <- c(.5, 2, 16)
x0 <- 25
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rwcc(sigma[1], x0, N)
qqplot(rw1$x,rw2$x,xlab = "r",ylab = "c++")
abline(0,1)

## -----------------------------------------------------------------------------
set.seed(1425)
ts<-microbenchmark(rw1=rw.Metropolis(sigma[2], x0, N),
                   rw2=rwcc(sigma[2], x0, N))
summary(ts)[,c(1,3,5,6)]

