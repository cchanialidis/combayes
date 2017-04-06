combayes
========

combayes implements Bayesian inference for COM-Poisson regression models using exact samplers. It also provides functions for sampling exactly from the COM-poisson distribution (using rejection sampling) and for evaluating exact bounds for the normalisation constant of the probability mass function of the COM-Poisson distribution.



Sampling from COM-Poisson distributions with different dispersion levels
---------------------

```r
# Choose sample size
n <- 100
# Sampling from an underdispersed COM-Poisson distribution
comp_under <- rcmpois(mu=10,nu=2,n=n)
# Sampling from a COM-Poisson distribution where nu=1 (i.e. Poisson distribution)
comp_poisson <- rcmpois(mu=10,nu=1,n=n)
# Sampling from an overdispersed COM-Poisson distribution
comp_over <- rcmpois(mu=10,nu=0.5,n=n)
# Check mean and variance for each one
distributions <- matrix(0,nrow = n,ncol=3)
distributions[,1]<- comp_under
distributions[,2]<- comp_poisson 
distributions[,3]<- comp_over
apply(distributions,2,mean)
apply(distributions,2,var)
```

Estimating the logarithm of the normalisation constant 
----------------

```r
logzcmpois(mu=10,nu=2)
logzcmpois(mu=10,nu=1)
logzcmpois(mu=10,nu=0.5)
```

Estimating the probability mass function 
-----------------------

```r
# Estimating p.m.f. of COM-Poisson distribution with different dispersion levels
 x <- 0:25
dcmpois(x, mu=10, nu=1)
dcmpois(x, mu=10, nu=0.5)
dcmpois(x, mu=10, nu=2)
matplot(x, cbind(dcmpois(x, mu=10, nu=1),
                 dcmpois(x, mu=10, nu=0.5),
                 dcmpois(x, mu=10, nu=2)), 
                 type="o", col=2:4, pch=16, ylab="p.m.f.")  
legend("topright", col=2:4, lty=1:3, 
                 c(expression(nu*"="*1),
                   expression(nu*"="*0.5),
                   expression(nu*"="*2)))
```


Bayesian COM-Poisson regression on the fertility data
-----------------


```r
# Load data from library Countr
library(Countr)
data(fertility)
# Standardise all non-binary covariates
fertility[,c(2,9,10)] <- scale(fertility[,c(2,9,10)],center=TRUE,scale=TRUE)
result <- cmpoisreg(y=fertility$Y, X=fertility[,-11], num_samples=1e4, burnin=1e3)
colMeans(result$posterior_beta)
colMeans(result$posterior_delta)
```

Bayesian COM-Poisson regression on the PhD publications data
----------------


```r
# Load data from library Rchoice
library(Rchoice)
data(Articles)
# Focusing only on the students with at least one publication
phdpublish <- subset(Articles, art>0)
phdpublish <- transform(phdpublish, art=art-1)
# Standardise all non-binary covariates
phdpublish <- cbind(phdpublish[,c(1,2,3)],scale(phdpublish[,-c(1,2,3)],center=TRUE,scale=TRUE))
result <- cmpoisreg(y=phdpublish$art, X=phdpublish[,2:6], num_samples=1e4, burnin=1e3,prior_var_beta=diag(6),prior_var_delta=diag(6))
colMeans(result$posterior_beta)
colMeans(result$posterior_delta)
```
