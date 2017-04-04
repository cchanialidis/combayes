combayes
========

combayes implements Bayesian inference for COM-Poisson regression models using exact samplers. It also provides functions for sampling exactly from the COM-poisson distribution (using rejection sampling) and for evaluating exact bounds for the normalisation constant of the probability mass function of the COM-Poisson distribution.

Fertility example
-----------------


```r
# Load data from library Countr
library(Countr)
data(fertility)
result <- cmpoisreg(y=fertility$Y, X=fertility[,c(-9,-11)], num_samples=1e4, burnin=1e3)
colMeans(result$posterior_beta)
colMeans(result$posterior_delta)
```

PhD publications
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
result <- cmpoisreg(y=phdpublish$art, X=phdpublish[,2:6], num_samples=1e4, burnin=1e3)
colMeans(result$posterior_beta)
colMeans(result$posterior_delta)
```
