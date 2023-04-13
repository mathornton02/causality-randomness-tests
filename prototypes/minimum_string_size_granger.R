# An R code that simulates the granger causality test using a granger caused
#  binary sequence and attempts to determine the minimum bitstring size required
#  in order to determine whether the string was actually granger causal or not. 
library(zoo) # for the rollapply function

generate_granger_causal_sequence <- function(size = 10, lags = 1, seed = 0xFACE){
    set.seed(seed)
    X <- rbinom(size+1, 1, 0.5)
    coefs <- rnorm(lags)
    ptced <- rollapply(X, lags, function(x) {1/(1+exp(-(sum(x*coefs))))})
    Y <- Vectorize(rbinom, vectorize.args = 'prob')(1,1,ptced)
    retvals <- list(X=X, Y=Y, coefs=coefs, ptc=ptced)
    return(retvals)
}

generate_granger_causal_sequence(size = 100, lags = 1, seed = 0xBEEF) -> gcs


# This simulation is for the case where X actually logistically granger causes Y
M <- 1000
sz <- 80000
err <- c()
pvals <- c()
set.seed(0xFACE)

for (i in 1:M){
    X <- rbinom(sz, 1, 0.5)
    coef <- rnorm(1)
    rbinomv <- Vectorize(rbinom, vectorize.args = "prob")
    ptc <- 1/(1+exp(-X*coef))
    rbinomv(1,1,ptc) -> Y 
    glm(Y~X, family=binomial(link = "logit")) -> modres 
    sqerr <- (coef(modres)[2] - coef)^2
    err <- c(err,sqerr)
    pvals <- c(pvals, coef(summary(modres))[,4][2])
    }

mean(err)
sum(pvals <= 0.05)/M

# This simulation is for the case where X does not logistically granger cause Y
M <- 1000
sz <- 10000
pvals <- c()
set.seed(0xFACE)

for (i in 1:M){
    X <- rbinom(sz, 1, 0.5)
    Y <- rbinom(sz, 1, 0.5)
    glm(Y~X, family=binomial(link = "logit")) -> modres 
    pvals <- c(pvals, coef(summary(modres))[,4][2])
    }
sum(pvals <= 0.05)/M
