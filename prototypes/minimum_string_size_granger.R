# An R code that simulates the granger causality test using a granger caused
#  binary sequence and attempts to determine the minimum bitstring size required
#  in order to determine whether the string was actually granger causal or not. 
library(zoo)

generate_granger_causal_sequence <- function(size = 10, lags = 1, seed = 0xFACE){
    set.seed(seed)
    X <- rbinom(size+1, 1, 0.5)
    coefs <- rnorm(lags)
    prob_true_caused <- rollapply(X, lags, function(x) {1/(1+exp(-(sum(x*coefs))))})
    Y <- Vectorize(rbinom, vectorize.args = 'prob')(1,1,prob_true_caused)
    retvals <- list(X, Y, coefs, prob_true_caused)
    print(retvals)
}

