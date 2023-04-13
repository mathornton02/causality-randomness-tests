#  Title: Granger Causality Test for Randomness 
#
#  Author: Micah A. Thornton 
# 
#  Date: 04-13-2023
#
#  Description: An R code for the granger causality randomness test. 
#  this is a test which takes an input rbg string and determines whether 
#  there are a greater than expected number of granger-causal relationships 
#  among the data. 
# 
#  Algorithm: 
# 
#  1. The user sets the maximum number of lags to check 
#  2. The user sets the size of bitstrings to break the overall string into
#  3. The bitstring is broken into the desired number of substrings, 
#  4. pairwise granger causality tests using the maximum lag are conducted 
#     a. For each string, we conduct logistic auto-regression with the desired 
#         number of lags 
#     b. including the significant lag values only, include the lagged values 
#         of each other sequence and determine whether any are significant. 
#     c. if any are significant consider a granger causal relationship to be 
#         discovered. 
#     d. count the number of granger causal relationships that are discovered 
#         among the data. 
#  5. Some spurious granger causal relationships are to be expected under the 
#     null hypothesis that the data are truly random, but if there is a greater 
#     than expected number of granger causal relationships discovered in the 
#     data then the test should reject the null hypothesis that the data are 
#     a truly random bit string. 

splitsize <- 10000
maxlags <- 10

# First in this script I will create an example run
X <- rbinom(100000,1,0.5)

splitstring <- function(string,splitsize){
    # The splitsize should be less than or equal to half the string length 
    #  so that there are at least two causal relationships investigated. 
    if (splitsize > length(string/2)){
        cat("Error: splitsize not less than or equal to half the string length")
        stop()
    }
    split_vector <- split(string, ceiling(seq_along(string)/splitsize))
    return(split_vector)
}

Xsubs <- splitstring(X,10000)

lag_matrix <- function(x, n) {
  # Check if n is a positive integer
  if (n <= 0 || n %% 1 != 0) {
    stop("n must be a positive integer")
  }
  # Check if the input sequence is long enough
  if (length(x) < n + 1) {
    stop("The input sequence must be longer than n")
  }
  # Initialize an empty matrix with appropriate dimensions
  result_matrix <- matrix(0, nrow = length(x) - n, ncol = n + 1)
  # Fill in the matrix with the trailing part and lagged values
  for (i in 1:(length(x) - n)) {
    result_matrix[i, ] <- x[(i + n):(i)]
  }
  return(result_matrix)
}

gcdat <- c()
for (subs in Xsubs){ 
    lag_matrix(subs,5) -> tm 
    colnames(tm) <- paste("Y", 1:ncol(tm), sep = "")
    tm <- data.frame(tm) 
    glm(Y1~0+.,data=tm,family=binomial(link="logit")) -> mod1res
    names(which(coef(summary(mod1res))[,4] <= 0.1)) -> Ykeep
    for (sx in Xsubs){
        if (identical(sx,subs)){
            gcdat <- c(gcdat,0)
            next
        } else {
            lag_matrix(sx,5) -> lsx
            if (length(Ykeep) > 0){
                cbind(tm[,'Y1'],lsx[,-1],tm[,Ykeep]) -> tm2
                colnames(tm2) <- c("Y1",paste("X",2:ncol(lsx),sep=""),Ykeep)
            } else {
                tm2 <- cbind(tm[,'Y1'],lsx[,-1])
                colnames(tm2) <- c("Y1", paste("X",2:ncol(lsx),sep=""))
            }
            tm2 <- data.frame(tm2)
            glm(Y1~0+., data=tm2, family=binomial(link="logit")) -> mod2res
            names(which(coef(summary(mod2res))[,4] <= 0.1)) -> XKeep
            if (length(setdiff(XKeep,Ykeep)) != 0){
                gcdat <- c(gcdat,1)
            } else {
                gcdat <- c(gcdat,0)
            }
        }
    }
}

matrix(gcdat,nrow = length(Xsubs), ncol = length(Xsubs))
