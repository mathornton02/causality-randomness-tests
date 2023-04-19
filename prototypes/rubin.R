#  Title: Rubin Causality Test for Randomness 
#
#  Author: Micah A. Thornton 
# 
#  Date: 04-17-2023
#
#  Description: An R code for the rubin causality randomness test. 
#  this is a test which takes an input rbg string and determines whether 
#  there are a greater than expected number of rubin-causal relationships 
#  among the data. 
# 
#  Algorithm: 
# 
#  1. The user sets the size of bitstrings to break the overall string into
#  2. The bitstring is broken into the desired number of substrings, 
#  3. bitstrings are assigned to treatment or control groups based on the first 
#     bit of the substring.
#  4. a number of bits at the beginning of the bitstring are choosen to be the 
#     outcomes, and the remaining set of bits are choosen to be the covariates.
#  5. propensity matching occurs based on the covariates, and treatment effects 
#     are estimated for each of the outcome variables. 
#  6. the test for the significance of a treatment effect should be non-rejecting 
#     for the majority of variables investigated, the proportion of treatment 
#     effects that are significant in other words should be close to zero, the 
#     higher the proportion of significant effects, the more likely the string 
#     was not originally random. 

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

Xsubs <- splitstring(X,500)

# deteremine propensity matching model for covariate bits 
require(rlist)
list.rbind(Xsubs) -> Xmat
effectbits <- 100
Xmat[,1:(ncol(Xmat)-effectbits)] -> Xmat_prop
Xmat_prop <- data.frame(Xmat_prop)
colnames(Xmat_prop) <- paste("X", 1:ncol(Xmat_prop), sep = "")
propensity_match_model <- glm(X1~.,data=Xmat_prop,family=binomial(link="logit"))
propensity_scores <- predict(propensity_match_model,data=Xmat_prop[2:ncol(Xmat_prop)], type="response")

# Get the matched untreated for each treated (or as many as possible)
treatment_vector <- Xmat_prop[,1]
propensity_scores <- cbind(1:length(propensity_scores),propensity_scores)
propensity_scores <- data.frame(propensity_scores) 
propensity_scores_treated <- propensity_scores[treatment_vector == 1,]
propensity_scores_control <- propensity_scores[treatment_vector == 0,]

propensity_distances <- matrix(0, nrow = nrow(propensity_scores_treated), ncol = nrow(propensity_scores_control))
for (i in 1:nrow(propensity_scores_treated)){
    for (j in 1:nrow(propensity_scores_control)){
        propensity_distances[i,j] <- sqrt((propensity_scores_treated[i,"propensity_scores"] - propensity_scores_control[j,"propensity_scores"])^2)
    }
}

kmm <- propensity_distances
colnames(kmm) <- propensity_scores_control[,1]
rownames(kmm) <- propensity_scores_treated[,1]
matches <- data.frame(treated = NULL, control = NULL)
for (i in 1:nrow(propensity_distances)){
    if (is.vector(kmm)){
        propensity_scores_control[unlist(apply(propensity_distances, 2, function(x) {
            sum(x==kmm)
        })) == nrow(propensity_distances),1] -> matchid
        rowid <- names(kmm)[length(kmm)]
        matches <- rbind(matches, data.frame(treated=rowid,control=matchid))
        kmm <- NULL
        break
    } else {
        rowid <- rownames(kmm)[i]
        matchid <- colnames(kmm)[which.min(kmm[i,])]
        matches <- rbind(matches, data.frame(treated=rowid,control=matchid))
        kmm <- kmm[,-which.min(kmm[i,])]
    }
}

# For each pair of matched treated and untreated compute the treatment effects 
#  and determine how many of them are statistically significant. 
distlist <- data.frame(ncol = effectbits)
for (i in 1:nrow(matches)){
    treatmatch <- matches[i,1]
    contrmatch <- matches[i,2]
    Xmat[c(treatmatch,contrmatch),(ncol(Xmat)-effectbits + 1):ncol(Xmat)] -> Xmat_prop
    if (i == 1){
        distlist <- Xmat_prop[1,] - Xmat_prop[2,]
        next
    }
    distlist <- rbind(distlist,(Xmat_prop[1,] - Xmat_prop[2,]))
}

apply(distlist,2, function(x) {
    t.test(x)$p.value
}) -> pvalueList


ks.test(pvalueList, "punif", min(pvalueList), mac(pvalueList))
require(spgs)
chisq.unif.test(pvalueList, bins = 5)
