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

