# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma", version = "3.8")
# BiocManager::install("multtest", version = "3.8")
# BiocManager::install("qvalue", version = "3.8")
# install.packages("cp4p")


# There are various estimators of pi0 implemented in R. All of them, which use the
# ecdf of the p-values, could be combined with our marginal bootstrap algorithm. 
# Some example are given below.

# The function calculate_Pi0(p, lambda) with p-value vector p and tuning paramter 
# lambda is used in the simulations of the paper. This function can be replaced
# with the following functions. The following methods only require p as input.


# Phipson (2013)
library(limma)
propTrueNull(p, method = "lfdr")


# vector containing the following methods:
# Storey and Tibshirani (2003), Storey et al. (2004), Jiang and Doerge (2008), Nettleton et al. (2006),
# Langaas et al. (2005), Pounds and Cheng (2006), Benjamini and Hochberg (1995), Wang et al. (2011)
library(cp4p)
estim.pi0(p)[[1]]


# For further details see the documentation of these functions.