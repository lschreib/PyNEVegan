"""
Ultimately this module will try to emulate the functionality of the adonis2 function of R's vegan package.

Specifically, it will perform a permutational multivariate analysis of variance (PERMANOVA) using distance matrices.

Disclaimer: The goal of the module is not to emulate the full functionality of the adonis2 function, but primarily to simply enable users to replicate my personal use case of the function. I typically call this function as follows:

library(permute)
library(vegan)
perm <- how(nperm = 5000)
result.adonis2 <- adonis2(rpoB_dist_subset ~ .,
                          data =  mapping_chosen[,c("Incubation_time",
                                                    "Temperature",
                                                    "Horizon")],
                          by = "margin",
                          permutations = perm)
result.adonis2
"""

# Okay lets work through this, what does the initial permute function do?
# In the adonis2 function the permutations argument is passed ot the getPermuteMatrix function like this:
# perm <- getPermuteMatrix(permutations, NROW(data), strata = strata)
# Let's see what happens with the argument in the function then ...

# `getPermuteMatrix` <-
#     function(perm, N,  strata = NULL)
# {
#     ## 'perm' is either a single number, a how() structure or a
#     ## permutation matrix
#     if (length(perm) == 1) {
#         perm <- how(nperm = perm)
#     }
#     ## apply 'strata', but only if possible: ignore silently other cases
#     if (!missing(strata) && !is.null(strata)) {
#         if (inherits(perm, "how") && is.null(getBlocks(perm)))
#             setBlocks(perm) <- strata
#     }
#     ## now 'perm' is either a how() or a matrix
#     if (inherits(perm, "how"))
#         perm <- shuffleSet(N, control = perm)
#     else { # matrix: check that it *strictly* integer
#         if(!is.integer(perm) && !all(perm == round(perm)))
#            stop("permutation matrix must be strictly integers: use round()")
#     }
#     ## now 'perm' is a matrix (or always was). If it is a plain
#     ## matrix, set minimal attributes for printing. This is a dirty
#     ## kluge: should be handled more cleanly.
#     if (is.null(attr(perm, "control")))
#         attr(perm, "control") <-
#             structure(list(within=list(type="supplied matrix"),
#                            nperm = nrow(perm)), class = "how")
#     perm
# }

# Observations: 
# 1) The function actually creates a 'how' object intermittently anyway, so there is no need to do this outside the adonis2 function seperately.
# 2) The function as actually doing more then simply shuffling the data, as it offers the option to stratifiy the permutations
# 3) The function doesn't actually shuffle the data, but rather creates a permutation matrix that can be used to shuffle the data later

# Let's try to implement this function in python, but in a simplified manner
import numpy as np
import pandas as pd

#Let's generate some data to test the function
np.random.seed(42)
dist_mat = np.random.rand(100, 100) # 100 samples
strata = np.random.choice([0, 1, 2], 100)

N = dist_mat.shape[0]
nperm = 999

def get_permutation_matrix(N, nperm = 999, strata = None):
    """
    The function will take three arguments:
    Args:
        N:      the number of samples in the dataset 
        nperm:  a single number specifying the number of permutations to do
        strata: an optional 1D numpy array that indicates strata memberships of the individual samples
    Returns:
        perm:   a permutation matrix that can be used to shuffle the data 
    """
    # First we will write the case where strata were supplied as an argument
    if strata is not None:
        if strata.ndim != 1 or strata.size != N:
            raise ValueError("Strata must be a 1D numpy array")
        unique_strata = np.unique(strata)
        perm_matrix = np.empty((nperm, N), dtype=int)
        # For each requested permutation:
        for i in range(nperm):
            # Initialize the permutation indices
            perm_indices = np.arange(N)
            # Go through each stratum and permute the indices within each stratum
            for stratum in unique_strata:
                stratum_indices = np.where(strata == stratum)[0]
                perm_indices[stratum_indices] = np.random.permutation(stratum_indices)
            perm_matrix[i] = perm_indices
    # If no strata were supplied
    else:
        perm_matrix = np.empty((nperm, N), dtype=int)
        for i in range(nperm):
            perm_matrix[i] = np.random.permutation(N)
    
    return perm_matrix
# Done let's export it to a seperate file