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
# Let's export it to a seperate file - Done

# Let's continue with the actual adonis2 function
# The adonis2 function is a wrapper around the adonis function, which is a wrapper around the betadisper function

# So let's do the betadisper function first
# The betadisper function calculates the distance of each sample to the centroid of the group
# Here is the function:
# `betadisper` <-
#     function(d, group, type = c("median","centroid"), bias.adjust=FALSE,
#              sqrt.dist = FALSE, add = FALSE)
# {
#     ## inline function for double centring. We used .C("dblcen", ...,
#     ## PACKAGE = "stats") which does not dublicate its argument, but
#     ## it was removed from R in r60360 | ripley | 2012-08-22 07:59:00
#     ## UTC (Wed, 22 Aug 2012) "more conversion to .Call, clean up".
#     dblcen <- function(x, na.rm = TRUE) {
#         cnt <- colMeans(x, na.rm = na.rm)
#         x <- sweep(x, 2L, cnt, check.margin = FALSE)
#         cnt <- rowMeans(x, na.rm = na.rm)
#         sweep(x, 1L, cnt, check.margin = FALSE)
#     }
#     ## inline function for spatial medians
#     spatialMed <- function(vectors, group, pos) {
#         axes <- seq_len(NCOL(vectors))
#         spMedPos <- ordimedian(vectors, group, choices = axes[pos])
#         spMedNeg <- ordimedian(vectors, group, choices = axes[!pos])
#         cbind(spMedPos, spMedNeg)
#     }
#     ## inline function for centroids
#     centroidFUN <- function(vec, group) {
#         cent <- apply(vec, 2,
#                       function(x, group) tapply(x, INDEX = group, FUN = mean),
#                       group = group)
#         if(!is.matrix(cent)) { ## if only 1 group, cent is vector
#             cent <- matrix(cent, nrow = 1,
#                            dimnames = list(as.character(levels(group)),
#                            paste0("Dim", seq_len(NCOL(vec)))))
#         }
#         cent
#     }
#     ## inline function for distance computation
#     Resids <- function(x, c) {
#         if(is.matrix(c))
#             d <- x - c
#         else
#             d <- sweep(x, 2, c)
#         rowSums(d^2)
#     }
#     ## Tolerance for zero Eigenvalues
#     TOL <- sqrt(.Machine$double.eps)
#     ## uses code from stats:::cmdscale by R Core Development Team
#     if(!inherits(d, "dist"))
#         stop("distances 'd' must be a 'dist' object")
#     ## Someone really tried to analyse correlation like object in range -1..+1
#     if (any(d < -TOL, na.rm = TRUE))
#         stop("dissimilarities 'd' must be non-negative")
#     ## adjust to avoid negative eigenvalues (if they disturb you)
#     if (sqrt.dist)
#         d <- sqrt(d)
#     if (is.logical(add) && isTRUE(add))
#         add <- "lingoes"
#     if (is.character(add)) {
#         add <- match.arg(add, c("lingoes", "cailliez"))
#         if (add == "lingoes") {
#             ac <- addLingoes(as.matrix(d))
#             d <- sqrt(d^2 + 2 * ac)
#         }
#         else if (add == "cailliez") {
#             ac <- addCailliez(as.matrix(d))
#             d <- d + ac
#         }
#     }
#     if(missing(type))
#         type <- "median"
#     type <- match.arg(type)
#     ## checks for groups - need to be a factor for later
#     group <- if(!is.factor(group)) {
#         as.factor(group)
#     } else { ## if already a factor, drop empty levels
#         droplevels(group, exclude = NA) # need exclude = NA under Rdevel r71113
#     }
#     n <- attr(d, "Size")
#     x <- matrix(0, ncol = n, nrow = n)
#     x[row(x) > col(x)] <- d^2
#     ## site labels
#     labs <- attr(d, "Labels")
#     ## remove NAs in group
#     if(any(gr.na <- is.na(group))) {
#         group <- group[!gr.na]
#         x <- x[!gr.na, !gr.na]
#         ## update n otherwise C call crashes
#         n <- n - sum(gr.na)
#         ## update labels
#         labs <- labs[!gr.na]
#         message("missing observations due to 'group' removed")
#     }
#     ## remove NA's in d
#     if(any(x.na <- apply(x, 1, function(x) any(is.na(x))))) {
#         x <- x[!x.na, !x.na]
#         group <- group[!x.na]
#         ## update n otherwise C call crashes
#         n <- n - sum(x.na)
#         ## update labels
#         labs <- labs[!x.na]
#         message("missing observations due to 'd' removed")
#     }
#     x <- x + t(x)
#     x <- dblcen(x)
#     e <- eigen(-x/2, symmetric = TRUE)
#     vectors <- e$vectors
#     eig <- e$values
#     ## Remove zero eigenvalues
#     eig <- eig[(want <- abs(eig) > max(TOL, TOL * eig[1L]))]
#     ## scale Eigenvectors
#     vectors <- vectors[, want, drop = FALSE] %*% diag(sqrt(abs(eig)),
#                                nrow = length(eig))
#     ## store which are the positive eigenvalues
#     pos <- eig > 0
#     ## group centroids in PCoA space
#     centroids <-
#         switch(type,
#                centroid = centroidFUN(vectors, group),
#                median = spatialMed(vectors, group, pos)
#                )
#     ## for each of the groups, calculate distance to centroid for
#     ## observation in the group
#     ## Uses in-line Resids function as we want LAD residuals for
#     ## median method, and LSQ residuals for centroid method
#     dist.pos <- Resids(vectors[, pos, drop=FALSE],
#                        centroids[group, pos, drop=FALSE])
#     dist.neg <- 0
#     if(any(!pos))
#         dist.neg <- Resids(vectors[, !pos, drop=FALSE],
#                            centroids[group, !pos, drop=FALSE])

#     ## zij are the distances of each point to its group centroid
#     if (any(dist.neg > dist.pos)) {
#         ## Negative squared distances give complex valued distances:
#         ## take only the real part (which is zero). Github issue #306.
#         warning("some squared distances are negative and changed to zero")
#         zij <- Re(sqrt(as.complex(dist.pos - dist.neg)))
#     } else {
#         zij <- sqrt(dist.pos - dist.neg)
#     }
#     if (bias.adjust) {
#         n.group <- as.vector(table(group))
#         zij <- zij*sqrt(n.group[group]/(n.group[group]-1))
#     }
#     ## pre-compute group mean distance to centroid/median for `print` method
#     grp.zij <- tapply(zij, group, "mean")
#     ## add in correct labels
#     if (any(want))
#         colnames(vectors) <- names(eig) <-
#             paste("PCoA", seq_along(eig), sep = "")
#     if(is.matrix(centroids))
#         colnames(centroids) <- names(eig)
#     else
#         names(centroids) <- names(eig)
#     rownames(vectors) <- names(zij) <- labs
#     retval <- list(eig = eig, vectors = vectors, distances = zij,
#                    group = group, centroids = centroids,
#                    group.distances = grp.zij, call = match.call())
#     class(retval) <- "betadisper"
#     attr(retval, "method") <- attr(d, "method")
#     attr(retval, "type") <- type
#     attr(retval, "bias.adjust") <- bias.adjust
#     retval
# }

# Okay let's do this step by step
# Let's generate a distance matrix
from scipy.spatial import distance
import numpy as np
np.random.seed(42)
data = np.random.rand(100, 10)
dist_mat = distance.pdist(data, 'euclidean')
dist_mat = distance.squareform(dist_mat)
# Let's generate a group vector
group = np.random.choice([0, 1, 2], 100)

# Now let's start with the integrated helper functions
import numpy as np

def double_center_matrix(x, na_rm = True):
    """
    The function will take a data matrix and double center it
    Args:
        x:      a 2D numpy array
        na_rm:  a boolean value indicating whether to remove NA values
    Returns:
        x:      a double centered 2D numpy array
    """
    # Calculate the column means
    if na_rm:
        col_means = np.nanmean(x, axis = 0)
    else:
        col_means = np.mean(x, axis = 0)
    # Center the columns
    x = x - col_means

    # Calculate row means
    if na_rm:
        row_means = np.nanmean(x, axis = 1)
    else:
        row_means = np.mean(x, axis = 1)
    # Center the rows
    x = x - row_means[:, np.newaxis]
    
    return x


# Now this one:
#     spatialMed <- function(vectors, group, pos) {
#         axes <- seq_len(NCOL(vectors))
#         spMedPos <- ordimedian(vectors, group, choices = axes[pos])
#         spMedNeg <- ordimedian(vectors, group, choices = axes[!pos])
#         cbind(spMedPos, spMedNeg)
#     }
def ordination_median(vectors, group, choices):
    """
    Calculate the spatial median for each group.
    Args:
        vectors: A 2D numpy array where each row is an observation and each column is a variable.
        group: A 1D numpy array indicating group membership for each observation.
        choices: A list of column indices to consider for the median calculation.
    Returns:
        A 2D numpy array where each row is the spatial median for a group.
    """
    unique_groups = np.unique(group)
    medians = []
    for grp in unique_groups:
        grp_vectors = vectors[group == grp][:, choices]
        median = np.median(grp_vectors, axis=0)
        medians.append(median)
    return np.array(medians)

def spatial_medians(vectors, group, pos):
    """
    The function will calculate the spatial medians of the group
    Args:
        vectors:    a 2D numpy array of eigenvectors
        group:      a 1D numpy array of group memberships
        pos:        a 1D boolean numpy array indicating positive eigenvalues
    Returns:
        spatial_medians: a 2D numpy array of spatial medians
    """
    axes = np.arange(vectors.shape[1])
    spMedPos = ordination_median(vectors, group, choices = axes[pos])
    spMedNeg = ordination_median(vectors, group, choices = axes[~pos])
    spatial_medians = np.concatenate((spMedPos, spMedNeg), axis = 1)
    
    return spatial_medians
