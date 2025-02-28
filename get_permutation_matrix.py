import numpy as np

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