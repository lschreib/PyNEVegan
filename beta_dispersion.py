import numpy as np

def double_center_matrix(x, na_rm = True):
    """
    The function will take a data matrix and double center it.
    Args:
        x:      a 2D numpy array of shape (n, m) where n is the number of rows and m is the number of columns
        na_rm:  a boolean value indicating whether to remove NA values. Defaults to True. If set to False and there are NA values in the input array, the result will contain NaNs.
    Returns:
        x:      a double centered 2D numpy array of shape (n, m)
    """
    # Check if input matrix is a 2D numpy array
    if not isinstance(x, np.ndarray) or x.ndim != 2:
        raise ValueError("Input matrix must be a 2D numpy array")
    else:
        # Create a copy of the input array to avoid modifying it in place
        x_centered = np.copy(x)
    
        # Calculate the column means
        if na_rm:
            col_means = np.nanmean(x_centered, axis = 0)
        else:
            col_means = np.mean(x_centered, axis = 0)
        # Center the columns
        x_centered = x_centered - col_means

        # Calculate row means
        if na_rm:
            row_means = np.nanmean(x_centered, axis = 1)
        else:
            row_means = np.mean(x_centered, axis = 1)
        # Center the rows
        x_centered = x_centered - row_means[:, np.newaxis] # the np.newaxis conversion is necessary to turn the 1D array into a column vector
    return x_centered

def ordination_median(vectors, group, choices):
    """
    Calculate the spatial median for each group considering specified dimensions.
    Args:
        vectors:    a 2D numpy array where each row is an observation and each column is a variable/dimension.
        group:      a 1D numpy array indicating group membership for each observation.
        choices:    a list of column indices (dimensions) to consider for the median calculation.
    Returns:
        A 2D numpy array where each row is the spatial median for a group based on the specified dimensions.
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
    Calculate the spatial medians for groups by separating positive and negative eigenvalues.
    Args:
        vectors:    a 2D numpy array of eigenvectors
        group:      a 1D numpy array indicating group membership for each observation
        pos:        a 1D boolean numpy array indicating positive eigenvalues
    Returns:
        spatial_medians: a 2D numpy array of spatial medians, concatenating medians calculated from positive and negative eigenvalues.
    """
    axes = np.arange(vectors.shape[1])
    spatial_medians_positive = ordination_median(vectors, group, choices = axes[pos])
    spatial_medians_negative = ordination_median(vectors, group, choices = axes[~pos])
    spatial_medians = np.concatenate((spatial_medians_positive, spatial_medians_negative), axis = 1)
    
    return spatial_medians
