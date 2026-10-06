"""

adapted from SuPreMo scoring.py

"""

from itertools import chain
import math

from matplotlib.collections import LineCollection
from scipy import stats
from scipy.sparse import coo_matrix
from skimage.metrics import structural_similarity as ssim
from skimage.transform import resize
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BINS = 448  # length of side of square matrix
DIAG_OFFSET = 2  # if the diagonal is offset by a number of bins:
input_map_size = 2 ** 20

# -------------------- HELPER FUNCTIONS ----------------- #

# Set lower triangle to nans
# Assert no nans

def remove_missing_points_flat(flat_a, flat_b):
    """
    Remove points in vector where either flat_a OR flat_b are NaN
    Input:
        flat_a: numpy vector
        flat_b: numpy vector
    Returns:
        flat_a: numpy vector with nan positions removed
        flat_b: numpy vector with nan positions removed
    """
    mask = np.logical_or(np.isnan(flat_a), np.isnan(flat_b))
    return flat_a[~mask], flat_b[~mask]

def fill_missing_points_map(map_a, map_b, fill=0):
    """
    Fills points in matrix where either map_a OR map_b are NaN
    Input:
        map_a: numpy array
        map_b: numpy array
        fill: value to fill missing points with (default: 0)
    Returns:
        map_a: numpy array with missing positions filled with specified value
        map_b: numpy array with missing positions filled with specified value
    """
    mask = np.logical_or(np.isnan(map_a), np.isnan(map_b))
    map_a_copy, map_b_copy = np.copy(map_a), np.copy(map_b)
    map_a_copy[mask] = fill
    map_b_copy[mask] = fill
    return map_a_copy, map_b_copy

def fill_tril(contact_map, fill=np.nan):
    """ Fill the lower triangle of a matrix with a specified value
    Input:
        map: n x n numpy array
        fill: optional value for what to fill the lower triangle with (default: nan)
    Returns:
        map_filled: n x n numpy array with lower triangle filled
    """
    map_filled = contact_map.copy()
    fill_indices = np.tril_indices(map_filled.shape[0])
    map_filled[fill_indices] = fill
    return map_filled

def spearman_1D(vector_a, vector_b):
    """
    Function to calculate the spearman correlation between two 1D arrays.

    Input:
        vector_a: 1D numpy array of length(n)
        vector_b: 1D numpy array of length(n)
    Returns:
        scalar value
    """

    vector_a, vector_b = remove_missing_points_flat(vector_a, vector_b)

    spearmanr_val, pval = stats.spearmanr(vector_a, vector_b)
    return spearmanr_val

def pearson_1D(vector_a, vector_b):
    """
    Function to calculate the pearson correlation between two 1D arrays.

    Input:
        vector_a: 1D numpy array of length(n)
        vector_b: 1D numpy array of length(n)
    Returns:
        scalar value
    """

    vector_a, vector_b = remove_missing_points_flat(vector_a, vector_b)

    pearsonr_val, pval = stats.pearsonr(vector_a, vector_b)
    return pearsonr_val

def mse_1D(vector_a, vector_b):
    """
    Function to calculate the mean squared error between two 1D arrays.

    Input:
        vector_a: 1D numpy array of length(n)
        vector_b: 1D numpy array of length(n)

    Returns:
        scalar value
    """

    vector_a, vector_b = remove_missing_points_flat(vector_a, vector_b)

    mse = np.mean(np.square(vector_a - vector_b))
    return mse

#### MSE #####
def mse(map_a, map_b):
    """
    Mean Squared Error
    Input:
        map_a: n x n numpy array
        map_b: n x n numpy array
    Output:
        scalar: MSE between flattened map_a and map_b
    """

    flat_a = map_a.reshape(-1)
    flat_b = map_b.reshape(-1)
    mse = mse_1D(flat_a, flat_b)

    return mse

# -------------------- BASIC METHODS -------------------- #

#### SPEARMAN'S RANK CORRELATION COEF ####
def spearman(map_a, map_b):
    """
    Spearman correlation between two maps.
    Input:
        map_a: n x n numpy array
        map_b: n x n numpy array
    Returns:
        scalar: spearmanr
    """

    flat_a = map_a.reshape(-1)
    flat_b = map_b.reshape(-1)

    spearmanr = spearman_1D(flat_a, flat_b)
    return spearmanr

#### Pearson correlation #####
def pearson(map_a, map_b):
    """
    Spearman correlation between two maps.
    Input:
        map_a: n x n numpy array
        map_b: n x n numpy array
    Returns:
        scalar: pearsonr
    """

    flat_a = map_a.reshape(-1)
    flat_b = map_b.reshape(-1)

    pearsonr = pearson_1D(flat_a, flat_b)
    return pearsonr