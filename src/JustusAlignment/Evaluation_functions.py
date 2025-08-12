##########################################################################################
# This script provides the evalutation functions for embedding networks.
# It also provides a set of functions that alter the datasets
# for more concise evaluations of the procedure
##########################################################################################
# Imports:
# External Imports:
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import matplotlib.cm as cm
import plotly.graph_objects as go
from scipy.stats import ortho_group
# Internal Imports:
from Main_functions import *
from Datasets_Plotting import * # plot_3d() and plot_2d() (+Embedding versions)
##########################################################################################
# Functions to calculate performance metrics:
##########################################################################################
# Function to calculate error of summed-mean-distances for corresponding points:
def calculate_mean_distance_error(f_net, g_net, X_tensor, Y_tensor):
    fX = f_net(X_tensor)
    gY = g_net(Y_tensor)
    # Standardize the data
    eps = 1e-8
    fX = (fX - fX.mean(dim=0)) / (fX.std(dim=0) + eps)
    gY = (gY - gY.mean(dim=0)) / (gY.std(dim=0) + eps)
    return torch.norm(fX - gY, dim=1).mean().item()
# Neighborassignment based metric:
# Assert you have (two) dataset from which you know the correct Neighborgraph,
# Reconstruct the KNN graph (unweighted/ weighted?) from the embedding space.
# Then compare the agreement (one could do this with different k's and observe
# the agreement for both graphs (based on counting metrics?). In the weighted
# case, it is also possible to generate a continuous error term).
###########################################################################################
# Functions to alter datassets in a structured manner:
###########################################################################################
# Function to increase random scattering of datapoints in their respective dimension:
import numpy as np
def scatter_points(X, p=1.0, seed=None):
    """
    Scatters each point in X in a random direction by a distance in [0, p].
    Parameters:
    - X : np.ndarray of shape (n, m)
        Original dataset.
    - p : float
        Maximum scattering distance.
    - seed : int or None
        Random seed for reproducibility.
    Returns:
    - X_noisy : np.ndarray of shape (n, m)
        Scattered dataset.
    """
    if seed is not None:
        np.random.seed(seed)
    n, m = X.shape
    # Step 1: Sample distances uniformly from [0, p]
    distances = np.random.uniform(0, p, size=n)  # shape: (n,)
    # Step 2: Sample random directions (unit vectors) in R^m
    random_directions = np.random.normal(0, 1, size=(n, m))  # Gaussian for isotropy
    norms = np.linalg.norm(random_directions, axis=1, keepdims=True)
    unit_directions = random_directions / norms  # shape: (n, m)
    # Step 3: Scale directions by distances
    displacements = unit_directions * distances[:, np.newaxis]  # broadcast multiply
    # Step 4: Add displacements to original data
    X_noisy = X + displacements
    return X_noisy
# Function to introduce randomly scattered outlier points in the dataset:
def add_spherical_outliers(X, n_outliers=10, radius=1.0, seed=None):
    """
    Adds outliers sampled uniformly from a hypersphere of specified radius around the dataset's center.
    Parameters:
    - X : np.ndarray of shape (n, m)
        Original dataset.
    - n_outliers : int
        Number of outliers to add.
    - radius : float
        Radius of the hypersphere.
    - seed : int or None
        Random seed for reproducibility.
    Returns:
    - X_augmented : np.ndarray of shape (n + n_outliers, m)
        Dataset with outliers appended.
    - labels : np.ndarray of shape (n + n_outliers,)
        0 for original data, 1 for outliers.
    """
    if seed is not None:
        np.random.seed(seed)
    n, m = X.shape
    center = np.mean(X, axis=0)
    # Step 1: Sample random directions (uniform on sphere)
    directions = np.random.normal(0, 1, size=(n_outliers, m))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    # Step 2: Sample radii with volume-correct distribution
    # This avoids clustering at the center
    uniform_radii = np.random.uniform(0, 1, size=n_outliers) ** (1 / m)
    scaled_radii = uniform_radii * radius
    # Step 3: Compute outlier positions
    outliers = center + directions * scaled_radii[:, np.newaxis]
    # Step 4: Append to dataset
    X_augmented = np.vstack([X, outliers])
    labels = np.concatenate([np.zeros(n), np.ones(n_outliers)])
    return X_augmented, labels
# Function to project a dataspace into high dimensional space:
# This can include isomorphic or non-isomorphic transformations
def project_to_high_dim(X, target_dim, mode='isometric', seed=None):
    """
    Projects data X from d to D dimensions bijectively.
    Parameters:
    - X: np.ndarray, shape (n_samples, d)
    - target_dim: int, target dimension D > d
    - mode: 'isometric' or 'random'
    - seed: int or None, for reproducibility
    Returns:
    - X_proj: np.ndarray, shape (n_samples, D)
    - inverse_fn: function to recover original X from X_proj
    """
    assert target_dim > X.shape[1], "Target dim must be higher than original dim"
    n, d = X.shape
    D = target_dim
    rng = np.random.default_rng(seed)
    # Step 1: Zero-pad
    X_pad = np.zeros((n, D))
    X_pad[:, :d] = X
    # Step 2: Apply orthogonal or random transform
    if mode == 'isometric':
        # Generate a D x D orthogonal matrix (preserves distances)
        Q = ortho_group.rvs(dim=D, random_state=rng)
    elif mode == 'random':
        # Random invertible (not necessarily orthogonal) matrix
        Q = rng.standard_normal((D, D))
        while np.linalg.matrix_rank(Q) < D:
            Q = rng.standard_normal((D, D))  # ensure invertible
    else:
        raise ValueError("Mode must be 'isometric' or 'random'")
    X_proj = X_pad @ Q  # shape: (n, D)
    return X_proj