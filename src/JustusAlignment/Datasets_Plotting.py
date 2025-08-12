# This script contains all the functions required to generate the test data for the paradigm tests. 
# Generally, all functions generate a Dataset and a corresponding colour-code which highligts particular regions,
# which can then be decerned before and after the alignment process  

################################################################################################################
# Imports 
##############################################################################################################

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from itertools import combinations
import torch

##############################################################################################################
# Generalized (rudiment) plotting Functions. 
# All Plotting functions take data and labels as input. Point size, colour and title are adaptible
##############################################################################################################
 
# 3D plotting function
def plot_3d(X, labels=None, point_size=5, title="3D Scatter Plot", cmap="Viridis"):
    """
    Parameters:
    - X : numpy.ndarray of shape (n_samples, 3)
        3D coordinates of the data points.
    - labels : array-like of shape (n_samples,), optional
        Values used for coloring the data points. If None, the Z-axis values are used.
    - point_size : int, optional (default=5)
        Size of the plotted points.
    - title : str, optional
        Title of the plot.
    - cmap : str, optional (default="Viridis")
        Colormap used to map labels to colors.

    Returns:
    - None (writes and opens an interactive HTML plot)
    """
    
    x, y, z = X[:, 0], X[:, 1], X[:, 2]

    fig = go.Figure(data=[go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers',
        marker=dict(
            size=point_size,
            color=labels if labels is not None else z,
            colorscale=cmap,
            opacity=0.8,
            colorbar=dict(title='Label') if labels is not None else None
        )
    )])

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z'
        ),
        margin=dict(l=0, r=0, b=0, t=40)
    )

    fig.write_html("3d_plot.html", auto_open=True)

# 2D plotting function: 
def plot_2d(X, labels=None, title="2D Scatter Plot", cmap='viridis', point_size=10):
    """
    Parameters:
    - X : numpy.ndarray of shape (n_samples, 2)
        2D coordinates of the data points.
    - labels : array-like of shape (n_samples,), optional
        Values used for coloring the data points.
    - title : str, optional
        Title of the plot.
    - cmap : str, optional (default='viridis')
        Colormap used for the labels.
    - point_size : int, optional (default=10)
        Size of the plotted points.

    Returns:
    - None (displays the plot)
    """
    X = np.asarray(X)
    if X.shape[1] != 2:
        raise ValueError("Input X must have exactly 2 columns for 2D plotting.")

    plt.figure(figsize=(6, 5))
    scatter = plt.scatter(X[:, 0], X[:, 1], c=labels, cmap=cmap, s=point_size)

    if labels is not None:
        plt.colorbar(scatter, label="Labels")

    plt.title(title)
    plt.axis("equal")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


# 2D plotting function for joint embeddings: 
def plot_embeddings_2d(X, Y, labels_X, f_net, g_net, title="Joint Embedding Plot"):
    device = next(f_net.parameters()).device

    # Compute embeddings
    X_t = torch.tensor(X, dtype=torch.float32, device=device)
    Y_t = torch.tensor(Y, dtype=torch.float32, device=device)
    fX = f_net(X_t).cpu().detach().numpy()
    gY = g_net(Y_t).cpu().detach().numpy()

    # Plot
    plt.figure(figsize=(7, 6))
    plt.title(title)

    plt.scatter(fX[:, 0], fX[:, 1], c=labels_X, cmap='viridis', label='f(X)', alpha=0.8)
    plt.scatter(gY[:, 0], gY[:, 1], c="black", cmap='plasma', label='g(Y)', alpha=1, marker='x')

    plt.xlabel("Dim 1")
    plt.ylabel("Dim 2")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# 3d plotting function for joint embeddings:

def plot_joint_embeddings_3d(X, Y, labels_X, f_net, g_net, point_size=5,
                             title="3D Joint Embedding", cmap_X="Viridis", cmap_Y="Plasma"):
    """
    Plots 3D embeddings of two datasets after transformation with neural networks.

    Parameters:
    - X, Y : np.ndarray, shape (n_samples, features)
        Input datasets.
    - labels_X : array-like, shape (n_samples,)
        Labels used for coloring one of the datasets (X).
    - f_net, g_net : torch.nn.Module
        Neural networks to transform X and Y respectively.
    - point_size : int
        Size of the scatter points.
    - title : str
        Title of the plot.
    - cmap_X, cmap_Y : str
        Colormaps for X and Y datasets.

    Returns:
    - None (saves and opens interactive HTML plot)
    """
    
    device = next(f_net.parameters()).device

    # Compute embeddings
    with torch.no_grad():
        fX = f_net(torch.tensor(X, dtype=torch.float32, device=device)).cpu().numpy()
        gY = g_net(torch.tensor(Y, dtype=torch.float32, device=device)).cpu().numpy()

    # Check dimensions
    assert fX.shape[1] >= 3 and gY.shape[1] >= 3, "Embeddings must have at least 3 dimensions"

    # Take first 3 dimensions for plotting
    fX_plot = fX[:, :3]
    gY_plot = gY[:, :3]

    # Scatter3d traces
    trace_fX = go.Scatter3d(
        x=fX_plot[:, 0], y=fX_plot[:, 1], z=fX_plot[:, 2],
        mode='markers',
        name='f(X)',
        marker=dict(
            size=point_size,
            color=labels_X,
            colorscale=cmap_X,
            opacity=0.8,
            symbol='circle',
            colorbar=dict(title='Labels f(X)')
        )
    )

    trace_gY = go.Scatter3d(
        x=gY_plot[:, 0], y=gY_plot[:, 1], z=gY_plot[:, 2],
        mode='markers',
        name='g(Y)',
        marker=dict(
            size=point_size,
            color="black", 
            colorscale=cmap_Y,
            symbol='cross',
            opacity=0.8,
            colorbar=dict(title='Labels g(Y)')
        )
    )

    # Layout and figure
    fig = go.Figure(data=[trace_fX, trace_gY])

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='Dim 1',
            yaxis_title='Dim 2',
            zaxis_title='Dim 3'
        ),
        legend=dict(x=0.8, y=0.9),
        margin=dict(l=0, r=0, b=0, t=40)
    )

    fig.write_html("joint_3d_plot.html", auto_open=True)



##############################################################################################################
# Data/manifold generation functions: 
# All Functions provide a Noise parameter and a size parameter that specifies the number of points sampled 
# on the manifold. Each function also contains specific parameter based on the manifold
##############################################################################################################

# Roman T 
def generate_roman_T(total_points=1500, noise_std=0.0):
    """
    Generate a Roman-style 'T' shape with area-proportional point density.

    Parameters:
    - total_points (int): Total number of points in the shape.
    - noise_std (float): Standard deviation of Gaussian noise added to the points.

    Returns:
    - T (ndarray): Array of shape (total_points, 2) containing the coordinates.
    - labels (ndarray): Array of shape (total_points,) containing part labels (0: top, 1: stem, 2: serif).
    """
    rng = np.random.default_rng(0)

    # Define dimensions
    top_width, top_height = 2.0, 0.3
    stem_width, stem_height = 0.3, 1.5
    serif_width, serif_height = 0.6, 0.15

    # Compute areas
    area_top = top_width * top_height
    area_stem = stem_width * stem_height
    area_serif = serif_width * serif_height
    total_area = area_top + area_stem + area_serif

    # Allocate points proportionally
    n_top = int(total_points * area_top / total_area)
    n_stem = int(total_points * area_stem / total_area)
    n_serif = total_points - n_top - n_stem  # fill remainder

    # Generate top bar
    x_top = rng.uniform(-top_width / 2, top_width / 2, n_top)
    y_top = rng.uniform(0.8, 0.8 + top_height, n_top)
    top = np.stack([x_top, y_top], axis=1)

    # Generate stem
    x_stem = rng.uniform(-stem_width / 2, stem_width / 2, n_stem)
    y_stem = rng.uniform(0.8 - stem_height, 0.8, n_stem)
    stem = np.stack([x_stem, y_stem], axis=1)

    # Generate serif
    x_serif = rng.uniform(-serif_width / 2, serif_width / 2, n_serif)
    y_serif = rng.uniform(0.8 - stem_height - serif_height, 0.8 - stem_height, n_serif)
    serif = np.stack([x_serif, y_serif], axis=1)

    # Combine and label
    T = np.vstack([top, stem, serif])
    labels = np.array(
        [0] * n_top +
        [1] * n_stem +
        [2] * n_serif
    )

    # Add Gaussian noise if specified
    if noise_std > 0.0:
        noise = rng.normal(0, noise_std, T.shape)
        T += noise

    return T, labels

# Generate Datahubs with varying specs: 

def generate_clusters_with_bridges_data_v2(
    hubs_config,
    connect=1.0,
    noise=0.1,
    noise_c=0.1,
    bridge_points=20,
    random_seed=42
):
    """
    Generiert Cluster mit optionalen Verbindungen (Brücken) dazwischen.
    
    Parameter:
    - hubs_config: Liste von Tuples ([x, y, z], n_points, variance)
    - connect: float zwischen 0.0 und 1.0, wie viele Verbindungen erzeugt werden
    - noise: Rauschen für die Clusterpunkte
    - noise_c: Rauschen für die Verbindungsbrücken
    - bridge_points: Anzahl Punkte pro Verbindung
    - random_seed: Seed für Reproduzierbarkeit
    
    Rückgabe:
    - data_points: np.array mit allen Punkten (N, dim)
    - colors: Liste mit Hex-Farben (Länge N)
    """
    np.random.seed(random_seed)

    # Extrahiere Cluster-Zentren, Punktanzahl und Varianz
    hub_centers = np.array([np.array(cfg[0]) for cfg in hubs_config])
    points_per_hub = [cfg[1] for cfg in hubs_config]
    variances = [cfg[2] for cfg in hubs_config]

    n_hubs, dim = hub_centers.shape

    # Clusterpunkte generieren
    hubs = []
    for i in range(n_hubs):
        pts = hub_centers[i] + variances[i] * np.random.randn(points_per_hub[i], dim)
        pts += noise * np.random.randn(*pts.shape)
        hubs.append(pts)

    # Brückenpaare auswählen
    possible_pairs = list(combinations(range(n_hubs), 2))
    num_connections = int(connect * len(possible_pairs))
    if num_connections > 0:
        chosen_indices = np.random.choice(len(possible_pairs), size=num_connections, replace=False)
        pairs = [possible_pairs[i] for i in chosen_indices]
    else:
        pairs = []

    # Brückenpunkte generieren
    bridges = []
    for i, j in pairs:
        start, end = hub_centers[i], hub_centers[j]
        bridge_pts = np.array([
            start + t * (end - start) + noise_c * np.random.randn(dim)
            for t in np.linspace(0, 1, bridge_points)
        ])
        bridges.append(bridge_pts)

    # Farben generieren
    cmap = plt.cm.get_cmap('tab10', n_hubs)
    def rgba_to_hex(rgba):
        return '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))
    hub_colors = [rgba_to_hex(cmap(i)) for i in range(n_hubs)]
    bridge_color = '#000000'

    # Punkte und Farben sammeln
    all_points = []
    all_colors = []

    for i, pts in enumerate(hubs):
        all_points.append(pts)
        all_colors.extend([hub_colors[i]] * pts.shape[0])

    for pts in bridges:
        all_points.append(pts)
        all_colors.extend([bridge_color] * pts.shape[0])

    data_points = np.vstack(all_points)
    return data_points, all_colors


##############################################################################################################
# Projecting into high D. space + Adding noise and distortions: 
##############################################################################################################

def project_to_high_dim_space(data, output_dim=300, noise_range=0.0, random_seed=42):
    """
    Projects input data into a higher-dimensional space with a random bijective linear transformation
    and optional noise scattering.

    Parameters:
    - data: np.ndarray of shape (n_samples, n_features)
    - output_dim: int, target dimensionality (must be >= original dimension)
    - noise_range: float, max distance of noise scatter (uniformly sampled in [0, noise_range])
    - random_seed: int, for reproducibility

    Returns:
    - projected_data: np.ndarray of shape (n_samples, output_dim)
    """
    np.random.seed(random_seed)
    n_samples, input_dim = data.shape
    
    if output_dim < input_dim:
        raise ValueError("Output dimension must be greater than or equal to input dimension.")
    
    # Create a full-rank random projection matrix (bijective linear transform)
    A = np.random.randn(output_dim, input_dim)
    while np.linalg.matrix_rank(A) < input_dim:
        A = np.random.randn(output_dim, input_dim)

    # Linear transformation
    projected = data @ A.T  # shape: (n_samples, output_dim)
    
    # Add noise
    if noise_range > 0:
        directions = np.random.randn(n_samples, output_dim)
        norms = np.linalg.norm(directions, axis=1, keepdims=True)
        directions /= norms  # Normalize to unit vectors
        distances = np.random.uniform(0, noise_range, size=(n_samples, 1))
        noise = directions * distances
        projected += noise

    return projected

# Note, the noise scatter term becomes increasingly negligable in high dimensions
# due to the curse of dimensionality...