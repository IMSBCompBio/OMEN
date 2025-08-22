# All functions required for integration: 
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np

def pairwise_distances(x):
    # x: (n, d)
    x_norm = (x**2).sum(dim=1).view(-1,1)
    dist = x_norm + x_norm.t() - 2.0 * x @ x.t()
    return torch.clamp(dist, min=0)

def knn_graph(x, k=10):
    # Returns adjacency (n,n) with exponential decay weights on kNN
    with torch.no_grad():
        dist = pairwise_distances(x)
        knn_idx = dist.topk(k+1, largest=False)[1][:,1:]  # exclude self
        n = x.size(0)
        W = torch.zeros(n, n, device=x.device)
        sigma = dist.median().item()**0.5 + 1e-8  # bandwidth heuristic
        for i in range(n):
            neighbors = knn_idx[i]
            for j in neighbors:
                W[i,j] = torch.exp(-dist[i,j]/(2*sigma**2))
        # row normalize
        W = W / W.sum(dim=1, keepdim=True)
    return W

# KNN-aware Cauchy kernel
def cauchy_kernel(x, gamma=1.0,neighborindex = None):
    # x: (n,d)
    dist = pairwise_distances(x)
    K = 1 / (1 + dist / gamma)
    K[neighborindex == 0] = 0 # set non-neighbors to 0
    K = K / K.sum(dim=1, keepdim=True)  # row normalize to diffusion aka sum=1
    return K

def kl_divergence(P, Q):
    # P,Q: (n,n) row stochastic transition matrices
    # D_KL(P || Q) = sum_i sum_j P_ij * log(P_ij / Q_ij)
    # Add eps for numerical stability
    eps = 1e-8
    P = torch.clamp(P, eps, 1)
    Q = torch.clamp(Q, eps, 1)
    return (P * (P.log() - Q.log())).sum(dim=1).mean()

def compute_mmd(x, y, sigma=1.0):
    # RBF kernel MMD, can swap to Cauchy kernel if desired
    xx = torch.cdist(x, x, p=2)
    yy = torch.cdist(y, y, p=2)
    xy = torch.cdist(x, y, p=2)
    Kxx = torch.exp(-xx**2/(2*sigma**2)).mean()
    Kyy = torch.exp(-yy**2/(2*sigma**2)).mean()
    Kxy = torch.exp(-xy**2/(2*sigma**2)).mean()
    return Kxx + Kyy - 2*Kxy

class EmbedNet(nn.Module):
    def __init__(self, input_dim, embed_dim=2):
        super().__init__()  
        self.net = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(),
            nn.Linear(64, embed_dim)
        )
    def forward(self, x):
        return self.net(x)
    
def train_joint_embedding(X, Y, epochs=1000, lr=1e-3, k=20, embed_dim=2, device='cuda', verbose=False, alpha=0.5):

    # Convert to pytorch-tensors datastructure for input matrix R^(n x f): 
    device = 'cuda' if torch.cuda.device_count() > 0 else 'cpu'
    X = torch.tensor(X, dtype=torch.float32, device=device)
    Y = torch.tensor(Y, dtype=torch.float32, device=device)

    # Construction of MLP, maps R^n -> R^2 for all feature vectors in X,Y 
    # Use DataParallel for multi-GPU support
    
    f_net = nn.DataParallel(EmbedNet(X.shape[1], embed_dim).to(device))
    g_net = nn.DataParallel(EmbedNet(Y.shape[1], embed_dim).to(device))
    optimizer = torch.optim.Adam(list(f_net.parameters()) + list(g_net.parameters()), lr=lr)

    # Adapt K for non-matching spot dimensions to define similar KNN radius: 
    k_P = round(np.sqrt(X.numel() / Y.numel()) * k)
    k_Q = round(np.sqrt(Y.numel() / X.numel()) * k)

    # Precompute diffusion transition matrices P and Q on X and Y
    P = knn_graph(X, k=k_P)
    Q = knn_graph(Y, k=k_Q)

    # Training iteration: 
    for epoch in range(epochs):
        optimizer.zero_grad()
        # Put Data matrices through network: 
        fX = f_net(X)
        gY = g_net(Y)
        # Diffusion on embeddings
        T_fx = cauchy_kernel(fX)
        T_gy = cauchy_kernel(gY)

        # Loss functions: 
        kl_loss = kl_divergence(P, T_fx) + kl_divergence(Q, T_gy)
        mmd_loss = compute_mmd(fX, gY, sigma=1.0)

        # Scale losses
        loss = alpha * kl_loss + (1-alpha) * mmd_loss

        # Use this value as backpropagation loss function: 
        loss.backward()
        optimizer.step()
        if (epoch+1) % 100 == 0 and verbose == True:
            print(f"Epoch {epoch+1}, Loss: {loss.item():.4f}, KL: {kl_loss.item()*alpha :.4f}, MMD: {mmd_loss.item()* (1-alpha):.4f}")
    return f_net, g_net, loss
    
# More efficient loop for the training process

def proceed_train_joint_embedding(X, Y, f_net, g_net, epochs=2000, lr=1e-3, k=20, embed_dim=2, device='cuda', verbose=False, alpha=0.5):

    # Convert to pytorch-tensors datastructure for input matrix R^(n x f): 
    X = torch.tensor(X, dtype=torch.float32, device=device)
    Y = torch.tensor(Y, dtype=torch.float32, device=device)


    optimizer = torch.optim.Adam(list(f_net.parameters()) + list(g_net.parameters()), lr=lr)

    # Adapt K for non-matching spot dimensions to define similar KNN radius: 
    k_P = round(np.sqrt(X.numel() / Y.numel()) * k)
    k_Q = round(np.sqrt(Y.numel() / X.numel()) * k)
    #Previous way: 
    #k_P = round(X.numel() / (X.numel() + Y.numel()) * k)
    #k_Q = round(Y.numel() / (X.numel() + Y.numel()) * k)

    # Precompute diffusion transition matrices P and Q on X and Y
    P = knn_graph(X, k=k_P)
    Q = knn_graph(Y, k=k_Q)

    # Wrap models for multi-GPU
    # Check if CUDA is available and set device accordingly
    device = 'cuda' if torch.cuda.device_count() > 0 else 'cpu'
    f_net = nn.DataParallel(f_net).to(device)
    g_net = nn.DataParallel(g_net).to(device)

    # Training iteration: 
    for epoch in range(epochs):
        optimizer.zero_grad()
        # Put Data matrices through network: 
        fX = f_net(X)
        gY = g_net(Y)
        # Diffusion on embeddings
        T_fx = cauchy_kernel(fX)
        T_gy = cauchy_kernel(gY)

        # Loss functions: 
        kl_loss = kl_divergence(P, T_fx) + kl_divergence(Q, T_gy)
        mmd_loss = compute_mmd(fX, gY, sigma=1.0)

        # Scale losses
        loss = alpha * kl_loss + (1-alpha) * mmd_loss

        # Use this value as backpropagation loss function: 
        loss.backward()
        optimizer.step()
        if (epoch+1) % 100 == 0 and verbose == True:
            print(f"Epoch {epoch+1}, Loss: {loss.item():.4f}, KL: {kl_loss.item()*alpha :.4f}, MMD: {mmd_loss.item()* (1-alpha):.4f}")
    return f_net, g_net, loss