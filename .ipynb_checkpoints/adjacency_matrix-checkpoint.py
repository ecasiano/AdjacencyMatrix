# Build the adjacency matrix of a 1D,2D,3D
# Bose-Hubbard Lattice

import numpy as np
import matplotlib.pyplot as plt

def build_adjacency_matrix(L,D,boundary_condition='pbc'):
    
    # Number of lattice sites
    M = L**D
    
    # Define the basis vectors
    a1_vec = np.array((1,0,0))
    a2_vec = np.array((0,1,0))
    a3_vec = np.array((0,0,1))
    
    # Norms of basis vectors
    a1 = np.linalg.norm(a1_vec)
    a2 = np.linalg.norm(a2_vec)
    a3 = np.linalg.norm(a3_vec)
   
    # Initialize array that will store the lattice vectors
    points = np.zeros(M,dtype=(float,D))
    
    # Build the lattice vectors
    ctr = 0 # iteration counter
    if D == 1:
        for i1 in range(L):
            points[i1] = i1*a1
    elif D == 2:
        for i1 in range(L):
            for i2 in range(L):
                points[ctr] = np.array((i1*a1,i2*a2))
                ctr += 1
    else: # D == 3
        for i1 in range(L):
            for i2 in range(L):
                for i3 in range(L):
                    points[ctr] = np.array((i1*a1,i2*a2,i3*a3))
                    ctr += 1
                    
    # Initialize adjacency matrix
    A = np.zeros((M,M))
    
    # Set Nearest-Neighbor (NN) distance
    r_NN = a1
    
    # Build the adjacency matrix by comparing internode distances
    for i in range(M):
        for j in range(i+1,M):
            
            if boundary_condition=='pbc':
                A[i,j] = (np.linalg.norm(points[i]-points[j]) <= r_NN \
                or np.linalg.norm(points[i]-points[j]) == L-1)
            elif boundary_condition=='obc':
                A[i,j] = (np.linalg.norm(points[i]-points[j]) <= r_NN)
            
            # Symmetric elements
            A[j,i] = A[i,j]
    
    return A

# Main
L=3
D=2
A=build_adjacency_matrix(L,D,'obc')