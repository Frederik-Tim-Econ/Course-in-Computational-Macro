import numpy as np
import math

# interpolation functions:
def binary_search(gridsize,x,xi):
    # gridsize: number of gridpoints
    # x: grid
    # xi: value
    
    # output variable
    # imin: left index of value xi 
    
    # Initialize upper bound for search
    Nx = gridsize
    
    # a. Check if in grid
    if xi <= x[0]: # Return zero if xi is smaller than first gridpoint
        return 0
    elif xi >= x[Nx-2]: # Return penultimate gridpoint if xi is larger than last gridpoint
        return Nx-2
    else:
    # b. Binary search for index closest to xi | imin<xi
        half = Nx//2 # Floor division to find midpoint
        
        imin = 0
        
        while half!=0: # While the midpoint exists
            imid = imin + half
            if x[imid] <= xi:
                imin = imid
            Nx = Nx - half
            half = Nx//2
        
    return imin      


def interp_linear_1d_scalar(grid,value,xi):
    """ raw 1D interpolation """
    # Input
    # grid:     Grid to interpolate over
    # value:    Array of values defined over grid
    # xi:       Point for interpolation
    
    # Output
    # yi:        Interpolated value

    # a. search
    ix = binary_search(grid.size,grid,xi) # Find closest index to the left
    
    # b. relative positive
    rel_x = (xi - grid[ix])/(grid[ix+1] - grid[ix]) # Compute relative distance between left and right point
    
    # c. interpolate
    yi = value[ix] + rel_x * (value[ix+1] - value[ix])
    
    return yi


def interp_linear_1d(grid,value,X):
    # Input
    # grid:     Grid to interpolate over
    # value:    Array of values defined over grid
    # X:        Vector of points for interpolation
    
    # Output
    # Y:        Vector of interpolated values
    
    Y = np.empty(X.size)

    for i in range(X.size):
        Y[i] = interp_linear_1d_scalar(grid,value,X[i])
    
    return Y



def gauss_hermite(S):
    # Input
    # S: Number of Gauss-Hermite nodes
    
    # Output
    # x: Location of Gauss-Hermite nodes
    # w: Gauss-Hermite weights
    
    # a. calculations
    i = np.arange(1,S)
    a = np.sqrt(i/2)
    CM = np.diag(a,1) + np.diag(a,-1)
    L,V = np.linalg.eig(CM)
    I = L.argsort()
    V = V[:,I].T

    # b. nodes and weights
    x = L[I]
    w = np.sqrt(math.pi)*V[:,0]**2

    return x,w   




