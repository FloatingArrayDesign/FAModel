
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import linalg

def fd_solver(n, N, h, D, t, fy, EI, Ha, Va, zlug, z0, k_secant):
    '''Solves the finite difference equations from 'py_analysis_1'. This function should be run iteratively for
    non-linear p-y curves by updating 'k_secant' using 'y'. A single iteration is sufficient if the p-y curves
    are linear.

    Parameters
    ----------
    n : int
        Number of elements
    N : int
        Total number of nodes (real + imaginary)
    h : float
        Element size (m)
    D : float
        Pile diameter (m)
    t : float
        Pile wall thickness (m)
    fy : float
        Yield strength of pile material (Pa)
    EI : float
        Flexural rigidity of the pile (Nm²)
    Ha : float       
        Horizontal load at pile lug elevation (N)
    Va : float          
        Vertical load at pile lug elevation (N)
    zlug : float
        Depth of padeye from pile head (m)
    z0 : float
        Mudline depth from pile head (m)
    k_secant : array
        Secant stiffness from p-y curves at each node

    Returns
    -------
    y : array
        Lateral displacement at each node
    Mi : float
        Maximum internal bending moment (Nm)
    Mp : float
        Plastic moment capacity of the pile (Nm)
    hinge_formed : bool
        Whether plastic hinge is formed
    hinge_location : int
        Index of the node with hinge formation
    '''  
    
    # Initialize and assemble matrix
    X = np.zeros((N, N))

    # (n+1) finite difference equations for (n+1) real nodes
    for i in range(0, n+1):
        X[i, i]   =  1.0
        X[i, i+1] = -4.0 + Va*h**2/EI
        X[i, i+2] =  6.0 - 2*Va*h**2/EI + k_secant[i+2]*h**4/EI
        X[i, i+3] = -4.0 + Va*h**2/EI
        X[i, i+4] =  1.0

    # Curvature at pile head
    X[n+1, 1] =  1.0
    X[n+1, 2] = -2.0
    X[n+1, 3] =  1.0

    # Shear at pile head
    X[n+2, 0] = -1.0
    X[n+2, 1] =  2.0 - Va*h**2/EI
    X[n+2, 2] =  0.0
    X[n+2, 3] = -2.0 + Va*h**2/EI
    X[n+2, 4] =  1.0

    # Curvature at pile tip
    X[n+3, -2] =  1.0
    X[n+3, -3] = -2.0
    X[n+3, -4] =  1.0

    # Shear at pile tip
    X[n+4, -1] =  1.0
    X[n+4, -2] = -2.0 + Va*h**2/EI
    X[n+4, -3] =  0.0
    X[n+4, -4] =  2.0 - Va*h**2/EI
    X[n+4, -5] = -1.0

    # Initialize vector q
    q = np.zeros(N)
    
    # Always apply shear
    # Index of the node where the horizontal load is applied (padeye)
    zlug_index = int(zlug/h)
    q[zlug_index] = 2*Ha*h**3
                 
    y = linalg.solve(EI*X, q)

    # Compute the plastic moment capacity Mp
    Zp = (1/6)*(D**3 - (D - 2*t)**3)  # Plastic section modulus for hollow pile (m3)
    Mp = Zp*fy                        # Plastic moment capacity (N/m)

    # Check for plastic hinge formation
    Mi, Mp, hinge_formed, hinge_location = plastic_hinge(y, h, EI, Mp)
    
    return y, Mi, Mp, hinge_formed, hinge_location
   
def plastic_hinge(y, h, EI, Mp):
    '''Check for plastic hinge formation along the pile.
    
    Parameters
    ----------
    y : array
        Lateral displacements at each node
    h : float
        Element size (m)
    EI : float
        Flexural rigidity of the pile (Nm²)
    Mp : float
        Plastic moment capacity (Nm)

    Returns
    -------
    Mi_max : float
        Maximum internal moment along the pile (Nm)
    Mp : float
        Plastic moment capacity (Nm)
    hinge_formed : bool
        True if plastic hinge is formed
    hinge_location : int
        Node index where hinge forms (if any)
    '''
    
    hinge_formed = False
    hinge_location = -1
    Mi_all = []
        
    # Loop through each internal node and compute the bending moment
    for i in range(1, len(y) - 1):
        Mi = EI * (y[i+1] - 2*y[i] + y[i-1])/h**2
        Mi_all.append(Mi)

        if abs(Mi) >= Mp and not hinge_formed:
            hinge_formed = True
            hinge_location = i

    Mi_max = max(Mi_all, key=abs) if Mi_all else 0.0

    return Mi_max, Mp, hinge_formed, hinge_location