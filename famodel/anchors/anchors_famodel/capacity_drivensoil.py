
import numpy as np
import matplotlib.pyplot as plt

###################################
#### Pile Geometry and Loading ####
###################################

def getCapacityDrivenSoil(profile, soil_type, L, D, zlug, V, H, t, fy):
    
    '''Models a laterally loaded pile using the p-y method. The solution for
    lateral displacements is obtained by solving the 4th order ODE, EI*d4y/dz4
    - F*d2y/dz2 + ky = 0 using the finite difference method.

    Assumes that EI remains constant with respect to curvature i.e. pile
    material remains in the elastic region.

    Input:
    -----
    profile     - A 2D array of depths (m) and corresponding undrained shear strength(Pa)
                  Eg: array([[z1,Su1],[z2,Su2],[z3,Su3]...])
                  Use small values for Su (eg: 0.001) instead of zeros to avoid divisions by zero but always start z at 0.0
                  Example of a valid data point at the mudline is [0.0, 0.001]
    soil_type   - Select soil condition, 'clay' or 'sand'
                  Assigns which p-y model to use, 'Matlock' or 'API'.
    L           - Length of pile (m)
    D           - Outer diameter of pile (m)
    zlug        - Load eccentricity above the mudline or depth to mudline relative to the pile head (m)
    V           - Axial force at pile head (N), vertically downwards is postive.
    H           - Force at pile head (N), shear causing clockwise rotation of pile is positive.
    M           - Moment at pile head (N*m), moments causing tension on left side of pile is positive.
    n           - Number of elements (50 by default)
    iterations  - Number of iterations to repeat calculation in order obtain convergence of 'y'
                  (A better approach is to iterate until a predefined tolerance is achieved but this requires additional
                  coding so, I will implement this later.)

    Output:
    ------
    y           - Lateral displacement at each node, length = n + 5, (n+1) real nodes and 4 imaginary nodes
    z           - Vector of node locations along pile
    '''

    # Extract optional keyword arguments

    ls = 'x'
    n = 25; iterations = 10; loc=2
    
    # Convert L and D to floating point numbers to avoid rounding errors
    L = float(L)
    D = float(D)
    # t = (6.35 + D*20)/1e3       # Pile wall thickness (m), API RP2A-WSD
    E = 200e9                   # Elastic modulus of pile material (Pa)
    fy = 350e6                  # Yield strength of pile material (Pa)
    
    # Pile geometry
    I = (np.pi/64.0)*(D**4 - (D - 2*t)**4)
    EI = E*I
    h  = L/n                    # Element size
    N  = (n + 1) + 4            # (n+1) Real + 4 Imaginary nodes

    # Array for displacements at nodes, including imaginary nodes.
    y = np.ones(N)*(0.01*D)     # An initial value of 0.01D was arbitrarily chosen

    # Initialize and assemble array/list of p-y curves at each real node
    z = np.zeros(N)
    py_funs  = []
    k_secant = np.zeros(N)

    for i in [0, 1]:            # Top two imaginary nodes
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0

    # Extract soil profile data
    if soil_type == 'clay':
        z0, f_Su, f_sigma_v_eff, f_gamma = clay_profile(profile)

    elif soil_type == 'sand':
        z0, f_phi, f_sigma_v_eff = sand_profile(profile)

    for i in range(2, n+3):    # Real nodes
        z[i] = (i - 2)*h      
        if soil_type == 'clay':
            Su, sigma_v_eff, gamma = f_Su(z[i]), f_sigma_v_eff(z[i]), f_gamma(z[i])
            py_funs.append(py_Matlock(z[i], D, zlug, Su, sigma_v_eff, gamma))
        elif soil_type == 'sand':
            phi, sigma_v_eff = f_phi(z[i]), f_sigma_v_eff(z[i])
            py_funs.append(py_API(z[i], D, zlug, phi, sigma_v_eff))
        
        k_secant[i] = py_funs[i](y[i])/y[i]

    for i in [n+3, n+4]:        # Bottom two imaginary nodes
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0


    y1 = np.linspace(-2.*D, 2.*D, 500)
    plt.plot(y1, py_funs[loc](y1))
    plt.xlabel('y (m)'), plt.ylabel('p (N/m)')
    plt.grid(True)

    for j in range(iterations):
        # if j == 0: print 'FD Solver started!'
        y, hinge_formed, hinge_location = fd_solver(n, N, h, D, t, fy, EI, V, H, zlug, k_secant)

        for i in range(2, n+3):
            k_secant[i] = py_funs[i](y[i])/y[i]

    print(f'y_max = {max(y):.3f} m')
    
    resultsDrivenSoil = {}
    resultsDrivenSoil['Lateral displacement'] = y[2]
    resultsDrivenSoil['Rotational displacement'] = np.rad2deg((y[2] - y[3])/h) 
    
    return y[2:-2], z[2:-2], hinge_formed, hinge_location

#################
#### Solvers ####
#################

def fd_solver(n, N, h, D, t, fy, EI, V, H, zlug, k_secant):
    '''Solves the finite difference equations from 'py_analysis_1'. This function should be run iteratively for
    non-linear p-y curves by updating 'k_secant' using 'y'. A single iteration is sufficient if the p-y curves
    are linear.

    Input:
    -----
    n         - Number of elements
    N         - Total number of nodes
    h         - Element size
    EI        - Flexural rigidity of pile
    V         - Axial force at pile head/zlug depth
    H         - Shear at pile head/zlug depth
    M         - Moment at pile head/zlug depth
    zlug      - Load eccentricity above the mudline or depth to mudline relative to the pile head (m)
    k_secant  - Secant stiffness from p-y curves

    Output:
    ------
    y_updated - Lateral displacement at each node
    '''
    
    from scipy import linalg
    
    # Identify the node corresponding to zlug
    zlug_index = int(zlug/h)  # Index for the node corresponding to zlug
    print(zlug)
    print(zlug_index)

    # Initialize and assemble matrix
    X = np.zeros((N, N))

    # (n+1) finite difference equations for (n+1) real nodes
    for i in range(0, n+1):
        X[i,i]   =  1.0
        X[i,i+1] = -4.0 + V*h**2/EI
        X[i,i+2] =  6.0 - 2*V*h**2/EI + k_secant[i+2]*h**4/EI
        X[i,i+3] = -4.0 + V*h**2/EI
        X[i,i+4] =  1.0

    # Curvature at pile head
    X[n+1,1]   =  1.0
    X[n+1,2]   = -2.0
    X[n+1,3]   =  1.0

    # Shear at pile head
    X[n+2,0]   = -1.0
    X[n+2,1]   =  2.0 - V*h**2/EI
    X[n+2,2]   =  0.0
    X[n+2,3]   = -2.0 + V*h**2/EI
    X[n+2,4]   =  1.0

    # Curvature at pile tip
    X[n+3,-2]   =  1.0
    X[n+3,-3]   = -2.0
    X[n+3,-4]   =  1.0

    # Shear at pile tip
    X[n+4,-1]   =   1.0
    X[n+4,-2]   =  -2.0 + V*h**2/EI
    X[n+4,-3]   =   0.0
    X[n+4,-4]   =   2.0 - V*h**2/EI
    X[n+4,-5]   =  -1.0

    # Initialize vector q
    q = np.zeros(N)
    M = H*abs(zlug)

    # Populate q with boundary conditions
    if zlug <= 0:
        q[-3] = 2*H*h**3      # Shear at pile head
        q[-4] = M*h**2        # Moment at pile head
    else:
        q[zlug_index]      = 2*H*h**3     # Shear at pile head
        #q[zlug_index + 1]  = M*h**2      # Moment at pile head
               
    y = linalg.solve(EI*X, q)

    # Compute the plastic moment capacity Mp
    Zp = (1/6)*(D**3 - (D - 2*t)**3)  # Plastic section modulus for hollow pile
    Mp = Zp*fy                        # Plastic moment capacity

    # Check for plastic hinge formation
    hinge_formed, hinge_location = plastic_hinge(y, h, EI, Mp)
    
    return y, hinge_formed, hinge_location

###############################
#### P-Y Curve Definitions ####
###############################

def py_Matlock(z, D, zlug, Su, sigma_v_eff, gamma):
    
    '''Returns an interp1d interpolation function which represents the Matlock (1970) p-y curve at the depth of interest.
    Important: Make sure to import the interp1 function by running 'from scipy.interpolate import interp1d' in the main program.

    Input:
    -----
    z            - Depth relative to pile head (m)
    D            - Pile diameter (m)
    zlug         - Load eccentricity above the mudline or depth to mudline relative to the pile head (m)
    Su           - Undrained shear strength (Pa)
    sigma_v_eff  - Effective vertical stress (Pa)
    gamma        - Effective unit weight of the soil (kN/m3)

    Output:
    ------
    Returns an interp1d interpolation function which represents the p-y curve at the depth of interest.
    'p' (N/m) and 'y' (m).
    '''

    from scipy.interpolate import interp1d
    
    z0 = 0
    
    # Strain at half the strength as defined by Matlock (1970). 
    # Typically ranges from 0.005 (stiff clay) to 0.02 (soft clay).
    epsilon_50 = 0.02          
    # p-y curve properties
    J = 0.5
    
    if zlug < 0:
        # Scenario 1: zlug is negative (above mudline)
        if (z - z0) < 0:
            # No p-y curve between z = 0 and zlug
            Nc = 0.0
            z_cr = 1.0  # Dummy value to avoid crashing
        else:
            # Calculate p-y curve below zlug
            Nc = 3.0 + sigma_v_eff/Su + J*(z - abs(zlug))/D
            if Nc > 9.0: 
                Nc = 9.0
            z_cr = 6.0*D/(gamma*D/Su + J)

    else:
        # Scenario 2: zlug is positive (below mudline)
        # Calculate p-y curve for the entire pile (all depths)
        Nc = 3.0 + sigma_v_eff/Su + J*(z - zlug)/D
        if Nc > 9.0: 
            Nc = 9.0
        z_cr = 6.0 * D/(gamma*D/Su + J)

    p_ult = Su*Nc*D
    y_50 = 2.5*epsilon_50*D

    # Normalized lateral displacement
    Y = np.concatenate((-np.logspace(3,-4,100),[0],np.logspace(-4,3,100)))

    # Normalized p-y curves
    P = 0.5*np.sign(Y)*abs(Y)**(1.0/3.0)  # sign(Y) and abs(Y) used since negative numbers cannot be raised to fractional powers
                                          # Expression equivalent to P = 0.5*Y**(1.0/3.0) for Y>=0
    for i in range(0,len(Y)):
        if P[i] > 1.0:    P[i] = 1.0
        elif P[i] < -1.0: P[i] = -1.0

    # Un-normallized p-y curves
    p = P*p_ult
    y = Y*y_50

    f = interp1d(y, p, kind='linear')   # Interpolation function for p-y curve

    # Plot of p-y curve and check if 'k' is calculated correctly
    plt.plot(y, p,'-')
    plt.xlabel('y (m)') 
    plt.ylabel('p (N/m)')
    plt.title('PY Curves - Matlock (1970)')
    plt.grid(True)
    plt.xlim([-2*D, 2*D])

    return f   # This is f (linear interpolation of y-p)
       
def py_API(z, D, zlug, phi, sigma_v_eff):
    
    '''Returns an interp1d interpolation function which represents the Matlock (1970) p-y curve at the depth of interest.

    Important: Make sure to import the interp1 function by running 'from scipy.interpolate import interp1d' in the main program.

    Input:
    -----
    z            - Depth relative to pile head (m)
    D            - Pile diameter (m)
    zlug         - Load eccentricity above the mudline or depth to mudline relative to the pile head (m)
    phi          - Internal friction angle (deg)
    sigma_v_eff  - Effectve vertical stress (Pa)

    Output:
    ------
    Returns an interp1d interpolation function which represents the p-y curve at the depth of interest.
    'p' (N/m) and 'y' (m).
    '''
    
    from scipy.interpolate import interp1d
       
    # Interpolate coefficients depending on the effective friction angle
    phi_ref = [ 20,   25,   30,   35,   40]
    C1_ref = [0.80, 1.25, 1.90, 3.00, 4.50]
    C2_ref = [1.60, 2.10, 2.60, 3.40, 4.30]
    C3_ref = [  10,   15,   30,   55,  105]
    
    C1 = np.interp(phi, phi_ref, C1_ref)
    C2 = np.interp(phi, phi_ref, C2_ref)
    C3 = np.interp(phi, phi_ref, C3_ref)
    
    if (z - zlug) < 0:
        # p-y curves for the virtual soil layer between the pile head and the mudline should have p=0
        p_ult = 0.0
    else:
        try:
            p_ult = min(C1*z + C2*D, C3*D)*sigma_v_eff
        except ZeroDivisionError:
            print("Division by zero! phi = 0.0 so z_cr cannot be calculated.")
    
    Dr = 0.75     # Relative density of the soil (assumed)
    k = 54.6*Dr**2 + 0.8*Dr + 1.8
    
    # Normalized lateral displacement
    N = 20
    y = np.concatenate((-np.logspace(3,-4,N),[0],np.logspace(-4,3,N)))
    A = 3 - 0.8*z/D
    epsilon = 1e-8
    p = A*p_ult*np.tanh(k*z*y/(A*p_ult + epsilon))

    f = interp1d(y, p, kind='linear')   # Interpolation function for p-y curve
    
    # Plot of p-y curve and check if 'k' is calculated correctly
    plt.plot(y, p,'-') 
    plt.xlabel('y (m)') 
    plt.ylabel('p (N/m)')
    plt.title('PY Curves - API (1993)')
    plt.grid(True)
    plt.xlim([-2*D, 2*D])
    plt.ylim([min(y), max(y)])  # Adjust x-axis limits to match y values
        
    return f  # This is f (linear interpolation of y-p)

########################
#### Plastic Hinge #####
########################

def plastic_hinge(y, h, EI, Mp):
    '''
    Check for plastic hinge formation along the pile.
    
    Parameters:
    ----------
    y : ndarray
        Lateral displacements at each node.
    h : float
        Element size (distance between nodes).
    EI : float
        Flexural rigidity of the pile (N*m²).
    Mp : float
        Plastic moment capacity of the pile section.
    
    Returns:
    -------
    hinge_formed : bool
        True if a plastic hinge forms, False otherwise.
    hinge_location : int
        Index of the node where the plastic hinge forms (if any).
    '''
    
    hinge_formed = False
    hinge_location = -1
    
    # Loop through each internal node and compute the bending moment
    for i in range(1, len(y) - 1):
        # Approximate the bending moment at node i
        Mi = EI*(y[i+1] - 2*y[i] + y[i-1])/h**2
        
        # Check if the moment exceeds the plastic moment capacity
        if Mi >= Mp:
            hinge_formed = True
            hinge_location = i
            break  # Stop at the first plastic hinge formation
    
    return hinge_formed, hinge_location


########################
#### Soil Profiles #####
########################

def clay_profile(profile):
    '''Define the clay profile used by the p-y analyzer. Outputs 'interp1d' functions containing Su and sigma'_v
    profiles to be used by the p-y curve functions.

    Input:
    -----
    profile      - A 2D tuple in the following format: ([Depth (m), Su (kPa), gamma (kN/m^3), py-model, model parameter])
                   The soil profile should be defined relative to the pile/tower head (i.e. point of lateral load application)
                   so that any load eccentricities can be taken into account. An example soil profile is shown below.
                   Eg: array([[z0,Su0,gamma0,  'Matlock', 0.02],
                              ...])

                    *The current program cannot define layers with different p-y models. But it will added in the future.

    plot_profile - Plot Su vs depth profile. Choose 'Yes' to plot.

    Output:
    ------
    z0            - Depth of mudline relative to the pile head (m)
    f_Su          - 'interp1d' function containing undrained shear strength profile (Pa)
    f_sigma_v_eff - 'interp1d' function containing effective vertical stress profile (Pa)
    f_gamma       - 'interp1d' function containing effective unit weight (kN/m3)
    '''

    from scipy.interpolate import interp1d

    # Depth of mudline relative to pile head
    z0 = profile[0,0].astype(float)

    # Extract data from profile array and zero strength virtual soil layer
    # from the pile head down to the mudline
    depth = np.concatenate([np.array([0,z0]), profile[:,0].astype(float)])  # m
    Su    = np.concatenate([np.array([0, 0]), profile[:,1].astype(float)])  # kPa
    gamma = np.concatenate([np.array([0, 0]), profile[:,2].astype(float)])  # kN/m3

    # Calculate sigma_v_eff at each depth
    sigma_v_eff = np.zeros(len(depth))

    for i in range(1, len(depth)):
        sigma_v_eff[i] = sigma_v_eff[i-1] + gamma[i-1]*(depth[i] - depth[i-1])

    # Define interpolation functions
    f_Su          = interp1d(depth, Su*1000, kind='linear')          # Pa
    f_sigma_v_eff = interp1d(depth, sigma_v_eff*1000, kind='linear') # Pa
    f_gamma       = interp1d(depth, gamma, kind='linear')            # kN/m3

    return z0, f_Su, f_sigma_v_eff, f_gamma

def sand_profile(profile):
    '''Define the sand profile used by the p-y analyzer. Outputs 'interp1d' functions containing Su and sigma'_v
    profiles to be used by the p-y curve functions.

    Input:
    -----
    profile      - A 2D tuple in the following format: ([Depth (m), Su (kPa), gamma (kN/m^3), py-model, model parameter])
                   The soil profile should be defined relative to the pile/tower head (i.e. point of lateral load application)
                   so that any load eccentricities can be taken into account. An example soil profile is shown below.
                   Eg: array([[z0,Su0,gamma0,  'API', 0.02],
                              ...])

                    *The current program cannot define layers with different p-y models. But it will added in the future.

    plot_profile - Plot Su vs depth profile. Choose 'Yes' to plot.

    Output:
    ------
    z0            - Depth of mudline relative to the pile head (m)
    f_phi         - 'interp1d' function containing effective friction angle (deg)
    f_sigma_v_eff - 'interp1d' function containing effective vertical stress profile (Pa)
    '''

    from scipy.interpolate import interp1d

    # Depth of mudline relative to pile head
    z0 = profile[0,0].astype(float)

    # Extract data from profile array and zero strength virtual soil layer
    # from the pile head down to the mudline
    depth = np.concatenate([np.array([0,z0]), profile[:,0].astype(float)])  # m
    phi   = np.concatenate([np.array([0, 0]), profile[:,1].astype(float)])  # deg
    gamma = np.concatenate([np.array([0, 0]), profile[:,2].astype(float)])  # kN/m^3
   
    # Calculate sigma_v_eff and static loading factor at each depth
    sigma_v_eff = np.zeros(len(depth))

    for i in range(1, len(depth)):  
        sigma_v_eff[i] = sigma_v_eff[i-1] + gamma[i-1]*(depth[i] - depth[i-1])

    # Define interpolation functions
    f_phi         = interp1d(depth, phi, kind='linear')              # deg
    f_sigma_v_eff = interp1d(depth, sigma_v_eff*1000, kind='linear') # Pa
    
    return z0, f_phi, f_sigma_v_eff

if __name__ == '__main__':
    
    #                   Depth   Su  γ_sub    p-y model   p-y parameters
    profile = np.array([[0.0,  250.0, 20., 'Name of p-y model', 0.02],
                        [75.0, 250.0, 20., 'Name of p-y model', 0.02]])

    z0, f_Su, f_σ_v_eff, f_gamma_sub = clay_profile(profile)
    
    # #                 Depth  phi  γ_sub    p-y model   p-y parameters
    #profile = np.array([[0.0,  38.0, 20., 'Name of p-y model', 0.02],
    #                    [75.0, 40.0, 22., 'Name of p-y model', 0.02]])

    #zlug, f_phi, f_σ_v_eff = sand_profile(profile)
    
    #Pile dimensions
    L = 25                      # Pile length (m)
    D = 1.5                     # Pile diameter (m)
    t = (6.35 + D*20)/1e3       # Pile wall thickness (m), API RP2A-WSD
    fy = 355e6
    zlug = 6*D                  # Lug depth (m)

    
    #Pile head loads
    H = 280e4    # Horizontal load on pile head (N)
    V = 0       # Vertical load on pile head (N)
    
    y, z, hinge_formed, hinge_location = getCapacityDrivenSoil(profile, soil_type='clay', L=L, D=D, zlug=zlug, V=V, H=H, t=t, fy=fy)
    
    y0 = np.zeros(len(z))
    #Plot deflection profile of pile
    fig, ax = plt.subplots(figsize=(3,5))    
    ax.plot(y0,z,'black')
    ax.plot(y,z,'r')
    ax.set_xlabel('Displacement [m]')
    ax.set_ylabel('Depth below pile head [m]')
    ax.set_ylim([L + 2,-2])
    ax.set_xlim([-0.1*D, 0.1*D])
    ax.grid(ls='--')