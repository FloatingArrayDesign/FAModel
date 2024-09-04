
import numpy as np
import matplotlib.pyplot as plt

###################################
#### Pile Geometry and Loading ####
###################################

def getCapacityDrivenSoil(profile, soil_type, L, D, zlug, V, H):
    
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
                - Select which p-y model to use, 'Matlock' or 'API'.
    L           - Length of pile         (m)
    D           - Outer diameter of pile (m)
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
    # ls = 'x'
    n = 10, iterations = 10
    
    # Convert L and D to floating point numbers to avoid rounding errors
    L = float(L)
    D = float(D)
    t = (6.35 + D*20)/1e3       # Pile wall thickness (m), API RP2A-WSD
    E = 200e9                   # Elastic modulus of pile material (Pa)

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
        zlug, f_Su, f_sigma_v_eff, f_gamma_sub = clay_profile(profile)
    elif soil_type == 'sand':
        zlug, f_phi, f_sigma_v_eff = sand_profile(profile)

    for i in range(2, n+3):     # Real nodes
        z[i] = (i - 2)*h      
        if soil_type == 'clay':
            Su, Su0, sigma_v_eff, gamma_sub = f_Su(z[i]), f_Su(zlug + 0.01), f_sigma_v_eff(z[i]), f_gamma_sub(z[i])
            py_funs.append(py_Matlock(z[i], D, Su, sigma_v_eff, gamma_sub, zlug=zlug, epsilon_50=epsilon_50))
        elif soil_type == 'sand':
            phi, sigma_v_eff = f_phi(z[i]), f_sigma_v_eff(z[i])
            py_funs.append(py_API(z[i], D, phi, sigma_v_eff, zlug=zlug))
        
        k_secant[i] = py_funs[i](y[i])/y[i]

    for i in [n+3, n+4]:       # Bottom two imaginary nodes
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0

    # Track k_secant and current displacements
    if convergence_tracker == 'Yes':
        y1 = np.linspace(-2.*D, 2.*D,500)
        plt.plot(y1, py_funs[loc](y1))
        plt.xlabel('y (m)'), plt.ylabel('p (N/m)')
        plt.grid(True)

    for j in range(iterations):
        # if j == 0: print 'FD Solver started!'
        y = fd_solver(n, N, h, EI, V, H, zlug, k_secant)

        if convergence_tracker == 'Yes':
            plt.plot(y[loc], k_secant[loc]*y[loc], ls)

        for i in range(2, n+3):
            k_secant[i] = py_funs[i](y[i])/y[i]

    if print_output == 'Yes':
        print(f'y_0 = {y[2]:.3f} m')
    
    resultsDrivenSoil = {}
    resultsDrivenSoil['Lateral displacement'] = y[2]
    resultsDrivenSoil['Rotational displacement'] = (y[2] - y[3])/h 
    
    return resultsDrivenSoil

#################
#### Solvers ####
#################

def fd_solver(n, N, h, EI, getV, H, zlug, k_secant):
    '''Solves the finite difference equations from 'py_analysis_1'. This function should be run iteratively for
    non-linear p-y curves by updating 'k_secant' using 'y'. A single iteration is sufficient if the p-y curves
    are linear.

    Input:
    -----
    n        - Number of elements
    N        - Total number of nodes
    h        - Element size
    EI       - Flexural rigidity of pile
    V        - Axial force at pile head
    H        - Shear at pile head
    M        - Moment at pile head
    k_secant - Secant stiffness from p-y curves

    Output:
    ------
    y_updated - Lateral displacement at each node
    '''

    from scipy import linalg
    
    M = H*zlug

    # Initialize and assemble matrix
    X = np.zeros((N,N))

    # (n+1) finite difference equations for (n+1) real nodes
    for i in range(0,n+1):
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

    # Populate q with boundary conditions
    q[-3] = 2*H*h**3     # Shear at pile head
    q[-4] = M*h**2       # Moment at pile head

    y = linalg.solve(EI*X, q)

    return y

###############################
#### P-Y Curve Definitions ####
###############################

def py_Matlock(z, D, zlug, Su, sigma_v_eff, gamma_sub):
    
    '''Returns an interp1d interpolation function which represents the Matlock (1970) p-y curve at the depth of interest.
    Important: Make sure to import the interp1 function by running 'from scipy.interpolate import interp1d' in the main program.

    Input:
    -----
    z            - Depth relative to pile head (m)
    D            - Pile diameter (m)
    zlug          - Load eccentricity above the mudline or depth to mudline relative to the pile head (m)
    Su           - Undrained shear strength (Pa)
    sigma_v_eff  - Effective vertical stress (Pa)
    gamma_sub    - Effective unit weight of the soil (kN/m3)

    Output:
    ------
    Returns an interp1d interpolation function which represents the p-y curve at the depth of interest.
    'p' (N/m) and 'y' (m).
    '''

    from scipy.interpolate import interp1d
    epsilon_50 = 0.02          # Strain at half the strength as defined by Matlock (1970). Typically ranges from 0.005 (stiff clay) to 0.02 (soft clay).
    
    # p-y curve properties
    J = 0.5

    if (z - zlug) < 0:
        # p-y curves for the virtual soil layer between the pile head and the mudline should have p=0
        Nc = 0.0
        z_cr = 1.0 # Dummy value to keep program from crashing

    else:
        try:
            Nc = 3.0 + sigma_v_eff/Su + J*(z - zlug)/D
            if Nc > 9.0: Nc = 9.0
            z_cr = 6.0*D /(gamma_sub*D/Su + J)  
        except ZeroDivisionError:
            print("Division by zero! Su = 0.0 so z_cr cannot be calculated.")

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

    if print_curves == 'Yes':
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
    zlug          - Load eccentricity above the mudline or depth to mudline relative to the pile head (m)
    phi          - Internal friction angle (deg)
    sigma_v_eff  - Effectve vertical stress (Pa)
    epsilon_50   - Strain at half the strength as defined by Matlock (1970).
                   Typically ranges from 0.005 (stiff clay) to 0.02 (soft clay).

    Optional argument:
    return_Np    - Returns the Np values that in addtion to the p-y curve. This option was added to visualize
                   the gapping effect. It should only be used when this function is used by itself.
                   DO NOT set it to 'Yes' for p-y analysis as the program will crash!

    Output:
    ------
    Returns an interp1d interpolation function which represents the p-y curve at the depth of interest.
    'p' (N/m) and 'y' (m).
    '''
    
    from scipy.interpolate import interp1d
       
    # Interpolate coefficients depending on the effective friction angle
    phi_ref = [ 20,   25,   30,   35,   40]
    C1_ref = [0.80, 1.25, 1.90, 3.00, 4.50]
    C2_ref = [1.60, 2.00, 2.60, 3.40, 4.30]
    C3_ref = [  15,   15,   30,   55,  100]
    
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
            print("Division by zero! Su = 0.0 so z_cr cannot be calculated.")
    
    Dr = 0.75     # Relative density of the soil (assumed)
    k = 54.6*Dr**2 + 0.8*Dr + 1.8
    
    # Normalized lateral displacement
    N = 20
    y = np.concatenate((-np.logspace(3,-4,N),[0],np.logspace(-4,3,N)))
    A = 3 - 0.8*z/D
    epsilon = 1e-8
    p = A*p_ult*np.tanh(k*z*y/(A*p_ult + epsilon))

    f = interp1d(y, p, kind='linear')   # Interpolation function for p-y curve
    
    if print_curves == 'Yes':
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
#### Soil Profiles #####
########################

def clay_profile(profile):
    '''Define the clay profile used by the p-y analyzer. Outputs 'interp1d' functions containing Su and sigma'_v
    profiles to be used by the p-y curve functions.

    Input:
    -----
    profile      - A 2D tuple in the following format: ([Depth (m), Su (kPa), gamma_sub (kN/m^3), py-model, model parameter])
                   The soil profile should be defined relative to the pile/tower head (i.e. point of lateral load application)
                   so that any load eccentricities can be taken into account. An example soil profile is shown below.
                   Eg: array([[z0,Su0,gamma_sub0,  'Matlock', 0.02],
                              ...])

                    *The current program cannot define layers with different p-y models. But it will added in the future.

    plot_profile - Plot Su vs depth profile. Choose 'Yes' to plot.

    Output:
    ------
    z0            - Depth of mudline relative to the pile head (m)
    f_Su          - 'interp1d' function containing undrained shear strength profile (Pa)
    f_sigma_v_eff - 'interp1d' function containing effective vertical stress profile (Pa)
    f_gamma_sub   - 'interp1d' function containing effective unit weight (kN/m3)
    '''

    from scipy.interpolate import interp1d

    # Depth of mudline relative to pile head
    z0 = profile[0,0].astype(float)

    # Extract data from soil_profile array and zero strength virtual soil layer
    # from the pile head down to the mudline
    depth     = np.concatenate([np.array([0,z0]), profile[:,0].astype(float)])  # m
    Su        = np.concatenate([np.array([0, 0]), profile[:,1].astype(float)])  # kPa
    gamma_sub = np.concatenate([np.array([0, 0]), profile[:,2].astype(float)])  # kN/m^3

    if plot_profile == 'Yes':
        # Plot Su vs z profile for confirmation
        plt.plot(Su, depth, '-', label=r'$S_u$')
        plt.legend(loc='lower left')
        plt.xlabel('Undrained shear strength (kPa)')
        plt.ylabel('Depth below the pile head (m)')
        plt.grid(True)

        # Plot mudline/ground surface
        plt.plot([-0.5*max(Su), max(Su)], [z0,z0], '--', color='brown')
        plt.text(-0.5*max(Su), 0.95*z0, 'Mudline', color='brown')

        ax = plt.gca()
        ax.invert_yaxis(), ax.xaxis.tick_top(), ax.xaxis.set_label_position('top')

    # Calculate sigma_v_eff at each depth
    sigma_v_eff = np.zeros(len(depth))

    for i in range(1, len(depth)):
        sigma_v_eff[i] = sigma_v_eff[i-1] + gamma_sub[i-1]*(depth[i] - depth[i-1])

    # Define interpolation functions
    f_Su          = interp1d(depth, Su*1000, kind='linear')          # Pa
    f_sigma_v_eff = interp1d(depth, sigma_v_eff*1000, kind='linear') # Pa
    f_gamma_sub   = interp1d(depth, gamma_sub,kind='linear')         # kN/m3

    return z0, f_Su, f_sigma_v_eff, f_gamma_sub

def sand_profile(profile):
    '''Define the sand profile used by the p-y analyzer. Outputs 'interp1d' functions containing Su and sigma'_v
    profiles to be used by the p-y curve functions.

    Input:
    -----
    profile      - A 2D tuple in the following format: ([Depth (m), Su (kPa), gamma_sub (kN/m^3), py-model, model parameter])
                   The soil profile should be defined relative to the pile/tower head (i.e. point of lateral load application)
                   so that any load eccentricities can be taken into account. An example soil profile is shown below.
                   Eg: array([[z0,Su0,gamma_sub0,  'API', 0.02],
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

    # Extract data from soil_profile array and zero strength virtual soil layer
    # from the pile head down to the mudline
    depth     = np.concatenate([np.array([0,z0]), profile[:,0].astype(float)])  # m
    phi       = np.concatenate([np.array([0, 0]), profile[:,1].astype(float)])  # deg
    gamma_sub = np.concatenate([np.array([0, 0]), profile[:,2].astype(float)])  # kN/m^3
   
    if plot_profile == 'Yes':
        # Plot friction angle vs z profile for confirmation
        plt.plot(phi, depth, '-', label=r'$S_u$')
        plt.legend(loc='lower left')
        plt.xlabel('Friction angle (deg)')
        plt.ylabel('Depth below the pile head (m)')
        plt.grid(True)

        # Plot mudline/ground surface
        plt.plot([-0.5*max(phi), max(phi)], [z0,z0], '--', color='brown')
        plt.text(-0.5*max(phi), 0.95*z0, 'Mudline', color='brown')

        ax = plt.gca()
        ax.invert_yaxis(), ax.xaxis.tick_top(), ax.xaxis.set_label_position('top')

    # Calculate sigma_v_eff and static loading factor at each depth
    sigma_v_eff = np.zeros(len(depth))

    for i in range(1, len(depth)):  
        sigma_v_eff[i] = sigma_v_eff[i-1] + gamma_sub[i-1]*(depth[i] - depth[i-1])

    # Define interpolation functions
    f_phi = interp1d(depth, phi, kind='linear')                      # deg
    f_sigma_v_eff = interp1d(depth, sigma_v_eff*1000, kind='linear') # Pa
    
    return z0, f_phi, f_sigma_v_eff