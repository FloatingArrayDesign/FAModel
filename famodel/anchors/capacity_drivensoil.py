# Copyright Asitha Senanayake 2020

import numpy as np
import matplotlib.pyplot as plt

###################################
#### Pile Geometry and Loading ####
###################################

def py_analysis(soil_profile, L, D, t, E, F,
                     V_0, M_0, V_n, M_n, n=50, iterations=10, 
                     soil_type='clay', print_output='Yes',
                     convergence_tracker='Yes', loc=2, **kwargs):
    '''Models a laterally loaded pile using the p-y method. The solution for
    lateral displacements is obtained by solving the 4th order ODE, EI*d4y/dz4
    - F*d2y/dz2 + ky = 0 using the finite difference method.

    Assumes that EI remains constant with respect to curvature i.e. pile
    material remains in the elastic region.

    Input:
    -----
    Su_profile  - A 2D array of depths (m) and corresponding undrained shear strength(Pa)
                  Eg: array([[z1,Su1],[z2,Su2],[z3,Su3]...])
                  Use small values for Su (eg: 0.001) instead of zeros to avoid divisions by zero but always start z at 0.0
                  Example of a valid data point at the mudline is [0.0, 0.001]

    L           - Length of pile         (m)
    D           - Outer diameter of pile (m)
    t           - Wall thickness of pile (m)
    E           - Elastic modulus of pile material (Pa)
    F           - Axial force at pile head (N), vertically downwards is postive.
    V_0, V_n    - Force at pile head/tip  (N),  shear causing clockwise rotation of pile is positive.
    M_0, M_n    - Moment at pile head/tip (N-m), moments causing tension on left side of pile is positive.
    n           - Number of elements (50 by default)
    iterations  - Number of iterations to repeat calculation in order obtain convergence of 'y'
                  (A better approach is to iterate until a predefined tolerance is achieved but this requires additional
                  coding so, I will implement this later.)
    soil        - Select soil condition, 'clay' or 'sand'
                - Select which p-y model to use, 'Matlock' or 'API'.

    Optional:
    convergence_tracker - Track how k_secant converges to actual p-y curve at a selected node
    loc                 - Node number at which k_secant to be tracked (2 to n+1)

    Output:
    ------
    y           - Lateral displacement at each node, length = n + 5, (n+1) real nodes and 4 imaginary nodes
    z           - Vector of node locations along pile
    '''

    # Extract optional keyword arguments
    epsilon_50, A, gapping, N_p_max = 0.02, 550, 'No', 12.0  # Default parameters if no **kwargs are defined
    custom_py, a, strain_f          = 'No', 0.0, 0.0
    ls, alpha                       = 'x', 0.0

    for arg in kwargs:
        if arg == 'epsilon_50':
            epsilon_50 = kwargs[arg]
        if arg == 'Gmax_Su_ratio':
            A = kwargs[arg]
        if arg == 'gapping':
            gapping = kwargs[arg]
        if arg == 'N_p_max':
            N_p_max = kwargs[arg]
        if arg == 'alpha':
            alpha = kwargs[arg]
        if arg == 'custom_py':
            custom_py = kwargs[arg]
        if arg == 'a':
            a = kwargs[arg]
        if arg == 'strain_f':
            strain_f = kwargs[arg]
        if arg == 'ls':
            ls = kwargs[arg]

    # Convert L and D to floating point numbers to avoid rounding errors
    L = float(L)
    D = float(D)

    # Pile geometry
    I  = np.pi*(D**4 - (D - 2*t)**4)/64.0
    EI = E*I
    h  = L/n              # Element size
    N  = (n + 1) + 4      # (n+1) Real + 4 Imaginary nodes

    # Array for displacements at nodes, including imaginary nodes.
    y = np.ones(N)*(0.01*D)   # An initial value of 0.01D was arbitrarily chosen

    # Initialize and assemble array/list of p-y curves at each real node
    z = np.zeros(N)
    py_funs  = []
    k_secant = np.zeros(N)

    for i in [0, 1]:        # Top two imaginary nodes
        z[i] = (i-2)*h
        py_funs.append(0)
        k_secant[i] = 0.0

    # Extract soil profile data
    if soil_type == 'clay':
        z_0, f_Su, f_sigma_v_eff, f_gamma_sub = soil_clay(soil_profile)
    elif soil_type == 'sand':
        z_0, f_phi, f_sigma_v_eff = soil_sand(soil_profile)

    for i in range(2, n + 3):  # Real nodes
        z[i] = (i - 2)*h      
        if soil_type == 'clay':
            Su, Su0, sigma_v_eff, gamma_sub = f_Su(z[i]), f_Su(z_0+0.01), f_sigma_v_eff(z[i]), f_gamma_sub(z[i])
            py_funs.append(py_Matlock(z[i], D, Su, sigma_v_eff, gamma_sub, z_0=z_0, epsilon_50=epsilon_50))
        elif soil_type == 'sand':
            phi, sigma_v_eff = f_phi(z[i]), f_sigma_v_eff(z[i])
            py_funs.append(py_API(z[i], D, phi, sigma_v_eff, z_0=z_0))
        
        k_secant[i] = py_funs[i](y[i])/y[i]

    for i in [n+3, n+4]:   # Bottom two imaginary nodes
        z[i] = (i-2)*h
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
        y = fd_solver(n, N, h, EI, F, V_0, V_n, M_0, M_n, k_secant)

        if convergence_tracker == 'Yes':
            plt.plot(y[loc], k_secant[loc]*y[loc], ls)

        for i in range(2, n+3):
            k_secant[i] = py_funs[i](y[i])/y[i]

    if print_output == 'Yes':
        print(f'y_0 = {y[2]:.3f} m')

    return y[2:-2], z[2:-2]

#################
#### Solvers ####
#################

def fd_solver(n,N,h,EI,F,V_0,V_n,M_0,M_n,k_secant):
    '''Solves the finite difference equations from 'py_analysis_1'. This function should be run iteratively for
    non-linear p-y curves by updating 'k_secant' using 'y'. A single iteration is sufficient if the p-y curves
    are linear.

    Input:
    -----
    n - Number of elements
    N - Total number of nodes
    h - Element size
    EI - Flexural rigidity of pile
    F  - Axial force at pile head
    V_0, V_n - Shear at pile head/tip
    M_0, M_n - Moment at pile head/tip
    k_secant - Secant stiffness from p-y curves

    Output:
    ------
    y_updated - Lateral displacement at each node
    '''

    from scipy import linalg

    # Initialize and assemble matrix
    X = np.zeros((N,N))

    # (n+1) finite difference equations for (n+1) real nodes
    for i in range(0,n+1):
        X[i,i]   =  1.0
        X[i,i+1] = -4.0 + F*h**2/EI
        X[i,i+2] =  6.0 - 2*F*h**2/EI + k_secant[i+2]*h**4/EI
        X[i,i+3] = -4.0 + F*h**2/EI
        X[i,i+4] =  1.0

    # Curvature at pile head
    X[n+1,1]   =  1.0
    X[n+1,2]   = -2.0
    X[n+1,3]   =  1.0

    # Shear at pile head
    X[n+2,0]   = -1.0
    X[n+2,1]   =  2.0 - F*h**2/EI
    X[n+2,2]   =  0.0
    X[n+2,3]   = -2.0 + F*h**2/EI
    X[n+2,4]   =  1.0

    # Curvature at pile tip
    X[n+3,-2]   =  1.0
    X[n+3,-3]   = -2.0
    X[n+3,-4]   =  1.0

    # Shear at pile tip
    X[n+4,-1]   =   1.0
    X[n+4,-2]   =  -2.0 + F*h**2/EI
    X[n+4,-3]   =   0.0
    X[n+4,-4]   =   2.0 - F*h**2/EI
    X[n+4,-5]   =  -1.0

    # X*y = q

    # Initialize vector q
    q = np.zeros(N)

    # Populate q with boundary conditions
    q[-1] = 2*V_n*h**3     # Shear at pile tip
    q[-2] = M_n*h**2       # Moment at pile tip
    q[-3] = 2*V_0*h**3     # Shear at pile head
    q[-4] = M_0*h**2       # Moment at pile head

    y = linalg.solve(EI*X, q)

    return y

###############################
#### P-Y Curve Definitions ####
###############################

def py_Matlock(z, D, Su, sigma_v_eff, gamma_sub, z_0=0.0, epsilon_50=0.02, print_curves='Yes'):
    
    '''Returns an interp1d interpolation function which represents the Matlock (1970) p-y curve at the depth of interest.

    Important: Make sure to import the interp1 function by running 'from scipy.interpolate import interp1d' in the main program.

    Input:
    -----
    z            - Depth relative to pile head (m)
    D            - Pile diameter (m)
    Su           - Undrained shear strength (Pa)
    sigma_v_eff  - Effective vertical stress (Pa)
    gamma_sub    - Effective unit weight of the soil (kN/m3)
    z_0          - Load eccentricity above the mudline or depth to mudline relative to the pile head (m)
    epsilon_50   - Strain at half the strength as defined by Matlock (1970).
                   Typically ranges from 0.005 (stiff clay) to 0.02 (soft clay).

    Output:
    ------
    Returns an interp1d interpolation function which represents the p-y curve at the depth of interest.
    'p' (N/m) and 'y' (m).
    '''

    from scipy.interpolate import interp1d

    # p-y curve properties
    J = 0.5

    if (z - z_0) < 0:
        # p-y curves for the virtual soil layer between the pile head and the mudline should have p=0
        Nc = 0.0
        z_cr = 1.0 # Dummy value to keep program from crashing

    else:
        try:
            Nc = 3.0 + sigma_v_eff/Su + J*(z - z_0)/D
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
       
def py_API(z, D, phi, sigma_v_eff, z_0=0.0, print_curves='Yes'):
    
    '''Returns an interp1d interpolation function which represents the Matlock (1970) p-y curve at the depth of interest.

    Important: Make sure to import the interp1 function by running 'from scipy.interpolate import interp1d' in the main program.

    Input:
    -----
    z            - Depth relative to pile head (m)
    D            - Pile diameter (m)
    Su           - Undrained shear strength (Pa)
    sigma_v_eff  - Effectve vertical stress (Pa)
    z_0          - Load eccentricity above the mudline or depth to mudline relative to the pile head (m)
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
    
    if (z - z_0) < 0:
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

def soil_clay(soil_profile, plot_profile='Yes'):
    '''Define the clay profile used by the p-y analyzer. Outputs 'interp1d' functions containing Su and sigma'_v
    profiles to be used by the p-y curve functions.

    Input:
    -----
    soil_profile - A 2D tuple in the following format: ([Depth (m), Su (kPa), gamma_sub (kN/m^3), py-model, model parameter])
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
    z0 = soil_profile[0,0].astype(float)

    # Extract data from soil_profile array and zero strength virtual soil layer
    # from the pile head down to the mudline
    depth     = np.concatenate([np.array([0,z0]), soil_profile[:,0].astype(float)])  # m
    Su        = np.concatenate([np.array([0, 0]), soil_profile[:,1].astype(float)])  # kPa
    gamma_sub = np.concatenate([np.array([0, 0]), soil_profile[:,2].astype(float)])  # kN/m^3

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
    f_Su = interp1d(depth, Su*1000, kind='linear')                   # Pa
    f_sigma_v_eff = interp1d(depth, sigma_v_eff*1000, kind='linear') # Pa
    f_gamma_sub = interp1d(depth, gamma_sub,kind='linear')           # kN/m3

    return z0, f_Su, f_sigma_v_eff, f_gamma_sub

def soil_sand(soil_profile, plot_profile='Yes'):
    '''Define the sand profile used by the p-y analyzer. Outputs 'interp1d' functions containing Su and sigma'_v
    profiles to be used by the p-y curve functions.

    Input:
    -----
    soil_profile - A 2D tuple in the following format: ([Depth (m), Su (kPa), gamma_sub (kN/m^3), py-model, model parameter])
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
    z0 = soil_profile[0,0].astype(float)

    # Extract data from soil_profile array and zero strength virtual soil layer
    # from the pile head down to the mudline
    depth     = np.concatenate([np.array([0,z0]), soil_profile[:,0].astype(float)])  # m
    phi       = np.concatenate([np.array([0, 0]), soil_profile[:,1].astype(float)])  # deg
    gamma_sub = np.concatenate([np.array([0, 0]), soil_profile[:,2].astype(float)])  # kN/m^3
   
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

if __name__ == '__main__':
    
    #                        Depth   Su  γ_sub       p-y model   p-y parameters
    soil_profile = np.array([[0.0,  10.0, 10., 'Name of p-y model', 0.02],
                             [5.0,  10.0, 10., 'Name of p-y model', 0.02],
                             [10.0, 15.0, 10., 'Name of p-y model', 0.02],
                             [25.0, 15.0, 10., 'Name of p-y model', 0.02]])

    z0, f_Su, f_σ_v_eff, f_gamma_sub = soil_clay(soil_profile)
    
    # #                        Depth  phi  γ_sub       p-y model   p-y parameters
    # soil_profile = np.array([[0.0,  30.0, 17., 'Name of p-y model', 0.02],
    #                          [5.0,  32.0, 18., 'Name of p-y model', 0.02],
    #                          [10.0, 35.0, 19., 'Name of p-y model', 0.02],
    #                          [35.0, 38.0, 20., 'Name of p-y model', 0.02]])

    # z0, f_phi, f_σ_v_eff = soil_sand(soil_profile)
    
    #Pile dimensions
    L = 25                 # Pile length (m)
    D = 1.5                # Pile diameter (m)
    t = (6.35 + D*20)/1e3  # Pile wall thickness (m), API RP2A-WSD
    E = 200e9              # Elastic modulus of steel (Pa)
    
    #Pile head loads
    V_0 = 50e4             # N,  horizontal load on pile head
    M_0 = 0                # Nm, moment on pile head
    F_0 = 0                # N,  vertical load on pile head
    
    y,z = py_analysis(soil_profile, L=L, D=D, t = t, E=E, F = F_0,
                         V_0=V_0, M_0=M_0, V_n=0.0, M_n=0.0, n=100)
    
    y0 = np.zeros(len(z))
    #Plot deflection profile of pile
    fig, ax = plt.subplots(figsize=(3,5))    
    ax.plot(y0,z,'black')
    ax.plot(y,z,'r')
    ax.set_xlabel('Displacement [m]')
    ax.set_ylabel('Depth below pile head [m]')
    ax.set_ylim([L + 2, -2])
    ax.set_xlim([-0.1*D, 0.1*D])
    ax.grid(ls='--')