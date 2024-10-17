
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
    t = (6.35 + D*20)/1e3       # Pile wall thickness (m), API RP2A-WSD
    E = 200e9                   # Elastic modulus of pile material (Pa)
    fy = 350e6                  # Yield strength of pile material (Pa)
    rhows = 66.90e3             # Submerged steel specific weight (N/m3)
    rhow = 10e3                 # Water specific weight (N/m3) 
    alpha = 0.6                 # Adhesion coefficient of clay (-)
    beta = 0.37                 # Shaft friction factor (-)
    
    # Pile geometry
    I = (np.pi/64.0)*(D**4 - (D - 2*t)**4)
    EI = E*I
    h  = L/n                    # Element size
    N  = (n + 1) + 4            # (n+1) Real + 4 Imaginary nodes
    
    # Outer and inner surface of the pile skirt
    def PileSurface(Len, Dia):
        Sp = np.pi*Dia*Len
        return Sp    
    # Dry and wet mass of the pile    
    def PileWeight(Len, Dia, tw, rho):
        Wp = ((np.pi/4)*((Dia**2 - (Dia - 2*tw)**2)*Len 
        return Wp 
    # Mass of the soil plug      
    def SoilWeight(Len, Dia, tw, gamma_soil): 
        Wsoil =(np.pi/4)*(Dia - 2*tw)**2*Len*gamma_soil
        return Wsoil

    # Array for displacements at nodes, including imaginary nodes.
    y = np.ones(N)*(0.01*D)     # An initial value of 0.01D was arbitrarily chosen

    # Initialize and assemble array/list of p-y curves at each real node
    z = np.zeros(N)
    py_funs  = []; PileShaft =[]
    k_secant = np.zeros(N)

    for i in [0, 1]:            # Top two imaginary nodes
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0

    # Extract soil profile data
    if soil_type == 'clay':
        z0, f_Su, f_sigma_v_eff, f_gamma, f_alpha = clay_profile(profile)

    elif soil_type == 'sand':
        z0, f_phi, f_sigma_v_eff, f_gamma, f_Dr, f_beta = sand_profile(profile)

    for i in range(2, n+3):    # Real nodes
        z[i] = (i - 2)*h      
        if soil_type == 'clay':
            Su, sigma_v_eff, gamma, alpha = f_Su(z[i]), f_sigma_v_eff(z[i]), f_gamma(z[i]), f_alpha(z[i])
            py_funs.append(py_Matlock(z[i], D, zlug, Su, sigma_v_eff, gamma, plot=plot))
            Vo = np.pi*D*alpha*Su*z[i]**2
            PileShaft.append(Vo)
            Vmax = PileWeight(L, D, t, rhows) + SoilWeight(L, D, t, gamma) + PileShaft[-1]
            
        elif soil_type == 'sand':
            phi, sigma_v_eff, gamma, Dr, beta = f_phi(z[i]), f_sigma_v_eff(z[i]), f_gamma(z[i]), f_Dr(z[i]), f_beta(z[i])
            py_funs.append(py_API(z[i], D, zlug, phi, sigma_v_eff, Dr, plot=plot))
            fs = beta*sigma_v_eff
            Vo = np.pi*D*fs*z[i]
            PileShaft.append(Vo)
            Vmax = PileWeight(L, D, t, rhows) + SoilWeight(L, D, t, gamma) + PileShaft[-1]
        
        k_secant[i] = py_funs[i](y[i])/y[i]

    for i in [n+3, n+4]:        # Bottom two imaginary nodes
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0

    y1 = np.linspace(-2.*D, 2.*D, 500)
    if plot:
        plt.plot(y1, py_funs[loc](y1))
        plt.xlabel('y (m)'), plt.ylabel('p (N/m)')
        plt.grid(True)

    for j in range(iterations):
        # if j == 0: print 'FD Solver started!'
        y, Mi, Mp, hinge_formed, hinge_location = fd_solver(n, N, h, D, t, fy, EI, V, H, zlug, k_secant)

        for i in range(2, n+3):
            k_secant[i] = py_funs[i](y[i])/y[i]


    resultsDrivenSoil = {}
    # Populate q with boundary conditions
    if zlug <= 0:      
        print(f'y_max = {max(y):.3f} m')
        print(f'rot_max = {np.rad2deg((y[2] - y[3])/h):.3f} deg')
               
        resultsDrivenSoil['Lateral displacement'] = max(y)
        resultsDrivenSoil['Rotational displacement'] = np.rad2deg((y[2] - y[3])/h)
        resultsDrivenSoil['Axial capacity'] = Vmax
        resultsDrivenSoil['Pile weight'] = PileWeight(L, D, t, (rhows + rhow))
   
    else:
        print(f'y_max = {max(y):.3f} m')
        print(f'Mi = {Mi:.3f} Nm')
        
        resultsDrivenSoil['Lateral displacement'] = max(y)
        resultsDrivenSoil['Bending moment'] = Mi
        resultsDrivenSoil['Plastic moment'] = Mp
        resultsDrivenSoil['Plastic hinge'] = hinge_formed
        resultsDrivenSoil['Hinge location'] = hinge_location
        resultsDrivenSoil['Axial capacity'] = Vmax
        resultsDrivenSoil['Pile weight'] = PileWeight(L, D, t, (rhows + rhow))
    
    return y[2:-2], z[2:-2], resultsDrivenSoil

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
    #M = H*abs(zlug)

    # Populate q with boundary conditions
    if zlug <= 0:
        q[-3] = 2*H*h**3      # Shear at pile head
        #q[-4] = M*h**2        # Moment at pile head
    else:
        q[zlug_index]      = 2*H*h**3     # Shear at pile head
        #q[zlug_index + 1]  = M*h**2      # Moment at pile head
               
    y = linalg.solve(EI*X, q)

    # Compute the plastic moment capacity Mp
    Zp = (1/6)*(D**3 - (D - 2*t)**3)  # Plastic section modulus for hollow pile (m3)
    Mp = Zp*fy                        # Plastic moment capacity (N/m)

    # Check for plastic hinge formation
    Mi, Mp, hinge_formed, hinge_location = plastic_hinge(y, h, EI, Mp)
    
    return y, Mi, Mp, hinge_formed, hinge_location

###############################
#### P-Y Curve Definitions ####
###############################

def py_Matlock(z, D, zlug, Su, sigma_v_eff, gamma, plot=True):
    
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
    if plot:
        plt.plot(y, p,'-')
        plt.xlabel('y (m)') 
        plt.ylabel('p (N/m)')
        plt.title('PY Curves - Matlock (1970)')
        plt.grid(True)
        plt.xlim([-2*D, 2*D])

    return f   # This is f (linear interpolation of y-p)
       
def py_API(z, D, zlug, phi, sigma_v_eff, Dr, plot=True):
    
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
    phi_ref = [  20,   25,   30,   35,   40]
    C1_ref  = [0.80, 1.25, 1.90, 3.00, 4.50]
    C2_ref  = [1.60, 2.10, 2.60, 3.40, 4.30]
    C3_ref  = [  10,   15,   30,   55,  105]
    
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
    
    # Dr = 0.75     # Relative density of the soil (assumed)
    k = (54.6*Dr**2 + 0.8*Dr + 1.8)*1e3
    
    # Normalized lateral displacement
    N = 20
    y = np.concatenate((-np.logspace(3,-4,N),[0],np.logspace(-4,3,N)))
    A = max(3 - 0.8*z/D, 0.9)
    ε = 1e-6
    p = A*p_ult*np.tanh(k*z*y/(A*p_ult + ε))

    f = interp1d(y, p, kind='linear')   # Interpolation function for p-y curve
    
    if plot:
        # Plot of p-y curve and check if 'k' is calculated correctly
        plt.plot(y, p,'-') 
        plt.xlabel('y (m)') 
        plt.ylabel('p (N/m)')
        plt.title('PY Curves - API (1993)')
        plt.grid(True)
        plt.xlim([-0.10*D, 0.10*D])
    # plt.ylim([min(y), max(y)])  # Adjust x-axis limits to match y values
        
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
    Mp : float
        Plastic moment of the pile section (Nm)
    hinge_formed : bool
        True if a plastic hinge forms, False otherwise.
    hinge_location : int
        Index of the node where the plastic hinge forms (if any).
    '''
    
    hinge_formed = False
    hinge_location = -1
    Mi = []
        
    # Loop through each internal node and compute the bending moment
    for i in range(1, len(y)-1):
        # Approximate the bending moment at node i
        Mint = EI*(y[i+1] - 2*y[i] + y[i-1])/h**2
        Mi.append(Mint)
        
        # Check if the moment exceeds the plastic moment capacity
        if Mint >= Mp:
            hinge_formed = True
            hinge_location = i
            break  # Stop at the first plastic hinge formation
    
    return max(Mi), Mp, hinge_formed, hinge_location


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
                   Eg: array([[z0, Su0, gamma0, 'Matlock', 0.02],
                              ...])

                   *The current program cannot define layers with different p-y models. But it will added in the future.

    Output:
    ------
    z0            - Depth of mudline relative to the pile head (m)
    f_Su          - 'interp1d' function containing undrained shear strength profile (Pa)
    f_sigma_v_eff - 'interp1d' function containing effective vertical stress profile (Pa)
    f_gamma       - 'interp1d' function containing effective unit weight (kN/m3)
    f_alpha       - Adhesion factor for clays
    '''

    from scipy.interpolate import interp1d

    # Depth of mudline relative to pile head
    z0 = float(profile[0][0])

    # Extract data from profile array and zero strength virtual soil layer
    # from the pile head down to the mudline
    depth = np.concatenate([np.array([z0]),np.array([row[0] for row in profile],dtype=float)]) # m
    Su    = np.concatenate([np.array([0]), np.array([row[1] for row in profile],dtype=float)]) # kPa
    gamma = np.concatenate([np.array([0]), np.array([row[2] for row in profile],dtype=float)]) # kN/m3

    # Calculate sigma_v_eff at each depth
    sigma_v_eff = np.zeros(len(depth))

    for i in range(1, len(depth)):
        sigma_v_eff[i] = sigma_v_eff[i-1] + gamma[i-1]*(depth[i] - depth[i-1])

    # Define interpolation functions
    f_Su          = interp1d(depth, Su*1000, kind='linear')           # Pa
    f_sigma_v_eff = interp1d(depth, sigma_v_eff*1000, kind='linear')  # Pa
    f_gamma       = interp1d(depth, gamma*1000, kind='linear')        # N/m3
    
    # Calculate f_psi and f_alpha at each depth (not as a scalar)
    f_psi = lambda z: f_Su(z) / f_sigma_v_eff(z)
    
    def calc_alpha(z):
        psi_val = f_psi(z)
        if psi_val <= 1.0:
            return min(0.5*psi_val**-0.50, 1)
        else:
            return min(0.5*psi_val**-0.25, 1)

    # Create an interpolated adhesion factor function
    f_alpha = lambda z: calc_alpha(z)

    return z0, f_Su, f_sigma_v_eff, f_gamma, f_alpha

def sand_profile(profile):
    '''Define the sand profile used by the p-y analyzer. Outputs 'interp1d' functions containing Su and sigma'_v
    profiles to be used by the p-y curve functions.

    Input:
    -----
    profile      - A 2D tuple in the following format: ([Depth (m), Su (kPa), gamma (kN/m^3), py-model, model parameter])
                   The soil profile should be defined relative to the pile/tower head (i.e. point of lateral load application)
                   so that any load eccentricities can be taken into account. An example soil profile is shown below.
                   Eg: array([[z0, phi, gamma0, Dr, 'API', 0.02],
                              ...])

                   *The current program cannot define layers with different p-y models. But it will added in the future.

    Output:
    ------
    z0            - Depth of mudline relative to the pile head (m)
    f_phi         - 'interp1d' function containing effective friction angle (deg)
    f_sigma_v_eff - 'interp1d' function containing effective vertical stress profile (Pa)
    f_gamma       - 'interp1d' function containing effective unit weight (N/m3)
    f_Dr          - Relative density of the soil (%)
    f_beta        - Adhesion factor for clays
    '''

    from scipy.interpolate import interp1d

    # Depth of mudline relative to pile head
    z0 = float(profile[0][0])

    # Extract data from profile array and zero strength virtual soil layer
    # from the pile head down to the mudline
    depth = np.concatenate([np.array([z0]),np.array([row[0] for row in profile],dtype=float)]) # m
    phi   = np.concatenate([np.array([0]), np.array([row[1] for row in profile],dtype=float)]) # deg
    gamma = np.concatenate([np.array([0]), np.array([row[2] for row in profile],dtype=float)]) # kN/m3
    Dr    = np.concatenate([np.array([0]), np.array([row[3] for row in profile],dtype=float)]) # %
   
    # Calculate sigma_v_eff and static loading factor at each depth
    sigma_v_eff = np.zeros(len(depth))

    for i in range(1, len(depth)):  
        sigma_v_eff[i] = sigma_v_eff[i-1] + gamma[i-1]*(depth[i] - depth[i-1])

    # Define interpolation functions
    f_phi         = interp1d(depth, phi, kind='linear')                    # deg
    f_sigma_v_eff = interp1d(depth, sigma_v_eff*1000, kind='linear')       # Pa
    f_gamma       = interp1d(depth, gamma*1000, kind='linear')             # N/m3
    f_Dr          = interp1d(depth, Dr, kind='linear')                     # %
    
    # Define beta as a function of Dr
    def calc_beta(Dr_val):
        if 35 <= Dr_val < 50:
            return 0.29
        elif 50 <= Dr_val < 65:
            return 0.37
        elif 65 <= Dr_val < 85:
            return 0.46
        elif Dr_val >= 85:
            return 0.56
        else:
            return 0  # Default or error value for very low Dr values
        
    # Apply beta calculation to Dr profile
    beta_values = np.array([calc_beta(Dr_val) for Dr_val in Dr])
    f_beta      = interp1d(depth, beta_values, kind='linear')  # Interpolated beta values
    
    return z0, f_phi, f_sigma_v_eff, f_gamma, f_Dr, f_beta

if __name__ == '__main__':

    profile = np.array([[ 0.0, 25, 8, 85, 'Name of p-y model'],
                        [75.0, 45, 5, 35, 'Name of p-y model']])
    
    L = 50
    D = 1
    zlug = 30*D
    
    H0 = 426000
    V0 = 15900

    H = 1e4; V = 1e4

    values_H =[]; values_V =[]      
    
    y, z, results = getCapacityDrivenSoil(profile, soil_type='sand', L=L, D=D, zlug=zlug, V=V, H=H)

    # while results['Lateral displacement']<= 0.05*D and results['Rotational displacement'] <= 0.25:
               
    #     y, z, results = getCapacityDrivenSoil(profile, soil_type='clay', L=L, D=D, zlug=zlug, V=V, H=H)
        
    #     H += 10000
                   
    # values_H.append(H); H_ratio = np.array(values_H)/H0       
   
    while results['Lateral displacement']<= 0.05*D and results['Bending moment'] <= results['Plastic moment']:
               
        y, z, results = getCapacityDrivenSoil(profile, soil_type='sand', L=L, D=D, zlug=zlug, V=V, H=H)
        
        H += 100000
                   
    values_H.append(H); H_ratio = np.array(values_H)/H0
    
    # y, z, results = getCapacityDrivenSoil(profile, soil_type='sand', L=L, D=D, zlug=zlug, V=V, H=H)
    
        
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
    fig.show()