import numpy as np
#from matplotlib.pyplot import plot, draw, show
#import matplotlib.backends
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import linalg
import inspect

###################################
#### Pile Geometry and Loading ####
###################################

def py_analysis(profile, L, D, t, E,  
                V, H, H_n, M, M_n, n=10, iterations=10,
                print_output='No', convergence_tracker='No', loc=2):
    '''
    Models a laterally loaded pile using the p-y method. The solution for
    lateral displacements is obtained by solving the 4th order ODE, 
    EI*d4y/dz4 - V*d2y/dz2 + ky = 0 using the finite difference method.
    EI*d4y/dz4 - V*d2y/dz2 + K*z*dy/dz + ky = 0 using the finite difference method.

    Assumes that EI remains constant with respect to curvature i.e. pile
    material remains in the elastic region.

    Input:
    -----
    Su_profile  - A 2D array of depths (m) and corresponding undrained shear strength(Pa)
                  Eg: array([[z1,Su1],[z2,Su2],[z3,Su3]...])
                  Use small values for Su (eg: 0.001) instead of zeros to avoid divisions 
                  by zero but always start z at 0.0
                  Example of a valid data point at the mudline is [0.0, 0.001]

    L           - Length of pile         (m)
    D           - Outer diameter of pile (m)
    t           - Wall thickness of pile (m)
    E           - Elastic modulus of pile material (Pa)
    F           - Axial force at pile head (N), vertically downwards is postive.
    V           - Force at pile head (N), shear causing clockwise rotation of pile is positive.
    M           - Moment at pile head (N*m), moments causing tension on left side of pile is positive.
    n           - Number of elements (50 by default)
    iterations  - Number of iterations to repeat calculation in order obtain convergence of 'y'
                  (A better approach is to iterate until a predefined tolerance is achieved but this requires additional
                  coding so, I will implement this later.)
    py_model    - Select which p-y model to use, 'Matlock' or 'Jeanjean_2009', 'MM-1', or 'Jeajean_etal_2017'.

    Optional:
    convergence_tracker - Track how k_secant converges to actual p-y curve at a selected node
    loc         - Node number at which k_secant to be tracked (2 to n+1)

    Output:
    ------
    y           - Lateral displacement at each node, length = n + 5, (n+1) real nodes and 4 imaginary nodes
    z           - Vector of node locations along pile
    '''
    
    # Extract optional keyword arguments
    # ls = 'x'

    # Resistance factor
    nhuc = 1; nhu = 0.3; gamma_f = 1.3
    delta_r = 0.08  # Mean roughness height [m]
    
    # Convert L and D to floating point numbers to avoid rounding errors
    L = float(L)
    D = float(D)

    # Pile geometry
    I = np.pi*(D**4 - (D - 2*t)**4)/64.0
    EI = E*I
    h = L/n           # Element size
    N = (n + 1) + 4   # (n+1) Real + 4 Imaginary nodes

    # Array for displacements at nodes, including imaginary nodes.
    y = np.ones(N)*(0.01*D)   # An initial value of 0.01D was arbitrarily chosen

    # Initialize and assemble array/list of p-y curves at each real node
    z = np.zeros(N)
    py_funs = []
    k_secant = np.zeros(N)
    DQ = []

    for i in [0,1]:           # Top two imaginary nodes
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0
   
    for i in range(2,n+3):  # Real nodes
        z[i] = (i - 2)*h
        # Extract rock profile data
        f_UCS, f_Em = rock_profile(profile)
        UCS, Em = f_UCS(z[i])/gamma_f, f_Em(z[i])/gamma_f
        py_funs.append(py_Reese(z[i], D, t, UCS, Em))
        k_secant[i] = py_funs[i](y[i])/y[i]
        SCR = nhuc*Em/(UCS*(1 + nhu))*delta_r/D
        alpha = 0.36*SCR - 0.0005
        fs = alpha*UCS
        Dq = np.pi*D*fs*z[i]
        DQ.append(Dq)

    for i in [n+3,n+4]:      # Bottom two imaginary nodes
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0

    # Track k_secant and current displacements
    if convergence_tracker == 'Yes':
        y1 = np.linspace(-2.*D, 2.*D, 500)
        plt.plot(y1, py_funs[loc](y1))
        plt.xlabel('y (m)'), plt.ylabel('p (N/m)')
        plt.grid(True)

    for j in range(iterations):
        # if j == 0: print 'FD Solver started!'
        y = fd_solver(n, N, h, EI, V, H, H_n, M, M_n, k_secant)

        if convergence_tracker == 'Yes':
            plt.plot(y[loc], k_secant[loc]*y[loc])

        for i in range(2, n+3):
            k_secant[i] = py_funs[i](y[i])/y[i]

    if print_output == 'Yes':
        print(f'y_0 = {y[2]:.3f} m')

    return y[2:-2], z[2:-2], DQ

#################
#### Solvers ####
#################

def fd_solver(n, N, h, EI, V, H, H_n, M, M_n, k_secant):
    '''
    Solves the finite difference equations from 'py_analysis_1'. This function should be run iteratively for
    non-linear p-y curves by updating 'k_secant' using 'y'. A single iteration is sufficient if the p-y curves
    are linear.

    Input:
    -----
    n        - Number of elements
    N        - Total number of nodes
    h        - Element size
    EI       - Flexural rigidity of pile
    V        - Axial force at pile head
    H        - Shear at pile head/tip
    M        - Moment at pile head/tip
    k_secant - Secant stiffness from p-y curves

    Output:
    ------
    y_updated - Lateral displacement at each node
    '''

    

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
    X[n+1,1] =  1.0
    X[n+1,2] = -2.0
    X[n+1,3] =  1.0

    # Shear at pile head
    X[n+2,0] = -1.0
    X[n+2,1] =  2.0 - V*h**2/EI
    X[n+2,2] =  0.0
    X[n+2,3] = -2.0 + V*h**2/EI
    X[n+2,4] =  1.0

    # Curvature at pile tip
    X[n+3,-2] =  1.0
    X[n+3,-3] = -2.0
    X[n+3,-4] =  1.0

    # Shear at pile tip
    X[n+4,-1] =  1.0
    X[n+4,-2] = -2.0 + V*h**2/EI
    X[n+4,-3] =  0.0
    X[n+4,-4] =  2.0 - V*h**2/EI
    X[n+4,-5] = -1.0

    # Initialize vector q
    q = np.zeros(N)

    # Populate q with boundary conditions
    q[-1] = 2*H_n*h**3    # Shear at pile tip
    q[-2] = M_n*h**2      # Moment at pile tip
    q[-3] = 2*H*h**3      # Shear at pile head
    q[-4] = M*h**2        # Moment at pile head

    y = linalg.solve(EI*X, q)

    return y

###############################
#### P-Y Curve Definitions ####
###############################

def py_Reese(z, D, t, UCS, Em, z_0=0.0, RQD=69, print_curves='Yes'):
    '''
    Returns an interp1d interpolation function which represents the Reese (1997) p-y curve at the depth of interest.

    Important: Make sure to import the interp1 function by running 'from scipy.interpolate import interp1d' in the main program.

    Input:
    -----
    z      - Depth relative to pile head (m)
    D      - Pile diameter (m)
    UCS    - Undrained shear strength (Pa)
    Em     - Effectve vertical stress (Pa)
    z_0    - Load eccentricity above the mudline or depth to mudline relative to the pile head (m)
    RQD    - Strain at half the strength as defined by Matlock (1970).
             Typically ranges from 0.005 (stiff clay) to 0.02 (soft clay).

    Output:
    ------
    Returns an interp1d interpolation function which represents the p-y curve at the depth of interest.
    'p' (N/m) and 'y' (m).
    '''

    #from scipy.interpolate import interp1d
    #global var_Reese
    
    Dref = 0.305; nhu = 0.3; E = 200e9
    I  = np.pi*(D**4 - (D - 2*t)**4)/64.0
    EI = E*I
    alpha = -0.00667*RQD + 1
    krm = 0.0005
    
    if z < 3*D:
        p_ur = alpha*UCS*D*(1 + 1.4*z/D)
        #kir = (100 +400*z/(3*D))
    else:
        p_ur = 5.2*alpha*UCS*D
        #kir = 500
        
    kir = (D/Dref)*2**(-2*nhu)*(EI/(Em*D**4))**0.284
    Kir = kir*Em     
    y_rm = krm*D
    y_a = (p_ur/(2*y_rm**0.25*Kir))**1.333    
   
    # Normalized lateral displacement
    N = 20
    y = np.concatenate((-np.logspace(5,-3,N),[0],np.logspace(-3,5,N)))
    
    p=[]; P=[];
    for i in range (len(y)):
        if abs(y[i]) < y_a: 
            P = np.sign(y[i])*Kir*y[i]
        elif abs(y[i]) > y_a:
            P = min((p_ur/2)*(abs(y[i])/y_rm)**0.25,p_ur)
        p.append(P)    
   
    p = np.array(p).squeeze()     
    for j in range(len(y)):
        if y[j] < 0:
            p[j] = -1*p[j]
        elif y[j] > 0:
            p[j] = p[j]                            
 
    #var_Reese = inspect.currentframe().f_locals          
    
    f = interp1d(y, p)   # Interpolation function for p-y curve
    
    if print_curves == 'Yes':           
        plt.plot(y, p) 
        plt.xlabel('y (m)') 
        plt.ylabel('p (kN/m)'),
        plt.title('PY Curves - Reese (1997)')
        plt.grid(True)
        plt.xlim([-0.03*D,0.03*D])
        plt.ylim([min(p),max(p)])     
        
    return f      # This is f (linear interpolation of y-p)
   
#######################
#### Rock Profile #####
#######################

def rock_profile(profile, plot_profile='No'):
    '''
    Define the (weak) rock profile used by the p-y analyzer. Outputs 'interp1d' functions containing 
    UCS and Em profiles to be used by the p-y curve functions.

    Input:
    -----
    profile - A 2D tuple in the following format: ([depth (m), UCS (MPa), Em (MPa), py-model])
              The soil profile should be defined relative to the pile/tower head (i.e. point of lateral load application)
              so that any load eccentricities can be taken into account. An example soil profile is shown below.
              Eg: array([[z0,UCS0,Em0, 'Reese'],
                         [z1,UCS1,Em1, 'Reese'],
                         [z2,UCS2,Em2, 'Reese'],
                          ...])
              *The current program cannot define layers with different p-y models. But it will added in the future.

    plot_profile - Plot Su vs depth profile. Choose 'Yes' to plot.

    Output:
    ------
    z0       - Depth of mudline relative to the pile head (m)
    f_UCS    - 'interp1d' function containing undrained shear strength profile (Pa)
    f_Em     - 'interp1d' function containing effective vertical stress profile (Pa)
    '''

    
    #global var_rock_profile

    # Depth of mudline relative to pile head
    z0 = profile[0,0].astype(float)

    # Extract data from soil_profile array and zero strength virtual soil layer
    # from the pile head down to the mudline
    depth = np.concatenate([np.array([z0]),profile[:,0].astype(float)]) # m
    UCS = np.concatenate([np.array([0]),profile[:,1].astype(float)])    # MPa
    Em = np.concatenate([np.array([0]),profile[:,2].astype(float)])     # MPa

    if plot_profile == 'Yes':
        # Plot UCS vs z profile for confirmation
        #fig2, ax2 = plt.subplots(1,1)
        plt.plot(UCS,depth,'-',label=r'$S_u$',color='blue')
        plt.legend(loc='lower left')
        plt.xlabel('Uncompressed confined strength (MPa)'), 
        plt.ylabel('Depth below the pile head (m)'), plt.grid(True)
        # Plot mudline/ground surface
        plt.plot([-0.5*max(UCS),max(UCS)],[z0,z0],'--',color='red')
        plt.text(-0.5*max(UCS),0.95*z0,'Mudline',color='red')
        ax = plt.gca(); ax.invert_yaxis(),
        ax.xaxis.tick_top(), ax.xaxis.set_label_position('top')

    # Define interpolation functions
    f_UCS = interp1d(depth, UCS*1e6, kind='linear') # Pa
    f_Em = interp1d(depth, Em*1e6, kind='linear')   # Pa
    
    #var_rock_profile = inspect.currentframe().f_locals

    return f_UCS, f_Em

if __name__ == '__main__':

    profile = np.array([[0.0,  7, 50, 'Name of p-y model'],
                        [1e10,  7, 50, 'Name of p-y model']])
    
    L = 10
    D = 1
    t = 0.05
    H = 3187638.5081103545
    V = 2975548.3405682378
    H = H
    V = V
    '''
    y,z,DQ = py_analysis(profile, L=L, D=D, t=t, E=200e9, 
                        V=V, H=H, H_n=0, M=0, M_n=0)
    print(y[0])
    '''
    ys = []
    ls = np.arange(2, 40, 1)
    ds = np.arange(1,20,1)
    for l in ls:
    #for d in ds:
        y,z,DQ = py_analysis(profile, L=l, D=D, t=t, E=200e9, 
                        V=V, H=H, H_n=0, M=0, M_n=0)
        ys.append(y[0])

    fig, ax = plt.subplots(1,1)
    ax.plot(ls, ys)
    ax.set_xlabel('Pile Length [m]')
    ax.set_ylabel('Deflection (y) [m]')
    ax.text(ls[0], min(ys), f'D={D} m; t={t*1000} mm; H={H/1e6:.1f} MN; V={V/1e6:.1f} MN')
    plt.show()
    

    
    #                   depth  UCS   Em        p-y model    
    profile = np.array([[0.0,  7.0, 150., 'Name of p-y model'],
                        [10.0, 7.0, 150., 'Name of p-y model'],
                        [18.0, 7.0, 150., 'Name of p-y model'],
                        [50.0, 7.0, 150., 'Name of p-y model']])

    #f_UCS, f_Em = rock_profile(profile,plot_profile='No')
       
    #Pile dimensions
    L = 5                  # Pile length (m)
    D = 1                  # Pile diameter (m)
    t = (6.35 + D*20)/1e3  # Pile wall thickness (m), API RP2A-WSD
    E = 200e9              # Elastic modulus of steel (Pa)
    
    #Pile head loads
    z0 = L/10        # Depth of the rock elevation from pile head (m)
    H = 4.4e6        # Horizontal load on pile head (N)
    V = 3.0e6        # Vertical load on pile head (N)
    M = H*z0         # Moment on pile head (N*m)
    H_n = 0          # Horizontal load on pile tip (N)
    M_n = H*z0       # Moment on pile tip (N*m)
    
    y,z,DQ = py_analysis(profile, L=L, D=D, t=t, E=E, 
                      V=V, H=H, H_n=H_n, M=M, M_n=M_n)
    
    Qz = np.sum(DQ)
    
    if y[0] < 0.1*D:
        print('The lateral deflection at pile head is below the maximum allowable')
    else:
        print('PILE NOT COMPLIANT with lateral loading criteria (ASSESS D, L or t)')
    
    if V < Qz:
        print('The vertical load of the pile is below the axial capacity')
    else:
        print('PILE NOT COMPLIANT with vertical capacity criteria (ASSESS D or L)')
        
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