# class for wind turbine

import numpy as np

from famodel.seabed.seabed_tools import interpFromGrid


class Turbine():
    '''
    Class for Turbine - used for holding design info and making lookup
    table of thrust force for different conditions. This would typically
    be an entry in a Project turbineTypes list (one per turbine type).
    '''
    
    def __init__(self, dd):
        '''
        Initialize turbine object based on dictionary from ontology or RAFT
        input file.
        
        Parameters
        ----------
        dd : dict
            Dictionary describing the design, in RAFT rotor format.
        '''
        
        # Design description dictionary for this Turbine
        self.dd = deepcopy(dd)

        # Dictionaries for addition information
        self.loads = {}
        self.reliability = {}
        self.cost = {}
    
    
    def makeRotor(self):
        '''
        Create a RAFT Rotor object for the turbine.
        '''
        
        from raft.raft_rotor import Rotor
        
        # Make RAFT Rotor based on dd with no frequencies and rotor index 0
        self.rotor = Rotor(self.dd, [], 0)
        
    
    def calcThrustForces(self):
        '''
        Compute thrust force vector in various cases.
        '''
        
        # grid of pitch and yaw angles to evaluate
        self.grid_pitches = np.arange(0, 5, 1)*np.pi/180
        self.grid_yaws    = np.arange(-40, 40, 5)*np.pi/180
        
        n_p = len(pitches)
        n_y = len(yaws)
        
        self.grid_f6 = np.zeros([n_y, n_p, 6])  # 6 DOF aero load about PRP [N, N-m]
        self.grid_p  = np.zeros([n_y, n_p])  # power [W]
        
        
        for i_p in range(n_p):
            for i_y in range(n_y):
        
                loads, derivs = self.rotor.runCCBlade(U0, ptfm_pitch=self.grid_pitches(i_p), 
                                             yaw_misalign=self.grid_yaws(i_y))
                
                # Rotation matrix from rotor axis orientation to wind inflow direction
                R_inflow = rotationMatrix(0, turbine_tilt, yaw_misalign)
                
                # Set up vectors in axis frame (Assuming CCBlade forces (but not 
                # moments) are relative to inflow direction.
                forces_inflow = np.array([loads["T"][0], loads["Y"][0], loads["Z" ][0]])
                moments_axis = np.array([loads["My"][0], loads["Q"][0], loads["Mz"][0]])
                forces_axis = np.matmul(R_inflow, forces_inflow)
                
                # Rotation matrix from FOWT orientation to rotor axis oriention 
                R_axis = rotationMatrix(0, np.arctan2(self.rotor.axis[2], self.rotor.axis[0]),
                                           np.arctan2(self.rotor.axis[1], self.rotor.axis[0]) ) 
                                           
                # Rotate forces and moments to be relative to FOWT orientations
                
                # >>> also need to translate to be about PRP rather than hub <<<
                
                self.grid_f6[i_y, i_p, :3] = np.matmul(R_axis, forces_axis)
                self.grid_f6[i_y, i_p, 3:] = np.matmul(R_axis, moments_axis)
                
                # Save power too
                self.grid_p[i_y, i_p] = loads["P"]
    
    
    def getForces(self, yaw, pitch=0)
        '''Compute aero force/moment vector at a given yaw misalignment angle
        and platform pitch angle if provided.
        
        Parameters
        ----------
        yaw : float
            yaw misalgment angle [rad].
        pitch : float
            platform pitch angle [rad].
        
        Returns
        -------
        f6 : array
            Vector of aero forces and moments [N, N-m].
        '''
        
        f6, _,_,_,_ = interpFromGrid(yaw, pitch, self.grid_yaw, 
                                        self.grid_pitch, self.grid_f6)
        
        return f6
        