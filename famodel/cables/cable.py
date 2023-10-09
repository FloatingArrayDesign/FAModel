# class for a power cable in a floating array


class Cable():
    def __init__(self, system, number, rtA, rtB, cableType):
        ''' '
        d : dict
            Dictionary of cable information (see ontology).
        project : CableSystem
            For saving a reference to the CableSystem this Cable belonds to
        '''
        
        self.system = system
        
        self.AttachA   =d['AttachA']
        self.AttachB   =d['AttachB']
        self.DynCableA =self.system.dynamicCableConfigs[d['DynCableA']]
        self.DynCableB =self.system.dynamicCableConfigs[d['DynCableB']]
        self.headingA  =np.radians(d['headingA'])
        self.headingB  =np.radians(d['headingB'])
        self.cableType =d['cableType']
        
        # static cable end points
        self.rA = np.zeros(3)
        self.rB = np.zeros(3)
        
        # cable routing info
        self.x = []
        self.y = []
        self.r = []
        
        self.L = 0  # total length (to be computed) [m]
        
        
    def initializeCableRoute(self):
        # Initially assume that the static portion of the cable goes straight
        # between the two ends of the dynamic cables.
        
        # x y coordinates of the FOWTs that the cable runs between
        rtA = self.system.coords[self.AttachA-1]
        rtB = self.system.coords[self.AttachB-1]
        
        # set cable end points
        self.rA[0] = rtA[0] + self.DynCableA.span*np.cos(self.headingA)  # x
        self.rA[1] = rtA[1] + self.DynCableA.span*np.sin(self.headingA)  # y
        self.rB[0] = rtB[0] + self.DynCableB.span*np.cos(self.headingB)  # x
        self.rB[1] = rtB[1] + self.DynCableB.span*np.sin(self.headingB)  # y
        
        # fill in x, y, r lists for just straight connection between end points
        self.x = np.array([self.rA[0], self.rB[0]])
        self.y = np.array([self.rA[1], self.rB[1]])
        self.r = np.array([0, 0])
    

    def resolveCableRoute(self):
        '''Takes established cable route points, considers seabed
        bathymetry and embedment depth, and computes points along
        the path along with its length.'''
        
        # coordinates for plotting
        self.xs = []
        self.ys = []
        
        # figure out x and y coordinates
        # add points along each stretch
        for i in range(len(self.x))-1:
            n = 3 # number of segments along a stretch
            np.append(self.xs, np.linspace(self.x[i], self.x[i+1], n)
            np.append(self.ys, np.linspace(self.y[i], self.y[i+1], n)
        # add last point
        np.append(self.xs, self.x[-1])
        np.append(self.ys, self.y[-1])
        # >>> Stein to add radii around points <<<
        
        # get z coordinates along seabed
        self.zs = self.system.projectAlongSeabed(self.xs, self.ys)
        
        # calculate cable length (discretized approach for now)
    
    
    def calcCableLength(self):
        '''Calculates a cable's length based on its routing.
        '''

        
        # figure out cable length considering end coordinates and path
        
        self.L = length
        
        return length
    
    
    def checkCableExclusions(self, cable):
        '''Checks whether a cable crosses over any exclusions
        or other out of bounds areas.
        '''

        # select cable
        
        # check its path against any exclusion areas or boundaries
        
        # make a list of any exclusion/nodes that it is too close to
        
        return score, list_of_violations
    