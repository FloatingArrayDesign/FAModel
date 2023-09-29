# class for a power cable in a floating array


class Cable():
    def __init__(self):
        ''' '''
        
        self.L = 0  # total length (to be computed) [m]
        
        
        
        # Initially assume that the static portion of the cable goes straight
        # between the two ends of the dynamic cables.
        


    def refreshCableRoute(self):
        '''Takes established cable route points, considers seabed
        bathymetry and embedment depth, and computes points along
        the path along with its length.'''
        
        self.x = []
        self.y = []
        self.z = []
        
    
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
    