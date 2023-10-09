# class for the collection of power cables of an array

from cable import Cable


class CableSystem():

    def __init__(self, coords, cableProps=None):
        '''Initialize a cable system.
        
        Parameters
        ----------
        coords : 2D array
            Table of x-y coordinates of every floating platform in the array. [m]
        cableProps
            similar to MoorPy System lineProps...        
        '''
        
        self.coords = coords   # list of turbine coordinates
        self.nPtfm = coords.shape[0]
        
        self.cableList = []
        self.connections = []   # list of [iA, iB] indices of platforms that are connected by cables
        self.nc = 0  # number of cables
        
        self.cableTypes = {}  # dictionary containing dictionaries of cable cross section type information
        self.dynamicCableConfigs = {}  # dictionary containing dynamic cable config information (key, name, span [m], props dictionary)
        
        # load cable property scaling coefficients for easy use when creating cable types
        self.cableProps = loadCableProps(cableProps)
        
    # for calculations that involve both cables and seabed, 
    # should they be done in the cable classes (referencing/storing
    # the higher-level site/seabed data), or should it be done at
    # a higher level?
    
    
    def setCableType(self, A, cable_type, source=None, name="", **kwargs):
        '''Add or update a CableSystem cableType using the dictionary-based method.

        Parameters
        ----------
        A : float
            conductor cross-sectional area [mm^2].
        cable_type : string
            string identifier of the cable class to be used (corresponds to 
            a YAML properties entry).
        source : dict or filename (optional)
            YAML file name or dictionary containing line property scaling 
            coefficients. If not provided, whatever has already been loaded 
            into the MoorPy system will be used.
        name : string (optional)
            Identifier for the line type (otherwise will be generated 
            automatically).

        Returns
        -------
        None.
        '''
 
        # compute the actual values for this line type
        if source==None:
            lineType = getCableProps(A, material, lineProps=self.cableProps, name=name, rho=self.rho, g=self.g)  
        else:
            lineType = getCableProps(A, material, source=source, name=name, rho=self.rho, g=self.g)  
        
        lineType.update(kwargs)                      # add any custom arguments provided in the call to the lineType's dictionary
        
        # add the dictionary to the System's lineTypes master dictionary
        if lineType['name'] in self.lineTypes:                                # if there is already a line type with this name
            self.lineTypes[lineType['name']].update(lineType)                 # update the existing dictionary values rather than overwriting with a new dictionary
        else:
            self.lineTypes[lineType['name']] = lineType                       # otherwise save a new entry

        return lineType                              # return the dictionary in case it's useful separately

    
    def addDynamicCableConfig(self, d, name):
        '''Add a dynamic cable configuration based on a dictionary description.'''
        
        self.dynamicCableConfigs[name] = dict(key=name, name=d['name'], 
                                              span=d['span'], props=d)
        
    
    def addCable(self, d):
        '''Convenience function to add a Cable to a cable system.

        Parameters
        ----------
        cableType : string or dict
            string identifier of cableType for this line already added to the system, or dict specifying custom line type.
        d : dict
            dictionary of cable info

        Returns
        -------
        None.
        '''
        
        
        pointA   =d['AttachA']  # index of FOWT it is attached to
        pointB   =d['AttachB']
        DynCableA =d['DynCableA']
        DynCableB =d['DynCableB']
        headingA  =d['headingA']
        headingB  =d['headingB']
        cableType =d['cableType']
        
        
        rtA = self.coords[pointA-1]
        rtB = self.coords[pointB-1]
        
        ...
        
        if not isinstance(cableType, dict):
            if cableType in self.cableTypes:
                cableType = self.cableTypes[cableType]             # in which case that entry will get passed to Line.init
            else:
                raise ValueError(f"The specified cableType name ({cableType}) does not correspond with any cableType stored in this CableSystem")
        
        
        
        # update system-level lists
        self.connections.append([pointA-1, pointB-1])
        self.nc += 1
        
        # create new Cable object
        self.cableList.append( Cable(self, len(self.cableList)+1, rtA, rtB, cableType )
        
        if pointA > 0:
            if pointA <= len(self.nPtfm):
                pass #self.pointList[pointA-1].attachLine(self.lineList[-1].number, 0)
            else:
                raise Exception(f"Provided pointA of {pointA} exceeds number of platforms.")
        if pointB > 0:
            if pointB <= len(self.nPtfm):
                pass #self.pointList[pointB-1].attachLine(self.lineList[-1].number, 1)
            else:
                raise Exception(f"Provided pointB of {pointB} exceeds number of platforms.")
        