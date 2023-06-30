"""Project class for FAModel, containing information and key methods for
the site information and design information that make up a project."""

import moorpy as mp
from anchors.anchor_capacity import anchorCapacity



class Project():
    '''
    The overall object that defines a floating array for analysis and design 
    purposes. Its main model component is the RAFT model but it also includes
    more general information such as seabed info, metocean info, and...
    
    Ideally, this class can function even if RAFT and MoorPy are not loaded,
    at least for some high-level site processing functions.
    
    '''
    
    def __init__(self, lat=0, lon=0, depth=100):
    
    
        #self.array = RAFT_model
        
        self.lattitude  = lat  # lattitude of site reference point [deg]
        self.longitude  = lon  # longitude of site reference point [deg]
        
        self.grid_x      = np.array([0])
        self.grid_y      = np.array([0])
        self.grid_depth  = np.array([[depth]])  # depth at each grid point
        
        self.seabed_type = 'clay'  # switch of which soil property set to use (TBD)
        
        # soil parameters at each grid point
        self.soil_class = [["none"]]       # soil classification name (TBD)
        self.soil_gamma = np.zeros((1,1))  # soil effective unit weight [kPa] (all soils)
        self.soil_Su0   = np.zeros((1,1))  # undrained shear strength at mudline [kPa] (clay soils)
        self.soil_K     = np.zeros((1,1))  # undrained shear strength gradient [kPa/m] (clay soils)
        self.soil_alpha = np.zeros((1,1))  # soil skin friction coefficient [-] (clay soils)
        self.soil_phi   = np.zeros((1,1))  # angle of internal friction [deg] (sand soils)

        # project boundaries
        self.boundary_Xs = np.zeros(0)
        self.boundary_Ys = np.zeros(0)

    #Component property set processing functions
    # ----- These functions set up component properties -----


    # ----- Site conditions processing functions -----

    def setGrid(self, xCoords, yCoords):
        '''
        Set up the rectangular grid over which site or seabed
        data will be saved and worked with. Directions x and y are 
        generally assumed to be aligned with the East and North 
        directions, respectively, at the array reference point.
        
        Parameters
        ----------        
        xCoords : float array
            x coordinates relative to array reference point [m]
        yCoords : float array
            x coordinates relative to array reference point [m]
        '''
        
        self.grid_x = np.array(xCoords)
        self.grid_y = np.array(yCoords)
        
        #TODO: add check for existing seabed data. If present, convert or raise warning <<<
        
        

    def loadBathymetry(self, filename):
        '''
        Load bathymetry information from an input file (format TBD), convert to
        a rectangular grid, and save the grid to the floating array object (TBD).
        
        Paramaters
        ----------
        filename : path
            path/name of file containing bathymetry data (format TBD)
        '''
        
        # load data from file
        Xs, Yx, Zs = processASC(filename, lat, lon)
        
        # interpolate onto grid defined by grid_x, grid_y
        
        # save in object
        
        
        # also save in RAFT, in its MoorPy System(s)

        pass
        
        
    def loadSoil(self, filename):
        '''
        Load geoetechnical information from an input file (format TBD), convert to
        a rectangular grid, and save the grid to the floating array object (TBD).
        
        The parameters contained in the input file should be:
        k0
        k
        unit weight
        friction angle
        soil class
        ...        
        
        Paramaters
        ----------
        filename : path
            path/name of file containing soil data (format TBD)
        '''
        
        # load data from file
        
        # interpolate onto grid defined by grid_x, grid_y
        
        # save in object
        
        pass
        

    def getSoilAtLocation(self, x, y):
        '''
        Interpolate soil properties at specified location from the soil
        properties grid and return a dictionary of soil properties that
        can be used in anchor capacity calculations.
        
        Parameters
        ----------        
        x : float
            x coordinate in array reference frame [m].        
        y : float
            y coordinate in array reference frame [m].

        Returns
        -------            
        soilProps : dictionary
            Dictionary of standard MoorPy soil properties.
        '''
        
        soilProps = {}
        

        if self.seabed_type == 'clay':
            
            soilProps['class'] = 'clay'
            soilProps['gamma'] = interp2d(x, y, self.seabed_x, self.seabed_y, self.soil_gamma)
            soilProps['Su0'  ] = interp2d(x, y, self.seabed_x, self.seabed_y, self.soil_Su0  )
            soilProps['k'    ] = interp2d(x, y, self.seabed_x, self.seabed_y, self.soil_k    )
            soilProps['alpha'] = interp2d(x, y, self.seabed_x, self.seabed_y, self.soil_alpha)
            soilProps['phi'  ] = None
        
        elif self.seabed_type == 'sand':
            soilProps['class'] = 'sand'
            soilProps['gamma'] = interp2d(x, y, self.seabed_x, self.seabed_y, self.soil_gamma)
            soilProps['Su0'  ] = None
            soilProps['k'    ] = None
            soilProps['alpha'] = None
            soilProps['phi'  ] = interp2d(x, y, self.seabed_x, self.seabed_y, self.soil_phi  )
            
            # note: for sand, can assume homogeneous angle of internal fricton
        else:
            raise ValueError(f"Unsupported seabed type '{self.seabed_type}'.")
            
        return soilProps
        

    # ----- Anchor capacity calculation functions -----


    def calcAnchorCapacity(self, anchor):
        '''Compute holding capacity of a given anchor based on the soil
        info at its position. The anchor object's anchor properties and
        location will be used to determine the holding capacity, which
        will be saved to the anchor object.
        
        Parameters
        ----------
        anchor : MoorPy Anchor object (derived from Point)
            The anchor object in question.
        '''

        # interpolate soil properties/class based on anchor position
        anchor.soilProps = getSoilAtLocation(anchor.r[0], anchor.r[1])
        
        # apply anchor capacity model
        ...(call to appropriate model)
            capacity, info = anchorCapacity(anchor, 
        
        # save anchor capacity to anchor
        anchor.capacity = capacity
        
        # also return it
        return capacity
    

'''
Other future items:
Cost calc functions
System Reliability/failure analysis functions
Full scenario visualization functions
Load case setup and constraint eval?
'''
