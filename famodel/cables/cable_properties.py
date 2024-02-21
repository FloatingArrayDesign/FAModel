# Functions for loading and using cable property scaling coefficients

def loadCableProps(source):
    '''Load a set of power cable property scaling coefficients from
    a specified YAML file or passed dictionary. Any coefficients not included
    will take a default value Return a dictionary containing the complete 
    cable property scaling coefficient set to use for any provided cable types.
    
    Parameters
    ----------
    source : dict or filename
        YAML file name or dictionary containing cable property scaling 
        coefficients.
    
    Returns
    -------
    dictionary
        CableProps dictionary listing each supported cable type and 
        subdictionaries of scaling coefficients for each.
    '''

    if type(source) is dict:
        source = source
        
    elif source is None or source=="default":
        import os
        dir = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(dir,"CableProps_default.yaml")) as file:
            source = yaml.load(file, Loader=yaml.FullLoader)
        
    elif type(source) is str:
        with open(source) as file:
            source = yaml.load(file, Loader=yaml.FullLoader)

    else:
        raise Exception("loadCableProps supplied with invalid source")

    if 'cableProps' in source:
        cableProps = source['cableProps']
    else:
        raise Exception("YAML file or dictionary must have a 'cableProps' field containing the data")

    
    output = dict()  # output dictionary of default values with loaded coefficients
    
    # combine loaded coefficients and default values into dictionary for each type
    for ctd, props in cableProps.items():  
        output[ctd] = {}
        output[ctd]['kV'         ] = getFromDict(props, 'kV'           , shape=0)
        output[ctd]['mass_d'     ] = getFromDict(props, 'mass_d'       , shape=3)
        output[ctd]['EA_d'       ] = getFromDict(props, 'EA_d'         , shape=3)
        output[ctd]['EI_d'       ] = getFromDict(props, 'EI_d'         , shape=3)
        output[ctd]['MBL_D_coefs'] = getFromDict(props, 'MBL_D_coefs'  , shape=3)
        output[ctd]['MBR_D_coefs'] = getFromDict(props, 'MBR_D_coefs'  , shape=2)
        output[ctd]['D_A_coefs'  ] = getFromDict(props, 'D_A_coefs'    , shape=2)
        output[ctd]['P_A_coefs'  ] = getFromDict(props, 'P_A_coefs'    , shape=2)        
        output[ctd]['cost_A'     ] = getFromDict(props, 'cost_A'       , shape=2)
    return output


def getCableProps(A, cable_type, cableProps=None, source=None, name="", rho=1025.0, g=9.81, **kwargs):
    '''Sets up a dictionary that represents a cable type based on the 
    specified conductor area and cable_type type. The returned dictionary can serve as
    a MoorPy line type. Data used for determining these properties is a cabeType
    dictionary data structure, created by loadCableProps. This data
    can be passed in via the cableProps parameter, or a new data set can be
    generated based on a YAML filename or dictionary passed in via the source 
    parameter. The cableProps dictionary should be error-checked at creation,
    so it is not error check in this function for efficiency.
        
    Parameters
    ----------
    A : float
        cable conductor cross-sectional area [mm^2].
    cable_type : string
        string identifier of the cable_type type be used.
    cableProps : dictionary
        A MoorPy cableProps dictionary data structure containing the property scaling coefficients.
    source : dict or filename (optional)
        YAML file name or dictionary containing cable property scaling coefficients
    name : any dict index (optional)
        Identifier for the line type (otherwise will be generated automatically).
    rho : float (optional)
        Water density used for computing apparent (wet) weight [kg/m^3].
    g : float (optional)
        Gravitational constant used for computing weight [m/s^2].
    '''
    
    if cableProps==None and source==None:
        raise Exception("Either cableProps or source keyword arguments must be provided")
    
    # deal with the source (is it a dictionary, or reading in a new yaml?)
    if not source==None:
        cableProps = loadLineProps(source)
        if not cableProps==None:
            print('Warning: both cableProps and source arguments were passed to getLineProps. cableProps will be ignored.')
        
    # raise an error if the cable_type isn't in the source dictionary
    if not cable_type in cableProps:
        raise ValueError(f'Specified cable cable_type, {cable_type}, is not in the database.')
    
    # calculate the relevant properties for this specific line type
    ctd = cableProps[cable_type]       # shorthand for the sub-dictionary of properties for the cable_type in question    
    
    d = ctd['D_A_coefs'][0]*np.sqrt(A) + ctd['D_A_coefs'][1]  # cable outer diameter [mm]
    
    power = A**ctd['P_A_coefs'][0] + ctd['P_A_coefs'][1]  # cable rated power capacity [MW]
    
    mass = ctd['mass_d2']*d**2  # linear density [kg/m]
    w = (mass - np.pi/4*(d/1000)**2 *rho)*g  # apparent (wet) weight per unit length [N/m]
    
    EA   = ctd['EA_d'][0]*d**2 + ctd['EA_d'][1]*d + ctd['EA_d'][2]  # axial stiffness [N]
    EI   = ctd['EI_d'][0]*d**2 + ctd['EI_d'][1]*d + ctd['EI_d'][2]  # bending stiffness [Nm]
    
    # minimum breaking load in tension [N]
    MBL  = ctd['MBL_D_coefs'][0]*d**2 + ctd['MBL_D_coefs'][1]*d + ctd['MBL_D_coefs'][2]
    
    # minimum bending radius [m]
    MBR  = ctd['MBR_D_coefs'][0]*d + ctd['MBR_D_coefs']
    
    cost = ctd['cost_A'][0]*A + ctd['cost_A'][1]  # cost per unit length of cable [USD/m]
    
    
    # Set up a main identifier for the cable type unless one is provided
    if name=="":
        typestring = f"{cable_type}{dnommm:.0f}"
    else:
        typestring = name
    
    notes = f"made with getCableProps"
    
    # save dictionary (diameter converted to m)
    lineType = dict(name=typestring, d=d/1000, m=mass, EA=EA, EI=EI, w=w,
                    MBL=MBL, MBR=MBR, A=A, power=power,
                    cost=cost, notes=notes, cable_type=cable_type)
    
    lineType.update(kwargs)   # add any custom arguments provided in the call to the lineType's dictionary
          
    return lineType


def loadBuoyancyModuleProps():
    '''Load buoyancy module property scaling coefficients from a YAML.'''
    pass
    
def getBuayancyModuleProps():
    '''Compute properties for a given buoyancy module type and size.'''
    pass
