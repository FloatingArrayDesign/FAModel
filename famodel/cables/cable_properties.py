# Functions for loading and using cable property scaling coefficients

import numpy as np
import matplotlib.pyplot as plt
import yaml
from moorpy.helpers import getFromDict  # should adjust to avoid MoorPy dependency


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

    if 'cable_types' in source:
        cableProps = source['cable_types']
    else:
        raise Exception("YAML file or dictionary must have a 'cable_types' field containing the data")

    
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
        output[ctd]['R_A_coefs'  ] = getFromDict(props, 'R_A_coefs'    , shape=1)        
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
        cableProps = loadCableProps(source)
        if not cableProps==None:
            print('Warning: both cableProps and source arguments were passed to getLineProps. cableProps will be ignored.')
        
    # raise an error if the cable_type isn't in the source dictionary
    if not cable_type in cableProps:
        raise ValueError(f'Specified cable cable_type, {cable_type}, is not in the database.')
    
    # calculate the relevant properties for this specific line type
    ctd = cableProps[cable_type]       # shorthand for the sub-dictionary of properties for the cable_type in question    
    
    d = ctd['D_A_coefs'][0]*np.sqrt(A) + ctd['D_A_coefs'][1]  # cable outer diameter [mm]
    
    
    mass = ctd['mass_d'][0]*d**2 + ctd['mass_d'][1]*d + ctd['mass_d'][2] # linear density [kg/m]
    w = (mass - np.pi/4*(d/1000)**2 *rho)*g  # apparent (wet) weight per unit length [N/m]
    
    EA   = ctd['EA_d'][0]*d**2 + ctd['EA_d'][1]*d + ctd['EA_d'][2]  # axial stiffness [N]
    EI   = ctd['EI_d'][0]*d**2 + ctd['EI_d'][1]*d + ctd['EI_d'][2]  # bending stiffness [Nm]
    
    # minimum breaking load in tension [N]
    MBL  = ctd['MBL_D_coefs'][0]*d**2 + ctd['MBL_D_coefs'][1]*d + ctd['MBL_D_coefs'][2]
    
    # minimum bending radius [m]
    MBR  = ctd['MBR_D_coefs'][0]*d + ctd['MBR_D_coefs'][1]
    
    cost = ctd['cost_A'][0]*A + ctd['cost_A'][1]  # cost per unit length of cable [USD/m]
    
    # Electrical properties
    power = ctd['P_A_coefs'][0] * A**(ctd['P_A_coefs'][1])  # cable rated power capacity [MW]
    resistance = ctd['R_A_coefs'][0]/A
    
    # Set up a main identifier for the cable type unless one is provided
    if name=="":
        typestring = f"{cable_type}{A:.0f}mm2"
    else:
        typestring = name
    
    notes = f"made with getCableProps"
    
    # save dictionary (diameter converted to m)
    lineType = dict(name=typestring, d=d/1000, m=mass, EA=EA, EI=EI, w=w,
                    MBL=MBL, MBR=MBR, A=A, power=power, resistance=resistance,
                    cost=cost, notes=notes)
    
    lineType.update(kwargs)   # add any custom arguments provided in the call to the lineType's dictionary
          
    return lineType


def loadBuoyProps(source):
    '''Load buoyancy module property scaling coefficients from a YAML.
    
    Parameters
    ----------
    source : dict or filename
        YAML file name or dictionary containing buoyancy module property 
        scaling coefficients.
    
    Returns
    -------
    dictionary
        BuoyProps dictionary listing each supported bouyancy module type and 
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

    if 'buoyancy_module_types' in source:
        buoyProps = source['buoyancy_module_types']
    else:
        raise Exception("YAML file or dictionary must have a 'buoyancy_module_types' field containing the data")
    
    output = dict()  # output dictionary of default values with loaded coefficients
    
    # combine loaded coefficients and default values into dictionary for each type
    for ctd, props in buoyProps.items():  
        output[ctd] = {}
        output[ctd]['L_BM'   ] = getFromDict(props, 'L_BM'    , shape=0)
        output[ctd]['D_BM'   ] = getFromDict(props, 'D_BM'    , shape=0)
        output[ctd]['density'] = getFromDict(props, 'density' , shape=0)
        output[ctd]['cost_B' ] = getFromDict(props, 'cost_B'  , shape=2)
    
    return output
    
    
def getBuoyProps(V, buoy_type, buoyProps=None, source=None, name="", rho=1025.0, g=9.81, **kwargs):
    '''Compute properties for a given buoyancy module type and size.
    
    Parameters
    ----------
    V : float
        buoyancy module volume [m^3].
    buoy_type : string
        string identifier of the buoy_type type be used.
    buoyProps : dictionary
        A MoorPy buoyProps dictionary data structure containing the property scaling coefficients.
    source : dict or filename (optional)
        YAML file name or dictionary containing cable property scaling coefficients
    name : any dict index (optional)
        Identifier for the line type (otherwise will be generated automatically).
    rho : float (optional)
        Water density used for computing apparent (wet) weight [kg/m^3].
    g : float (optional)
        Gravitational constant used for computing weight [m/s^2].
    '''
    
    if buoyProps==None and source==None:
        raise Exception("Either buoyProps or source keyword arguments must be provided")
    
    # deal with the source (is it a dictionary, or reading in a new yaml?)
    if not source==None:
        buoyProps = loadLineProps(source)
        if not buoyProps==None:
            print('Warning: both buoyProps and source arguments were passed to getLineProps. buoyProps will be ignored.')
        
    # raise an error if the buoy_type isn't in the source dictionary
    if not buoy_type in buoyProps:
        raise ValueError(f'Specified cable buoy_type, {buoy_type}, is not in the database.')
    
    # calculate the relevant properties for this specific line type
    ctd = buoyProps[buoy_type]       # shorthand for the sub-dictionary of properties for the buoy_type in question    
    
    d = ctd['D_BM']*np.cbrt(V)  # buoyancy module outer diameter [m]
    l = ctd['L_BM']*np.cbrt(V)  # buoyancy module length [m]
    
    B = V
    
    mass = ctd['density']*B  # mass [kg]
    w = (mass - B*rho)*g  # apparent (wet) weight [N]
    
    cost = ctd['cost_B']*B   # cost per buoy
    
    # Set up a main identifier for the cable type unless one is provided
    if name=="":
        typestring = f"{buoy_type}{V:.0f}m3"
    else:
        typestring = name
    
    notes = f"made with getBuoyProps"
    
    # save dictionary (diameter converted to m)
    buoyType = dict(name=typestring, d=d/1000, m=mass, EA=EA, EI=EI, w=w,
                    MBL=MBL, MBR=MBR, A=A, power=power, resistance=resistance,
                    cost=cost, notes=notes)
    
    buoyType.update(kwargs)   # add any custom arguments provided in the call to the buoyType's dictionary
          
    return buoyType


if __name__ == '__main__':
    
    cableProps = loadCableProps('CableProps_default.yaml')
    
    As = [95,120, 150, 185, 200, 300, 400, 500, 630, 800]
    
    cableTypes33 = []
    cableTypes66 = []
    cableTypes132 = []
    
    for A in As:
        cableTypes33.append(getCableProps(A, 'dynamic_cable_33', cableProps=cableProps))
        cableTypes66.append(getCableProps(A, 'dynamic_cable_66', cableProps=cableProps))
        cableTypes132.append(getCableProps(A, 'dynamic_cable_132', cableProps=cableProps))
    
    # Print a table of values for a given cable type
    print(list(cableTypes66[0].keys()))
    for types in cableTypes66:
        print(list(types.values()))
    
    # Plot the values
    def plotProps(x, typeLists, labels, parameter):
        
        fig, ax = plt.subplots(1,1)
        
        for i in range(len(typeLists)):
            ax.plot(x, [ type_j[parameter] for type_j in typeLists[i] ], '*-', label = labels[i])
        
        ax.set_ylabel(parameter)
        ax.set_xlabel(r'A (mm$^2$)')
        ax.set_ylim([0, ax.get_ylim()[1]])
        ax.legend()
    
    plotProps(As, [cableTypes33, cableTypes66, cableTypes132], ['33 kV', '66 kV', '132 kV'], 'd')
    plotProps(As, [cableTypes33, cableTypes66, cableTypes132], ['33 kV', '66 kV', '132 kV'], 'm')
    plotProps(As, [cableTypes33, cableTypes66, cableTypes132], ['33 kV', '66 kV', '132 kV'], 'EA')
    plotProps(As, [cableTypes33, cableTypes66, cableTypes132], ['33 kV', '66 kV', '132 kV'], 'EI')
    plotProps(As, [cableTypes33, cableTypes66, cableTypes132], ['33 kV', '66 kV', '132 kV'], 'MBL')
    plotProps(As, [cableTypes33, cableTypes66, cableTypes132], ['33 kV', '66 kV', '132 kV'], 'MBR')
    plotProps(As, [cableTypes33, cableTypes66, cableTypes132], ['33 kV', '66 kV', '132 kV'], 'power')
    plotProps(As, [cableTypes33, cableTypes66, cableTypes132], ['33 kV', '66 kV', '132 kV'], 'resistance')
    
    plt.show()
    
    
    
    
    
    
    
    
    
    