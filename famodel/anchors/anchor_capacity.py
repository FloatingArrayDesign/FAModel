"""The anchor capacity calculation 'switchboard' that holds generic
anchor capacity functions and calls the specific calculation functions
from other modules depending on the soil and anchor information."""

import matplotlib.pyplot as plt
import numpy as np

import moorpy.MoorProps as mprop 

from anchors.capacity_plate import getCapacityPlate
from anchors.capacity_suction import getCapacitySuction





def anchorCapacity(anchor, soil, display=0):
    '''Calculate anchor holding capacity based on specified anchor and soil
    information.
    
    Parameters
    ----------
    anchor : dictionary
        anchor description
    soil : dictionary 
        soil description. Can be a keyword ([_/soft/medium/hard] clay, or sand)
        for the level 1 model, or a soilProps dict for the level 2 model.
    model_level : int
        1 or 2.

    Returns
    -------
    UHC: float
        required anchor ultimate holding capacity [kN]
    info: dict
        dictionary with additional information depending on the model.
    '''

    
    if model_level == 1:   # soil keyword indicates level 1 models
    
    
        # calls level 1 anchor capacity function, with anchor/soil types and default assumptions
        uhc, mass, info = mprop.getAnchorMass(uhc_mode=True, mass_int=anchor['mass'], 
                                              anchor=anchor['type'], soil_type=soil['class'], 
                                              method='static', display=0)
        
        #fx, fz = anchor_curves.anchorCapacity(0, 0, 0, anchor=anchor['type'],
        #                                      soil_type=soil['class'], display=display)


    elif model_level==2:  # dict indicates a soilProps dictionary
    
        # >>> we probably need anchor details too then ...
        
        
        # For now the anchor properties get checked in this function 
        # but in the future they coudl be moved to the individual functions.
        
        if anchor['type'] == 'DEA':
            # make curves from 
            pass
        
        elif anchor['type'] == 'SCA':
        
            L     = getFromDict(anchor, 'length')
            D     = getFromDict(anchor, 'diameter', default=L/6)
            thick = getFromDict(anchor, 'thickness', default=L/100)
            F_ang = np.degrees(np.atan2(Fz, Fx))  # load inclination angle [deg]
            
            if soil['class'] == 'clay':
            
                gamma = getFromDict(soil, 'gamma', default=4.7)
                Su0   = getFromDict(soil, 'So0'  , default=2.39)
                k     = getFromDict(soil, 'k'    , default=1.41)
                alpha = getFromDict(soil, 'alpha', default=0.7)
                SF = 2
            
                results = getCapacitySuction(L, L_D_aspect=L/D, D_t_aspect=D/thick, 
                                    A_angle=F_ang, Su0=Su0, k=k, 
                                    Alpha=alpha, gamma=gamma, J=1/SF)
            
            elif soil['class'] == 'sand':
            
                gamma = getFromDict(soil, 'gamma', default=9.0)
                phi   = getFromDict(soil, 'phi'  , default=30)            
                results = getCapacitySuction(L, L_D_aspect=L/D, D_t_aspect=D/thick, 
                                    A_angle=F_ang, gamma=gamma, Phi=phi)
            
            else:
                #raise Exception(f"soil class '{soil.class}' is not supported.")
                pass
                
            
        elif anchor['type'] == 'VLA':
        
            # same plate capacity calc as SEPLA for now - will in future consider angle
        
            A     = getFromDict(anchor, 'area')
            thick = getFromDict(anchor, 'thickness', default=np.sqrt(L)/40)
            H     = getFromDict(anchor, 'embedment')  # embedment depth [m]
            
            if soil['class'] == 'clay':
            
                gamma = getFromDict(soil, 'gamma', default=4.7)
                Su0   = getFromDict(soil, 'So0'  , default=2.39)
                k     = getFromDict(soil, 'k'    , default=1.41)
                
                results = getCapacityPlate(A, B_t_aspect=np.sqrt(L)/thick, 
                                           Hs=H, Bita=30, Los=0.05,
                                           gamma=gamma, So0=So0, k=k)
            else:
                raise Exception("Only clay soil is supported for this anchor type.")
         
        
        elif anchor['type'] == 'SEPLA':
            
            A     = getFromDict(anchor, 'area')
            thick = getFromDict(anchor, 'thickness', default=np.sqrt(L)/40)
            H     = getFromDict(anchor, 'embedment')  # embedment depth [m]
            
            if soil['class'] == 'clay':
            
                gamma = getFromDict(soil, 'gamma', default=4.7)
                Su0   = getFromDict(soil, 'So0'  , default=2.39)
                k     = getFromDict(soil, 'k'    , default=1.41)
                
                results = getCapacityPlate(A, B_t_aspect=np.sqrt(L)/thick, 
                                           Hs=H, Bita=30, Los=0.05,
                                           gamma=gamma, So0=So0, k=k)
            else:
                raise Exception("Only clay soil is supported for this anchor type.")

        else:
            raise Exception(f"Anchor type '{anchor.type}' is not yet supported in hte intermediate anchor model set")
    

        
        print(f"UHC input: fx:{fx} fz:{fz} -- Mass: {mass}, Cost: {cost}")
        info["UHC input"] = fx,fz           #[kN]
        info["Capacity_sf"] = capacity_sf   #[kN]
        info["Mass"] = mass                 #[mT]
        info["Cost"] = cost                 #[$/mT]
        #info["Length"] = L
        info["Area"] = area    

    else:
        raise Exception("Model level must be 1 or 2")


    return capacity, info

