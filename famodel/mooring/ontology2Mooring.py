
import yaml
import numpy as np
import matplotlib.pyplot as plt
import moorpy as mp
from moorpy.MoorProps import getLineProps
from moorpy.subsystem import Subsystem
from famodel import Mooring

from copy import deepcopy

#load yaml Ontology file
with open('mooringOntology.yaml','r') as stream:
    try:
        ont = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

#reorganize data . . .


#create manual subsystem


#find out how many different mooring configurations are used in each system
lineconfig = []  # Name of mooring line configuration of every mooring line listed in the yaml
msNum = []       # Number of the mooring system (1...) corresponding to every mooring line listed in the yaml
lNum = []        # Each mooring line's row number in the mooring system table it appears in

# Loop through each mooring system name and contents
for key, val in ont['mooring_systems'].items():  

    # Parse mooring system table into a list of dictionaries
    ms_info = [dict(zip(val['keys'], row))  for row in val['data']]
    
    for j in range(len(ms_info)):  # go through each row of the mooring system info
        lineconfig.append(ms_info[j]['MooringConfigID'])
        msNum.append(i+1)#mooring system number
        lNum.append(j)#row number in the mooring system table (essentially mooring line number)
        
        if ms_info[j]['lengthAdjust'] > 0:
            print('Length adjust is listed as non-zero, but currently there is no functionality for adjusting length. Will be updated soon.')


mlist = []  # List to hold created Mooring objects

#loop through all lines
for i,lc in enumerate(lineconfig):
    #set up dictionary of information on the mooring configurations
    m_config = {'line_types':{},'rAnchor':{},'zAnchor':{},'rFair':{},'zFair':{},'EndPositions':{}}
    
    #loop through all line sections of each line
    for j in range(0,len(ont['mooring_line_configs'][lc]['sections'])):
        yamltype = ont['mooring_line_configs'][lc]['sections'][j]['type']#get the name of the line section type from the yaml file
        ltype = str(j)+'_'+yamltype #assign unique name for line section w/increasing # to keep linetype order correct (you could have 2 sections w/same line type so need unique name)
        
        #set up a sub-dictionary that will contain all the information on the line section type
        m_config['line_types'][ltype] = {}
        
        if not yamltype in ont['mooring_line_types']:
            print('Mooring line type ',yamltype,' listed in mooring_line_configs is not found in mooring_line_types')
            
        #set properties for line type
        
        ont_ltype = ont['mooring_line_types'][yamltype]  # short handle to condense the next 10 lines
        
        m_config['line_types'][ltype]['d_nom'] =  ont['mooring_line_types'][yamltype]['d_nom']
        m_config['line_types'][ltype]['material'] = ont['mooring_line_types'][yamltype]['material']
        m_config['line_types'][ltype]['d_vol'] = float(ont['mooring_line_types'][yamltype]['d_vol'])
        m_config['line_types'][ltype]['m'] = float(ont['mooring_line_types'][yamltype]['m'])
        m_config['line_types'][ltype]['EA'] = float(ont['mooring_line_types'][yamltype]['EA'])
        
        #add dynamic stretching if there is any
        if 'EAd' in ont['mooring_line_types'][yamltype]: 
            m_config['line_types'][ltype]['EAd'] = float(ont['mooring_line_types'][yamltype]['EAd'])
            m_config['line_types'][ltype]['EAd_Lm'] = float(ont['mooring_line_types'][yamltype]['EAd_Lm'])
        #set more properties
        m_config['line_types'][ltype]['MBL'] = float(ont['mooring_line_types'][yamltype]['MBL'])
        m_config['line_types'][ltype]['cost'] = float(ont['mooring_line_types'][yamltype]['cost'])
        m_config['line_types'][ltype]['length'] = float(ont['mooring_line_configs'][lc]['sections'][j]['length'])
    
    #set general information on the whole line (not just a section/line type)
    heading = ont['mooring_systems']['ms'+str(msNum[i])]['data'][lNum[i]][1]  # (to make later code lines shorter and clearer)
    
    m_config['rAnchor'] = ont['mooring_line_configs'][lc]['anchoring_radius']
    m_config['zAnchor'] = -ont['site']['general']['water_depth']
    m_config['zFair'] = ont['mooring_line_configs'][lc]['fairlead_depth']
    m_config['rFair'] = ont['mooring_line_configs'][lc]['fairlead_radius']
    m_config['EndPositions']['endA'] = [np.cos(np.radians(ont['mooring_systems']['ms'+str(msNum[i])]['data'][lNum[i]][1]))*ont['mooring_line_configs'][lc]['anchoring_radius'],np.sin(np.radians(ont['mooring_systems']['ms'+str(msNum[i])]['data'][lNum[i]][1]))*ont['mooring_line_configs'][lc]['anchoring_radius'],-float(ont['site']['general']['water_depth'])]
    m_config['EndPositions']['endB'] = [np.cos(np.radians(ont['mooring_systems']['ms'+str(msNum[i])]['data'][lNum[i]][1]))*ont['mooring_line_configs'][lc]['fairlead_radius'],np.sin(np.radians(ont['mooring_systems']['ms'+str(msNum[i])]['data'][lNum[i]][1]))*ont['mooring_line_configs'][lc]['fairlead_radius'], float(ont['mooring_line_configs'][lc]['fairlead_depth'])]

    #create mooring class instance for this line
    mlist.append(Mooring(dd = m_config, rA = m_config['EndPositions']['endA'], rB = m_config['EndPositions']['endB'], rad_anch =  m_config['rAnchor'], rad_fair = m_config['rFair'], z_anch = -ont['site']['general']['water_depth'], z_fair = m_config['zFair']))
