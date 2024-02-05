
import yaml
import numpy as np
import matplotlib.pyplot as plt
from famodel import Mooring

#load yaml Ontology file
with open('mooringOntology.yaml','r') as stream:
    try:
        ont = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

#reorganize data . . .


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
        msNum.append(j+1)#mooring system number
        lNum.append(j)#row number in the mooring system table (essentially mooring line number)
        
        if ms_info[j]['lengthAdjust'] > 0:
            print('Length adjust is listed as non-zero, but currently there is no functionality for adjusting length. Will be updated soon.')


mlist = []  # List to hold created Mooring objects


#loop through all lines
for i,lc in enumerate(lineconfig):
    #set up dictionary of information on the mooring configurations
    m_config = {'sections':{},'connectors':{},'rAnchor':{},'zAnchor':{},'rFair':{},'zFair':{},'EndPositions':{}}
    
    #loop through all line sections of each line
    for j in range(0,len(ont['mooring_line_configs'][lc]['sections'])):
        yamltype = ont['mooring_line_configs'][lc]['sections'][j]['type']#get the name of the line section type from the yaml file
        name = str(j)+'_'+yamltype #assign unique name for line section w/increasing # to keep linetype order correct (you could have 2 sections w/same line type so need unique name)
        
        #set up a sub-dictionary that will contain all the information on the line section type
        m_config['sections'][name] = {}
        m_config['sections'][name]['type'] = {}
        
        if not yamltype in ont['mooring_line_types']:
            print('Mooring line type ',yamltype,' listed in mooring_line_configs is not found in mooring_line_types')
            
        #set properties for line type
        
        ont_ltype = ont['mooring_line_types'][yamltype]  # short handle to condense the next 10 lines
        
        m_config['sections'][name]['type']['name'] = name
        m_config['sections'][name]['type']['d_nom'] =  ont_ltype['d_nom']
        m_config['sections'][name]['type']['material'] = ont_ltype['material']
        m_config['sections'][name]['type']['d_vol'] = float(ont_ltype['d_vol'])
        m_config['sections'][name]['type']['m'] = float(ont_ltype['m'])
        m_config['sections'][name]['type']['EA'] = float(ont_ltype['EA'])
        #need to calculate the submerged weight of the line (not currently available in ontology yaml file)
        m_config['sections'][name]['type']['w'] = (ont_ltype['m']-np.pi/4*ont_ltype['d_vol']**2*1025)*9.81
        
        #add dynamic stretching if there is any
        if 'EAd' in ont['mooring_line_types'][yamltype]: 
            m_config['sections'][name]['type']['EAd'] = float(ont_ltype['EAd'])
            m_config['sections'][name]['type']['EAd_Lm'] = float(ont_ltype['EAd_Lm'])
            
        #set more line section properties
        m_config['sections'][name]['type']['MBL'] = float(ont_ltype['MBL'])
        m_config['sections'][name]['type']['cost'] = float(ont_ltype['cost'])
        m_config['sections'][name]['length'] = float(ont['mooring_line_configs'][lc]['sections'][j]['length'])
    
    #set general information on the whole line (not just a section/line type)
    heading = ont['mooring_systems']['ms'+str(msNum[i])]['data'][lNum[i]][1]  # (to make later code lines shorter and clearer)
    
    m_config['rAnchor'] = ont['mooring_line_configs'][lc]['anchoring_radius']
    m_config['zAnchor'] = -ont['site']['general']['water_depth']
    m_config['zFair'] = ont['mooring_line_configs'][lc]['fairlead_depth']
    m_config['rFair'] = ont['mooring_line_configs'][lc]['fairlead_radius']
    m_config['EndPositions']['endA'] = [np.cos(np.radians(heading))*ont['mooring_line_configs'][lc]['anchoring_radius'],np.sin(np.radians(heading))*ont['mooring_line_configs'][lc]['anchoring_radius'],-float(ont['site']['general']['water_depth'])]
    m_config['EndPositions']['endB'] = [np.cos(np.radians(heading))*ont['mooring_line_configs'][lc]['fairlead_radius'],np.sin(np.radians(heading))*ont['mooring_line_configs'][lc]['fairlead_radius'], float(ont['mooring_line_configs'][lc]['fairlead_depth'])]

    #create mooring class instance for this line
    mlist.append(Mooring(dd = m_config, rA = m_config['EndPositions']['endA'], rB = m_config['EndPositions']['endB'], rad_anch =  m_config['rAnchor'], rad_fair = m_config['rFair'], z_anch = -ont['site']['general']['water_depth'], z_fair = m_config['zFair']))

############## Attempt at creating an array in a moorpy system ###########################
#CURRENTLY DOES NOT WORK WHEN THERE IS MORE THAN ONE PLATFORM/BODY

#platform locations
PFlocs = np.array([[-1200, 0, 1200],[-1200, -1200, -1200]])

import moorpy as mp
ms = mp.System(depth=-m_config['zAnchor'])

for i in range(0,len(PFlocs[0])):#loop through for each body/platform    #int(len(mlist)/3)):
    #add a new body at the correct location
    r6 = [PFlocs[0][i],PFlocs[1][i],0,0,0,0]
    ms.addBody(0, r6,m=1.784e7,v=20206,rM=100,AWP=1011)
    #loop through each line onthe body (for now assuming 3 lines per body)
    for j in range(0,3):
        #create subsystem
        mlist[i*3+j].createSubsystem()
        #set the location of subsystem to make next few lines shorter
        ssloc = mlist[i*3+j].subsystem
        
        #add start and end points adjusted to include location offsets, attach subsystem
        ms.addPoint(1,[ssloc.rA[0]+PFlocs[0][i],ssloc.rA[1]+PFlocs[1][i],ssloc.rA[2]])#anchor
        #add subsystem to mp system line list
        ms.lineList.append(ssloc)
        ssloc.number = i*3+j
        ms.pointList[-1].attachLine(i*3+j,0)
        ms.addPoint(1,[ssloc.rB[0]+PFlocs[0][i],ssloc.rB[1]+PFlocs[1][i],ssloc.rB[2]])#fairlead
        ms.pointList[-1].attachLine(i*3+j,1)
        #attach body to fairlead line
        ms.bodyList[i].attachPoint(len(ms.pointList),ssloc.rB)#attach to fairlead

#initialize, solve equilibrium, and plot the system       
ms.initialize()
ms.solveEquilibrium()
fig,ax = ms.plot()


##########Testing out making a moorpy system in the platform class##########
    # import famodel.platform.platform as pf
    # PF = pf.Platform(r=[0,0],heading=0,mooring_headings=[60,180,300])
    # pf.Platform.mooringSystem(PF, rotateBool=0,mList=mlist[0:3])
    # PF2 = pf.Platform(r=[1200,0],heading=0,mooring_headings=[60,180,300])
    # pf.Platform.mooringSystem(PF2, rotateBool=0,mList=mlist[3:6])
    # PF3 = pf.Platform(r=[0,1200],heading=0,mooring_headings=[60,180,300])
    # pf.Platform.mooringSystem(PF3, rotateBool=0,mList=mlist[6:9])

