# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 14:30:13 2024

@author: lsirkis
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 14:25:53 2024

@author: lsirkis
"""
import yaml
import numpy as np
import matplotlib.pyplot as plt
import moorpy as mp
from moorpy.MoorProps import getLineProps
from moorpy.subsystem import Subsystem
from copy import deepcopy

#load yaml Ontology file
with open('C:/Users/LSIRKIS/Documents/Python Scripts/Ontology/mooringOntology.yaml','r') as stream:
    try:
        ont = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

#reorganize data . . .


#create manual subsystem

#find out how many different mooring configurations are used in each system
lineconfig = set()
jj=len(ont['mooring_systems'])
for i in range(0,jj):
    for j in range(1,len(ont['mooring_systems']['ms'+str(i+1)]['data'])):
        lineconfig.add(ont['mooring_systems']['ms'+str(i+1)]['data'][j][0])
        if ont['mooring_systems']['ms'+str(i+1)]['data'][j][3]>0:
            print('Length adjust is listed as non-zero, but currently there is no functionality for adjusting length. Will be updated soon.')

ss0 = []
ss1 = []
ss2 = []

for i,lc in enumerate(lineconfig):
    lengths = []
    types = []
    ltype = []
    ss0.append(Subsystem(depth=ont['site']['general']['water_depth'], spacing=ont['mooring_line_configs'][lc]['anchoring_radius'], rBfair=[ont['mooring_line_configs']['taut-poly_1']['fairlead_radius'],0,ont['mooring_line_configs'][lc]['fairlead_depth']]))
    for j in range(0,len(ont['mooring_line_configs'][lc]['sections'])):
        ltype.append(ont['mooring_line_configs'][lc]['sections'][j]['type'])
        types.append(str(j)+'_'+ltype[-1])
        if not ltype[-1] in ont['mooring_line_types']:
            print('Mooring line type listed in mooring_line_configs is not found in mooring_line_types')
        ss0[i].setLineType(ont['mooring_line_types'][ltype[-1]]['d_nom']*1000,ont['mooring_line_types'][ltype[-1]]['material'],name=types[-1])
        ss0[i].lineTypes[types[-1]]['d_vol'] = float(ont['mooring_line_types'][ltype[-1]]['d_vol']) 
        ss0[i].lineTypes[types[-1]]['m'] = float(ont['mooring_line_types'][ltype[-1]]['m'])
        ss0[i].lineTypes[types[-1]]['EA'] = float(ont['mooring_line_types'][ltype[-1]]['EA'])
        if 'EAd' in ont['mooring_line_types'][ltype[-1]]: 
            ss0[i].lineTypes[types[-1]]['EAd'] = float(ont['mooring_line_types'][ltype[-1]]['EAd'])
            ss0[i].lineTypes[types[-1]]['EAd_Lm'] = float(ont['mooring_line_types'][ltype[-1]]['EAd_Lm'])
        ss0[i].lineTypes[types[-1]]['MBL'] = float(ont['mooring_line_types'][ltype[-1]]['MBL'])
        ss0[i].lineTypes[types[-1]]['cost'] = float(ont['mooring_line_types'][ltype[-1]]['cost'])
        lengths.append(ont['mooring_line_configs'][lc]['sections'][j]['length'])
    # set up the lines and points and stuff
    ss0[i].makeGeneric(lengths, types)
    ss0[i].setEndPosition([np.cos(np.radians(ont['mooring_systems']['ms'+str(i+1)]['data'][0][1]))*ont['mooring_line_configs'][lc]['anchoring_radius'],np.sin(np.radians(ont['mooring_systems']['ms'+str(i+1)]['data'][0][1]))*ont['mooring_line_configs'][lc]['anchoring_radius'],-ont['site']['general']['water_depth']], endB=0)
    ss0[i].setEndPosition([np.cos(np.radians(ont['mooring_systems']['ms'+str(i+1)]['data'][0][1]))*ont['mooring_line_configs'][lc]['fairlead_radius'],np.sin(np.radians(ont['mooring_systems']['ms'+str(i+1)]['data'][0][1]))*ont['mooring_line_configs'][lc]['fairlead_radius'], ont['mooring_line_configs'][lc]['fairlead_depth']], endB=1)
    ss0[i].staticSolve()
    
    # 3D plot in the Subsystem's local frame
    fig, ax = ss0[i].plot()
    # 2D view of the same
    ss0[i].plot2d()
    #Line-like plot (in global reference frame)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
    ss1.append(deepcopy(ss0[i]))
    ss1[i].setEndPosition([np.cos(np.radians(ont['mooring_systems']['ms'+str(i+1)]['data'][1][1]))*ont['mooring_line_configs'][lc]['anchoring_radius'],np.sin(np.radians(ont['mooring_systems']['ms'+str(i+1)]['data'][1][1]))*ont['mooring_line_configs'][lc]['anchoring_radius'],-ont['site']['general']['water_depth']], endB=0)
    ss1[i].setEndPosition([np.cos(np.radians(ont['mooring_systems']['ms'+str(i+1)]['data'][1][1]))*ont['mooring_line_configs'][lc]['fairlead_radius'],np.sin(np.radians(ont['mooring_systems']['ms'+str(i+1)]['data'][1][1]))*ont['mooring_line_configs'][lc]['fairlead_radius'],ont['mooring_line_configs'][lc]['fairlead_depth']], endB=1)
    ss1[i].staticSolve()
    ss2.append(deepcopy(ss0[i]))
    ss2[i].setEndPosition([np.cos(np.radians(ont['mooring_systems']['ms'+str(i+1)]['data'][2][1]))*ont['mooring_line_configs'][lc]['anchoring_radius'],np.sin(np.radians(ont['mooring_systems']['ms'+str(i+1)]['data'][2][1]))*ont['mooring_line_configs'][lc]['anchoring_radius'],-ont['site']['general']['water_depth']], endB=0)
    ss2[i].setEndPosition([np.cos(np.radians(ont['mooring_systems']['ms'+str(i+1)]['data'][2][1]))*ont['mooring_line_configs'][lc]['fairlead_radius'],np.sin(np.radians(ont['mooring_systems']['ms'+str(i+1)]['data'][2][1]))*ont['mooring_line_configs'][lc]['fairlead_radius'],ont['mooring_line_configs'][lc]['fairlead_depth']], endB=1)
    ss2[i].staticSolve()

    ss0[i].drawLine(0, ax, color='r')
    ss1[i].drawLine(0, ax, color='b')
    ss2[i].drawLine(0, ax, color='g')
