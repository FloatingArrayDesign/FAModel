# ssBool=1
# import yaml
# import numpy as np
# import matplotlib.pyplot as plt
# from famodel import Mooring

from famodel.project import Project
import moorpy as mp
import numpy as np
# create project class instance from yaml file
#Array = Project(file='mooringOntology.yaml')
###########No marine growth####################################
Array = Project(file='C:/Users/LSIRKIS/Documents/FAModel/famodel/mooring/mooringOntology.yaml')
settings = {}
settings["linelabels"] = True
settings["pointlabels"] = True                          
Array.ms.plot( **settings)
# model = Array.array
# # model.analyzeUnloaded()
# # model.solveEigen()
# model.analyzeCases(display=1)
# model.plotResponses()
# model.plot()

# model = Array.array
# # model.analyzeUnloaded()
# # model.solveEigen()
# model.analyzeCases(display=1)
# model.plotResponses()
# model.plot()

#############NON-FLIPPED TURBINE EXAMPLE (clumped weight reduced to 20,000)##########
# mgDict = {'th':[[0.0,-600,-100],[0.05,-100,-40],[0.1,-40,0]],'rho':1325}
# #mgDict = {'th':[[0.0,-200,-100],[0.05,-100,-40],[0.1,-40,0]],'rho':[1330,1330,1330]}
# Array.getMarineGrowth(mgDict,tol=3)
# settings = {}
# settings["linelabels"] = True
# settings["pointlabels"] = True                          
# Array.ms.plot( **settings)

print('Running RAFT')
model = Array.array
# model.analyzeUnloaded()
# model.solveEigen()
model.analyzeCases(display=1)
model.plotResponses()
model.plot()
#####################FLIPPED TURBINE EXAMPLE (clumped weight is original 80,000)#################################
# Array1 = Project(file='../OntologySample600m_FLIPPED_TURBINE.yaml')
# mgDict = {'th':[[0.0,-600,-100],[0.05,-100,-40],[0.1,-40,0]],'rho':1325}
# #mgDict = {'th':[[0.0,-200,-100],[0.05,-100,-40],[0.1,-40,0]],'rho':[1330,1330,1330]}
# Array1.getMarineGrowth(mgDict,tol=3)
# Array1.ms.plot( **settings)
# model = Array1.array
# # model.analyzeUnloaded()
# # model.solveEigen()
# model.analyzeCases(display=1)
# model.plotResponses()
# model.plot()
###############################
# # mgDict = {'changeDepths':{'depth':[-90,-50,0],'th':[.005,.1,.3]},'depthTol':3}

# # Array.getMoorPyArray(mgDict=mgDict,plt=1)

# # model2 = Array.array
# # # model.analyzeUnloaded()
# # # model.solveEigen()
# # model2.analyzeCases(display=1)
# # model2.plotResponses()

# ms00 = mp.System(file='C:/Users/LSIRKIS/Downloads/SharedMooring2 (1) 1.dat')
# print(len(ms00.bodyList))
# for body in ms00.bodyList:
#     body.m = 19911423.956678286
#     body.v = 19480.104108645974
#     body.rCG = np.array([ 1.49820657e-15,  1.49820657e-15, -2.54122031e+00])
#     body.AWP = 446.69520543229874
#     body.rM = np.array([2.24104273e-15, 1.49402849e-15, 1.19971829e+01])
#     body.type = -1
 
# ms00.bodyList[1].setPosition([1600,0,0,0,0,0])
 
# ms00.initialize()
 
# ms00.solveEquilibrium()
# model.ms= ms00
# for i in range(0,len(ms00.bodyList)):
#     model.fowtList[i].body = ms00.bodyList[i]
# model.analyzeCases(display=1)
# model.plotResponses()

# ms = mp.System(depth=600)
# ms.addBody(-1,[0,0,0,0,0,0],m=19911423.956678286,v=19480.104108645974,rCG = np.array([ 1.49820657e-15,  1.49820657e-15, -2.54122031e+00]),AWP = 446.69520543229874,rM = np.array([2.24104273e-15, 1.49402849e-15, 1.19971829e+01]))
# ms.addBody(-1,[1600,0,0,0,0,0],m=19911423.956678286,v=19480.104108645974,rCG = np.array([ 1.49820657e-15,  1.49820657e-15, -2.54122031e+00]),AWP = 446.69520543229874,rM = np.array([2.24104273e-15, 1.49402849e-15, 1.19971829e+01]))

# ss = mp.Subsystem(depth=600,span=1557.57,rBfair=[1557.57,0,-20])

# lengths = []
# types = []
# for i in range(0,3):
#     lengths.append(ms00.lineList[i].L)
#     types.append(ms00.lineList[i].type['name'])
#     ss.lineTypes[types[-1]] = ms00.lineTypes['rope']

# lengths[1] = 1131.37
# ss.makeGeneric(lengths,types,suspended=1)
# ss.setEndPosition([42.43,0,-20],endB=0)
# ss.setEndPosition([1557.57,0,-20],endB=1)
# ss.pointList[1].m = 80000
# ss.pointList[1].v = 0.0
# ss.pointList[2].m = 80000
# ss.pointList[2].v = 0.0

# ms.lineList.append(ss)
# ss.number = len(ms.lineList)
# ms.addPoint(1,ss.rB)
# ms.pointList[-1].attachLine(ss.number,1)
# ms.bodyList[1].attachPoint(len(ms.pointList),[ss.rB[0]-1600,ss.rB[1],ss.rB[2]])
# ms.addPoint(1,ss.rA)
# ms.pointList[-1].attachLine(ss.number,1)
# ms.bodyList[0].attachPoint(len(ms.pointList),ss.rA)

# ms.addLineType('rope',ms00.lineTypes['rope']['d_vol'],ms00.lineTypes['rope']['m'],ms00.lineTypes['rope']['EA'],name='rope')

# j=4
# bdnum = 0
# for i in range(3,6):
#     if i>4:
#         bdnum = 1
#     ms.addLine(ms00.lineList[i].L,ms00.lineList[i].type['name'],ms00.lineList[i].nNodes)
#     ms.addPoint(1,ms00.pointList[j].r)
#     j=j+1
#     ms.pointList[-1].attachLine(len(ms.lineList),1)
#     if bdnum == 0:
#         ms.bodyList[1].attachPoint(len(ms.pointList),ms.pointList[-1].r)
#     else:
#         ms.bodyList[1].attachPoint(len(ms.pointList),[ms.pointList[-1].r[0]-1600,ms.pointList[-1].r[1],ms.pointList[-1].r[2]])
#     ms.addPoint(1,ms00.pointList[j].r)
#     ms.pointList[-1].attachLine(len(ms.lineList),1)
#     j=j+1

# ss.staticSolve()

# for body in ms00.bodyList:
#     body.m = 19911423.956678286
#     body.v = 19480.104108645974
#     body.rCG = np.array([ 1.49820657e-15,  1.49820657e-15, -2.54122031e+00])
#     body.AWP = 446.69520543229874
#     body.rM = np.array([2.24104273e-15, 1.49402849e-15, 1.19971829e+01])
 
# ms00.bodyList[1].setPosition([1600,0,0,0,0,0])
 
# ms00.initialize()
 
# ms00.solveEquilibrium()
# fig, ax = ms00.plot(color='green')



# ###########ms0 = Array.getMoorPyArray(plt = 1)
# model = Array.array
# # model.analyzeUnloaded()
# # model.solveEigen()
# model.analyzeCases(display=1)
# model.plotResponses()

#ms = Array.getMoorPyArray(mgDict = mgDict, plt = 1)
#Array.plot3d(draw_boundary=False,boundary_on_bath=False)

# # create a shared anchor (for visualization purposes)
# import moorpy as mp
# import numpy as np

# ms1 = mp.System(depth=200)
# mp.helpers.loadLineProps(None)
# mp.helpers.getLineProps(180,'polyester',source='C:/Users/LSIRKIS/Documents/MoorPy/moorpy/MoorProps_default.yaml')
# # ms = mp.System(file="C:/Users/LSIRKIS/Downloads/SharedMooring4.dat", depth=600)


# mss = mp.System(depth=600)

# angles = np.radians([60,300])
# rAnchor = 1600
# zFair = -21
# rFair = 20
# lineLength = 1800

# mss.setLineType(120, material='chain', name='chain')
# # Add a free, body at [0,0,0] to the system (including some properties to make it hydrostatically stiff)
# mss.addBody(0, np.zeros(6), m=1e6, v=1e3, rM=100, AWP=1e3)

# for i, angle in enumerate(angles):
#     mss.addPoint(1, [rAnchor*np.cos(angle), rAnchor*np.sin(angle), -600])
#     mss.addPoint(1, [  rFair*np.cos(angle),   rFair*np.sin(angle),  zFair])
    
#     # attach the fairlead Point to the Body (so it's fixed to the Body rather than the ground)
#     mss.bodyList[0].attachPoint(2*i+2, [rFair*np.cos(angle), rFair*np.sin(angle), zFair])     
    
#     # add a Line going between the anchor and fairlead Points
#     mss.addLine(lineLength, 'chain', pointA=2*i+1, pointB=2*i+2)

# i=len(mss.pointList)

# yy = 2*rAnchor*np.sin(np.radians(60))
# # Add a free, body at [0,0,0] to the system (including some properties to make it hydrostatically stiff)
# mss.addBody(0, [0,yy,0,0,0,0], m=1e6, v=1e3, rM=100, AWP=1e3)
# angles = np.radians([60,180,300])
# for j, angle in enumerate(angles[0:2]):
#     mss.addPoint(1, [rAnchor*np.cos(angle), yy+rAnchor*np.sin(angle), -600])
#     mss.addPoint(1, [rFair*np.cos(angle),   yy+rFair*np.sin(angle),  zFair])
    
#     # attach the fairlead Point to the Body (so it's fixed to the Body rather than the ground)
#     mss.bodyList[1].attachPoint(i+2*j+2, [rFair*np.cos(angle), rFair*np.sin(angle), zFair])     
    
#     # add a Line going between the anchor and fairlead Points
#     mss.addLine(lineLength, 'chain', pointA=i+2*j+1, pointB=i+2*j+2)

# mss.addPoint(1, [ rFair*np.cos(angles[2]), yy+rFair*np.sin(angles[2]),  zFair])
# mss.bodyList[1].attachPoint(len(mss.pointList), [mss.pointList[-1].r[0],mss.pointList[-1].r[1]-yy,mss.pointList[-1].r[2]])
# mss.addLine(lineLength, 'chain', pointA=1, pointB=len(mss.pointList))

# # now create a shared mooring
# mss.setLineType(190, material='polyester',name='poly')
# mss.addBody(0, np.array([-1600,0,0,0,0,0]), m=1e6,v=1e3,rM=100,AWP=1e3)
# angs = np.radians([120,240])
# for k, angle in enumerate(angs):
#     mss.addPoint(1, [-1600+rAnchor*np.cos(angle), rAnchor*np.sin(angle), -600])
#     mss.addPoint(1, [-1600+rFair*np.cos(angle),   rFair*np.sin(angle),  zFair])
    
#     # attach the fairlead Point to the Body (so it's fixed to the Body rather than the ground)
#     mss.bodyList[2].attachPoint(len(mss.pointList), [rFair*np.cos(angle), rFair*np.sin(angle), zFair])     
    
#     # add a Line going between the anchor and fairlead Points
#     mss.addLine(lineLength, 'chain', pointA=len(mss.pointList)-1, pointB=len(mss.pointList))


# mss.addPoint(1, [rFair-1600, 0, zFair]) # point B
# # attach point B to body
# mss.bodyList[2].attachPoint(len(mss.pointList), [mss.pointList[-1].r[0]+1600,mss.pointList[-1].r[1],mss.pointList[-1].r[2]])
# mss.addPoint(1,[rFair,0,zFair]) # point A
# if ssBool: # shared line with subsystem
#     from moorpy.subsystem import Subsystem

#     # create subsystem
#     ss=Subsystem(depth=600, spacing=rAnchor, rBfair=[rFair-1600,0,zFair])
#     lengths = []
#     types = []

#     # get length and type
#     lengths.append(1650)
#     types.append('poly')
#     # add subsystem line type
#     ss.lineTypes['poly'] = mss.lineTypes['poly']


#     # make the lines and set the points 
#     ss.makeGeneric(lengths,types,suspended=1)
#     ss.setEndPosition([rFair,0,zFair],endB=0) # end A
#     ss.setEndPosition([rFair-1600,0,zFair],endB=1) # end B

#     # solve the system
#    # ss.staticSolve()

#     mss.lineList.append(ss)
#     ss.number = len(mss.lineList)
#     mss.pointList[-2].attachLine(ss.number,1) # attach end B to line
#     mss.pointList[-1].attachLine(ss.number,0) # attach end A to line

#     # attach point A to body
#     mss.bodyList[0].attachPoint(len(mss.pointList),[mss.pointList[-1].r[0],mss.pointList[-1].r[1],mss.pointList[-1].r[2]])

# else: # shared line without subystem 
#     mss.addLine(1650,'poly',pointA=len(mss.pointList),pointB=len(mss.pointList)-1)
#     # attach point A to body
#     mss.bodyList[0].attachPoint(len(mss.pointList), [mss.pointList[-1].r[0], mss.pointList[-1].r[1], mss.pointList[-1].r[2]])

# # mss.addLine(lineLength, 'poly', pointA=4, pointB=len(mss.pointList))

# mss.initialize()                                             # make sure everything's connected

# mss.solveEquilibrium()                                       # equilibrate
# fig2, ax2 = mss.plot()

# #fill in project class instance information with yaml file
# #Array.load('mooringOntology.yaml')

# # #load yaml Ontology file
# # with open('mooringOntology.yaml','r') as stream:
# #     try:
# #         ont = yaml.safe_load(stream)
# #     except yaml.YAMLError as exc:
# #         print(exc)

# # #reorganize data . . .


# # #find out how many different mooring configurations are used in each system
# # lineconfig = []  # Name of mooring line configuration of every mooring line listed in the yaml
# # msNum = []       # Number of the mooring system (1...) corresponding to every mooring line listed in the yaml
# # lNum = []        # Each mooring line's row number in the mooring system table it appears in
# # count = 0        #keeping track of the times through the outside loop
# # # Loop through each mooring system name and contents
# # for key, val in ont['mooring_systems'].items():  
# #     count = count + 1
# #     # Parse mooring system table into a list of dictionaries
# #     ms_info = [dict(zip(val['keys'], row))  for row in val['data']]
    
# #     for j in range(len(ms_info)):  # go through each row of the mooring system info
# #         lineconfig.append(ms_info[j]['MooringConfigID'])
# #         msNum.append(count)#mooring system number
# #         lNum.append(j)#row number in the mooring system table (essentially mooring line number)
        
# #         if ms_info[j]['lengthAdjust'] > 0:
# #             print('Length adjust is listed as non-zero, but currently there is no functionality for adjusting length. Will be updated soon.')


# # mlist = []  # List to hold created Mooring objects

# # # mSystems = {}
# # # for k, v in ont['mooring_systems'].items():
# # #     mSystems[k] = v
# # # lineConfigs = {}
# # # for k, v in ont['mooring_line_configs'].items():
# # #     lineConfigs[k] = v
# # #     for j in range(0,len(lineConfigs[k]['sections'])):
# # #         for key,val in lineConfigs[k]['sections'][j].items():
# # #             if not val['type'] in mSystems:
# # #                 print('NOOOOOO')

# # #loop through all lines
# # for i,lc in enumerate(lineconfig):
# #     #set up dictionary of information on the mooring configurations
# #     m_config = {'sections':{},'connectors':{},'rAnchor':{},'zAnchor':{},'rFair':{},'zFair':{},'EndPositions':{}}
    
# #     #loop through all line sections of each line
# #     for j in range(0,len(ont['mooring_line_configs'][lc]['sections'])):
# #         yamltype = ont['mooring_line_configs'][lc]['sections'][j]['type']#get the name of the line section type from the yaml file
# #         name = str(j)+'_'+yamltype #assign unique name for line section w/increasing # to keep linetype order correct (you could have 2 sections w/same line type so need unique name)
        
# #         #set up a sub-dictionary that will contain all the information on the line section type
# #         m_config['sections'][name] = {}
# #         m_config['sections'][name]['type'] = {}
        
# #         if not yamltype in ont['mooring_line_types']:
# #             print('Mooring line type ',yamltype,' listed in mooring_line_configs is not found in mooring_line_types')
            
# #         #set properties for line type
        
# #         ont_ltype = ont['mooring_line_types'][yamltype]  # short handle to condense the next 10 lines
        
# #         m_config['sections'][name]['type']['name'] = name
# #         m_config['sections'][name]['type']['d_nom'] =  ont_ltype['d_nom']
# #         m_config['sections'][name]['type']['material'] = ont_ltype['material']
# #         m_config['sections'][name]['type']['d_vol'] = float(ont_ltype['d_vol'])
# #         m_config['sections'][name]['type']['m'] = float(ont_ltype['m'])
# #         m_config['sections'][name]['type']['EA'] = float(ont_ltype['EA'])
# #         #need to calculate the submerged weight of the line (not currently available in ontology yaml file)
# #         m_config['sections'][name]['type']['w'] = (ont_ltype['m']-np.pi/4*ont_ltype['d_vol']**2*1025)*9.81
        
# #         #add dynamic stretching if there is any
# #         if 'EAd' in ont_ltype: 
# #             m_config['sections'][name]['type']['EAd'] = float(ont_ltype['EAd'])
# #             m_config['sections'][name]['type']['EAd_Lm'] = float(ont_ltype['EAd_Lm'])
            
# #         #set more line section properties
# #         m_config['sections'][name]['type']['MBL'] = float(ont_ltype['MBL'])
# #         m_config['sections'][name]['type']['cost'] = float(ont_ltype['cost'])
# #         m_config['sections'][name]['length'] = float(ont['mooring_line_configs'][lc]['sections'][j]['length'])
    
# #     #set general information on the whole line (not just a section/line type)
# #     heading = ont['mooring_systems']['ms'+str(msNum[i])]['data'][lNum[i]][1]  # (to make later code lines shorter and clearer)
    
# #     m_config['rAnchor'] = ont['mooring_line_configs'][lc]['anchoring_radius']
# #     m_config['zAnchor'] = -ont['site']['general']['water_depth']
# #     m_config['zFair'] = ont['mooring_line_configs'][lc]['fairlead_depth']
# #     m_config['rFair'] = ont['mooring_line_configs'][lc]['fairlead_radius']
# #     m_config['EndPositions']['endA'] = [np.cos(np.radians(heading))*ont['mooring_line_configs'][lc]['anchoring_radius'],np.sin(np.radians(heading))*ont['mooring_line_configs'][lc]['anchoring_radius'],-float(ont['site']['general']['water_depth'])]
# #     m_config['EndPositions']['endB'] = [np.cos(np.radians(heading))*ont['mooring_line_configs'][lc]['fairlead_radius'],np.sin(np.radians(heading))*ont['mooring_line_configs'][lc]['fairlead_radius'], float(ont['mooring_line_configs'][lc]['fairlead_depth'])]

# #     #create mooring class instance for this line
# #     mlist.append(Mooring(dd = m_config, rA = m_config['EndPositions']['endA'], rB = m_config['EndPositions']['endB'], rad_anch =  m_config['rAnchor'], rad_fair = m_config['rFair'], z_anch = -ont['site']['general']['water_depth'], z_fair = m_config['zFair']))



# # #platform locations (in future would pull from YAML file)
# # PFlocs = []
# # PF = [ont['array']['data'][0][3],ont['array']['data'][0][4]]
# # for i in range(len(ont['array']['data'])):
# #     PFlocs.append(np.array([int(ont['array']['data'][i][3]),int(ont['array']['data'][i][4])]))#np.array([[-1200,-1200],[0,-1200],[1200,-1200],[-1200,0],[0,0],[1200,0],[-1200,1200],[0,1200],[1200,1200]])

# # from famodel.platform.platform import Platform

# # PFlist = []#list of platform class instances

# # import moorpy as mp
# # ms = mp.System(depth=-m_config['zAnchor'])

# # import time
# # time1=time.time()
# # counter = 0 #keep track of number of times through the loop (since # of lines per turbine might not be consistent)
# # for i in range(0,len(PFlocs)):#loop through for each body/platform    #int(len(mlist)/3)):
# #     #create a platform class instance at the correct location
# #     PFlist.append(Platform(r=PFlocs[i]))
# #     #get number of mooring lines for this platform from YAML file
# #     nLines = len(ont['mooring_systems']['ms'+str(i+1)]['data'])
# #     #add a moorpy body at the correct location
# #     r6 = [PFlocs[i][0],PFlocs[i][1],0,0,0,0]
# #     ms.addBody(0,r6,m=1.784e7,v=20206,rM=100,AWP=1011)
# #     #loop through each line on the body (for now assuming 3 lines per body)
# #     for j in range(0,nLines):
# #         #create subsystem
# #         mlist[counter].createSubsystem()
# #         #set the location of subsystem to make next few lines shorter
# #         ssloc = mlist[counter].subsystem
        
# #         #add start and end points adjusted to include location offsets, attach subsystem
# #         ms.addPoint(1,[ssloc.rA[0]+PFlocs[i][0],ssloc.rA[1]+PFlocs[i][1],ssloc.rA[2]])#anchor
# #         #add subsystem to mp system line list
# #         ms.lineList.append(ssloc)
# #         ssloc.number = counter+1
# #         ms.pointList[-1].attachLine(counter+1,0)
# #         ms.addPoint(1,[ssloc.rB[0]+PFlocs[i][0],ssloc.rB[1]+PFlocs[i][1],ssloc.rB[2]])#fairlead
# #         ms.pointList[-1].attachLine(counter+1,1)
# #         #attach body to fairlead line
# #         ms.bodyList[i].attachPoint(len(ms.pointList),ssloc.rB)#attach to fairlead
# #         #increase counter
# #         counter = counter + 1
# # time2=time.time()        
# # #initialize, solve equilibrium, and plot the system       
# # ms.initialize()
# # time3=time.time()
# # ms.solveEquilibrium()
# # fig,ax = ms.plot()

# # time4=time.time()

# # print('Loop took ',time2-time1)
# # print('intialize ',time3-time2)
# # print('solveEquilibrium took: ',time4-time3)
# # ##########Testing out making a moorpy system in the platform class##########
# #     # import famodel.platform.platform as pf
# #     # PF = pf.Platform(r=[0,0],heading=0,mooring_headings=[60,180,300])
# #     # pf.Platform.mooringSystem(PF, rotateBool=0,mList=mlist[0:3])
# #     # PF2 = pf.Platform(r=[1200,0],heading=0,mooring_headings=[60,180,300])
# #     # pf.Platform.mooringSystem(PF2, rotateBool=0,mList=mlist[3:6])
# #     # PF3 = pf.Platform(r=[0,1200],heading=0,mooring_headings=[60,180,300])
# #     # pf.Platform.mooringSystem(PF3, rotateBool=0,mList=mlist[6:9])

