# class for a mooring line

import numpy as np
from copy import deepcopy
from moorpy.subsystem import Subsystem
from moorpy import helpers
from famodel.mooring.connector import Connector, Section
from famodel.famodel_base import Edge, Node
from famodel.helpers import calc_midpoint
from famodel.platform.fairlead import Fairlead

class Mooring(Edge):
    '''
    Class for a floating array mooring line (anchored or shared).
    The idea is to have a common object that can be used for 2D
    layout optimization as well as (optionally) 3D design.
    Work in progress. Eventually will inherit from Edge.
    '''
    
    def __init__(self, dd=None, subsystem=None, anchor=None,
                 rho=1025, g=9.81,id=None,shared=0):
        '''
        Parameters
        ----------
        dd: dictionary
            Design dictionary that contains all information on a mooring line needed to create a MoorPy subsystem
            Layout: {
                     subcomponents: # always starts and ends with connectors even if connector dict is blank
                         {
                             0   
                                 { # connector
                                 type: {m, v, CdA}
                                 }
                             1
                                 { # section
                                  type:
                                      {
                                         name, d_nom, material, d_vol, m, EA, EAd, EAd_Lm, MBL, cost, weight
                                      }
                                  L # length in [m]
                                 }
                         }
                     connectors:
                     [
                        ...
                     ]
                     span
                     zAnchor
                     z_fair
                     rad_fair
                    }
        
        '''    
        Edge.__init__(self, id)  # initialize Edge base class
        # Design description dictionary for this Mooring
        self.dd = dd
        
        # MoorPy subsystem that corresponds to the mooring line
        self.ss = subsystem
        self.ss_mod = None
        # workaround for users who are making mooring objects based on pre-existing subsystems
        if self.ss and not self.dd:
            self.dd = {}
            if not 'zAnchor' in self.dd:
                self.dd['zAnchor'] = -self.ss.depth
            if not 'span' in self.dd:
                self.dd['span'] = self.ss.span
            if not 'rad_fair' in self.dd:
                self.dd['rad_fair'] = self.ss.rad_fair
            if not 'z_fair' in self.dd:
                self.dd['z_fair'] = self.ss.z_fair
                
            self.dd['subcomponents'] = []
            # find the starting point index
            for lp in self.ss.pointList:
                if not any(lp.attachedEndB):
                    # this is the starting point at end A - add this point first
                    self.dd['subcomponents'].append({'CdA':lp.CdA, 'm':lp.m, 'v':lp.v})
                    break
            # now iterate through and add the line and point B
            for s in range(len(self.ss.lineList)):
                # find what entry in the attached list is for a line attached at end A
                endA = [lp.attachedEndB[i] for i in lp.attachedEndB if i==0][0]
                # save the line number attached to the point at its end A
                line_num = lp.attached[endA]
                # pull out the line section and save it to subcomponents
                ls = self.ss.lineList[line_num-1]
                self.dd['subcomponents'].append({'type':ls.type, 'L':ls.L})
                # go through the point list again and pull out the point attached to end B of the line
                for lb in self.ss.pointList:
                    if line_num in lb.attached and lb != lp:
                        lp = lb
                        # save end B point to subcomponents
                        self.dd['subcomponents'].append({'CdA':lp.CdA, 'm':lp.m, 'v':lp.v})
            

        self.parallels = False  # True if there are any parallel sections in the mooring
                

        # let's turn the dd into something that holds subdict objects of connectors and sections
        if self.dd:
            # >>> list of sections ?  And optional of section types (e,g, chian, poly) 
            # and dict of scaling props (would be used by linedesign) ?
            # self.n_sec = len(self.dd['sections'])

            self.i_con = []
            self.i_sec = []
                            
            # # Turn what's in dd and turn it into Sections and Connectors
            # con_i = 0
            # for i, con in enumerate(self.dd['connectors']):
            #     if isinstance(self.dd['connectors'][i],list):
            #         for j,subcon in enumerate(self.dd['connectors'][i]):
            #             if isinstance(self.dd['connectors'][i][j],list):
            #                 if con and 'type' in con:
            #                     Cid = con['type']+str(con_i)
            #                 else:
            #                     Cid = None
            #                 self.addConnector(con, con_i, id=Cid, insert=False)
            #                 con_i += 1
            #     else:            
            #         if con and 'type' in con:
            #             Cid = con['type']+str(con_i)
            #         else:
            #             Cid = None
            #         self.addConnector(con, con_i, id=Cid, insert=False)
            #         con_i += 1
            # sub_i = 0
            # for i, sec in enumerate(self.dd['sections']):
            #     if isinstance(self.dd['sections'][i],list):
            #         for j,subcon in enumerate(self.dd['sections'][i]):
            #             if isinstance(self.dd['sections'][i][j],list):
            #                 self.addSection(sec['L'],sec['type'],sub_i, insert=False)
            #                 sub_i += 1
            #     else:
            #         self.addSection(sec['L'],sec['type'],sub_i, insert=False)
            #         sub_i += 1
                
            # convert subcomponents list into actual objects
            self.convertSubcomponents(self.dd['subcomponents'])
            
            # connect subcomponents
            self.addSubcomponents(self.dd['subcomponents'])
            
            # point dd['subcomponents'] list to self.subcomponents
            self.dd['subcomponents'] = self.subcomponents

        
        
        # relative positions (variables could be renamed)
        self.rad_anch = self.dd['span'] + self.dd['rad_fair']
        self.rad_fair = self.dd['rad_fair']
        self.z_anch   = self.dd['zAnchor']
        self.z_fair   = self.dd['z_fair']
        self.span     = self.dd['span']
        
        # end point absolute coordinates, to be set later
        self.rA = np.array([-self.rad_anch, 0, self.z_anch])
        self.rB = np.array([-self.rad_fair, 0, self.z_fair])
        self.heading = 270  # compass heading from B to A [deg]
        
        self.adjuster = None  # custom function that can adjust the mooring
        
        self.shared = shared # int for if the mooring line is a shared line (1) or part of a shared line (2)
        self.symmetric = False # boolean for if the mooring line is a symmetric shared line
        
        # relevant site info
        self.rho = rho
        self.g = g
        
        # Dictionaries for additional information
        self.envelopes = {}  # 2D motion envelope, buffers, etc.
        self.loads = {}
        self.safety_factors = {} # calculated safety factor
        self.safety_factors_required = {} # minimum allowable safety factor 
        self.reliability = {}
        self.cost = {}
        self.failure_probability = {}
        self.env_impact = {
            "disturbedSeabedArea": 0
        }

        self.raftResults = {}
    
    def update(self, dd=None):
        '''Update the Mooring object based on the current state of the design
        dictionary (self.dd) or, if passed in, a new design dictionary (dd).
        '''
        
        if not dd == None:  # if dd passed in
            self.dd.update(dd)  # move contents of dd into Mooring.dd
            self.convertSubcomponents(dd['subcomponents'])
            self.addSubcomponents(self.dd['subcomponents'])
            
        # Update section lengths and types
        for i in range(len(self.i_sec)):
            sec = self.getSubcomponent(self.i_sec[i])
            
            if self.ss:
                self.ss.lineList[i].setL(sec[i]['L'])
                self.ss.lineTypes[i] = sec[i]['type']
        
        #TODO: update any other things (connectors, positions...)
        
        
    def setSectionLength(self, L, i):
        '''Sets length of section, including in the subsystem if there is
        one.'''
        sec = self.getSubcomponent(self.i_sec[i])
        sec['L'] = L  # set length in dd (which is also Section/subcomponent)
        
        if self.ss:  # is Subsystem exists, adjust length there too
            self.ss.lineList[i].setL(L)
    
    """
    def setSectionDiameter(self, d, i):
        '''Adjusts diameter of section (including in the subsystem if there is
        one, since it should point to the lineType dict in the Mooring.'''
        
        # just adjust the dict? ss.lineTypes should already point to it
        self.sectionType[i].update( setLineType( ... d))
    """
    
    def setSectionType(self, lineType, i):
        '''Sets lineType of section, including in the subsystem 
        if there is one.'''
        sec = self.getSubcomponent(self.i_sec[i])
        # set type dict in dd (which is also Section/subcomponent)
        sec['type'] = lineType  
        
        if self.ss:  # is Subsystem exists, adjust length there too
            self.ss.lineTypes[i] = lineType
    
    
    
    
    def reposition(self, r_center=None, heading=None, project=None, 
                   degrees=False, rad_fair=[], z_fair=[], adjust=True, **kwargs):
        '''Adjusts mooring position based on changed platform location or
        heading. It can call a custom "adjuster" function if one is
        provided. Otherwise it will just update the end positions.
        
        Parameters
        ----------
        r_center : list or nested list
            The x, y, z coordinates of the platform(s) (undisplaced) [m]. If shared mooring, must be a list of lists, with the
            coordinates for each platform connection. In this case, end A platform connection is the first entry and end B 
            platform connection is the second entry. If not given, r_center will be populated from what self is attached to.
        heading : float
            The absolute heading compass direction of the mooring line from end B
            [deg or rad] depending on degrees parameter (True or False). Must account for the platform heading as well.
        project : FAModel Project, optional
            A Project-type object for site-specific information used in custom
            mooring line adjustment functions (mooring.adjuster).
        rad_fair : list, optional
            fairlead radius of node connected on each end of the mooring line (list should be length 2)
            If not provided, the fairlead radius will be determined from the attached nodes' listed
            fairlead radius (or, if it's an anchor, 0)
        z_fair : list, optional
            fairlead depth (relative to platform depth) of node connected on each end of mooring line (list should be length 2)
            If not provided, the fairlead depth will be determined from the attached nodes' listed fairlead depths (or, if its an anchor, 0)
        **kwargs : dict
            Additional arguments passed through to the designer function.
        '''
        
        # Adjust heading if provided
        if not heading == None:
            if degrees:
                self.heading = heading
            else:
                self.heading = np.degrees(heading)

        phi = np.radians(90-self.heading) # heading in x-y radian convention [rad] (CONVERTS FROM COMPASS HEADINGS!)
        # heading 2D unit vector
        u = np.array([np.cos(phi), np.sin(phi)])
        
        if r_center is not None:
            if self.shared == 1:
                r_centerA = np.array(r_center)[0]
                r_centerB = np.array(r_center)[1]
            else:
                r_centerB = np.array(r_center)
        else:
            r_centerA = self.attached_to[0].r
            r_centerB = self.attached_to[1].r
            
        # check if there are fairlead objects attached to end connectors        
        fairs = True if len(self.subcons_B[0].attachments)>1 else False
        # if there is no fairlead object, use traditional method to determine new fairlead location and set it, otherwise end B should be set already
        if not fairs:
            # create fairlead radius list for end A and end B if needed
            if not rad_fair:
                rad_fair = [self.attached_to[x].rFair if (hasattr(self.attached_to[x],'rFair') and self.attached_to[x].rFair) else 0 for x in range(2)]
            # create fairlead depth list for end A and end B if needed
            if not z_fair:
                z_fair = [self.attached_to[x].zFair if (hasattr(self.attached_to[x],'zFair') and self.attached_to[x].zFair) else 0 for x in range(2)]
                
            # Set the updated end B location
            self.setEndPosition(np.hstack([r_centerB[:2] + rad_fair[1]*u, z_fair[1] + r_centerB[2]]), 'b')
            
        
        # Run custom function to update the mooring design (and anchor position)
        if self.adjuster and adjust:
            
            #if i_line is not defined, assumed segment 0 will be adjusted
            if not hasattr(self,'i_line'):
                self.i_line = 0
       
            
            if hasattr(self,'slope'):
                self.adjuster(self, method = 'pretension', r=r_centerB, project=project, target = self.target, i_line = self.i_line, slope = self.slope)
            
            else:
                
                #move anchor based on set spacing then adjust line length
                xy_loc = self.rB[:2] + self.span*u #r_centerB[:2] + (self.span + rad_fair[1])*u
                if project:
                    self.dd['zAnchor'] = -project.getDepthAtLocation(xy_loc[0],xy_loc[1])
                    self.z_anch = self.dd['zAnchor']
                else:
                    print('Warning: depth of mooring line, anchor, and subsystem must be updated manually.')

                self.setEndPosition(np.hstack([self.rB[:2] + self.span*u, self.z_anch]), 'a', sink=True)                
                self.adjuster(self, method = 'horizontal', r=r_centerB, project=project, target = self.target, i_line = self.i_line)
            
        elif self.shared == 1: # set position of end A at platform end A if no fairlead objects
            if not len(self.subcons_A[0].attachments) > 1:
                self.setEndPosition(np.hstack([r_centerA[:2] - rad_fair[0]*u, z_fair[0] + r_centerA[2]]),'a')
        
        else: # otherwise just set the anchor position based on a set spacing (NEED TO UPDATE THE ANCHOR DEPTH AFTER!)
            if not fairs:          
                xy_loc = self.rB[:2] + self.span*u #r_centerB[:2] + (self.span + rad_fair[1])*u
            else:  
                xy_loc = calc_midpoint([sub.r[:2] for sub in self.subcons_B]) + self.span*u
            if project:
                self.dd['zAnchor'] = -project.getDepthAtLocation(xy_loc[0],xy_loc[1])
                self.z_anch = self.dd['zAnchor']
            else:
                print('Warning: depth of mooring line, anchor, and subsystem must be updated manually.')
            self.setEndPosition(np.hstack([xy_loc, self.z_anch]), 'a', sink=True)

        # Update the mooring profile given the repositioned ends
        if self.ss:
            self.ss.staticSolve()
            
        for i,att in enumerate(self.attached_to):
            iend = self.rA if i == 0 else self.rB
            if type(att).__name__ in 'Anchor':
                # this is an anchor, move anchor location & get soil (if possible)
                if project:
                    project.updateAnchor(att) 
                else: 
                    att.r = iend
                if att.mpAnchor:
                    att.mpAnchor.r = att.r
                    
        # reposition the subcomponents           
        self.positionSubcomponents()
            
    
    
    def setEndPosition(self, r, end, sink=False):
        '''Set the position of an end of the mooring.
        
        Parameters
        ----------
        r : list
            Cordinates to set the end at [m].
        end
            Which end of the edge is being positioned, 'a' or 'b'.
        sink : bool
            If true, and if there is a subsystem, the z position will be on the seabed.
        '''
        
        if end in ['a', 'A', 0]:
            self.rA = np.array(r)
            
            if self.ss:
                if self.rA[2]<-self.ss.depth:
                    self.ss.depth = -self.rA[2]
                self.ss.setEndPosition(self.rA, False, sink=sink)
            
        elif end in ['b', 'B', 1]:
            self.rB = np.array(r)
            
            if self.ss:
                if self.rB[2]<-self.ss.depth:
                    self.ss.depth = -self.rB[2]
                self.ss.setEndPosition(self.rB, True, sink=sink)
                
        else:
            raise Exception('End A or B must be specified with either the letter, 0/1, or False/True.')
    
    
    def getCost(self,from_ss=True):
        '''
        Obtain cost of the Mooring object either from the subsystem costs (if from_ss = True)
        or from the design dictionary (which will calculate line and connector costs)

        Parameters
        ----------
        from_ss : bool, optional
            Describe whether to get costs from the subsystem (True) or design dict (False). The default is True.

        Returns
        -------
        Total cost (float)
            Returns the sum of all costs in the design dictionary

        '''
        
        # mooring line cost
        line_cost = 0
        conn_cost = 0
        if self.ss and from_ss:
            for line in self.ss.lineList:
                try:
                    line_cost += line.getCost()
                except:
                    line_cost += 0
                    print('Could not find line cost for',self.id)
            for point in self.ss.pointList:
                try:
                    if self.loads:
                        ccost,self['MBL'],_ = point.getCost_and_MBL(peak_tension=max([val for val in self.loads.values() if isinstance(val,float)]))
                        conn_cost += ccost
                    elif point.cost:
                        conn_cost += point.cost
                except:
                    pass
        else:
            for sub in self.subcomponents:
                if isinstance(sub,Section):
                    if 'cost' in sub['type']:
                        line_cost += sub['type']['cost']*sub['L']
                elif isinstance(sub,Connector):                  
                    if self.loads:
                        conn_cost += sub.getCost(peak_tension=max([val for val in self.loads.values() if isinstance(val,float)]))
                    elif 'cost' in sub:
                        conn_cost += sub.cost
                        
            self.cost['connector'] = conn_cost
        
        self.cost['line'] = line_cost

        
        # sum up the costs in the dictionary and return
        return sum(self.cost.values())
    
    def updateTensions(self, DAF=1):
        ''' Gets tensions from subsystem and updates the max tensions dictionary if it is larger than a previous tension
        '''
        Ts = []
        Tc = []
        # get tensions for each section
        for sec in self.sections():
            if not 'Tmax' in sec.loads:
                sec.loads['Tmax'] = 0
            Tmax = max([abs(sec.mpLine.TA), abs(sec.mpLine.TB)])
            if Tmax*DAF > sec.loads['Tmax']:
                sec.loads['Tmax'] = deepcopy(Tmax)*DAF
            Ts.append(sec.loads['Tmax'])
        for conn in self.connectors():
            if not 'Tmax' in conn.loads:
                conn.loads['Tmax'] = 0
            Tmax = np.linalg.norm(conn.mpConn.getForces())
            if Tmax*DAF > conn.loads['Tmax']:
                conn.loads['Tmax'] = deepcopy(Tmax)*DAF
            Tc.append(conn.loads['Tmax'])
            
        return max(Ts)
    
    def updateSafetyFactors(self,key='tension',load='Tmax', prop='MBL', 
                            sections=True, connectors=True, info={}):
        """Update safety factors for desired factor type, load type, and property
        
        Parameters
        ---------
        key: str/int, optional
            safety_factor dictionary key. Default is 'tension'. 
        load: str, optional
            key in loads dictionary. Default is 'Tmax'
        prop: str, optional
            key in line type dictionary. Default is 'MBL'
        info: str, optional
            information string to add to safety_factors dictionary
            
        Returns
        -------
        Minimum safety factor for the given key across all sections in the mooring line
        """
        
        # get safety factors for each section
        if sections:
            for sec in self.sections():
                if prop in sec['type']:
                    sec.safety_factors[key] = sec['type'][prop]/sec.loads[load]
                sec.safety_factors['info'] = info
        if connectors:
            for con in self.connectors():
                if 'type' in con and prop in con['type']:
                    con.safety_factors[key] = con['type'][prop]/con.loads[load]
                sec.safety_factors['info'] = info

            
            
    
    
    def createSubsystem(self, case=0, pristine=True, dd=None, ms=None):
        ''' Create a subsystem for a line configuration from the design dictionary
        
        Parameters
        ----------
        pristine : bool, optional
            0/False: modified line (marine growth, corrosion, etc) 1/True: pristine line (default)
        case : int, optional
            Selector shared/suspended cases:
                - 0 (default): end A is on the seabed
                - 1: assembly is suspended and end A is at another floating system
                - 2: the assembly is suspended and assumed symmetric, end A is the midpoint
        dd : dict, optional
            Dictionary describing the design
        ms : MoorPy System, optional
            MoorPy system this subsystem is a part of. Necessary if
        '''
        
        # set design dictionary as self.dd if none given, same with connectorList
        if not dd:
            dd = self.dd
        
        # get list of sections and connectors, send in dd in case it is not from self.dd
        secs = self.sections(dd)
        conns = self.connectors(dd)
        
        if self.parallels:  # make parts of a MoorPy system
            
            if not ms:
                raise Exception('A MoorPy system (ms) must be provided for a Mooring with parallel/bridle parts.')
            # Make Points
            # check where the mooring line is connected to anchor, don't make that connector
            for i in [0,-1]: # no bridles for anchors, so shouldn't be multiple anchor connectors
                # if attached to anchor, remove connector from list
                if hasattr(self.attached_to[i], 'mpAnchor'):
                    conns.pop(i)
            # make connector points
            for i,con in enumerate(conns):
                # >>> leah had some checks here that I didn't understand <<<
                con.makeMoorPyConnector(ms)
            # Make Lines
            for sec in secs:
                sec.makeMoorPyLine(ms) # this also will connect the Lines to Points
        
        else:
            ss=Subsystem(mooringSys=ms, depth=-dd['zAnchor'], rho=self.rho, g=self.g, 
                          span=dd['span'], rad_fair=self.rad_fair,
                          z_fair=self.z_fair)#, bathymetry=dict(x=project.grid_x, y=project.grid_y, depth=project.grid_depth))    # don't necessarily need to import anymore
            
            lengths = []
            types = []
            # run through each line section and collect the length and type
            for sec in secs:
                lengths.append(sec['L'])
                types.append(sec['type']) # list of type names
            
            
            # make the lines and set the points 
            ss.makeGeneric(lengths, types, 
                connectors=[conns[ic+1] for ic in range(len(conns)-2)], 
                suspended=case)
            ss.setEndPosition(self.rA,endB=0)
            ss.setEndPosition(self.rB,endB=1)
            
            for i,sec in enumerate(self.sections(dd)):
                sec.mpLine = ss.lineList[i]
            for i,con in enumerate(self.connectors(dd)):
                con.mpConn = ss.pointList[i]
            
            # add in connector info to subsystem points
            if case == 0: # has an anchor - need to ignore connection for first point because anchor is a point itself so can't have a point attached to a point
                startNum = 1
            else: # no anchor - need to include all connections
                startNum = 0 

            for i in range(startNum,len(ss.pointList)): 
                conn = conns[i]                             
                conn.mpConn = ss.pointList[i]
                conn.mpConn.CdA = conns[i]['CdA']
                conn.getProps()
                
            # solve the system
            ss.initialize()
            ss.staticSolve()
            
            
            # add to the parent mooring system if applicable
            if ms:
                ms.lineList.append(ss)
                ss.number = len(ms.lineList)
            
            # save ss to the correct Mooring variable
            if pristine:
                # save to ss
                self.ss = ss
                return(self.ss)
            else:
                # save to modified ss (may have marine growth, corrosion, etc)
                self.ss_mod = ss
                return(self.ss_mod)
    
    
    def positionSubcomponents(self):
        '''Puts any subcomponent connectors/nodes along the mooring in 
        approximate positions relative to the endpoints based on the 
        section lengths.'''
        
        
        # Tabulate the section lengths
        L = []
        n_serial_nodes = 0  # number of serial nodes, including first and last
                
        # ----- First pass, going through each section in series -----
        
        # Figure out lengths
        for item in self.subcomponents:
        
            if isinstance(item, list):  # indicates there are parallel sections here
                pLtot = []  # total length of each parallel string
                for j, parallel in enumerate(item):  # go through each parallel string
                    if isinstance(parallel, list):  # if it's a concatenation of multiple things
                        pLtot.append(0)
                        # go through each item along the parallel path
                        for subitem in parallel:  
                            if isinstance(subitem, Edge):
                                pLtot[j] += subitem['L'] # add the L of each edge
                    else:
                        raise Exception("Unsupported situation ... parallel subitems must be lists")
                
                L.append(min(pLtot))  # save minimum parallel string length
                
            elif isinstance(item, Node):
                n_serial_nodes += 1
                
            elif isinstance(item, Edge):
                L.append(item['L'])  # save length of section
        
        
        # Position nodes along main serial string between rA and rB
        Lsum = np.cumsum(np.array(L))
        j = 0  # index of node along serial string (at A is 0)

        for i, item in enumerate(self.subcomponents):
            if isinstance(item, list) or isinstance(item, Edge):
                j = j+1  # note that we're moving a certain length along the string
            
            # if it's a node, but not the first or last one
            elif isinstance(item, Node) and i > 0 and i < len(self.subcomponents)-1:
                r = self.rA + (self.rB-self.rA)*Lsum[j-1]/Lsum[-1]
                item.setPosition(r)
        
        
        
        # ----- Second pass, to position any nodes that are along parallel sections -----
        for i, item in enumerate(self.subcomponents):
            
            if isinstance(item, list):  # indicates there are parallel sections here
                for j, parallel in enumerate(item):  # go through each parallel string
            
                    # --- go through each item along the parallel path ---
                    # Note: this part repeats some logic, could be replaced with recursive fn
                    L = []
                    n_serial_nodes = 0
                    
                    for subitem in parallel:  
                        if isinstance(item, Node):
                            n_serial_nodes += 1
                            
                        elif isinstance(subitem, Edge):
                            L.append(subitem['L'])  # save length of section
    
                    # --- Figure out the end points of this paralle string ---
                    
                    # if this parallel is on the first section, then it's a bridle at A
                    if i == 0:  
                        if isinstance(parallel[0], Edge):  # if first object is an Edge
                            rA = parallel[0].rA
                        else:
                            rA = parallel[0].r
                    else:
                        rA = self.subcomponents[i-1].r
                    
                    # if this parallel is on the last section, then it's a bridle at B
                    if i == len(self.subcomponents)-1:  
                        if isinstance(parallel[-1], Edge):  # if last object is an Edge
                            rB = parallel[-1].rB
                        else:
                            rB = parallel[-1].r
                    else:
                        rB = self.subcomponents[i+1].r
                    
                    # --- Do the positioning ---
                    Lsum = np.cumsum(np.array(L))
                    
                    for subitem in parallel:
                        if isinstance(subitem, Edge):
                            j = j+1  # note that we're moving a certain length along the string
                        
                        # if it's a node, but not the first or last one
                        elif isinstance(item, Node):
                            if j > 0 and j < n_serial_nodes-1:
                                r = rA + (rB-rA)*Lsum[j]/Lsum[-1]
                                item.setPosition(r)
                            else:
                                print('end of parallel')
                                breakpoint()

    def mirror(self, create_subsystem=True):
        '''
        Mirrors a half design to create a full line (symmetry assumption).
        Works with self.sections() and self.connectors() instead of self.dd.
        '''

        # detach existing subcomponents
        for i in self.i_sec:
            sec = self.getSubcomponent(i)
            sec.detachFrom(end='a')
            sec.detachFrom(end='b')

        # copy base sections and connectors
        sections = deepcopy(self.sections())
        connectors = deepcopy(self.connectors())

        # check middle connector (connector[0] assumed at line center)
        if not connectors[0] or connectors[0].get('m', 0) == 0:
            doubleL = True   # double the middle section length instead
        else:
            doubleL = False  # double connector mass/volume instead

        # prepare mirrored versions
        addSections = deepcopy(sections)
        addConns = deepcopy(connectors)

        # handle the middle connector
        if doubleL:
            # drop first section + connector, double middle section length
            addSections.pop(0)
            addConns.pop(0)
            connectors.pop(0)  # drop first original connector
            sections[0]['L'] *= 2
        else:
            # drop first mirrored connector, double last connector’s props
            addConns.pop(0)
            connectors[0]['m'] *= 2
            connectors[0]['v'] *= 2

        # reverse and extend lists
        sections.reverse()
        connectors.reverse()
        sections.extend(addSections)
        connectors.extend(addConns)

        # update Mooring design dictionary with these sections and connectors



        # -----------------------------
        # build alternating subcomponent list
        # [C0, S0, C1, S1, ..., Sn, Cn+1], where n is the number of sections.
        # -----------------------------
        subs_list = []
        for i in range(len(sections)):
            subs_list.append(connectors[i])
            subs_list.append(sections[i])
        subs_list.append(connectors[-1])

        # convert dicts into Section/Connector objects
        self.dd['subcomponents'] = subs_list
        self.subcomponents = []
        self.i_con = []  # reset indices
        self.i_sec = []        
        self.convertSubcomponents(subs_list)        
        self.addSubcomponents(self.dd['subcomponents'])
        # create subsystem if asked
        if create_subsystem:
            self.createSubsystem(case=1)


    """
    # rough method ideas...maybe not necessary or can use existing dict methods
    def ddApply():
        '''Applies the current contents of self.dd to update all the parts
        of the Mooring object, including the subsystem if enabled.'''
    
    def ddLoad(self, new_dd):
        '''Takes in a new self.dd then applies it?'''
        
        self.dd = new_dd
        self.ddApply()
    
    def ddSave(self):
        '''Takes the current design info stored in the self attributes
        and writes the values into self.dd. Also returns it, for use in
        writing output yamls...'''
        
        ... save things into self.dd...
        
        return self.dd
    """
    
    def addMarineGrowth(self, mgDict):
        '''Re-creates sections part of design dictionary to account for marine 
        growth on the subystem, then calls createSubsystem() to recreate the line

        Parameters
        ----------
        mgDict : dictionary
            Provides marine growth thicknesses and the associated depth ranges
            {
                th : list with 3 values in each entry - thickness, range lower z-cutoff, range higher z-cutoff [m]
                    *In order from sea floor to surface*
                    example, if depth is 200 m: - [0.00,-200,-100]
                                                - [0.05,-100,-80]
                                                - [0.10,-80,-40]
                                                - [0.20,-40,0]
                rho : list of densities for each thickness, or one density for all thicknesses, [kg/m^3] (optional - default is 1325 kg/m^3)
                }

        Returns
        -------
        changePoints : list
            List of point indices in the moorpy subsystem that are at the changeDepth
        changeDepths : list
            List of cutoff depths the changePoints should be located at

        '''
        # set location of reference mooring object
        # if project: # use pristine line
        #     oldLine = project.mooringListPristine[idx]
        # else: # use current mooring object
        #     oldLine = self
        if self.parallels:
            raise Exception('addMarineGrowth not set up to work with parallels at this time')
        # create a reference subsystem if it doesn't already exist
        if not self.ss:
            self.createSubsystem()     
        oldLine = self.ss
        # set up variables
        LTypes = [] # list of line types for new lines (types listed are from reference object)
        Lengths = [] # lengths of each section for new line
        Mats = [] # materials list for new line        
        connList = [] # new list of connectors (need to add empty connector objects in between changeDepths)
        LThick = [] # list of mg thicknesses for new lines
        ln_raw = [] # list of line lengths from rA to current split in line (needed to determine length of new sections when there are multiple splits in one line section)
        # set up variables needed to check before/after of changeDepths
        changePoints = []
        changeDepths = [] # index of list that has the corresponding changeDepth
        
        # set first connector
        connList.append(self.connectors()[0])
        # go through each line section
        for i in range(0,len(oldLine.lineList)):
            slthick = [] # mg thicknesses for the section (if rA is above rB, needs to be flipped before being added to full subsystem list LThick)
            slength = [] # line lengths for the section (if rA is above rB, needs to be flipped before being added to full subsystem list)
            schangeDepth = [] # changeDepths for the section (if rA is above rB, needs to be flipped before being added to full subsystem list)
            # set reference subsystem line section location
            ssLine = oldLine.lineList[i]
            # add line material, type to list
            Mats.append(ssLine.type['material'])
            LTypes.append(ssLine.type['name'])
                       
            # check whether rA is above rB (can happen for sections of shared lines)
            if ssLine.rA[2]>ssLine.rB[2]: # set low and high point locations accordingly
                low = ssLine.rB
                high = ssLine.rA
                flip = 1
            else:
                low = ssLine.rA
                high = ssLine.rB
                flip = 0

            th = mgDict['th'] # set location for ease of coding
            # look up what thickness this line section starts at (if lowest point is not on the sea floor, first segment will have a thickness other than the sea floor thickness)
            rAth = 0 # exit while loop when we find thickness at low
            count1 = 0 # counter
            while rAth==0 and count1 <= len(th):
                if flip:
                    if high[2] <= th[count1][2]:
                        LThick.append(th[count1][0])
                        rAth = 1 # exit while loop
                else:
                    if low[2] <= th[count1][2]:
                        LThick.append(th[count1][0])
                        rAth = 1 # exit while loop
                count1 = count1 + 1 # iterate counter
                
            # determine if this line section will be split
            for j in range(0,len(th)): # go through all changeDepths
                if flip:
                    rs = 2
                else:
                    rs = 1
                if th[j][rs]>low[2] and th[j][rs]<high[2]:
                    # line section will be split - add line type, mg thickness, and material to list
                    LTypes.append(ssLine.type['name'])
                    slthick.append(th[j][0])
                    Mats.append(ssLine.type['material'])
                    # add an empty connector object to list for split location
                    connList.append(Connector(str(len(connList))+'_empty'))
                    changePoints.append(len(connList)-1)
                    schangeDepth.append([j,rs])
                    
                    # get length of line between each node
                    lenseg = ssLine.L/ssLine.nNodes
                    
                    old_line = ssLine.getLineCoords(Time=0) # get the coordinates of the line
                    #find length of each new section by finding node at changeDepth
                    for k in range(0, ssLine.nNodes-1): # go through each node in the line
                        if flip: # need to check the node ahead is <= the changeDepth to see which node is split
                            if old_line[2][k+1]<=th[j][rs] and old_line[2][k]>th[j][rs]:
                                nodeD = k+1 # nodeD is closest node below changeDepth
                                xp = old_line[2][::-1] # flip because np.interp doesn't take
                                yp = old_line[1][::-1]
                                fp = old_line[0][::-1]
                        else:
                            if old_line[2][k]<=th[j][rs] and old_line[2][k+1]>th[j][rs]:
                                nodeD = k # the node right below changeDepth depth
                                xp = old_line[2][:]
                                yp = old_line[1][:]
                                fp = old_line[0][:]
                    
                    # interpolate to find x & y coordinates at chosen depth (since node might not be exactly at the right depth)
                    xChange = float(np.interp(th[j][rs], xp, fp))
                    yChange = float(np.interp(th[j][rs], xp, yp))
                    
                    # get the "raw length" of the new lower line (from lower end of section to split point) - if there is multiple splits in one line section this raw length may not be the actual length of the new line
                    if flip: # node numbers start at end A (top) so need to subtract from total line length
                        # raw length = total line length - nodeD*(segment length) + 3d pythagorean theorem (to find length from nodeD to actual cutoff location)
                        ln_raw.append(ssLine.L - lenseg*nodeD + np.sqrt((xChange-old_line[0][nodeD])**2 + (yChange-old_line[1][nodeD])**2 + (th[j][rs]-old_line[2][nodeD])**2))
                    else:
                        # raw length = nodeD*(segment length) + 3d pythagorean theorem 
                        ln_raw.append(lenseg*nodeD + np.sqrt((xChange-old_line[0][nodeD])**2 + (yChange-old_line[1][nodeD])**2 + (th[j][rs]-old_line[2][nodeD])**2))
                    
                    
                    if len(slthick)>1: # line has multiple cuts (upper cut sections have to count the length only from previous nodeD)
                        slength.append(float(ln_raw[-1]-ln_raw[-2]))
                        
                    else: # first split (raw length is actual length)
                        slength.append(float(ln_raw[-1]))
                
            if slthick: # add the length of the top line (from last split to upper end of section) if there was a split in the line
                slength.append(float(ssLine.L-ln_raw[-1]))
                # if rA above rB, reverse the order of the section-level lists (we gathered info from lowest depths up, but this line segment travels from top to bottom)
                if flip:
                    slength.reverse()
                    slthick.reverse()
                    schangeDepth.reverse()
                # Append section-level lists to the subsystem-level lists
                Lengths.extend(slength)
                LThick.extend(slthick)
                changeDepths.extend(schangeDepth)
            else: # line section was not split, add full line length
                Lengths.append(ssLine.L)
                
            # add connector at end of section to list
            connList.append(self.connectors()[i+1])
                
        # Set up list variables for pristine line info
        EA = []
        m = []
        d_ve_old = []
        cd = []
        cdAx = []
                                            
        # create arrays
        d_nom_old = np.zeros((len(LTypes)))        
        ve_nom_adjust = np.zeros((len(LTypes)))
        mu_mg = np.zeros((len(LTypes)))
        rho_mg = np.ones((len(LTypes)))*1325
        # adjust rho value if alternative provided
        if 'rho' in mgDict:
            if not type(mgDict['rho']) is list:
                # just one density given for all marine growth on the line
                rho_mg = rho_mg*mgDict['rho']/1325
            else: # density given for each thickness of marine growth
                for i in range(0,len(rho_mg)):
                    # look up what thickness number this rho is related to
                    for j in range(0,len(th)):
                        # compare thickness to th list
                        if LThick == th[j][0]:
                            # assign rho_mg based on the rho_mg of the thickness
                            rho_mg[i] = mgDict['rho'][j]                   
                
    
        nd = [] # list of dictionaries for new design dictionary sections part
        
        for j,ltyp in enumerate(LTypes):
            st =  deepcopy(oldLine.lineTypes)
            for k in st:
                if st[k]['name'] == ltyp:
                    ltyp = k
            # add in information for each line type without marine growth
            EA.append(st[ltyp]['EA'])
            m.append(st[ltyp]['m'])
            d_ve_old.append(st[ltyp]['d_vol'])
            # new dictionary for this line type
            nd.append({'type':{}, 'L':{}}) # new design dictionary
            ndt = nd[j]['type']
            
            # load in line props from MoorProps
            opt = helpers.loadLineProps(None)
            
            if 'd_nom' in st[ltyp]:
                d_nom_old[j] = st[ltyp]['d_nom']
                # get ratio between ve and nom diameter normally
                ve_nom_adjust[j] = d_ve_old[j]/d_nom_old[j]
            elif Mats[j] in opt:
                # get ratio between ve and nom diameter from MoorProps yaml                
                ve_nom_adjust[j] = opt[Mats[j]]['dvol_dnom']
            # get cd and cdAx if given, or assign to default value
            if Mats[j] in opt and not 'Cd' in st[ltyp]:
                cd.append(opt[Mats[j]]['Cd'])
            elif 'Cd' in st[ltyp]:
                cd.append(st[ltyp]['Cd'])
            else:
                #print('No Cd given in line type and material not found in MoorProps yaml. Default Cd of 1 will be used.')
                cd.append(2)
            if Mats[j] in opt and 'CdAx' in opt[Mats[j]] and not 'CdAx' in st[ltyp]:
                cdAx.append(opt[Mats[j]]['CdAx'])
            elif 'CdAx' in st[ltyp]:
                cdAx.append(st[ltyp]['CdAx'])
            else:
                #print('No CdAx given in line type and material not found in MoorProps yaml. Default CdAx of 0.5 will be used.')
                cdAx.append(0.5)
            
            if LThick[j] == 0:
                nd[j]['type'] = deepcopy(st[ltyp])
                nd[j]['type']['name'] = j
            else:
                # get mu for material
                if Mats[j] == 'chain' or Mats[j] == 'chain_studlink':
                    mu_mg[j] = 2
                else:
                    mu_mg[j] = 1
                
                # re-form dictionaries with marine growth values            
                # calculate nominal diameter
                d_nom_old[j] = d_ve_old[j]/ve_nom_adjust[j] # m
                
                # calculate new line diameter that includes marine growth
                ndt['d_nom'] = float(d_nom_old[j]+2*LThick[j]) #m
                
                # calculate the new mass per meter including marine growth
                growthMass = np.pi/4*(ndt['d_nom']**2-d_nom_old[j]**2)*rho_mg[j]*mu_mg[j] # marine growth mass
                ndt['m'] =  float(growthMass + m[j]) # kg/m (total mass)
                
                # calculate the submerged weight per meter including marine growth
                ndt['w'] = float(growthMass*(1-self.rho/rho_mg[j])*self.g + (m[j]-np.pi/4*d_ve_old[j]**2*self.rho)*self.g) # N/m
                
                # calculate new volume-equivalent diameter (cannot use regular chain/polyester conversion because of added marine growth)
                ndt['d_vol'] = np.sqrt(4*((ndt['m']*self.g-ndt['w'])/self.rho/self.g)/np.pi)
                
                # calculate new increased drag coefficient from marine growth
                # convert cd to cd for nominal diameter, then multiply by inverse of new ve_nom_adjust (ratio of d_nom with mg to d_ve with mg) to return to cd for volume equivalent diameter
                ndt['Cd'] = float(cd[j]*ve_nom_adjust[j]*(ndt['d_nom']/ndt['d_vol']))
                ndt['CdAx'] = float(cdAx[j]*ve_nom_adjust[j]*(ndt['d_nom']/ndt['d_vol']))
                
                # add line details to dictionary
                ndt['material'] = Mats[j]
                ndt['name'] = str(j)
                if 'MBL' in oldLine.lineTypes[ltyp]:
                    ndt['MBL'] = oldLine.lineTypes[ltyp]['MBL']
                if 'cost' in oldLine.lineTypes[ltyp]:
                    ndt['cost'] = oldLine.lineTypes[ltyp]['cost']
                ndt['EA'] = EA[j]
                if 'EAd' in oldLine.lineTypes[ltyp]:
                    ndt['EAd'] = oldLine.lineTypes[ltyp]['EAd']
            # add lengths                 
            nd[j]['L'] = Lengths[j]
        
        # # overwrite design dictionary with new dictionary
        # self.dd['sections'] = nd
        # # overwrite connectorList with new connectorList
        # self.connectorList = connList
        
        # fill out rest of new design dictionary
        nd1 = deepcopy(self.dd)
        nd1['subcomponents'] = [None]*(len(nd)*2+1)
        for i in range(len(nd)):
            nd1['subcomponents'][2*i] = Connector('C'+str(i),**connList[i])
            nd1['subcomponents'][2*i+1] = Section('S'+str(i),**nd[i])
        nd1['subcomponents'][2*i+2] = Connector('C'+str(i),**connList[i+1])
        
        # call createSubsystem() to make moorpy subsystem with marine growth
        if self.shared:
            self.createSubsystem(case=1,dd=nd1,pristine=0)
        else:
            self.createSubsystem(dd=nd1,pristine=0)
            
        return(changeDepths,changePoints)


    def addCorrosion(self, lineProps, corrosion_mm=None):
        '''
        Calculates MBL of chain line with corrosion included

        Parameters
        ----------
        corrosion_mm : float, optional
            amount of corrosion in mm. If the value is not given, the corrosion rate from the lineProps dictionary will be used with a design life of 28 years.

        '''
        design_life = 28 # years  - changeable later if needed
        from moorpy.helpers import getLineProps

        if self.ss:
            for i, line in enumerate(self.ss.lineList):
                # check if the line type has a corrosion property in its MoorProps instead of checking for material name.
                mat = line.type['material']
                if mat not in lineProps:
                    raise ValueError(f'Line material {mat} not found in lineProps dictionary.')
                else:
                    if lineProps[mat].get('corrosion_rate', False):
                        if not corrosion_mm:
                            corrosion_mm = lineProps[mat]['corrosion_rate'] * design_life  # total corrosion over design life
                        corrosion_m = corrosion_mm / 1000  # convert to m
                        factor = ((line.type['d_nom'] - corrosion_m) / line.type['d_nom'])**2
                        line.type['MBL'] *= factor  # update MBL with corrosion factored in
                        # TODO: should we not also update d_nom to reflect corrosion?
                        
        # Old method
        # for i in self.i_sec:
        #     sec = self.getSubcomponent(i)
        #     if sec['type']['material']=='chain':
        #         MBL_cor = sec['type']['MBL']*( (sec['type']['d_nom']-(corrosion_mm/1000))/sec['type']['d_nom'] )**2  # corroded MBL
        #     else:
        #         MBL_cor = sec['type']['MBL']
        #     sec['type']['MBL'] = MBL_cor

    def addCreep(self, lineProps, creep_percent=None):
        '''
        Elongates the polyester lines (if exists) in the mooring by a certain creep percentage

        Parameters
        ----------
        lineProps : dict
            Dictionary of line properties from MoorProps yaml
        creep_percent : float, optionals
            Percentage of creep elongation to add to polyester lines. If not given, the creep rate from the lineProps dictionary will be used with a design life of 28.
        '''
        design_life = 28 # years  - changeable later if needed

        from moorpy.helpers import getLineProps
        if self.ss:
            for i, line in enumerate(self.ss.lineList):
                # check if the line type has a creep property in its MoorProps instead of checking for material name.
                mat = line.type['material']
                if mat not in lineProps:
                    raise ValueError(f'Line material {mat} not found in lineProps dictionary.')
                else:
                    if lineProps[mat].get('creep_rate', False):
                        if not creep_percent:
                            creep_percent = lineProps[mat]['creep_rate'] * design_life  # total creep over design life
                        L_creep = line.L * (1 + creep_percent)
                        self.setSectionLength(L_creep, i)
                        # Change the diameter size to account for creep thinning
                    
                        d_nom_creep = line.type['d_nom'] / np.sqrt(1 + creep_percent)
                        lineType_creep = getLineProps(d_nom_creep*1e3, mat, lineProps)  # convert to mm for getLineProps
                        line.type = lineType_creep  # update line type with new diameter [not sure if that's what we need to do.]
        else:
            raise ValueError('Mooring subsystem must be created before adding creep.')
    
    def switchStiffnessBounds(self, lineProps, lower=False):
        '''
        Switches the line stiffnesses to either the lower or upper bounds defined in the MoorProps yaml

        Parameters
        ----------
        lineProps : dict
            Dictionary of line properties from MoorProps yaml
        
        lower : bool, optional
            Whether to switch to lower bound (True) or upper bound (False). The default is False.
        '''

        suffix = '_LB' if lower else '_UB'

        from moorpy.helpers import getLineProps
        if self.ss:
            for i, line in enumerate(self.ss.lineList):
                mat = line.type['material']
                mat_suffix = mat + suffix
                if mat_suffix in lineProps:
                    lineType = getLineProps(line.type['d_nom']*1e3, mat_suffix, lineProps) # convert to mm for getLineProps
                    line.type = lineType  # update line type with new stiffness

    def getEnvelope(self,ang_spacing=45,SFs=True):
        '''Computes the motion envelope of the Mooring based on the watch 
        circle(s) of what it's attached to. If those aren't already 
        calculated, this method will call the relevant getWatchCircle method.
        
        Parameters
        ----------
        ang_spacing : int
            Spacing between angles to calculate for watch circle
        SFs : bool
            Whether or not to calculate safety factors etc for the lines in the watch circle function
        
        Returns
        x: list of x-coordinates of motion envelope
        y: list of y-coordinates of motion envelope
        -------
        '''
        
        from famodel.platform.platform import Platform
        import shapely as sh
        
        pts = []  # list of watch circle points (on either end)
        
        # Loop through ends A and B
        for i in range(2):
            att = self.attached_to[i] # the object the end is attached to
            
            # Only consider a watch circle if it's attached to a platform
            if isinstance(att, Platform):
                
                # If no watch circle saved, compute it 
                if not 'mean' in att.envelopes:  
                    att.getWatchCircle(ang_spacing=ang_spacing,SFs=SFs)
                
                # Add each vertex of the watch circle to the list
                for x, y in zip(att.envelopes['mean']['x'], 
                                att.envelopes['mean']['y']):
                    
                    pts.append((x, y))
            
            # Otherwise add a single point (e.g. for an anchor)
            else:
                pts.append((att.r[0], att.r[1]))
        
        # Convert to shapely polygon, then get convex hull
        poly = sh.Polygon(pts)  # make polygon object
        hull = sh.convex_hull(poly)  # convex hull (wraps around points)
        
        # Convert to np array and save in object envelope
        x = np.array(hull.exterior.xy[0])
        y = np.array(hull.exterior.xy[1])
        self.envelopes['mean'] = dict(x=x, y=y, shape=hull)
        
        return(x,y)


    def addSection(self, section_length, section_type, index, id=None, insert=True):
        '''
        Add a section to the design
        
        Parameters
        ----------
        section_length : float
            Length of new section in [m]
        section_type : dict
            Dictionary of section properties
        index : list
            New index of section in the mooring design dictionary sections list
            List of length 1 or 3 depending on if part of a subsection or not
        id : str/int, optional
            Id of section
        '''
        if not id:
            if insert:
                for i in self.i_sec:
                    # update ids of subcomponents after this in series
                    # first check if the i_sec index i is the same level as index to add in
                    # and the final entry in i is greater than the index to add in
                    if len(i)==len(index) and i[-1]>index[-1]:
                        # check if all indices within the index list are less 
                        # than or equal to the i_con index, i, list
                        if np.all([i[j]>=index[j] for j in range(len(i))]) and i[-1]>index[-1]:
                            sec = self.getSubcomponent(i)
                            sec.id = '_'.join(['S',*[str(j) for j in i]])
            # make the id start with S and add each component of the index separated by _
            id='_'.join(['S',*[str(j) for j in index]])
        newsection_dd = {'type':section_type,'L':section_length}
        # create section object
        newsec = Section(id,**newsection_dd)
        if insert:
            if len(index)==1:
                self.dd['subcomponents'].insert(index[0], 
                                                newsec)
            elif len(index)==2:
                self.dd['subcomponents'][index[0]][index[1]].insert(0, 
                                                                    newsec)
            elif len(index)==3:
                self.dd['subcomponents'][index[0]][index[1]].insert(index[2], 
                                                                    newsec)
            else:
                raise Exception('Length of index must be 1 or 3')
        else:
            if len(index)==1:
                self.dd['subcomponents'][index[0]] =  newsec
            elif len(index)==2:
                self.dd['subcomponents'][index[0]][index[1]][0] = newsec
            elif len(index)==3:
                self.dd['subcomponents'][index[0]][index[1]][index[2]] = newsec
            else:
                raise Exception('Length of index must be 1 or 3')
        
        
        return(newsec)
    
    def addConnector(self, conn_dd, index, id=None, insert=True):
        '''
        Add a connector to the design

        Parameters
        ----------
        conn_dd : dict
            Connector design dictionary
        index : int
            New index of connector in the mooring design dictionary subcomponents list
        id : str or int, optional
            ID of new connector
        insert : bool, optional
            Controls whether to insert a connector in the list or replace an entry with a connector

        Returns
        -------
        Connector object
            New connector object added to design dictionary

        '''
        if not id:
            if insert:
                for i in self.i_con:
                    # update ids of subcomponents after this in series
                    # first check if the i_con index i is the same level as index to add in
                    # and the final entry in i is greater than the index to add in
                    if len(i)==len(index) and i[-1]>index[-1]:
                        # check if all indices within the index list are less 
                        # than or equal to the i_con index, i, list
                        if np.all([i[j]>=index[j] for j in range(len(i))]):
                            conn = self.getSubcomponent(i)
                            conn.id = '_'.join(['C',*[str(j) for j in i]])
            # make the id start with C and add each component of the index separated by _
            id = '_'.join(['C',*[str(j) for j in index]])
        # create connector object
        newconn = Connector(id, **conn_dd)
        # insert it in self.dd['subcompoents'] list or replace entry in list as needed
        if insert:
            if len(index)==1:
                self.dd['subcomponents'].insert(index[0], 
                                                newconn)
            elif len(index)==2:
                self.dd['subcomponents'][index[0]][index[1]][0].insert(0,
                                                                       newconn)
            elif len(index)==3:
                self.dd['subcomponents'][index[0]][index[1]].insert(index[2], 
                                                                    newconn)
            else:
                raise Exception('Length of index must be 1 or 3')
        else:
            if len(index)==1:
                self.dd['subcomponents'][index[0]] =  newconn
            elif len(index)==2:
                self.dd['subcomponents'][index[0]][index[1]][0] = newconn
            elif len(index)==3:
                self.dd['subcomponents'][index[0]][index[1]][index[2]] = newconn
            else:
                raise Exception('Length of index must be 1 or 3')
        
        return(newconn)
    
    # def connectSubcomponents(self, subcons=None):
        
    #     # first disconnect any current subcomponents
    #     for ii in self.i_sec:
    #         self.subcomponents[ii].detachFrom('A')
    #         self.subcomponents[ii].detachFrom('B')
            
    #     # # detach end connectors from platforms/anchors just in case
    #     # if len(self.subcomponents)>0:
    #     #     endattsA = [att['obj'] for att in self.subcomponents[0].attachments.values()]
    #     #     endattsB = [att['obj'] for att in self.subcomponents[-1].attachments.values()]
    #     #     for att in endattsA:
    #     #         self.subcomponents[0].detach(att)
    #     #     for att in endattsB:
    #     #         self.subcomponents[-1].detach(att)
            
    #     # Now connect the new set of subcomponents and store them in self(Edge).subcomponents!
    #     if subcons is None:
    #         subcons = []  # temporary list of node-edge-node... to pass to the function
    #         for i in range(self.n_sec):
    #             subcons.append(self.dd['connectors'][i])
    #             subcons.append(self.dd['sections'][i])
    #         subcons.append(self.dd['connectors'][-1])
    #     self.addSubcomponents(subcons)  # Edge method to connect and store em

        
    def convertSubcomponents(self, subs_list):
        '''Create section and connector objects from the subcomponents dicts.
        '''
        # go through each entry in subcomponents list
        for i,sub in enumerate(subs_list):
            # if this entry is a list, go through each entry in that
            if isinstance(sub,list):
                self.parallels = True  # flag there is at least one parallel section
                for j,subsub in enumerate(sub):
                    # if this is a list (3rd level), make sections and connectors from entries
                    if isinstance(subsub, list):
                        for k, subsubsub in enumerate(subsub):
                            if 'L' in subsubsub:
                                id = '_'.join(['S',*[str(l) for l in [i,j,k]]])
                                # this is a section
                                subs_list[i][j][k] = self.addSection(subsubsub['L'], 
                                                                     subsubsub['type'],
                                                                     [i,j,k],
                                                                     id=id,
                                                                     insert=False)
                                self.i_sec.append([i, j, k])
                            else:
                                # this should be a connector (no length provided)
                                id = '_'.join(['C',*[str(l) for l in [i,j,k]]])
                                subs_list[i][j][k] = self.addConnector(subsubsub,
                                                                       [i,j,k],
                                                                       id=id,
                                                                       insert=False)
                                self.i_con.append([i, j, k])
                    else:
                        raise Exception('subcomponent list entry must be length 1 or 3')
            elif 'L' in sub:
                # this is a section
                id = 'S'+str(i)
                subs_list[i] = self.addSection(sub['L'], 
                                               sub['type'], 
                                               [i], 
                                               id=id,
                                               insert=False)
                self.i_sec.append([i])
            else:
                # this is a connector
                id = 'C'+str(i)
                subs_list[i] = self.addConnector(sub, [i], id=id, insert=False)
                self.i_con.append([i])
                
    def sections(self, dd=None):
        '''
        returns list of sections in the mooring
        '''
        secs = []
        # allow option to input dict of subcomponents and pull sections from that
        if dd:
            for sub in dd['subcomponents']:
                if 'L' in sub:                    
                    secs.append(sub)
                elif isinstance(sub, list):
                    for subsub in sub:
                        if isinstance(subsub, list):
                            for sss in subsub:
                                if 'L' in sss:
                                    secs.append(sss)
                        elif 'L' in sss:
                            secs.append(subsub)
        else:
            for i in self.i_sec:
                secs.append(self.getSubcomponent(i))
            
        return secs
            
    def connectors(self, dd=None):
        '''
        returns list of connectors in the mooring
        '''
        conns = []
        # allow option to input dict of subcomponents and pull sections from that
        if dd:
            for sub in dd['subcomponents']:
                if not 'L' in sub and isinstance(sub, dict):                    
                    conns.append(sub)
                elif isinstance(sub, list):
                    for subsub in sub:
                        if isinstance(subsub, list):
                            for sss in subsub:
                                if not 'L' in sss and isinstance(sss, dict):
                                    conns.append(sss)
                        elif not 'L' in sss and isinstance(sss, dict):
                            conns.append(subsub)
        else:
            for i in self.i_con:
                conns.append(self.getSubcomponent(i))
            
        return conns
    
    def fairleads(self, end):
        '''
        returns list of fairleads connected to the mooring
        '''
        fairs = []
        if end in [1, 'b', 'B']:
            for sub in self.subcons_B:
                fairs = [att['obj'] for att in sub.attachments.values() if isinstance(att['obj'],Fairlead)]
        if end in [0, 'a', 'A']:
            for sub in self.subcons_A:
                fairs = [att['obj'] for att in sub.attachments.values() if isinstance(att['obj'],Fairlead)]
            
        return fairs
            
    
    # def convertSubcomponents(self,subs_list, level=0, index=[0]):
    #     ind = index
    #     for i,sub in enumerate(subs_list):
    #         if isinstance(sub, list):
    #             lvl = level+1
    #             if len
    #             ind.append(0)
    #             self.convertSubcomponents(sub,level=lvl, index=ind)
    #         elif 'L' in sub:
    #             # this is a section
    #             id = '_'.join([str(j) for j in ind])
    #             self.addSection(sub['L'], sub['type'], ind, id=id)
    #             ind[level] += 1
    #         else:
    #             self.addConnector(sub, ind)
    #             ind[level] += 1
        
