# class for a power cable in a floating array

import numpy as np
from copy import deepcopy
from moorpy.subsystem import Subsystem
from moorpy import helpers
import shapely as sh

from famodel.cables.dynamic_cable import DynamicCable
from famodel.cables.static_cable import StaticCable
from famodel.cables.components import Joint, Jtube
from famodel.famodel_base import Edge


class Cable(Edge):
    '''Class for an entire subsea cable that transmits power between two 
    turbines or between a turbine and a substation. This can include a
    dynamic cable section (at each end), a static cable section, and
    subsea joints that connect them. Possible scenarios:
    DynamicCable - Joint - StaticCable - Joint - DynamicCable
    DynamicCable - Joint - StaticCable - Joint
    DynamicCable - Joint
                   Joint - StaticCable - Joint - DynamicCable
                                         Joint - DynamicCable
    DynamicCable (fully suspended)
    '''
    
    def __init__(self, id, d=None):
        ''' '
        d : dict
            Dictionary of cable information (see ontology).
        '''
        
        Edge.__init__(self, id)  # initialize Edge base class
        
        self.dd = {'joints':[],'cables':[]}  # save design dictionary
        self.n_sec = len(d['cables'])
        if 'upstream_turb_count' in d:
            self.upstream_turb_count = d['upstream_turb_count'] # number of upstream turbines from the cable
        
        # Turn what's in dd and turn it into Sections and Connectors (if there is more than one section)
        if len(d['cables'])>1:
            for i, joi in enumerate(d['joints']):
                if joi and 'type' in joi:
                    Jid = id+'_'+d['joints'][i]['type']+str(i)
                else:
                    Jid = id+'_'+str(i)
                self.dd['joints'].append(Joint(Jid, **d['joints'][i]))
            
            for i, sec in enumerate(d['cables']):
                Cid = id+'_'+sec['cable_type']['name']+str(i)
                if sec['type'].upper() == 'STATIC':
                    if 'routing' in sec:
                        d['cables'][i]['routing'] = sec['routing']
                    if 'burial' in sec:
                        d['cables'][i]['burial'] = sec['burial']
                    self.dd['cables'].append(StaticCable(Cid, dd=d['cables'][i], **d['cables'][i]))
                else:
                    if 'routing' in sec:
                        d['cables'][i]['routing'] = sec['routing']
                    self.dd['cables'].append(DynamicCable(Cid, dd=d['cables'][i], **d['cables'][i]))
            
            # Connect them and store them in self(Edge).subcomponents!
            subcons = []  # list of node-edge-node... to pass to the function)
            for i in range(0,self.n_sec-1):
                subcons.append(self.dd['cables'][i])
                subcons.append(self.dd['joints'][i])
            subcons.append(self.dd['cables'][-1])
            self.addSubcomponents(subcons)  # Edge method to connect and store em
            # Indices of connectors and sections in self.subcomponents list
            self.i_con = list(range(0, 2*self.n_sec+1, 2))
            self.i_sec = list(range(1, 2*self.n_sec+1, 2))
        
        else:
            # just create the singular cable object as a dynamic cable
            Cid = id+'_'+d['cables'][0]['cable_type']['name']+str(0)
            self.dd['cables'].append(DynamicCable(Cid, dd=d['cables'][0],**d['cables'][0]))
            self.addSubcomponents([self.dd['cables'][0]])
        
        # add overall routing to design dictionary
        if 'routing' in d:
            self.dd['routing'] = d['routing']
            
        # cost dictionary
        self.cost = {}
        # add any connector costs to cost dictionary - TO BE MOVED INTO DYNAMIC CABLE APPENDAGES SECTION 
        if 'connector_cost' in d:
            self.cost['connector_cost'] = d['connector_cost']

        '''
        self.system = system
        
        self.AttachA   =d['AttachA']
        self.AttachB   =d['AttachB']
        self.DynCableA =self.system.dynamicCableConfigs[d['DynCableA']]
        self.DynCableB =self.system.dynamicCableConfigs[d['DynCableB']]
        self.headingA  =np.radians(d['headingA'])
        self.headingB  =np.radians(d['headingB'])
        self.cableType =d['cableType']
        '''
        
        # static cable end points
        self.rA = np.zeros(3)
        self.rB = np.zeros(3)
        
        # cable routing info
        self.x = []
        self.y = []
        self.r = []
        
        # get cable length
        self.getL()
            
        
        # failure probability
        self.failure_probability = {}
        
    def reposition(self,headings=None,project=None,rad_fair=[]):
        '''
        Repositions the cable based on headings of the end subcomponents and the locations of attached nodes

        Parameters
        ----------
        headings : list, optional
            List of absolute compass headings associated with the platform/substation attached
            to each end of the cable. The default is None.
        project : FAModel project object, optional
            FAModel project object associated with this cable, only used if 
            the end points of cable sections in the middle need to be set as well.
        rad_fair : list, optional
            fairlead radius of node connected on each end of the cable (list should be length 2)
            If not provided, the fairlead radius will be determined from the attached nodes' listed
            fairlead radius

        Returns
        -------
        None.

        '''
        # reposition cable and set end points for the first and last cable sections (or the dynamic cable for a suspended cable)
        if not headings:
            # convert headings to unit circle
            headingA = np.pi/2 - self.subcomponents[0].headingA
            headingB = np.pi/2 - self.subcomponents[-1].headingB
        else:
            # convert headings to unit circle
            headingA = np.pi/2 - headings[0]
            self.subcomponents[0].headingA = headings[0]
            headingB = np.pi/2 - headings[1]
            self.subcomponents[-1].headingB = headings[1]
            
        if not isinstance(self.subcomponents[0].attached_to[0], Jtube):
            if not rad_fair:
                rf = self.attached_to[0].rFair if self.attached_to[0] else 0
            else:
                if rad_fair[0] == None:
                    rf = self.attached_to[0].rFair if self.attached_to[0] else 0
                else:
                    rf = rad_fair[0]

            # calculate fairlead locations
            Aloc = [self.attached_to[0].r[0]+np.cos(headingA)*rf, 
                    self.attached_to[0].r[1]+np.sin(headingA)*rf, 
                    self.attached_to[0].zFair+self.attached_to[0].r[2]]
            self.subcomponents[0].rA = Aloc; self.rA = Aloc
        if not isinstance(self.subcomponents[-1].attached_to[-1], Jtube):
            if not rad_fair:
                rf = self.attached_to[1].rFair if self.attached_to[1] else 0
            else:
                if rad_fair[1] == None:
                    rf = self.attached_to[1].rFair if self.attached_to[1] else 0
                else:
                    rf = rad_fair[1]

            Bloc = [self.attached_to[1].r[0]+np.cos(headingB)*rf, 
                    self.attached_to[1].r[1]+np.sin(headingB)*rf, 
                    self.attached_to[1].zFair+self.attached_to[1].r[2]]
            self.subcomponents[-1].rB = Bloc; self.rB = Bloc
        
        if project:
            # set end points of subcomponents
            lensub = len(self.subcomponents)
            # if there's more than 1 subcomponent, there will be joints
            # and will need to set the locs of inner subcomponents
            if lensub > 1:
                for i,sub in enumerate(self.subcomponents):                 
                    if i == 0:
                        # get depth at location of rB
                        xy = [sub.rA[0]+np.cos(headingA)*sub.span,
                              sub.rA[1]+np.sin(headingA)*sub.span]
                        z = project.getDepthAtLocation(xy[0],xy[1])
                        # set the end B of the first subsection
                        sub.rB = [xy[0],xy[1],-z]
                        sub.z_anch = -z
                        sub.depth = z
                        # set joint
                        self.subcomponents[i+1]['r'] = sub.rB
                        # set rA of next cable section
                        self.subcomponents[i+2].rA = self.subcomponents[i].rB
                    if i == lensub-1:
                        xy = [sub.rB[0]+np.cos(headingB)*sub.span,
                              sub.rB[1]+np.sin(headingB)*sub.span]
                        z = project.getDepthAtLocation(xy[0],xy[1])
                        # set the end A of the last subsection
                        sub.rA = [xy[0],xy[1],-z]
                        # update z_anch and depth of the subsection
                        sub.z_anch = -z
                        sub.depth = z
                        # set joint
                        self.subcomponents[i-1]['r'] = sub.rA
                        # set the rB of the previous cable section (likely a static cable)
                        self.subcomponents[i-2].rB = self.subcomponents[i].rA
    
    
    def getCost(self):
        '''
        Get the total cost of the cable (includes all cable subcomponent costs)
        
        Returns
        -------
        cost: total cost of the cable
        '''
        
        cost = 0
        for sub in self.subcomponents:
            if isinstance(sub,Joint):
                if 'cost' in sub:
                    cost += sub['cost']
            else:
                cost += sub.getCost()
            # elif 'cost' in sub.dd['cable_type']:
            #     cost += sub.dd['cable_type']['cost']*sub.dd['L']
        if self.cost:
            if 'connector_cost' in self.cost:
                cost += 2*self.cost['connector_cost']
                
        
        return(cost)
    
    def getL(self):
        '''
        Get the total length of the cable (includes any routing for static cables)

        Returns
        -------
        None.

        '''
        L = 0
        for cab in self.dd['cables']:
            if isinstance(cab,StaticCable):
                cab.getLength()
            L += cab.L
            
        self.L = L
        
    def makeLine(self,buff_rad=20,include_dc=True):
        '''Make a 2D shapely linestring of the cable'''
        coords = []
        for sub in self.subcomponents:
            if isinstance(sub,Joint):
                coords.append(sub['r'][:2])
            elif isinstance(sub,DynamicCable) and include_dc:
                coords.append(sub.rA[:2])
                coords.append(sub.rB[:2])
                pass
            elif isinstance(sub, StaticCable):
                if len(sub.x)>0:
                    for i in range(len(sub.x)):
                        coords.append([sub.x[i],sub.y[i]])
        line = sh.LineString(coords)

        return(line)
    
    def updateSpan(self,newSpan):
        '''
        Change the lengths of subcomponents based on the new total span
        For dynamic-static-dynamic configurations, just increase static span
        For suspended cable, increase lengths of non-buoyant sections

        Parameters
        ---------
        newSpan : float
            New span of full cable object from platform to platform or platform to substation
        
        Returns
        -------
        None.

        '''
        # determine old span 
        oldSpan = 0
        for sub in self.subcomponents:
            if not isinstance(sub,Joint):
                oldSpan += sub.span
        # determine difference between old and new span
        spanDiff = newSpan - oldSpan
        
        jointCount = 0
        # determine number of subcomponents
        if len(self.subcomponents) > 1:
            # this is not a suspended cable
            for i,sub in enumerate(self.subcomponents):
                if isinstance(sub,StaticCable):
                    sub.span = sub.span + spanDiff
                    sub.L = sub.L + spanDiff
                
                # elif isinstance(sub,Joint):
                #     jointCount += 1
                #     # move all joints after the first joint (since the dynamic cable span is unchanged, first joint loc is unchanged too)
                #     if jointCount > 1:
                #         self.estJointLoc(i)
        else:
            # this is a suspended cable
            sub = self.subcomponents[0]
            # update span
            sub.span = sub.span + spanDiff
            sub.L = sub.L + spanDiff
            # determine number of buoyancy sections
            nb = len(sub.dd['buoyancy_sections'])
            if sub.dd['buoyancy_sections'][0]['L_mid'] < sub.dd['buoyancy_sections'][0]['N_modules']*sub.dd['buoyancy_sections'][0]['spacing']/2:
                # starts with a buoyancy section -> # of regular sections = # of buoyancy sections
                addS = 0
            else:
                # starts with a regular section -> # of reg sections = # of buoyancy sections + 1
                addS = 1
            
            nrs = nb + addS
            
            # add length to each regular section 
            addL = spanDiff/nrs
            # if abs((oldSpan-newSpan)/newSpan) > 0.1:
            #     # add buoyancy as well
            #     addL = spanDiff/(nrs+len(sub.dd['buoyancy_sections']))
            #     for i,bs in enumerate(sub.dd['buoyancy_sections']):
            #         # determine number of add'l buoyancy modules
            for i,bs in enumerate(sub.dd['buoyancy_sections']):
                bs['L_mid'] = bs['L_mid'] + addL*(i+addS)
                


        