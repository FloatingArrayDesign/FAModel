import moorpy as mp # type: ignore
import fadesign.MoorSolve as msolve
import numpy as np
import scipy
import matplotlib as mpl
import matplotlib.pyplot as plt
import fadesign.conceptual.graph_helpers as gh
# Old shared moorings linear modeling code from 2021 / updated in 2024 by Rudy Alkarem

def getUnitAndLength( rA, rB ):

    dr = rB-rA
    l = np.linalg.norm(dr)
    u = dr/l
    
    return u, l


class LinearSystem():
    '''2D array representation with linear mooring properties for optimization.
    
    
    some conventions/flags:
    - shared (formerly profile)
    - net:  False for normal/shared line;  True for >2 interconnected lines
    
    
    '''
    
    
    def __init__(self, coords, intraMat,  nPtfm, interMats=[], 
        interCoords=None, inits=None, profileMap=None, intersectZ=None, 
        rFair=0, zFair=0, depth=600., fmax=1e6, xmax=40.0, plots=0, nonlin=1.0, 
        center=True, old_mode=True):
        '''Formerly makeSimpleSystem in Array.py, this sets up a LinearSystem 
        from a coordinates list and adjacency matrix.

        Parameters
        ----------
        coords : 2D list
            intra-cell coordinates
        intraMat : 2D array
            Adjacency matrix for intra-cell lines (within)...
        nPtfm : int
            Number of intra-cell platforms or turbines            
        interMats : list of 2D arrays
            Adjacency matrix for inter-cell lines (between)...
        interCoords : list of 2D lists
            inter-cell spacing (coordinates of center of neighboring cell w.r.t. the center of the unit cell)
        inits: dictionary
            initial tension and stiffness values and hybrid line-related characteristics.
            In the "old_mode", this requires mooring groups to be provided with the following:
            - tenA                                                    # [N]  Initial tension(Anchored)
            - tenS                                                    # [N]  Initial tension(Shared)
            - klA                                                     # [N/m] Initial stiffness (Anchored)
            - klS                                                     # [N/m] Initial stiffness (Shared)
            - tenTen                                                  # [N] Initial tendon tension (Hybrid)
            - w       
            In the newer more general mode, this requires mooring groups to be provided with:
            - kl
            - kt
            - 
        profileMap: 2D list
            allows the user to manually define the profiles of each mooring group (0: anchored, 1: shared, 2: hybrid)
        intersectZ: 2D list
            depth of midpoint for all anchors (if the z value is not zero, it is a hybrid system)
        rFair : float
            Radius of fairleads from platform centerline (just a uniform value for now..)
        zFair : float
            Z elevation of fairleads (just a uniform value for now..)
        depth : float, optional
            Depth of the water in the system. The default is 600
        nonlin : float
            A multiplier to compensate for nonlinearitize in a real mooring system. 1: linear, >1: add some margin from the watch circles
        center : bool
            If true, the system will be shifted so that its center is at 0,0
        old_mode : bool
            If true (default) functions the old way with assumed weight-tension
            -stiffness relations. If False, uses a more generic approach.
        '''
        
        self.old_mode = bool(old_mode)
        
        # key array layout info
        self.coords = np.array(coords)
        self.intraMat = np.array(intraMat)
        self.interMats = np.array(interMats)
        self.nPtfm = nPtfm 
        if np.any(intersectZ==None):
            self.intersectZ = np.zeros(len(coords))
        else:
            self.intersectZ = intersectZ
        self.interCoords = interCoords
        if inits:
            self.inits = inits
        
        self.profileMap = profileMap
        # lists of objects to be created
        #self.bodyList = []              # one for each FOWT
        #self.pointList = []             # for each anchor and also each attachment point (on the FOWT bodies)
        self.mooringGroups = []          # a list of dictionaries for the linear properties of each mooring group in the array 
        
        # arrays for each actual mooring line (not groups), i.e. each anchor or shared line
        self.u = []                     # list of line unit vectors
        self.l = []                     # list of line horiztonal spacings
        self.group = []                  # line group index
        self.endA = []                  # id of what the line is connected to (this corresponds to the row or column index in the adjacency matrix)
        self.endB = []
        self.rA = []            # coordinates of the connected lines at end A
        self.rB = []            # coordinates of the connected lines at end B
        self.angA = []                  # offset angles (deg about yaw axis) of fairlead attachment on platform [deg]
        self.angB = []
        self.boundary = []             # a boolean to check if the line is a boundary (inter-cell) line
        # parameters
        self.depth = depth
        self.fmax = fmax
        
        self.xmax = xmax                # max offset (watch circle radius) constraint
        self.nonlin = nonlin
        self.angles = np.arange(0,360, 15)  # wind angles to consider
        
        self.anchor_cost_factor = 1.0    # this factor is used to scale the cost of anchor lines (whereas shared line costs aren't scaled). Can be adjusted externally.
        
        
        # shift the coordinates of all the bodies and anchors so that the center of the plot is ideally at (0,0)
        if center:
            cx = np.mean([self.coords[i,0] for i in range(len(self.coords))])
            cy = np.mean([self.coords[i,1] for i in range(len(self.coords))])
            #mss.transform(trans=[-cx, -cy])
            self.coords = self.coords - [cx,cy]

        self.center = [np.mean([self.coords[i,0] for i in range(len(self.coords))]), 
                       np.mean([self.coords[i,1] for i in range(len(self.coords))])]
        
        # also add some generic line types (which map to LineDesigns), 
        # one for each line grouping as defined in the entries of intraMat
        maxGroupNumber = np.max(intraMat)
        netGroup = 0
        if interMats:
            for interMat in interMats:
                maxGroupNumber = np.max([maxGroupNumber, np.max(interMat)])
        
        if self.old_mode:  # make mooring groups using the old approach (compatible with conceptDesign)
            
            for i in range(maxGroupNumber):               # loop up to maximum number in adjacency matrix (one for each line grouping)
                #if profileMap:  # MH: not needed <<<
                #    shared = self.profileMap[i]
                #else:
                if interMats:
                    shared = i+1 in intraMat[:nPtfm, :nPtfm] or np.any([i+1 in interMat for interMat in interMats])
                else:
                    shared = i+1 in intraMat[:nPtfm, :nPtfm]
                
                if inits:
                
                    if not np.all(self.intersectZ==0) and i+1 in intraMat[:nPtfm, nPtfm:] and shared==1:
                        percent_droop = 100 - self.intersectZ[np.where(i+1==intraMat[:nPtfm, :])[1][0]]/self.depth * 100  # Check the depth of that anchor point to determine the percent_droop, if any
                        intersectDeg = np.sum(intraMat[self.intersectZ > 0, :] > 0, axis=1)[netGroup]
                        net = True
                        tendON = True
                        netGroup += 1
                        if percent_droop <= 20:   # if the line has low drop, a tendon is too expensive
                            tendON = False
                    else:
                        net = False
                        percent_droop = 50
                    
                    if shared==0:  # anchored
                        self.mooringGroups.append(dict(type=i+1, ten=self.inits['tenA'], kl=self.inits['klA'], w=self.inits['w'], shared=shared, net=False, cost=1))
                    elif shared==1:  # shared or hybrid (net)
                        if net:
                            self.mooringGroups.append(dict(type=i+1, ten=self.inits['tenS'], kl=self.inits['klS'], w=self.inits['w'], shared=shared, percent_droop=percent_droop, net=net, tendON=tendON, tenTen=self.inits['tenTen'], intersectDeg=intersectDeg, cost=1))    
                        else:
                            self.mooringGroups.append(dict(type=i+1, ten=self.inits['tenS'], kl=self.inits['klS'], w=self.inits['w'], shared=shared, percent_droop=percent_droop, net=net, cost=1))
                else:
                    self.mooringGroups.append(dict(type=i+1, kl=100, ten=1000, w=1500, tenTen=1000, shared=shared, net=False, cost=1))  
        
        else:  # new more general/simple approach for mooring groups
            # note: inits must be provided as an input in this case (defining each mooring group)
                
            # figure out shared from place in intraMat and interMats...
            groupNumbers = []
            profileMap = []  # says whether each group is anchored or shared line
            n = self.intraMat.shape[0]
            for i in range(n):
                for j in range(n):
                    a = self.intraMat[i,j]
                    # if the entry is nonzero and not already stored, add it
                    if a > 0:
                        if not a in groupNumbers:
                            groupNumbers.append(a)
                            if i < self.nPtfm and j < self.nPtfm:
                                profileMap.append(1)  # flag as shared line
                                shared = 1
                            elif i > self.nPtfm and j > self.nPtfm:
                                raise Exception("This appears to be an anchor-anchor mooring!")
                            else:
                                profileMap.append(0)  # flag as anchored line
                                shared = 0
                                
                            # set up the mooring Group
                            self.mooringGroups.append(dict(type=a, shared=shared, net=False, cost=1))
                            self.mooringGroups[-1].update(self.inits[a-1])  # add in input info on the mooring group
                            
            
            # now go through interMat (intercell matrices) and do similarly
            for interMat in self.interMats:
                n = interMat.shape[0]
                for i in range(n):
                    for j in range(n):
                        a = interMat[i,j]
                        # if the entry is nonzero and not already stored, add it
                        if a > 0:
                            if not a in groupNumbers:
                                groupNumbers.append(a)
                                if i < self.nPtfm and j < self.nPtfm:
                                    profileMap.append(1)  # flag as shared line
                                    shared = 1
                                elif i > self.nPtfm and j > self.nPtfm:
                                    raise Exception("This appears to be an anchor-anchor mooring!")
                                else:
                                    profileMap.append(0)  # flag as anchored line
                                    shared = 0
            
                                # set up the mooring Group
                                self.mooringGroups.append(dict(type=a, shared=shared, net=False, cost=1))
                                self.mooringGroups[-1].update(self.inits[a-1])  # add in input info on the mooring group
                                
            # not sure what to do about nets in the above!!
        
        
        # make lines using adjacency matrix (and add corresponding points if they're on a platform)
        #linenum = 1
        for iA in range(len(coords)):
            for iB in range(iA):
                k = int(intraMat[iA,iB]) # the number (if positive) indicates the lineDesign type or grouping
                if k > 0:
                    dr = self.coords[iB,:] - self.coords[iA,:]
                    l  = np.linalg.norm(dr)
                    self.l.append(l)           # store length   <<<<<<<< need to subtract fairlead radii... ?
                    self.u.append(np.round(dr/l, 2))        # store unit vector
                    self.group.append(k)        # store lineDesign type (starts at 1, need to subtract 1 for index)
        
                    self.endA.append(iA)       # end A attachment object index
                    self.endB.append(iB)
                    self.rA.append(self.coords[iA, :])
                    self.rB.append(self.coords[iB, :])
                    self.angA.append(0.0)  # to be handled later <<<<<<<
                    self.angB.append(0.0)
                    self.boundary.append(False)
                    # fill in ratios in the corresponding line design for convenience - eventually need to check/enforce consistency <<<<
                    mg = self.mooringGroups[k-1]
                    shared = mg['shared']
                    
                    if self.old_mode:
                        mg['ten__w'] = mg['ten']/mg['w'] 
                        mg['kl__w' ] = mg['kl' ]/mg['w'] 
                        mg['kt__kl'] = mg['ten']/l/mg['kl']         # kt/ k = ten/l/k
                    else:
                        pass  # <<< nothing needed here??
                        
                    mg['l'] = l  # store a copy so it's acccessible in the mooringGroups list                    
                    mg['dl_max'] =  xmax  # maximum extension of this mooring group (initial value equal to xmax)
                    mg['dl_min'] = -xmax  # maximum compression of this mooring group (must be negative, initial value equal to -xmax)
                    
        # Add inter-cell lines:
        if interMats:
            for interMat, interCoord in zip(interMats, interCoords):
                interMat = np.array(interMat)
                interCoord = np.array(interCoord)
                for iA in range(self.nPtfm):
                    for iB in range(iA):        
                        if iA < self.nPtfm and iB < self.nPtfm:
                            b = int(interMat[iA, iB])
                            if b > 0:
                                intercenter = self.center + interCoord
                                rotMat = np.array([[np.cos(np.pi), -np.sin(np.pi)], 
                                                [np.sin(np.pi),  np.cos(np.pi)]])
                                intercenterp = np.matmul(rotMat, intercenter)
                                intercoordsiB = intercenter + (self.coords[iB, :] - self.center)
                                intercoordsiA = intercenter + (self.coords[iA, :] - self.center)
                                intercoordsiBp = intercenterp + (self.coords[iB, :] - self.center)
                                intercoordsiAp = intercenterp + (self.coords[iA, :] - self.center)
                                # compute both distances and choose the smallest one 
                                drBA = intercoordsiB - self.coords[iA,:]
                                lBA = np.linalg.norm(drBA)
                                drBAp = intercoordsiBp - self.coords[iA,:]
                                lBAp = np.linalg.norm(drBAp)
                                if lBAp > lBA:
                                    dr = drBA
                                    # not that it matters but for consistency, choose the right A and B:
                                    self.endA.append(iA)
                                    self.endB.append(iB)   
                                    self.rA.append(self.coords[iA, :])                         
                                    self.rB.append(intercoordsiB)

                                    self.endA.append(iB)
                                    self.endB.append(iA)   
                                    self.rA.append(self.coords[iB, :])                         
                                    self.rB.append(intercoordsiAp)
                                else:
                                    dr = drBAp
                                    self.endA.append(iA)
                                    self.endB.append(iB)
                                    self.rA.append(self.coords[iA, :])                         
                                    self.rB.append(intercoordsiBp)       

                                    self.endA.append(iB)
                                    self.endB.append(iA)
                                    self.rA.append(self.coords[iB, :])                         
                                    self.rB.append(intercoordsiA)                         
                                l = np.linalg.norm(dr)
                                self.l.append(l)
                                self.l.append(l)
                                self.u.append(dr/l)
                                self.u.append(-dr/l)
                                self.group.append(b)
                                self.group.append(b)
                                self.angA.append(0.0)
                                self.angB.append(0.0)
                                self.angA.append(0.0)
                                self.angB.append(0.0)                            
                                self.boundary.append(True)
                                self.boundary.append(True)
                                # not sure about this part:
                                mg = self.mooringGroups[b-1]
                                shared = mg['shared']
                                if self.old_mode:
                                    mg['ten__w'] = mg['ten']/mg['w'] 
                                    mg['kl__w' ] = mg['kl' ]/mg['w'] 
                                    mg['kt__kl'] = mg['ten']/l/mg['kl']         # kt/ k = ten/l/k
                                else:
                                    pass  # <<< nothing needed here?
                                mg['l'] = l
                                mg['dl_max'] =  xmax  
                                mg['dl`_min'] = -xmax  
                    
            
            print("end of intermat setup?")

        self.nLines = len(self.l)
        
        # should make some error checks for consistent properties (spacing, shared) in mooring groups <<<
        
        # now also make the structure matrix        
        self.StructureMatrix = np.zeros([2*self.nPtfm, self.nLines])     # rows: DOFs; columns: lines
        
        for j in range(self.nLines):  
            if self.endA[j] < self.nPtfm: # only if not an anchor
                self.StructureMatrix[self.endA[j]*2    , j] =  self.u[j][0]
                self.StructureMatrix[self.endA[j]*2 + 1, j] =  self.u[j][1]
                
            if self.endB[j] < self.nPtfm: # only if not an anchor
                self.StructureMatrix[self.endB[j]*2    , j] = -self.u[j][0]
                self.StructureMatrix[self.endB[j]*2 + 1, j] = -self.u[j][1]
    
        
        # remember, you may want to  call calibrate (or similar) to set up better
        # values for ten__w, kl__w, and kt__kl for each mooring object assuming a continuous catenary line
        
    
    def preprocess(self, plots=0, display=0):
        '''Initializes things...
        
        Does all the things that can be done once the lineDesign characteristics are set (f/l/k and f/w)
        
        the mooring system objects to their initial positions if applicable?
        
        contents of former getTensionMatrix and getKnobs are now here'''
        
        # ensure the points (and line ends) are in the right positions
        #for point in self.pointList:
        #    point.setPosition(Point.r)
        
        # update line properties so that initial values are in place
        #self.calibrate_kt_over_k(plots=plots)
        
        # fill in the initial line stiffnesses and generate the system stiffness matrix so it's in place
        #self.updateStiffnesses(np.ones(self.nLines))  
                
        # draw initial mooring system if desired
        #if plots==1:
        #    self.plot(title="Mooring system at initialization")

    
        #def getTensionMatrix(self):
        '''Tension matrix defines the contribution of each line group's weight to each line's tension
        
        Essentially it is just a mapping from each line group's weight to each individual line's tension.
        There is only one nonzero entry per row - i.e. each line's tension is just based on a single group's stiffness.
        This seems simple and maybe doesn't need to be a matrix??
        '''
        
        self.TensionMatrix = np.zeros([self.nLines, len(self.mooringGroups)])
    
        for j in range(self.nLines):
            
            i = self.group[j]-1
            
            if self.old_mode:
                self.TensionMatrix[j, i] = self.mooringGroups[i]['ten__w']
            else:  # NEW - USING TENSIONS DIRECTLY RATHER THAN WEIGHT RATIOS
                self.TensionMatrix[j, i] = self.mooringGroups[i]['ten']
                
            if self.boundary[j]:
                self.TensionMatrix[j, i] *= 0.5  #MH: this seems suspect <<<
        
        
        
        #def getKnobs(self):
        '''based on structure and tension matrices, calculatesd self.Knobs_k, which is used by c_to_k when optimizing stiffness.'''
    
        # Null space of Structure Matrix
        N1 = scipy.linalg.null_space(self.StructureMatrix)#, rcond = 0.0001)
        
        # null space of N1 augmented with tension matrix 
        N2 = scipy.linalg.null_space(np.hstack([N1, -self.TensionMatrix])) #, rcond = 0.0001)
        #N2 = scipy.linalg.null_space(np.append(N1, -self.TensionMatrix,1))#, rcond = 0.0001)
        
        # nullspace matrix containing basis vectors of valid line weight solutions for equilibrium given line groupings 
        # (we skip the top part of the vectors--the coefficients for feasible tension combinations--because tensions can be calculated directly from weights)
        self.wBasis = N2[-len(self.mooringGroups):]
        
        # self.getSystemStiffness()
        # print(self.SystemStiffness)

        # check to make sure there is at least on design variable
        if np.prod(self.wBasis.shape)==0:   
            raise Exception("No knobs available") 
        
        # normalize each weight basis vector and flip signs on any that are mostly negative
        self.nBasis = self.wBasis.shape[1]       # store the number of vectors
        #print(self.nKnobs)
        
        for i in range(self.nBasis):
            self.wBasis[:,i] = self.wBasis[:,i] / np.linalg.norm(self.wBasis[:,i]) * np.sign(np.sum(self.wBasis[:,i]))
            #self.Knobs_k[:,i] = self.Knobs_k[:,i] / np.linalg.norm(self.Knobs_k[:,i])
        
        
        #Create Initial guess for q (the knob values multiplied by the weight basis vectors)
        self.q0 = np.zeros(self.nBasis)+100

        # lower any initial knob values for basis vectors that contain negatives to avoid any negative weights to start with
        #for j in range(len(self.mooringGroups)):
        #    for i in range(self.nBasis):  
        #        wtemp = np.sum(self.wBasis[:,i]*self.q0
        #        if any(self.wBasis[:,i] < 0):
        #            self.q0[i] = 0.0
        # raise the knob of the all-positive basis vector if needed to make all weights positive
        for i in range(self.nBasis):  
            w = np.matmul(self.wBasis, self.q0)               # initial weights per unit length of each line group
            if any(w < 0) and all(self.wBasis[:,i] > 0):       # if there's a negative line weight and this basis vector is all positive
                
                q_i_add = np.max(-w/self.wBasis[:,i])
                if display > 0:
                    print(f' to resolve negative initial line tension ({w}), increasing q0[{i}] by {2*q_i_add}')
                
                self.q0[i] += 2*q_i_add
                
        
        

    def getSystemStiffness(self):
        '''Calculates the stiffness matrix of the SimpleSystem, based on current positions and LineType info'''
    
        #If were to generalize so each point has 3 or more degrees of freedom, replace each 2 with a 3
        #If we were to further generalize to consider points and bodies, we have to put more carefull thought into the indexing
        
    
        # self.SystemStiffness = np.zeros([2*len(self.coords), 2*len(self.coords)])
        # MH: Changed back to nPtfm x nPtfm. Not sure why anchors were included. Maybe for hybrid? <<<
        self.SystemStiffness = np.zeros([2*self.nPtfm, 2*self.nPtfm])
        
        for j in range(self.nLines): 
        
            # first get the line's stiffness matrix (xA, yA, xB, yB)
            cc = self.u[j][0]*self.u[j][0]  # cos(theta)*cos(theta)
            ss = self.u[j][1]*self.u[j][1]
            cs = self.u[j][0]*self.u[j][1]
            
            # Find Transformation Matrices:
            transMat_inline        = np.array([      [ cc, cs,-cc,-cs],
                                                     [ cs, ss,-cs,-ss],
                                                     [-cc,-cs, cc, cs],
                                                     [-cs,-ss, cs, ss]])

            transMat_perpendicular = np.array([      [ ss,-cs,-ss, cs ],
                                                     [-cs, cc, cs,-cc ],
                                                     [-ss, cs, ss,-cs ],
                                                     [ cs,-cc,-cs, cc ]])
            # Lookup inline and perpendicular stiffness values for this line type (assumes a certain line spacing, etc.)            
            mg = self.mooringGroups[self.group[j]-1]
            if self.old_mode:
                kl = mg['kl__w']*mg['w']            
                kt = mg['kt__kl']*kl
                mg['kl'] = kl
            else:  # the new mode uses stiffness values that are already provided
                kl = mg['kl']            
                kt = mg['kt']
                
            # Multiply stiffness values by transformation matrix
            K_inline = kl * transMat_inline
                
            #Force in y direction from displacement in y direction caused by tension in x direction   
            K_perpendicular = kt * transMat_perpendicular
            
            # Note: Force in x direction from displacement in y direction caused by tension in x direction is neglected as second-order
            K_sum = K_inline + K_perpendicular

            # now apply to the appropriate spots in the system stiffness matrix        
            iA = self.endA[j]
            iB = self.endB[j]
            
            # MH: re-adding the old logic here >>>
            if self.endA[j]<self.nPtfm:
                self.SystemStiffness[iA*2:iA*2+2, iA*2:iA*2+2] += K_sum[:2,:2]
                
            if self.endB[j]<self.nPtfm:
                self.SystemStiffness[iB*2:iB*2+2, iB*2:iB*2+2] += K_sum[2:,2:]
                
            if self.endA[j]<self.nPtfm and self.endB[j]<self.nPtfm:
                self.SystemStiffness[iA*2:iA*2+2, iB*2:iB*2+2] += K_sum[:2,2:]
                self.SystemStiffness[iB*2:iB*2+2, iA*2:iA*2+2] += K_sum[2:,:2]
            
            '''
            # >>> MH: this part may be for hybrid shared moorings >>>
            boundaryLineCounti = np.sum(self.boundary[:j])
            if boundaryLineCounti % 2 == 0:  # only count the stiffness of the boundary line once
                self.SystemStiffness[iA*2:iA*2+2, iA*2:iA*2+2] += K_sum[:2,:2]
                self.SystemStiffness[iB*2:iB*2+2, iB*2:iB*2+2] += K_sum[2:,2:]
                self.SystemStiffness[iA*2:iA*2+2, iB*2:iB*2+2] += K_sum[:2,2:]
                self.SystemStiffness[iB*2:iB*2+2, iA*2:iA*2+2] += K_sum[2:,:2]                 
                if iA >= self.nPtfm:
                    tau = mg.get('tenTen', 1000)
                    tau__L = tau/self.intersectZ[iA]  # intentionally will be set to inf if it's an anchor
                    # Only apply if there is a tendon
                    if mg.get('net', False) and not mg.get('tendON', False):
                        tau__L = 0
                    
                    self.SystemStiffness[iA*2:iA*2+2, iA*2:iA*2+2] += (np.array([[1, 0],[0, 1]]) * tau__L) 
            '''
        """
        >>> MH: maybe this was a clever approach to remove anchor row/columns? >>>
        # remove any rows and columns in the stiffness matrix with infinity values
        finite_mask = np.isfinite(self.SystemStiffness).all(axis=1)
        self.SystemStiffness = self.SystemStiffness[np.ix_(finite_mask, finite_mask)]
        self.nAnch = int(np.sum(finite_mask==False)/2)

        >>> RA: This works for non-hybrid designs but need to think of a new way to include 
        net/hybrid designs.
        """
        self.nAnch = len(self.coords) - self.nPtfm  #MH: temporarily filling this in here <<<
        
        # self.SystemStiffness[np.abs(self.SystemStiffness) < 1e-5] = 0
        # calculate inverse of stiffness matrix to avoid solving multiple times later
        self.K_inverse = np.linalg.inv(self.SystemStiffness)
    

    def get_x(self, f, theta=0, heading=0, removeHybrid=True):
        '''Get displacement in all dof using linear model. Nonlinear factor included.
        This assumes system in equilibrium when no external force is applied.
        
        theta is wind directions in radians, heading is the wind direction in degrees
        '''
        
        #if watch_circles > 0:
        if not hasattr(self, 'K_inverse'):
            raise Exception("In a Linear System, getSystemStiffness must be called before calling get_x.")
        
        
        if theta==0:
            theta = np.radians(heading)
        
        if np.isscalar(f):                  # thrust force and direction
            F = [f*np.cos(theta), f*np.sin(theta)]*(len(self.coords)-self.nAnch)
            F[2*self.nPtfm:] = [0, 0]*(len(self.coords) - self.nPtfm - self.nAnch)
            
        elif len(f)==2:                  # x and y force components to be applied to each turbine
            F = [f[0], f[1]]*(len(self.coords)-self.nAnch)
            F[2*self.nPtfm:] = [0, 0]*(len(self.coords) - self.nPtfm - self.nAnch)
        
        elif len(f)==2*(len(self.coords)-self.nAnch):          # full vector of forces
            F = f
        
        else:
            raise ValueError("Invalid format of f provided to get_x")


        #Use linear algebra to solve for x vector (solve for offsets based on stiffness matrix and force vector)  Nonlinear factor included here.   
        xi = np.matmul(self.K_inverse, F)*self.nonlin        
        
        # also figure out peak tensions etc here? <<<<
        
        
        if removeHybrid:
            # remove hybrid and anchor points
            xi = xi[:2*self.nPtfm]
        return xi
     
        
        
    def windsweep(self, f=0):
        ''' gets offsets and changes in line spacing across wind directions. 
        '''
        
        if f == 0:
            f=self.fmax
            
        self.xi_sweep = np.zeros([len(self.angles), 2*self.nPtfm])    # holds displacement vector (x and y of each platform) for each wind direction
        self.dl_sweep = np.zeros([len(self.angles), self.nLines])      # change in each line's spacing for each wind direction
        
        for i,angle in enumerate(self.angles):
                                
            xi = self.get_x(f, heading=angle)                            # Get the offsets in each DOF
            self.xi_sweep[i,:] = xi            
            
            for il in range(self.nLines):                                   # loop through line indices
                dl = 0.0
                iA = self.endA[il]
                iB = self.endB[il]
                if iA < self.nPtfm:                                     # if this end is attached to a platform
                    dl += np.sum( xi[2*iA : 2*iA+2] * self.u[il])       # calculate extension as dot product of displacement and line unit vector
                if iB < self.nPtfm:                                     # if this end is attached to a platform
                    dl -= np.sum( xi[2*iB : 2*iB+2] * self.u[il])       # calculate extension as -dot product of displacement and line unit vector
            
                self.dl_sweep[i, il] = dl
        
        # also compute watch circle area (by summation of areas of triangles)
        self.areas = np.zeros(self.nPtfm)
        for i in range(self.nPtfm):
        
            for j in range(-1, len(self.angles)-1):                                
                self.areas[i] += self.xi_sweep[j,2*i] * (self.xi_sweep[j+1,2*i+1] - self.xi_sweep[j-1,2*i+1]) * 0.5   
        
        return
        
    def getCost(self):
        '''Calculate the cost function of the line for the spring model'''
        #Assume cost is proportional to line length and mass
                    
        #self.line_cost = [self.l[i]*self.mooringGroups[self.group[i]-1]['w'] for i in range(self.nLines)] <<< this was double counting
        self.line_cost = []
         
        for i in range(self.nLines):
            
            self.line_cost.append(self.l[i]*self.mooringGroups[self.group[i]-1]['w']) # Cost Function for each line [kg m]
            
            if self.mooringGroups[self.group[i]-1]['shared'] != 1:                    # If it's an anchor mooring
                self.line_cost[-1] = self.line_cost[-1]*self.anchor_cost_factor                 # include an anchorcost factor 
            elif self.mooringGroups[self.group[i]-1]['net'] and self.mooringGroups[self.group[i]-1]['tendON']:  # if there's an anchor in a hybrid system
                self.line_cost[-1] = self.line_cost[-1]*self.anchor_cost_factor/self.mooringGroups[self.group[i]-1]['intersectDeg']                 # include a shared anchorcost factor 
            # sloppily store the cost in the mooring group as well for use in vizualization
            self.mooringGroups[self.group[i]-1]['cost'] = self.line_cost[-1]  
            
        self.cost = np.sum(np.array(self.line_cost))   # save total cost
        
        return(self.cost)
           
           
    def getWeight(self):
        '''Calculate the total system mooring line (wet) weight'''
                    
        return sum([self.l[i]*self.mooringGroups[self.group[i]-1]['w'] for i in range(self.nLines)])
            
        
    def optimize(self, display=0):
        '''solve the cheapeast system for a given input force, wind angle, and maximum displacement'''
        
        if not self.old_mode:
            raise Exception("LinearSystem.optimize only works when old_mode = True")
        
        if display > 1:
            print(f'Beginning LinearSystem optimization.')
            
            print('weight basis vectors are:')
            print(self.wBasis)
        
        def dopt_fun(q):
            '''evaluation function for the dopt solver. This function includes both 
            design variables and constraints. This function inputs q which is an array of knob values.
            '''
            
            # ----- calculate line weight from q, and upate line types and system stiffness matrix ------------                
            w = np.matmul(self.wBasis, q)                       # weights per unit length of each line group
            
            for i, mg in enumerate(self.mooringGroups):
                mg['w']   = w[i]                                # update wet weight per unit length [N/m]
              
            self.getSystemStiffness()   # need this for evaluating constraints
            
            # ---- objective function - mooring system mass or cost -----
            f = self.getCost()
            
            
            # ----- constraint values - margin from max offsets, and line weights (must be positive) -----
            '''
            Finds how much a certain stiffness will undershoot the maximum design displacement
            This function returns a list of undershoots for each degree of freedom.             
            '''
            
            self.windsweep()   # update offset and line extension numbers for each wind direction
            
            self.offsets = np.zeros([len(self.angles), self.nPtfm])  # store offset distance (radius) for each direction
            
            for i,angle in enumerate(self.angles):  # Calculate the hypotenuse of each platform's offset                
                self.offsets[i,:] = [np.linalg.norm(self.xi_sweep[i, 2*k:2*k+2]) for k in range(self.nPtfm)]  
               
            peak_offsets = np.max(self.offsets, axis=0)             # get largest displacement of each platform (over the range of wind directions)
            
            offset_margin = [self.xmax - peak_offset for peak_offset in peak_offsets]         # substract maximum from resultant values to get undershoot
            
            g = np.hstack([offset_margin, w.tolist()])
            
            
            # return the objective and constraint values
            return f, g, w, 0, 0, False
            
                
        dX_last = np.zeros(self.nBasis)+1e5
        
        
        # call the optimizer (set display=2 for a desecription of what's going on)
        #q, f, res = msolve.dopt(dopt_fun, q0, tol=0.001, maxIter=70, a_max=1.5, dX_last=dX_last, display=max(display-1,0))
        #q, f, res = msolve.dopt(dopt_fun, self.q0, tol=0.002, stepfac=100, maxIter=100, a_max=1.5, dX_last=dX_last, display=display-1)
        q, f, res = msolve.dopt2(dopt_fun, self.q0, tol=0.002, stepfac=100, maxIter=100, a_max=1.5, dX_last=dX_last, display=display-1)
            
        
        # check the results - and retry up to 4 times
        for i in range(4):
            if display > 0:  print(f" Message from dopt: {res['message']}")

            f, g, w, _, _, _ = dopt_fun(q)     # get latest objective and constraint values
            
            #if display>0:                
            #    if res['success'] == False:
            #        print('LinearSystem Mooring optimization was UNSUCCESSFUL: '+res['message'])                
            #    else:
            #        print(f"LinearSystem Mooring optimization was successful after {res['iter']} iterations.")

            # check for overly stiff or soft solution (and the rerun with better starting points)
            if self.nPtfm==1:
                offset_margin = g[0]  # <<< can this be simplified?
            else:
                offset_margin = np.min(g[:-len(q)-1])   # this is the closest any watch circle comes to the limit
                
            if offset_margin > 0.05*self.xmax:  # if the closest it gets to the target watch circles is more than 5% short
                if display > 0: print(f' LinearSystem optimization attempt {i} detected overly small watch circles (largest is {offset_margin:5.1f} m from the limit).')                        
                if display > 1: print(' Retrying the optimization with lighter starting points (q0)')
                self.q0 = 0.3*self.q0
                                
                
            elif offset_margin < -0.1*self.xmax:      # if it overshoots the target watch circles by more than 10%
                if display > 0: print(f' LinearSystem optimization attempt {i} detected extreme watch circles (largest is {-offset_margin:5.1f} m over the limit).')                        
                if display > 1: print(' Retrying the optimization with heavier starting points (q0)') 
                self.q0 = 10.0*self.q0
                
            else:  # otherwise, call it succsessful            
                if display>0:  print(f" LinearSystem optimization attempt {i} was successful after {res['iter']} iterations.")
                break
            
            # this is where we rerun dopt with the modified settings
            q, f, res = msolve.dopt(dopt_fun, self.q0, tol=0.002, stepfac=100, maxIter=100, a_max=1.5, dX_last=dX_last, display=display-1)
                
        
        if display>1:
            
            if res['success'] == False:
                print('Final LinearSystem Mooring optimization was UNSUCCESSFUL: '+res['message'])                
            else:
                print(f"Final LinearSystem Mooring optimization was successful after {res['iter']} iterations.")
        

        
        # plotting
            
        if ((res['success'] == False and display >0) or display > 1) and self.nPtfm>1:  # plot the optimization if it's likely desired
            
            n = len(q)
            fig, ax = plt.subplots(n+3, 1, sharex=True)
            Xs = np.array(res["Xs"])
            Fs = np.array(res["Fs"])
            Gs = np.array(res["Gs"])
            iter        = res["iter"]
            
            for i in range(n):
                ax[i].plot(Xs[:iter+1,i])
                ax[i].set_ylabel(f"q{i}")
            
            ax[n].plot(Xs[:iter+1,n:])
            ax[n].set_ylabel("weights")
                
            ax[n+1].plot(Fs[:iter+1])
            ax[n+1].set_ylabel("cost")
            
            #m = len(self.Knobs_k)
            ax[n+2].plot(Gs[:iter+1,:-n-1])
            #ax[n+3].plot(Gs[:iter+1,-n-1:])
            ax[n+2].set_ylabel("g offsets")
            #ax[n+3].set_ylabel("g weights")
            ax[n+2].set_xlabel("iteration")
            
            #breakpoint()
            if display > 1:
                plt.show()
                breakpoint()
             
        
        #For debugging purposes:
        self.res = res
        self.q = q
        
        
        # get max extension of each line group's spacing and store it in the mooring group   
        dl_max = np.max(self.dl_sweep, axis=0)
        dl_min = np.min(self.dl_sweep, axis=0)        
        #print((" group:  "+"".join([" {:6d}"]*self.nLines)).format(*self.group ))
        #print((" dl_max: "+"".join([" {:6.1f}"]*self.nLines)).format(*dl_max.tolist() ))
        #print((" dl_min: "+"".join([" {:6.1f}"]*self.nLines)).format(*dl_min.tolist() ))
        for i, mg in enumerate(self.mooringGroups):
            mg['dl_max'] = np.mean(dl_max[[j for j, k in enumerate(self.group) if k==i+1]])    # take the mean from any lines in this mooring group (group i+1)
            mg['dl_min'] = np.mean(dl_min[[j for j, k in enumerate(self.group) if k==i+1]])  
                
        # update each mooring group's weight and tension values
        for i, mg in enumerate(self.mooringGroups):
            mg['w']   = w[i]                                # update wet weight per unit length [N/m]
            mg['ten'] = mg['ten__w']*mg['w']                # update line tension [N]
            
            if np.round(w[i], 3) == 0:
                w[i] = 0.0
            
            if w[i] < 0:
                raise ValueError("breakpoint due to negative weight")
          
        self.getSystemStiffness()   # need this for evaluating constraints
        
        return q
        
    
    
    def optimize2(self, display=0):
        '''NEW: Figure out what mooringGroup stiffnesses will achieve the 
        desired watch circles, for a given input force, wind angle, and 
        maximum displacement'''
        
        if self.old_mode:
            raise Exception("LinearSystem.optimize2 only works when old_mode = False")
        
        if display > 1:
            print(f'Beginning LinearSystem optimization2.')
            
            print('tension basis vectors are:')
            print(self.wBasis)
        
        
        self.iter=0 # reset iteration counter
        
        def dopt_fun(kls):
            '''evaluation function for the dopt solver. This function includes 
            both design variables and constraints. This function inputs kls, 
            inline stiffness values.
            '''
            
            # ----- upate line types and system stiffness matrix ------------                
            
            for i, mg in enumerate(self.mooringGroups):
                mg['kl']   = kls[i]  # update inline stiffness [N/m]
              
            self.getSystemStiffness()   # need this for evaluating constraints
            
            # ---- objective function - mooring system mass or cost -----
            # approximate cost as product of line length and stiffness
            line_stiffness_costs = [self.l[i]*self.mooringGroups[self.group[i]-1]['kl'] for i in range(self.nLines)]
            f = sum(line_stiffness_costs)
            
            
            # ----- constraint values - margin from max offsets, and line weights (must be positive) -----
            '''
            Finds how much a certain stiffness will undershoot the maximum design displacement
            This function returns a list of undershoots for each degree of freedom.             
            '''
            
            self.windsweep()   # update offset and line extension numbers for each wind direction
            
            # store offset distance (radius) for each direction
            self.offsets = np.zeros([len(self.angles), self.nPtfm])  
            for i,angle in enumerate(self.angles):            
                self.offsets[i,:] = [np.linalg.norm(self.xi_sweep[i, 2*k:2*k+2]) for k in range(self.nPtfm)]  
               
            peak_offsets = np.max(self.offsets, axis=0)             # get largest displacement of each platform (over the range of wind directions)
            
            offset_margin = [self.xmax - peak_offset for peak_offset in peak_offsets]         # substract maximum from resultant values to get undershoot
            
            g = np.hstack([offset_margin, kls.tolist()])  # constraints are offset and positive stiffness
            
            self.iter = self.iter + 1 # update counter (this isn't actually iterations, it's function calls)
            if display > 3 and self.iter%20 == 0:
                sys.plot2d(watch_circles=4, line_val="stiffness")
                plt.show()
            
            
            # return the objective and constraint values
            return f, g, [], 0, 0, False
            
        n = len(self.mooringGroups)
        kls0 = np.array([mg['kl'] for mg in self.mooringGroups])  # starting point
        dX_last = np.zeros(n) + 0.01*np.mean(kls0)  # step size
        
        
        # call the optimizer (set display=2 for a desecription of what's going on)
        #q, f, res = msolve.dopt(dopt_fun, q0, tol=0.001, maxIter=70, a_max=1.5, dX_last=dX_last, display=max(display-1,0))
        #q, f, res = msolve.dopt(dopt_fun, self.q0, tol=0.002, stepfac=100, maxIter=100, a_max=1.5, dX_last=dX_last, display=display-1)
        kls, f, res = msolve.dopt2(dopt_fun, kls0, tol=0.002, stepfac=20, maxIter=40, a_max=1.4, dX_last=dX_last, display=display-1)
            
        
        # check the results - and retry up to 4 times
        for i in range(4):
            if display > 0:  print(f" Message from dopt: {res['message']}")

            f, g, _, _, _, _ = dopt_fun(kls)     # get latest objective and constraint values
            
            #if display>0:                
            #    if res['success'] == False:
            #        print('LinearSystem Mooring optimization was UNSUCCESSFUL: '+res['message'])                
            #    else:
            #        print(f"LinearSystem Mooring optimization was successful after {res['iter']} iterations.")

            # check for overly stiff or soft solution (and the rerun with better starting points)
            if self.nPtfm==1:
                offset_margin = g[0]  # <<< can this be simplified?
            else:
                offset_margin = np.min(g[:-len(kls)-1])   # this is the closest any watch circle comes to the limit
                
            if offset_margin > 0.05*self.xmax:  # if the closest it gets to the target watch circles is more than 5% short
                if display > 0: print(f' LinearSystem optimization attempt {i} detected overly small watch circles (largest is {offset_margin:5.1f} m from the limit).')                        
                if display > 1: print(' Retrying the optimization with lighter starting points (q0)')
                self.q0 = 0.3*self.q0
                                
                
            elif offset_margin < -0.1*self.xmax:      # if it overshoots the target watch circles by more than 10%
                if display > 0: print(f' LinearSystem optimization attempt {i} detected extreme watch circles (largest is {-offset_margin:5.1f} m over the limit).')                        
                if display > 1: print(' Retrying the optimization with heavier starting points (q0)') 
                self.q0 = 10.0*self.q0
                
            else:  # otherwise, call it succsessful            
                if display>0:  print(f" LinearSystem optimization attempt {i} was successful after {res['iter']} iterations.")
                break
            
            # this is where we rerun dopt with the modified settings
            q, f, res = msolve.dopt(dopt_fun, kls0, tol=0.002, stepfac=100, maxIter=100, a_max=1.5, dX_last=dX_last, display=display-1)
                
        
        if display>1:
            
            if res['success'] == False:
                print('Final LinearSystem Mooring optimization was UNSUCCESSFUL: '+res['message'])                
            else:
                print(f"Final LinearSystem Mooring optimization was successful after {res['iter']} iterations.")
        

        
        # plotting
            
        if True: #((res['success'] == False and display >0) or display > 1) and self.nPtfm>1:  # plot the optimization if it's likely desired
            
            n = len(kls)
            fig, ax = plt.subplots(n+3, 1, sharex=True)
            Xs = np.array(res["Xs"])
            Fs = np.array(res["Fs"])
            Gs = np.array(res["Gs"])
            iter        = res["iter"]
            
            for i in range(n):
                ax[i].plot(Xs[:iter+1,i])
                ax[i].set_ylabel(f"kls{i}")
            
            ax[n].plot(Xs[:iter+1,n:])
            ax[n].set_ylabel("weights")
                
            ax[n+1].plot(Fs[:iter+1])
            ax[n+1].set_ylabel("cost")
            
            #m = len(self.Knobs_k)
            ax[n+2].plot(Gs[:iter+1,:-n-1])
            #ax[n+3].plot(Gs[:iter+1,-n-1:])
            ax[n+2].set_ylabel("g offsets")
            #ax[n+3].set_ylabel("g weights")
            ax[n+2].set_xlabel("iteration")
            
             
        
        #For debugging purposes:
        self.res = res
        self.kls = kls 
        
        
        # get max extension of each line group's spacing and store it in the mooring group   
        dl_max = np.max(self.dl_sweep, axis=0)
        dl_min = np.min(self.dl_sweep, axis=0)        
        #print((" group:  "+"".join([" {:6d}"]*self.nLines)).format(*self.group ))
        #print((" dl_max: "+"".join([" {:6.1f}"]*self.nLines)).format(*dl_max.tolist() ))
        #print((" dl_min: "+"".join([" {:6.1f}"]*self.nLines)).format(*dl_min.tolist() ))
        for i, mg in enumerate(self.mooringGroups):
            mg['dl_max'] = np.mean(dl_max[[j for j, k in enumerate(self.group) if k==i+1]])    # take the mean from any lines in this mooring group (group i+1)
            mg['dl_min'] = np.mean(dl_min[[j for j, k in enumerate(self.group) if k==i+1]])  
                
        # update each mooring group's kl
        for i, mg in enumerate(self.mooringGroups):
            mg['kl']   = kls[i] 
            
            if np.round(kls[i], 3) == 0:
                kls[i] = 0.0
            
            if kls[i] < 0:
                raise ValueError("breakpoint due to negative kl")
          
        self.getSystemStiffness()   # need this for evaluating constraints
        
        return kls
    
    
    def plot2d(self, ax=None, **kwargs):
        '''Plots 2d view of simple system mooring configuration, optionally including additional properties
        
        Parameters
        ----------
        ax : matplotlib axes
            The axes to draw the plot on. A new figure is created and returned if this is not provided.
        show_lines : string
            Specifies whether to show lines: none, anch, shared, all (default)
        watch_circles: float
            Specifies whether to draw watch circles (and at what scale, >0) or not to draw them (0)
        line_val : string
            Specifies what to show for the lines: uniform, two, groups (default), stiffness, cost, weight, tension
        colormap : int or string
            Specifies what colormap to use.
        colorbar : int
            0 - don't draw one, 1 - draw as normal, 2 - draw on seperate axes specified in kwarg cbax.
        colorscale : string
            Specify linear or log (for logarithmic)
        cbax : plt.Axes
            Only used if colorbar=2
        show_axes : bool
            Whether to show the axes of the figure or not (hide them).
        labels : string
            Whether to label lines (l), points (p or t for turbine, a for anchor), etc. '' means no labels.
        title : string
            Text to add above the figure (otherwise default text will be shown).
        line_color 
        line_style    
        line_width    
        '''
        
        
        #plt.ion()  #Turn on interactive mode

        # initialize some plotting settings
        n = self.nLines
        
        # some optional argument processing and setting default values if not supplied
        
        line_val      = kwargs.get("line_val"     , "groups"     )  # get the input value, or use "groups" as default
        show_lines    = kwargs.get("show_lines"   , "all"        )  # get the input value, or use "groups" as default
        watch_circles = kwargs.get("watch_circles", 0            )  # 
        colormap      = kwargs.get("colormap"     , 0            )  # 
        colorbar      = kwargs.get("colorbar"     , 1            )  # 
        colorscale    = kwargs.get("colorscale"   , "linear"     )  # 
        show_axes     = kwargs.get("show_axes"    , True         )  # 
        labels        = kwargs.get("labels"       , ''           )  # 
        title         = kwargs.get("title"        , []           )  # 
        figsize       = kwargs.get("figsize"      , (5,5)        )  # 
        wea           = kwargs.get("wea"          , None         )  #
        #center        = kwargs.get("center"       , 1            )  # turns on and off whether the plot is centered or not
        
        # receive or use default uniform line color/style/width (may be overriden by non-uniform color coding options in line_val)
        colors = [kwargs.get("line_color", 'black')]*n
        styles = [kwargs.get("line_style", 'solid')]*n
        thicks = [kwargs.get("line_width",       2)]*n
        
        
        
        # set up colormap
        if colormap == 0 or colormap == "rainbow" or colormap == "jet":
            #Create Rainbow colormap (I still incorrectly use 'jet' sometimes when I want a rainbow colormap, so I will keep support for that keyworkd)
            cmap = mpl.cm.rainbow
        
        elif colormap == 1 or colormap == 'aut':  #Create autumn colormap
            cmap = mpl.cm.autumn
        
        else:
            raise ValueError("invalide colormap input provided to plot2d.")
             

        # set whether colormap will be linear of logarithmic
        if colorscale == "linear":
            normalizer = mpl.colors.Normalize
        elif colorscale == "log" or colorscale == "logarithmic":
            normalizer = mpl.colors.LogNorm
        else:
            raise ValueError("colorscale must be 'linear' or 'log'.")
        
        # set up color map bounds if provided
        if "val_lim" in kwargs:
        
            def getLineColors(values):
                norm = normalizer(vmin=kwargs["val_lim"][0], vmax=kwargs["val_lim"][1])     #Scale to min and max values
                s_m = mpl.cm.ScalarMappable(norm=norm, cmap = cmap)                    # create Scalar Mappable for colormapping
                return s_m.to_rgba(values), s_m
                
        else:
              
            def getLineColors(values):
                if min(values) == max(values):
                    norm = normalizer(vmin=0, vmax=max(values))    # if only one value, start scale at zero
                else: 
                    norm = normalizer(vmin=min(values), vmax=max(values))    # set scaling to min and max values
                s_m = mpl.cm.ScalarMappable(norm = norm, cmap = cmap)                    # create Scalar Mappable for colormapping
                return s_m.to_rgba(values), s_m
            
        '''
        # get arrays of all the values of interest up-front
        line_k    = [0.001*self.lineTypes[key].k         for key in self.lineTypes]
        line_c    = [0.001*self.lineTypes[key].cost      for key in self.lineTypes]  # note: this is when the cost parameter has been repurposed from $/m to $/line
        line_m    = [      self.lineTypes[key].mlin      for key in self.lineTypes]
        line_w    = [      self.lineTypes[key].w         for key in self.lineTypes]
        line_t    = [0.001*self.lineTypes[key].t         for key in self.lineTypes]
        line_kt_k = [      self.lineTypes[key].kt_over_k for key in self.lineTypes]
        line_MBL =  [0.001*self.lineTypes[key].MBL       for key in self.lineTypes]
        line_MSF =  [ MBL/t for MBL,t in zip(line_MBL,line_t)]  # safety factor
        '''
        
        line_k    = [self.mooringGroups[self.group[i]-1]['kl' ]/1e3    for i in range(n)]   # stiffness in kN/m
        line_t    = [self.mooringGroups[self.group[i]-1]['ten']/1e6    for i in range(n)]   # tension in MN
        line_w    = [self.mooringGroups[self.group[i]-1]['w']          for i in range(n)]   # wet weight per length in N/m
        line_m    = [self.mooringGroups[self.group[i]-1]['w']/9.81     for i in range(n)]   # wet weight per length in kg/m
        #line_t_w  = [self.mooringGroups[self.group[i]-1]['ten__w']     for i in range(n)]      # can add this in later
        #line_k_w  = [self.mooringGroups[self.group[i]-1]['kl__w']      for i in range(n)]      # can add this in later
        if self.old_mode:
            line_kt_k = [self.mooringGroups[self.group[i]-1]['kt__kl']     for i in range(n)]
        line_cost = [self.mooringGroups[self.group[i]-1]['cost']       for i in range(n)]
        
        
        clist = ['tab:blue','tab:cyan','tab:green','tab:olive','tab:brown','tab:purple',
                     'tab:red','tab:orange','tab:blue','tab:pink','tab:gray']
        
        
        # set up line data display - Detetermine which line variable we are using            
        if line_val == 'uniform':   # all lines drawn black and solid
            pass
            
        # elif line_val == 'two':     # distinguishes shared vs. anchor lines
        #     for i in range(n):
        #         if self.mooringGroups[i]['shared']:
        #             colors[i] = "blue"
        #             styles[i] = "solid"
        #         else:
        #             colors[i] = "black"
        #             styles[i] = "dashed"
        
        elif line_val == 'groups':
            for i in range(n):
                ii = self.group[i]-1
                colors[i] = clist[ii]
        
        elif line_val == 'shared':
            for i in range(n):
                ii = self.group[i]-1
                if self.mooringGroups[ii]['shared']:
                    colors[i] = 'tab:cyan'
                else:
                    colors[i] = 'tab:pink'
            colorbar = 0
            
        elif line_val == 'stiffness':
            colors, s_m = getLineColors(line_k)              # get colors corresponding to each line type
            colorbar_label = 'Effective stiffness (kN/m)'
            line_var = 'k'


        elif line_val == 'weight':
            colors, s_m = getLineColors(line_w)
            colorbar_label = 'Wet weight (N/m)'
            line_var = 'weight'
            
        elif line_val == 'mass':
            colors, s_m = getLineColors(line_m)
            colorbar_label = 'Wet mass (kg/m)'
            line_var = 'weight'
            
        elif line_val == 'tension':
            colors, s_m = getLineColors(line_t)
            colorbar_label = 'Horizontal tension (MN)'
            line_var = 'T'
            
        elif line_val == 'kt_over_k':
            colors, s_m = getLineColors(line_kt_k)
            colorbar_label = 'Line kt/k (-)'
            line_var = 'kt/k'
        
        elif line_val == 'cost':
            colors, s_m = getLineColors(line_cost)
            colorbar_label = 'Line cost [?]'
            line_var = 'cost'
            
            
        else:
            raise ValueError('Incorrect line_val given')
        
        
        # set up axes        
        if ax == None:                    # if no axes passed in, make a new figure
            fig, ax = plt.subplots(1,1, figsize=figsize, constrained_layout=True)
        else:
            fig = ax.get_figure()         # otherwise plot on the axes passed in

        
        
        # # plot each mooring line, colored differently for each line type 
        for i in range(n):
            
            # shousner: I don't understand how the j var found an integer
            #j = int(Line.type[4:])-1     # index of LineType
            ii = self.group[i]-1
            
            shared = self.mooringGroups[ii]['shared']==1
            
            rA = self.rA[i]
            rB = self.rB[i]
            


            if not (show_lines=="none" or (show_lines=="anch" and shared) or (show_lines=="shared" and not shared)):
                if self.boundary[i]:
                    l, = ax.plot([rA[0], rB[0]],[rA[1], rB[1]], color=colors[i], linestyle='--', lw=thicks[i])
                else:
                    l, = ax.plot([rA[0], rB[0]],[rA[1], rB[1]], color=colors[i], linestyle=styles[i], lw=thicks[i])
                if 'l' in labels:
                    coord = 0.5*(rA + rB)  # position label at midpoint between line ends
                    ax.text(coord[0], coord[1], f"{i+1}", bbox=dict(facecolor='none', edgecolor='k'))                      
                
                
        # display colorbar   
        if not line_val in ["uniform", "two", "groups"]:
            if colorbar==2:
                if 'cbax' in kwargs:
                    #if isinstance(colorbar, plt.Axes):       # if an axes has been passed in via colorbar
                    plt.gcf().colorbar(s_m, label=colorbar_label, ax=kwargs["cbax"], shrink=0.4, aspect=12)  # put the colorbar on that axes
                else:
                    raise ValueError("An axes to put the colorbar beside must be provided (as 'cbax') when colorbar=2")
                    
            elif colorbar == 1:                           # make a regular colorbar on the current axes                
                cax = plt.gca().inset_axes([1.1, 0, 0.05, 1])
                plt.gcf().colorbar(s_m, label=colorbar_label, cax=cax)
            elif colorbar == 0:                           # don't make a colorbar
                pass         
            else:
                raise ValueError("Unrecognized entry for colorbar when calling plot2d.")
            
            
        #plot each platform and anchor
        #for i in range(self.coords.shape()[0]):    # loop through each platform or anchor 
        for i in range(len(self.intraMat)):
            
            # platform
            if i < self.nPtfm:
                ax.plot(self.coords[i,0], self.coords[i,1], 'ko', markersize = 6)
                
                # plot watch circles if requested
                if watch_circles > 0:
                    if not hasattr(self, 'xi_sweep'):
                        raise Exception("In a Linear System, windsweep must be called before trying to plot watch circles.")
                    
                    scale = watch_circles
                    
                    center_x = self.coords[i,0]
                    center_y = self.coords[i,1]
                    
                    # plot calculated displacement envelopes
                    #disps_x = self.xi_sweep[:,2*i] * scale 
                    #disps_y = self.xi_sweep[:,2*i+1] * scale                     
                    #ax.plot(center_x + disps_x, center_y + disps_y,'k',lw=1.5, alpha = 0.6)
                    
                    watch_circle_coords = np.column_stack([self.xi_sweep[:,2*i  ]*scale + center_x, 
                                                           self.xi_sweep[:,2*i+1]*scale + center_y])
                    
                    # ax.add_patch(mpl.patches.Polygon(watch_circle_coords, lw=1, ec=[(self.depth - np.max(self.intersectZ))/self.depth,0,0,1.0], fc=[0,0,0,0]))
                    ax.add_patch(mpl.patches.Polygon(watch_circle_coords, lw=1, ec=[0,0,0,0.6], fc=[0,0,0,0]))
                    
                    # Plot the boundaries
                    r = self.xmax * scale 
                    
                    thetas = np.linspace(0, 2 * np.pi, 201)
                    xs, ys = (np.array(()),np.array(()))
                    for theta in thetas:
                        xs = np.append(xs,r * np.cos(theta) + center_x)
                        ys = np.append(ys,r * np.sin(theta) + center_y)
                    ax.plot(xs,ys,'r--', lw=1, alpha = 0.5)
                                        
                
                if 't' in labels:
                    coord = np.array([self.coords[i,0], self.coords[i,1],0]) + np.array([250, 150,0])
                    ax.text(coord[0], coord[1], f"T{i+1}", fontweight='bold')#, bbox=dict(facecolor='none', edgecolor='k', boxstyle='circle,pad=0.3'))
                elif 'p' in labels:
                    coord = np.array([self.coords[i,0], self.coords[i,1],0]) + np.array([200, 200,0])
                    ax.text(coord[0], coord[1], str(i+1))#, bbox=dict(facecolor='none', edgecolor='k', boxstyle='circle,pad=0.3'))

            # anchor
            elif i >= self.nPtfm and self.intersectZ[i] > 0:
                ax.plot(self.coords[i,0], self.coords[i,1], 'ko', markersize=6, mfc='cyan')   
                if 'h' in labels:
                    coord = np.array([self.coords[i,0], self.coords[i,1],0]) + np.array([200, 200,0])
                    ax.text(coord[0], coord[1], "Hbrd"+str(i+1-self.nPtfm), bbox=dict(facecolor='none', edgecolor='c', boxstyle='circle,pad=0.3'))                
            else:
                ax.plot(self.coords[i,0], self.coords[i,1], 'ko', markersize=6, mfc='none')   
                if 'a' in labels:
                    coord = np.array([self.coords[i,0], self.coords[i,1],0]) + np.array([200, 200,0])
                    ax.text(coord[0], coord[1], "Anch"+str(i+1-self.nPtfm), bbox=dict(facecolor='none', edgecolor='k', boxstyle='circle,pad=0.3'))
        
        
        
           
        #Uncomment to hard code labels
        #plt.legend(shadow=True, loc="upper left")  <<< still need to sort out legend
        if line_val == 'groups':
        
            from matplotlib.lines import Line2D
            handles = []
            for i in range(len(self.mooringGroups)):
                handles.append(Line2D([0], [0], label=f'Group {i}', color=clist[i]))
        
            plt.legend(handles=handles)
        
        ax.set_aspect('equal')
        
        if not show_axes:
            ax.axis('off')
        
        if len(title) > 0:
            plt.title(title)

        # if this made a new figure, return its handles
        #if axes == None:

        # Plot the WEA boundaries if given:
        if wea:
            x, y = wea.exterior.xy
            plt.plot(x, y, color='green', linestyle='--')

        return fig, ax
        

    
        
    
    #note: method analyzeWind(self) made some nice plots with binning offsets by direction according to severity
    #      See code prior to March 20 for this capability.
        
    
    
    def eigenAnalysis(self, plot=0, M=1e6, deg=0):
        '''
        deg
            first desired direction of turbine 1 for organizing eigenmodes. Default 0 (deg)
        '''
        
        v1 = [[ np.cos(np.radians(deg))], [np.sin(np.radians(deg))]]
        v2 = [[-np.sin(np.radians(deg))], [np.cos(np.radians(deg))]]
        
        
        #Take code and ideas from patrick and run eigen analysis
        
        #Define and Populate Mass Matrix
        if np.isscalar(M):
            self.MassMatrix = np.zeros((2*self.nPtfm, 2*self.nPtfm))
            np.fill_diagonal(self.MassMatrix, M)
        else: # if it's a full matrix
            if M.shape != self.SystemStiffness.shape:
                S = self.SystemStiffness
                raise ValueError(f'The mass matrix is of size {M.shape[0]}x{M.shape[1]} and needs to be of size {S.shape[0]}x{S.shape[1]}')
            else:
                self.MassMatrix = M
            
        
        #Calculate eigenvalues and eigenvectors  (note: eigenvectors or mode shapes are *columns*)
        self.eigenvalues, self.eigenvectors = np.linalg.eig(np.matmul(np.linalg.inv(self.MassMatrix), self.SystemStiffness))
        
        #Calculate Natrual Frequency
        self.nat_freqs = np.real(np.sqrt(self.eigenvalues))
        
        #Find Indicies to sort from smallest to largest
        sort_indices = np.argsort(self.nat_freqs)
        
        
        #Use sort_indices to sort natrual frequency, eigenvalues and eigenvectors 
        self.nat_freqs = np.array([self.nat_freqs[i] for i in sort_indices])   
        self.eigenvalues = np.array([self.eigenvalues[i] for i in sort_indices])     
        self.eigenvectors = np.transpose([self.eigenvectors[:,i] for i in sort_indices]) 
        self.periods = np.pi * 2 / self.nat_freqs
        
        #Round periods to 5 decimals
        self.periods = np.round(self.periods,5)
        
        #Pretty plots 
        #Loop through each eigen vector
        #Re-orient eigenvector pairs to look nice 
        for period in set(self.periods): 
            count = np.count_nonzero(self.periods == period)
            
            #If there are duplicates, re-order them so that they are orthogonal 
            if count == 2:
                #print('re-normalizing modes for period {}'.format(period))
                
                #Get indicies
                ind = [i for i in range(len(self.periods)) if self.periods[i] == period]
                eigs = np.empty([len(self.periods),0])
                for i in ind:
                    #eigs is a nxc matrix where n is the number of DOF and c is the period count
                    eigs = np.column_stack((eigs,self.eigenvectors[:,i]))
                    #print(eigs)
                    
                #Make orthogonal
                eigs = scipy.linalg.orth(eigs)
                
                #A elegant bit of linear algebra used to get desired directions
                #desired directions
                dir1 = np.array(v1) 
                dir2 = np.array(v2) 
                dirs = np.append(dir1,dir2,axis = 1)
                
                #current directions
                current = eigs[:2,:2]
                
                #get the weights needed
                #weights1 = np.matmul(np.linalg.inv(current),dir1)
                #weights2 = np.matmul(np.linalg.inv(current),dir2) 
                weights = np.matmul(np.linalg.inv(current),dirs) 
                
                #trasform eigs using weights
                eigs = np.matmul(eigs,weights)
                
                #Update variables 
                for i in ind:  
                    self.eigenvectors[:,i] = eigs[:,0]
                    eigs = np.delete(eigs, 0, axis = 1)
        
        
        #Plot Things
        if plot == 1:
            
            def closestDivisors(n):
                a = round(np.sqrt(n))
                while n%a > 0: a -= 1
                return int(a),int(n//a)
            
            rows, cols = closestDivisors(len(self.eigenvalues))
            
            fig,ax = plt.subplots(rows,cols)
            #Loop through each eigen values
            for ind in range(len(self.eigenvalues)):
                # np.unravel_index() allows linear indexing of a 2D array, like in matlab
                if len(np.shape(ax)) == 2:  
                    plt_ind = np.unravel_index(ind,[rows,cols],'F')
                else:
                    plt_ind = ind
                
                self.eigenPlot(ind, ax=ax[plt_ind])
                
                '''
                #size eigenvector #TODO: Change this based on the size of the plot
                eigenvector = self.eigenvectors[:,ind] * 1500
                
                
                #Loop through each point
                for i in range(self.nPtfm):
                    r = np.array([self.coords[i,0], self.coords[i,1]])
                    
                    ax[plt_ind].plot(r[0], r[1], 'ko', markersize = 2)
                    
                    ax[plt_ind].quiver(np.array(r[0]),
                              np.array(r[1]),
                              np.array(eigenvector[2*i]),
                              np.array(eigenvector[2*i+1])
                              ,units='xy', scale=1)
                
                ax[plt_ind].set_aspect('equal')
                ax[plt_ind].set_xticks([])
                ax[plt_ind].set_yticks([])
                ax[plt_ind].set_axis_off()
                ax[plt_ind].set_xlim(ax[plt_ind].get_xlim()[0] - 1000, ax[plt_ind].get_xlim()[1] + 1000)
                ax[plt_ind].set_ylim(ax[plt_ind].get_ylim()[0] - 1000, ax[plt_ind].get_ylim()[1] + 1000)
                ax[plt_ind].set_title('T = {:.3f}s'.format(self.periods[ind]))
                '''
        #TODO
        #Nice way to make them perpindicular
        #Animation of them rotating
        
        #Collective Mode
        #is there a slick way to add the two together 
        #anti-collective mode
        
        #kx collective 
        #kt collective
       
        #kx anti-collective
        #kt anti-collective
            
        
        #1. Singlue turbine - 2 modes, same period 
        #2. 2 turbines - 4 modes, all permuatataions of kt/kx, col/anti-col
            # no two modes have the same peiord. 
        #3. 4 turbine square - 8 modes. 2 copies of all perumatations of kt/kx
        #4  3 turbines triangle - 6 modes. Things get weird because all motion
            # includes combinations of kt and kx
        #5 6 turbine hexagon - 12 modes. Theres 120 deg symmetry so I expect
            # atleast 3 copies of all permutations. I expect these permutations
            # to look similar to the triangle
        #6 7 turbine hexagon - 14 modes. There are 120 deg symmetry so I again 
            # expect 3 copies of all permutations. Things get weird because 14
            # is not divisable by 3, so I don't quite know where that leads 
    
    
    # new method to plot any given eigenmode
    def eigenPlot(self, ind, ax=None, period=True, figsize=(5,5), length=800, color='k', line_width=4, head_size=3):
        '''Plot an eigenmode of a Linear System.  i is the mode index. eigenAnalysis must be called first.'''

        if ax == None:
            fig, ax = plt.subplots(1,1, figsize=figsize)
        else:
            fig = ax.get_figure()

        # get largest length of an eigenvector horizontal motion vector
        maxLength = max(np.hypot(self.eigenvectors[0::2,ind], self.eigenvectors[1::2,ind])) 

        # scale eigenvector to desired length
        eigenvector = self.eigenvectors[:,ind] * length/maxLength
        
        #Loop through each point
        for i in range(self.nPtfm):
            r = np.array([self.coords[i,0], self.coords[i,1]])
            
            ax.plot(r[0], r[1], 'ko', markersize=2)
            
            ax.quiver(np.array(r[0]),
                               np.array(r[1]),
                               np.array(eigenvector[2*i]),
                               np.array(eigenvector[2*i+1]), 
                               units='dots', width=line_width, color=color, zorder=10,
                               headwidth=head_size, headlength=head_size, headaxislength=head_size,
                               angles='xy', scale_units='xy', scale=1)
                               #units='xy', scale=1)
        
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        #ax.set_xlim(ax.get_xlim()[0] - length, ax.get_xlim()[1] + length)
        #ax.set_ylim(ax.get_ylim()[0] - length, ax.get_ylim()[1] + length)
        if period:
            ax.set_title('T = {:.3f}s'.format(self.periods[ind]))
    
        return fig, ax
    
    # :::::::::::: methods below here to eventually be moved to separate code :::::::::::::::
    
    def calibrate(self, percent_droop=50, percent_drag=60, plots=0): 
        
        
        def laylength_eval(X, args):  
            '''Function to be solved for lay length target'''
            
            # Step 1. break out design variables and arguments into nice names
            L = X[0]
            [Xf,Zf,EA,W] = args
             
            # Step 2. do the evaluation (this may change mutable things in args)
            (Fx1, Fy1, Fx2, Fy2, info) = mp.catenary(Xf, Zf, L, EA, W)
        
            # Step 3. group the outputs into objective function value and others
            Y = info["LBot"]                    # objective function
            oths = dict(message="none")         # other outputs - returned as dict for easy use
            
            return np.array([Y]), oths, False
            
        
        def laylength_step(X, args, Y, oths, Ytarget, err, tols, iter, maxIter):
            '''Stepping functions for achieving lay length target'''
            
            L = X[0]
            [Xf,Zf,EA,W] = args
            LBot = Y[0]
                
            if LBot <= 0:    # if no seabed contact, increase line length by 10% of spacing
                dL = 0.1*Xf
                
            else:            # get numerical derivative
                deltaL = 2*tols[0]                                                  # step size
                (Fx1, Fy1, Fx2, Fy2, info) = mp.catenary(Xf, Zf, L+deltaL, EA, W)  # evaluate LBot in perturbed case
                LBot2 = info["LBot"] 
                dLBot_dL = (LBot2-LBot)/deltaL                                  # derivative
            
                # adjust as per Netwon's method
                dL = -err[0]/dLBot_dL
                
            return np.array([dL])   # returns dX (step to make)
        
        
        def droop_eval(X, args):  
            '''Function to be solved for shared droop target'''
            
            # Step 1. break out design variables and arguments into nice names
            L = X[0]
            [Xf,Zf,EA,W,cb] = args
             
            # Step 2. do the evaluation (this may change mutable things in args)
            (Fx1, Fy1, Fx2, Fy2, info) = mp.catenary(Xf, Zf, L, EA, W, cb)
        
            # Step 3. group the outputs into objective function value and others
            Y = info["Zextreme"]                # objective function
            oths = dict(message="none")         # other outputs - returned as dict for easy use
            
            return np.array([Y]), oths, False
            
            
        def droop_step(X, args, Y, oths, Ytarget, err, tols, iter, maxIter):
            '''Stepping functions for achieving shared droop target'''
            
            L = X[0]
            [Xf,Zf,EA,W,cb] = args
            Zmin = Y[0]
                
            if Zmin >= -tols[0]:    # if nearly no droop at all (in which case derivative will be near zero), add length
                dL = 0.1*Xf
                
            else:            # get numerical derivative
                deltaL = 2*tols[0]                                                     # step size
                (Fx1, Fy1, Fx2, Fy2, info) = mp.catenary(Xf, Zf, L+deltaL, EA, W, cb)  # evaluate droop in perturbed case
                Zmin2 = info["Zextreme"] 
                dZmin_dL = (Zmin2-Zmin)/deltaL                                      # derivative
            
                # adjust as per Netwon's method
                dL = -err[0]/dZmin_dL
                
            return np.array([dL])   # returns dX (step to make)
        
        
        # initialize 3D coordinates (can probably go in init)
        coords = np.zeros([len(self.coords),3])
        for j in range(len(self.coords)):
            if j < self.nPtfm:
                coords[j,:] = np.array([self.coords[j][0], self.coords[j][1], 0])
            else:
                coords[j,:] = np.array([self.coords[j][0], self.coords[j][1], -self.depth+self.intersectZ[j]])
        
        # Just need to get an initial Fx to send to LineDesign -> assume it's a simple catenary for now
        
        # Loop through each mooring line and update its properties
        for ii in range(np.max(self.intraMat)):
            mg = self.mooringGroups[ii]             # the mooring group shortcut
            i  = self.group.index(ii+1)             # the index of endA/endB where this mooring line object occurs first
            
            rA = coords[self.endA[i]]               # starting coordinates of the line
            rB = coords[self.endB[i]]               # ending coordinates of the line
            
            # initialize line parameters
            Xf = np.linalg.norm((rA - rB)[0:2])     # horizontal distance (a.k.a. L_xy)
            Zf = np.linalg.norm((rA - rB)[2  ])     # vertical distance (aka depth)
            L = 1.2*np.hypot(Xf, Zf)                # unstretched line length (design variable)
            EA = 1232572089.6                       # EA value of 120mm chain
            W = 2456.820077481978                   # W value of 120mm chain
            cb = -self.depth                        # the distance down from end A to the seabed
            
            
            # if anchored, adjust line length to have line on seabed for percent_drag of spacing
            if mg['shared']==0: # anchored line
                X0         = [L]
                LBotTarget = [percent_drag/100*Xf]                  # target lay length is percent_drag of horizontal anchor spacing
                args       = [Xf,Zf,EA,W]                           # the other (constant) arguments needed by catenary   
                X, Y, info = msolve.dsolve2(laylength_eval, X0, Ytarget=LBotTarget, step_func=laylength_step, args=args, maxIter=20)
            
                # set line length to the solved value
                L = X[0]
                # Call catenary function with resized line length
                (Fx1, Fy1, Fx2, Fy2, info) = mp.catenary(Xf, Zf, L, EA, W, plots=plots)
            
            # if shared, adjust the line length to have the line droop for percent_droop of spacing
            elif mg['shared']==1 and not mg['net']: # shared line
                X0          = [L]
                DroopTarget = [-percent_droop/100*self.depth - rB[2]]  # target droop elevation relative to fairlead (from percent_droop of depth)
                args        = [Xf,Zf,EA,W,cb]                       # the other (constant) arguments needed by catenary   
                X, Y, info  = msolve.dsolve2(droop_eval, X0, Ytarget=DroopTarget, step_func=droop_step, args=args, maxIter=20)
                
                # set line length to the solved value
                L = X[0]

                # Call catenary function with resized line length                
                (Fx1, Fy1, Fx2, Fy2, info) = mp.catenary(Xf, Zf, L, EA, W, cb, plots = plots)
            
            elif mg['net'] and not mg.get('tendON', False): # hybrid line (without tendon)
                Xf *= 2
                Zf *= 0
                X0  = [2*L]
                DroopTarget = [rA[2]]
                args = [Xf,Zf,EA,W,cb]
                X, Y, info  = msolve.dsolve2(droop_eval, X0, Ytarget=DroopTarget, step_func=droop_step, args=args, maxIter=20)
                
                # set line length to the solved value
                L = X[0]

                # Call catenary function with resized line length                
                (Fx1, Fy1, Fx2, Fy2, info) = mp.catenary(Xf, Zf, L, EA, W, cb, plots = plots)

            elif mg['net'] and mg.get('tendON', False): # hybrid line (with tendon) - assume shared with zero droop
                X0          = [L]
                DroopTarget = [0]  # target droop elevation relative to fairlead (from percent_droop of depth)==
                args        = [Xf,Zf,EA,W,cb]                       # the other (constant) arguments needed by catenary   
                X, Y, info  = msolve.dsolve2(droop_eval, X0, Ytarget=DroopTarget, step_func=droop_step, args=args, maxIter=20)
                
                # set line length to the solved value
                L = X[0]

                # Call catenary function with resized line length                
                (Fx1, Fy1, Fx2, Fy2, info) = mp.catenary(Xf, Zf, L, EA, W, cb, plots = plots)
            
            Fx = np.abs(info['HF'])                     # horizontal tension component at fairlead [N]
            Kx = info['stiffnessB'][0,0]                # effective horizontal stiffness at fairlead [N/m]
            kt_over_k = Fx / Kx / Xf                    # kt/Kx = Fx/L_xy / Kx
            
            
            if plots == 2:
                plt.title('Catenary Line Profile for Mooring Line: {}'.format(ii+1))
                
                print('Force for Mooring Line {}:   {}'.format(ii, Fx))
                print('Stiffness for Mooring Line {}:   {}'.format(ii, Kx))
                print('\n{} '.format('kt/k for Mooring Line {}:   {}'.format(ii, kt_over_k)))
            
            
            # Update mooringGroups dictionary values with 
            mg['ten__w'] = Fx/W
            mg['kl__w'] = Kx/W
            mg['kt__kl'] = kt_over_k
            mg['L'] = L     # <<<< this is a shortcut that should be done outside of LinearSystem in future
    
    
    
    def makeMoorPySystem(self):
        ''' sets up a *very simple* MoorPy system using points for the FOWTS, and no bodies '''
        
        ms = mp.System(depth=self.depth)
        
        
        # make free points for platforms (assumed at z=0) and fixed points for anchors 
        for i in range(len(self.coords)):
            if i < self.nPtfm:
                #ms.addPoint(0, np.hstack([self.coords[i][:2], 0]), m=1e9, v=2e9*ms.rho)   # add buoyancy as double the mass - this should make it equilibrate at z=0
                ms.addPoint(0, np.hstack([self.coords[i][:2], 0]), DOFs=[0,1])   # specify as free to move in x and y only
            else:
                ms.addPoint(1, np.hstack([self.coords[i][:2], -self.depth]))

        
        # also add some generic line types, one for each line grouping as defined in the entries of intraMat
        for i in range(np.max(self.intraMat)):    
            shared = i+1 in self.intraMat[:self.nPtfm, :self.nPtfm]    # detect if it's a shared line (True) or not (False, assumed anchored)
            massden = self.mooringGroups[i]['w']/9.81
            
            ms.lineTypes[f"type{i+1}"] = mp.LineType(f"type{i+1}", 0.0, massden, 1.0e15) #, shared=shared) 
        
        
        # make lines using adjacency matrix
        linenum = 1
        for i in range(len(self.coords)):
            for j in range(i):
                k = np.int(self.intraMat[i,j])   # the entry in the intraMat corresponds to the line type number (starting from 1)
                if k > 0:
                    ml = self.mooringGroups[k-1]
                    ms.addLine(ml['L'], f'type{k}')
                    ms.pointList[i].attachLine(linenum, 1)
                    ms.pointList[j].attachLine(linenum, 0)
                    linenum = linenum + 1
        '''  this should be done to coords if it happens anywhere    
        if self.center:
            cx = np.mean([point.r[0] for point in ms.pointList])
            cy = np.mean([point.r[1] for point in ms.pointList])
            ms.transform(trans=[-cx, -cy])
        '''
       
        ms.initialize()
        
        return ms
            
                
    def getYawStiffness(self, rf):
        yawStiffness = np.zeros(self.nPtfm)
        for i in range(self.nPtfm):
            connectedLines = np.where(np.abs(self.StructureMatrix[i*2:i*2+1, :]) > 0)[1]
            for j in connectedLines:
                if self.boundary[j]:
                    lineYawStiffness = self.mooringGroups[self.group[j] - 1]['ten']/(2*self.l[j]) * rf**2    
                else:
                    lineYawStiffness = self.mooringGroups[self.group[j] - 1]['ten']/self.l[j] * rf**2
                yawStiffness[i] += lineYawStiffness
        return yawStiffness

       
    def removeRedundantGroups(self):
        '''
        this method: 
        1) removes any zero weight groups (lower than 1% of mean weight)
        2) merge mooring groups that are similar to one another (mooring groups that their w are 5% different from the w range), and 
        3) reformulates the system and its mooring groups
        '''
        # calculate the mean weight for the LinearSystem:
        wSum  = sum(mg['w'] for mg in self.mooringGroups)
        wMax  = max(mg['w'] for mg in self.mooringGroups)
        wMin  = min(mg['w'] for mg in self.mooringGroups)
        wMean = wSum / len(self.mooringGroups)
        
        # find the indices of mooring groups to keep and the ones to remove
        remove_indices = [i for i, mg in enumerate(self.mooringGroups) if mg['w'] < 0.01 * wMean]
        keep_indices = list(set(range(len(self.mooringGroups))) - set(remove_indices))
        
        # Update mooringGroups and profileMap
        self.mooringGroups = [self.mooringGroups[i] for i in keep_indices]
        self.profileMap = [self.profileMap[i] for i in keep_indices]

        # create a mapping from old indices to new indices
        index_mapping = {old: new for new, old in enumerate(keep_indices)}        

        # remove lines that belong to redundant mooring groups
        new_l, new_u, new_endA, new_endB, new_rA, new_rB, new_angA, new_angB, new_boundary, new_group = [], [], [], [], [], [], [], [], [], []
        
        for i in range(len(self.l)):
            if self.group[i] - 1 in keep_indices:  # subtract 1 because group starts at 1, not 0
                new_group.append(index_mapping[self.group[i] - 1] + 1)
                new_l.append(self.l[i])
                new_u.append(self.u[i])
                new_endA.append(self.endA[i])
                new_endB.append(self.endB[i])
                new_rA.append(self.rA[i])
                new_rB.append(self.rB[i])
                new_angA.append(self.angA[i])
                new_angB.append(self.angB[i])
                new_boundary.append(self.boundary[i])

        # remove from intra-cell adjacency matrix
        for remIdx in remove_indices:
            for i, mg in enumerate(self.group):
                if mg == remIdx + 1:
                    remA = self.endA[i]
                    remB = self.endB[i]
                    self.intraMat[remA, remB], self.intraMat[remB, remA] = 0, 0

        unique_intra_groups = sorted(np.unique(self.intraMat[self.intraMat > 0]))
        intra_group_mapping = {old: new for new, old in enumerate(unique_intra_groups, start=1)}
        for old, new in intra_group_mapping.items():
            self.intraMat[self.intraMat == old] = new        
        # update the properties with the new filtered lists
        self.l = new_l
        self.u = new_u
        self.endA = new_endA
        self.endB = new_endB
        self.rA = new_rA
        self.rB = new_rB
        self.angA = new_angA
        self.angB = new_angB
        self.boundary = new_boundary
        self.group = new_group
        self.nLines = len(self.l)
        self.StructureMatrix = np.zeros([2*self.nPtfm, self.nLines])     # rows: DOFs; columns: lines
        
        # merge mooring Groups
        removeIndex = []
        for i, mg1 in enumerate(self.mooringGroups[:-1]):
            for j in range(i+1, len(self.mooringGroups)):
                mg2 = self.mooringGroups[j]
                # Check following conditions:
                #   if the difference in weight is minimial,
                #   if they both have the same shared map,
                #   and if they have the same length:                
                con1 = np.abs(mg2['w'] - mg1['w'])/(wMax - wMin) < 0.05
                con2 = self.profileMap[i]==self.profileMap[i+1]
                con3 = np.round(mg1['l'], 2)==np.round(mg2['l'], 2)
                if con1 and con2 and con3:  
                    self.mooringGroups[i]['w'] = (mg1['w'] + mg2['w']) / 2
                    removeIndex.append(j)
                    self.group = [i+1 if g==j+1 else g for g in self.group]
                    idx1, idx2 = np.where(self.intraMat==j+1)
                    self.intraMat[idx1, idx2] = i+1
        
        if removeIndex:
            for idx in sorted(np.unique(removeIndex), reverse=True):
                del self.mooringGroups[idx]
                del self.profileMap[idx]
            
            unique_groups = sorted(np.unique(self.group))
            unique_intra_groups = sorted(np.unique(self.intraMat[self.intraMat > 0]))
            group_mapping = {old: new for new, old in enumerate(unique_groups, start=1)}
            intra_group_mapping = {old: new for new, old in enumerate(unique_intra_groups, start=1)}
            self.group = [group_mapping[g] for g in self.group]
            
            for old, new in intra_group_mapping.items():
                self.intraMat[self.intraMat == old] = new
        
        for j in range(self.nLines):  
            if self.endA[j] < self.nPtfm: # only if not an anchor
                self.StructureMatrix[self.endA[j]*2    , j] =  self.u[j][0]
                self.StructureMatrix[self.endA[j]*2 + 1, j] =  self.u[j][1]
                
            if self.endB[j] < self.nPtfm: # only if not an anchor
                self.StructureMatrix[self.endB[j]*2    , j] = -self.u[j][0]
                self.StructureMatrix[self.endB[j]*2 + 1, j] = -self.u[j][1]
        
        # Check if any row/column in the intraMat is empty, delete it, and delete the corresponding self.coords:
        empty_rows = np.where(~self.intraMat.any(axis=1))[0]
        empty_cols = np.where(~self.intraMat.any(axis=0))[0]
        empty_indices = np.unique(np.concatenate((empty_rows, empty_cols)))
        if len(empty_indices) > 0:
            # Remove empty rows and columns from intraMat
            self.intraMat = np.delete(self.intraMat, empty_indices, axis=0)
            self.intraMat = np.delete(self.intraMat, empty_indices, axis=1)
            self.coords = np.delete(self.coords, empty_indices, axis=0)
            self.intersectZ  = np.delete(self.intersectZ, empty_indices, axis=0)

        # Reassign the 'type' of each mooring group after deletion
        for i, group in enumerate(self.mooringGroups):
            group['type'] = i + 1  # Reassign types from 1 to len(mooringGroups)

        self.preprocess()
        self.optimize()

# ------- test script

if __name__ == '__main__':
    
    import Array as array
    
    from moorpy.helpers import printVec, printMat
    
    # specify the array layouts and their parameters
    T     = 2000.
    A     = 1200.
    depth =  600.
    
    '''
    # ----- old examples -----
    
    #coords, intraMat, nPtfm, name = array.layout_pair_4_anchs(T, A, deg=120)
    #coords, intraMat, nPtfm, name = array.layout_triangle_3_anchs(T, A)
    coords, intraMat, nPtfm, name = array.layout_1_square_8_anchs(T, A)
    
    sys = LinearSystem(coords, intraMat, nPtfm, depth=600., fmax=1e6, 
        xmax=0.1*min(T,A))
    
    
    
    sys.preprocess()
    
    q= sys.optimize()
    print(q)
    sys.plot2d(watch_circles=1, line_val="stiffness")
    
    #sys.eigenAnalysis(plot=1)
    
    #sys.updateDesign()
    
    
    # ----- newer more advanced example -----
    '''
    print("New LinearSystem example")
    # def __init__(self, coords, intraMat,  nPtfm, interMats=None, 
    #     interCoords=None, inits=None, profileMap=None, intersectZ=None, 
    #     rFair=0, zFair=0, depth=600., fmax=1e6, xmax=40.0, plots=0, nonlin=1.0, 
    #     center=True, old_mode=True):
    
    #coords, intraMat, nPtfm, name = array.Grid3x3(T, A)
    #coords, intraMat, nPtfm, name = array.Fat_Hexagon(T, A, fathexagontype='min_linetypes')
    coords, intraMat, nPtfm, name = array.Square(T, A, type='water-strider')
    
    
    mooringGroupDict = [
        dict(w=1500, ten=100000, kl=10000, kt=50, shared=False),
        #dict(w=1500, ten=100000, kl=10000, kt=50, shared=True),
        dict(w=1500, ten=100000, kl=10000, kt=50, shared=True)] #,
    '''
        dict(w=1500, ten=100000, kl=10000, kt=50, shared=True),
        dict(w=1500, ten=100000, kl=10000, kt=50, shared=True),
        dict(w=1500, ten=100000, kl=10000, kt=50, shared=True),
        dict(w=1500, ten=100000, kl=10000, kt=50, shared=True),
        dict(w=1500, ten=100000, kl=10000, kt=50, shared=True),
        dict(w=1500, ten=100000, kl=10000, kt=50, shared=True),
        dict(w=1500, ten=100000, kl=10000, kt=50, shared=True),
        dict(w=1500, ten=100000, kl=10000, kt=50, shared=True),
        dict(w=1500, ten=100000, kl=10000, kt=50, shared=True)]
    '''
    
    sys = LinearSystem(coords, intraMat, nPtfm, depth=600., fmax=1e6, 
        xmax=0.1*min(T,A), inits=mooringGroupDict, old_mode=False)
    
    
    sys.preprocess()
    
    
    sys.getSystemStiffness()
    sys.windsweep() # figure out watch circles
    
    sys.plot2d(watch_circles=4, line_val="stiffness")
    
    
    # now try a stiffness optimization
    
    sys.optimize2(display=2)
    sys.plot2d(watch_circles=4, line_val="stiffness")
    
    '''
    
    print("New LinearSystem example - with inter-array shared lines")
    # def __init__(self, coords, intraMat,  nPtfm, interMats=None, 
    #     interCoords=None, inits=None, profileMap=None, intersectZ=None, 
    #     rFair=0, zFair=0, depth=600., fmax=1e6, xmax=40.0, plots=0, nonlin=1.0, 
    #     center=True, old_mode=True):
    
    coords, intraMat, nPtfm, name = array.Grid3x3(T, A)
    
    # remove anchored lines on side E-W turbines (turbine 3 and 5)
    intraMat[20,3] = 0
    intraMat[14,5] = 0
    
    # make interMats
    interMats = []
    a = np.zeros([9,9])
    a[5,3] = 3   # connect turbines 3 and 5 with a shared line
    print(a)
    interMats.append(np.array(a))
    
    # specify a lateral pattern spaced 6 km apart so shared lines all have same length
    interCoords = [[600,0]]  
    
    mooringGroupDict = [
        dict(w=1500, ten=100000, kl=10000, kt=500, shared=False),
        dict(w=1500, ten=100000, kl=10000, kt=500, shared=True),
        dict(w=1500, ten=100000, kl=10000, kt=500, shared=True)]
    
    
    sys = LinearSystem(coords, intraMat, nPtfm, depth=600., fmax=1e6, 
        xmax=0.1*min(T,A), 
        interMats = interMats, interCoords = interCoords,
        inits=mooringGroupDict, old_mode=False)
    
    sys.getSystemStiffness()
    sys.windsweep() # figure out watch circles
    
    sys.plot2d(watch_circles=1, line_val="stiffness")
    '''
    
    # ----- example with inter-array shared lines! -----
    
    
    '''
    
    import fadesign.conceptual.Cell as cell
    #coords, intraMat, nPtfm, interMats, interCoords = cell.Grid3x3(T, A) # no inter shared lines!
    #coords, intraMat, nPtfm, interMats, interCoords = cell.honeycombPattern(T) # 
    coords, intraMat, nPtfm, interMats, interCoords = cell.grid() # 
    
    breakpoint()
    
    mooringGroupDict = [
        dict(w=1500, ten=100000, kl=10000, kt=50, shared=False),
        dict(w=1500, ten=100000, kl=10000, kt=50, shared=True),
        dict(w=1500, ten=100000, kl=10000, kt=50, shared=True)]
    
    
    sys = LinearSystem(coords, intraMat, nPtfm, interCoords=interCoords, 
        depth=600., fmax=1e6, 
        xmax=0.1*min(T,A), inits=mooringGroupDict, old_mode=False)
    
    sys.getSystemStiffness()
    sys.windsweep() # figure out watch circles
    
    sys.plot2d(watch_circles=1, line_val="stiffness")
    
    '''
    
    
    plt.show()