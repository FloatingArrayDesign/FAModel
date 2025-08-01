import numpy as np
from famodel.helpers import calc_midpoint  
'''
famodel_base contains base classes that can be used for various classes
in the floating array model that relate to entities that act like either
nodes (represented by a single point in space) or edges (that run between
two nodes). The base Node and Edge classes and their methods provide the
functionality for representing 2D positions and relationships between
objects that can be used by many other classes.

Both edges and nodes keep registers of their attachments. This results in
duplicate information, meaning consistency is critical. For consistency 
with MoorPy and MoorDyn, actual attachments or detachments are controlled 
by a Node method. It then also calls the subordinate Edge or Node method. 
The registries are as follows:
'attachments' lists other objects that are bound to the present Node object.
'attached_to' states what object the present Edge end or Node is bound to. 
Each of the above is a dict of dicts, for now, where the key is the id of
the attached object, and the value is a dictionary of more info such as the
object reference, type, and any relative offset...

Linear assemblies of edges and nodes can be grouped into a higher-level
edge object. This allows mooring line sections and connectors to be 
part of an overall mooring line object, which is itself an edge. These 
relationships are tracked with subcomponents lists in the 
higher-level Edge object, and part_of entries in the lower-level node or edge
objects. 

NEW:
Edge subcomponents can now be in any arrangement
The subcomponent(s) at the end of an edge can attach to the node
the edge is attached to, or a sub-node of the node. 
If the super-edge end is attached to a node, 
the sub-objects could be automatically attached to the node, or
left to be manually attached <<<<< which one??
If a super-edge is disconnected, the sub-objects should be disconnected, right?

Is there some overriding logic of the above??

How do we identify/track the sub-objects?
How do we identify/track the attacheble end sub-objets?
Dicts???


The current approach is not top-down. If you connect higher-level objects,
the lower level objects aren't automatically connected. Instead, if you connect
lower objects (including specifying particular positions), the higher-level
connections are also made.


Nodes can attach to nodes as either subortinately or equally...


when attached equally, they reference each other (mirror) in self.attachments
but not in self.attached_to.  The latter is used only for subordinate
connections (i.e. to a higher node, and only to one).
Mutual/equal node attachments are done with the join method, and undone with
the separate method.


When a node or edge of a higher edge is attached to a node of a higher node,
the following rules apply:
The higher level objects can be attached regardless of sub-object attachment.
If any sub-objects are attached, the higher level objects must be attached.
- Detaching the higher level objects must detach the sub objects in the process.
- When all the sub objects are detached, it would be convenient to detach the higher objects.

Can an end of a higher level edge attach to multiple nodes? 
(using a list, corresponding to multiple sub objects at the end)


'''

class Node():
    '''A base class for objects that have a single spatial reference point and
    can be interconnected by edge-type objects.

    Its attachments dictionary has the following format for each entry:

    id  : { 'ref'   :  ,  #
            'r_rel' :  ,  #
            'type'  :     # Node or Edge
    
    
    For grouped/multilevel edges, connections are stored at the highest level,
    so 'ref' would be to the highest level edge object.
    '''

    def __init__(self, id):
        
        self.id = id  # id number or string, used as the key when attached to things
        
        self.attachments = {}  # dictionary listing attached edges or nodes
        # (key is the id, value is attachment object ref and other info)
        
        self.attached_to = None  # whether this object is subordinately bound to another object
        
        self.part_of = None  # whether this object is part of an Edge group
    
        # position/orientation variables
        self.r = np.zeros(2)  # position [m]
        self.theta = 0  # heading [rad] CCW+
        self.R = np.eye(2)  # rotation matrix
    
        # installation dictionary [checks for installation status]
        self.inst = {'mobilized': False,
                     'installed': False}

        
    def isAttached(self, object, end=None):
        '''Check if something is attached to this node. This works for both
        attached edges (the end can be specified), or nodes that are joined
        with this node.
        
        Parameters
        ----------
        object
            The Node or Edge object being attached to this one.
        end
            If an Edge is being considered, which end of it, 'a' or 'b' (optional).
        '''
        
        if object.id in self.attachments:  # See if it's attached
        
            if isinstance(object, Node):  # if it's a node, it's simple
                return True
                
            elif isinstance(object, Edge):
                if end == None:  # if end not specified, it's simple
                    return True
                else:  # otherwise check which end
                    if endToIndex(self.attachments[object.id]['end']) == endToIndex(end):
                        return True  # the end in question is attached
                    else:
                        return False
            else:
                raise Exception('Provided object is not an Edge or Node.')  
        else:
            return False
    
    
    def join(self, object):
        '''Join another node ot this node, in a mutual way.
        This could be multiple connectors within a higher level edge,
        or one connector in a higher level edge (potentially connecting
        to a connector in a higher level node).
        '''
        
        if not isinstance(object, Node):
            raise Exception('Provided object is not a Node.')  
        
        
        # Make sure they're not already attached
        if object.id in self.attachments:
            raise Exception(f"Object {object.id} is already attached to {self.id}")
        if self.id in object.attachments:
            raise Exception(f"{self.id} is already attached to {object.id}")
        
        
        # make sure there isn't some incompatibility in joining these nodes?
        if isinstance(self.part_of, Edge) and isinstance(object.part_of, Edge):
            if not self.part_of == object.part_of:
                raise Exception("Cannot join two nodes that are each part of a different edge")

        # do the mutual joining
        self.attachments[object.id] = dict(obj=object, id=object.id, 
                                           r_rel=np.array([0,0]), type='node')
        
        object.attachments[self.id] = dict(obj=self, id=self.id, 
                                           r_rel=np.array([0,0]), type='node')
        
        # Register the attachment in higher level objects if applicable
        if isinstance(self.part_of, Edge) and isinstance(object.part_of, Edge):
            raise Exception("This attachment would directly connect two higher-level edges to each other, which is not allowed.")
        
        elif isinstance(self.part_of, Edge) and isinstance(object.attached_to, Node):
            end = self.part_of.findEnd(self)
            object.attached_to.attach(self.part_of, end=end, 
                                      r_rel=object.attached_to.attachments[object.id]['r_rel'])
            #self.part_of._attach_to(object.part_of, end=end)
            
        elif isinstance(self.attached_to, Node) and isinstance(object.part_of, Edge):
            end = object.part_of.findEnd(object)
            self.attached_to.attach(object.part_of, end=end, 
                                    r_rel=self.attached_to.attachments[self.id]['r_rel'])
            
        elif isinstance(self.attached_to, Node) and isinstance(object.attached_to, Node):
            raise Exception("This would attach two higher-level nodes, which is not supported.")
    
    
    def isJoined(self):
        '''Check if this node is joined to anything else.'''
        
        for att in self.attachments.values():  # Look through everything attached
            object = att['obj']
            if isinstance(object, Node):  # Only another node could be joined
                if self.id in object.attachments:
                    return True
                    
        # If we've gotten this far, it's not joined with anything
        return False
    
    
    def separate(self, object):
        '''Opposite of join'''
        
        if not isinstance(object, Node):
            raise Exception('Provided object is not a Node.')  
        
        # Make sure they're already attached
        if not object.id in self.attachments:
            raise Exception(f"Object {object.id} is not attached to {self.id}")
        if not self.id in object.attachments:
            raise Exception(f"{self.id} is not attached to {object.id}")
        
        # do the mutual separating
        del self.attachments[object.id]
        del object.attachments[self.id]
        
        # Register the separation in higher level objects if applicable
        if isinstance(self.part_of, Edge) and isinstance(object.attached_to, Node):
            end = self.part_of.findEnd(self)
            object.attached_to.dettach(self.part_of, end=end)
            
        elif isinstance(self.attached_to, Node) and isinstance(object.part_of, Edge):
            end = object.part_of.findEnd(object)
            self.attached_to.detach(object.part_of, end=end)
        
    
    def attach(self, object, r_rel=[0,0], end=None):
        '''Attach something to this node.
        Attachments are noted with an entry in self.attachments
        consisting of a dictionary with ref, type, and end entries.
        
        Parameters
        ----------
        object
            The Node or Edge object being attached to this one.
        r_rel : list
            x and y coordinates of the attachment point [m].
        end
            If an Edge is being attached, which end of it, 'a' or 'b'.
        '''
        
        #object_parent
        #self_parent
        
        '''
        # find top-level edge end if applicable
        if isinstance(object, Node):
            object2, i_end = object.getTopLevelEdge()
        elif isinstance(object, Edge):
            i_end = endToIndex(end, estr='when attaching an edge to a node.')
            object2, _ = object.getTopLevelEdge(i_end)  
        else:
            raise Exception('Provided object is not an Edge or Node.')  
        '''
        # Make sure it's not already attached (note this doesn't distinguish end A/B)
        if object.id in self.attachments:
            # for bridles, the mooring will already be attached to platform
            # for second bridle section
            # need to calculate new r_rel that is average of end points
            if isinstance(object, Edge):
                # pull out relative dist of each point on end to self
                r_rel = self.calculate_r_rel(object,end=end)
                self.attachments[object.id]['r_rel'] = r_rel
                # update end position
                Node.setPosition(self, r=self.r,theta=self.theta)
                # don't need to attach, already attached- just return
                return
            else:   
                raise Exception(f"Object {object.id} is already attached to {self.id}")
        
        
        # Attach the object
        if isinstance(object, Node):  # (object is a node)
        
            if object.attached_to:  # object is already attached to something
                raise Exception("The object being attached is already attached to a higher node - it needs to be detached first.")
        
            self.attachments[object.id] = dict(obj=object, id=object.id, 
                                            r_rel=np.array(r_rel), type='node')
            object._attach_to(self)  # tell it it's attached to this Node
            
        elif isinstance(object, Edge):  # (object is an edge)
            i_end = endToIndex(end, estr='when attaching an edge to a node.')
            
            if object.attached_to[i_end]:   # object is already attached to something
                raise Exception("The object being attached is already attached to a higher node - it needs to be detached first.")
            
            self.attachments[object.id] = dict(obj=object, id=object.id, 
                                            r_rel=np.array(r_rel), type='edge', 
                                            end=i_end)
                                            
            object._attach_to(self, i_end)  # tell it it's attached to this Node
        
        else:
            raise Exception('Unrecognized object type')
        
        # See about attaching higher-level objects (new) (note: r_rel will be neglected at higher level)
        
        if isinstance(object.part_of, Edge):  # attached object is part of an edge
        
            # figure out which end of the edge object corresponds to
            if object in object.part_of.subcons_A and object in object.part_of.subcons_B:
                # there is only one subcomponent, keep end that was passed in
                i_end = end
            elif object in object.part_of.subcons_A:
                end = 0
            elif object in object.part_of.subcons_B:
                end = 1
            else:
                end = -1  # object isn't at the end of the higher level edge so do nothing
                if not self.part_of == object.part_of:  
                    raise Exception("Cannot attach two non-end subcomponents of different edges.")
            
            if self.part_of == None and end > -1:  # this node isn't part of an edge
                if self.attached_to:  # if self is part of a higher node
                    # attach higher edge to higher node
                    self.attached_to.attach(object.part_of, end=end)
                else:
                    # attach higher edge to this node
                    self.attach(object.part_of, r_rel=r_rel, end=end)
                
            else:   # if self is part of a higher level edge
                raise Exception("This attachment would directly connect two higher-level edges to each other, which is not allowed.")
        '''
        elif isinstance(object.attached_to, Node):  # attached object is attached to a higher node
            
            raise Exception("The object being attached is part of a higher node - this operation is not supported.")
            
            if self.part_of:  # self node is part of a higher level edge
            
                # figure out which end of the edge object corresponds to
                if self in self.part_of.subcons_A:
                    end = 0
                elif self in self.part_of.subcons_B:
                    end = 1
                
                # attach higher node and edge
                object.attached_to.attach(self.part_of, end=end)
                
                # attach higher edge to this node
                self.attach(object.part_of, r_rel=r_rel, end=end)
                
            elif if self.attached_to:  # if self is part of a higher node
                Exception("This attachment would directly connect two higher-level nodes to each other, which is not allowed.")
            
            else:  # self has nothing higher, so attach the object's higher level thing to self as well
                self.attach(object.attached_to, r_rel=r_rel, end=end) XXXX
        '''
    
    
    
    def detach(self, object, end=None):
        '''Detach the specified object from this node.
        Will also detach any attached sub-objects
        
        Parameters
        ----------
        object
            The object to detach (its ID should be in the node's attachments dict).
        end
            If an Edge is being detached, which end of it, 'a' or 'b'.
        '''
        
        # Make sure it's attached before trying to detach it
        if not object.id in self.attachments:
            raise Exception(f"Object {object.id} is not attached to {self.id}")
            # this exception could be optionally disabled
        
        # Handle attachment type and the end if applicable, and record
        # the detachment in the subordinate object.
        if isinstance(object, Node):
            object._detach_from()  # tell the attached object to record that it's detached
        
        elif isinstance(object, Edge):
            if end:
                i_end = endToIndex(end, estr='when detaching an edge from a node.')
            else: # otherwise figure out which end is attached
                i_end = self.attachments[object.id]['end']
                
            # Detach this edge from self 
            # This will also detach end subcomponents of the edge from anything.
            object._detach_from(i_end)
        
        else:
            raise Exception('Provided object is not an Edge or Node.')
    
        # Remove it from the attachment registry
        del self.attachments[object.id]
    

    def _attach_to(self, object):
        '''Internal method to update the Node's attached_to registry when
        requested by a node. With an error check.
        
        Parameters
        ----------
        object
            The Node object to attach to.
        '''

        # Make sure it's not already attached to something else
        if self.attached_to:  # True if populated, False if empty dict
            # self.detach() # commented out because I don't think this would work for shared anchors
            # could optionally have a warning or error here
            print(f'Warning: node {self.id} is attached to 2 objects')
            breakpoint()
            print("fix up this scenario in the code somewhere...")
                
        # Add it to the attached_to registry
        self.attached_to = object

        
    def _detach_from(self):
        '''Internal method to update the Node's attached_to registry when
        requested by a node to detach.
        '''
        self.attached_to = None
        

    """
    def getTopLevelObject(self):
        '''If this node is part of a higher-level object, and the request
        corresponds to an end of that higher-level group, return the higher object,
        otherwise return this same object. Can be recursive. 
        A similar method exists for edges.'''
        
        if self.part_of:  # if part of something higher
            supe = self.part_of  # shorthand for the super edge
            if supe.subcomponents[0] == self:  # if this node is at end A of supe
                return supe.getTopLevelEdge(0)  # return supe, and which end
    >>        elif supe.subcomponents[-1] == self:  # if this node is at end B of supe
                return supe.getTopLevelEdge(1)
        else:
            return self, -1  # if not part of something bigger, just return self
    """
    
    def setPosition(self, r, theta=0, force=False):
        '''Set the position of the node, as well as any attached objects.
        
        Parameters
        ----------
        r : list
            x and y coordinates to position the node at [m].
        theta, float (optional)
            The heading of the object [rad].
        force : bool (optional)
            When false (default) it will not allow movement of subordinate objects.
        '''
        
        # Don't allow this if this is part of another object
        if self.attached_to and not force:
            raise Exception("Can't setPosition of an object that's attached to a higher object unless force=True.")
        
        # Store updated position and orientation
        if len(r) > len(self.r): # default r is 2D, but can be adjusted to 3D
            self.r = np.array(r)
        else: # if just a portion of r is being adjusted, only change up to length of initial r
            self.r[:len(r)] = r

        self.theta = theta
        
        # Get rotation matrix...
        if len(self.r) == 2:
            self.R = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
        elif len(self.r) == 3:
            if np.isscalar(theta):
                angles = np.array([ 0, 0, theta])
            elif len(theta)==1:
                angles = np.array([0, 0, theta[0]])
            elif len(theta)==3:
                angles = np.array(theta)
            else:
                raise Exception("theta needs to be length 1 or 3.")
        
            self.R = rotationMatrix(*angles)
        else:
            raise Exception("Length of r must be 2 or 3.")
        
        # Update the position of any attach objects
        for att in self.attachments.values():
        
            if len(self.r) == len(att['r_rel']):
                pass  # all good
            elif len(self.r) == 3 and len(att['r_rel']) == 2:  # pad r_rel with z=0 if needed
                att['r_rel'] = np.hstack([att['r_rel'], [0]])
            else:
                raise Exception("r and attachment r_rel values have mismatched dimensions.")
        
            # Compute the attachment's position
            r_att = self.r + np.matmul(self.R, att['r_rel'])
            
            # set position of any attached node that isn't subordinate to another node
            # (prevents infinite loop of setPositioning for nodes)
            if isinstance(att['obj'], Node):
                if not isinstance(att['obj'].attached_to, Node) or att['obj'].attached_to == self:
                    att['obj'].setPosition(r_att, theta=theta, force=True)
                
            elif isinstance(att['obj'], Edge):
                att['obj'].setEndPosition(r_att, att['end'])
                
    def calculate_r_rel(self,object, end=None):
        '''Calculate the relative distance between node and object
        based on the combined relative distances of the subordinate/
        subcomponent nodes connecting them'''
        if isinstance(object,Edge):
            # pull out subcomponent(s) attached to self at the correct end
            end = endToIndex(end) # find end
            subs = object.subcons_A if end==0 else object.subcons_B           
                
            # go through all end subcomponents of edge at the correct end
            rel_locs = [] # relative location list (in case multiple points at end)
            for sub in subs:    
                # first check if subordinate/subcomponent joined
                att = [att for att in sub.attachments.values() if att['obj'].attached_to==self]
                if len(att)>0:
                    # find attachment of sub that is subordinately connected to self (Node)
                    att = att[0] # just need 1st entry
                    r_rel_att_self = self.attachments[att['id']]['r_rel']
                    r_rel_att_sub = att['obj'].attachments[sub.id]['r_rel']
                    # r_rel of sub to self is r_rel of attachment to self + r_rel of sub to attachment
                    if len(r_rel_att_self) < 3: # pad as needed
                        r_rel_att_self = np.hstack([r_rel_att_self,[0]])
                    if len(r_rel_att_sub) < 3: # pad as needed
                        r_rel_att_sub = np.hstack([r_rel_att_sub,[0]])
                    rel_locs.append(r_rel_att_self + r_rel_att_sub)
                # otherwise, check if directly connected
                elif self.isAttached(object):
                        # if no subordinate/subcomponent connection, should be 
                        # only 1 attachment point at this end
                        return(self.attachments[object.id]['r_rel'])
                else:
                    raise Exception(f'Cannot determine how {self.id} and {object.id} are connected')
            return calc_midpoint(rel_locs)
        elif isinstance(object, Node):
            # node to node - check if 2 subordinates connected
            att = [att for att in object.attachments.values() if self.isAttached(att['obj'])]
            if len(att)>0:
                att = att[0] # just need 1st entry
                # get relative distance of subordinately attached nodes
                r_rel_att_self = self.attachments[att['id']]['r_rel']
                r_rel_att_obj = object.attachments[att['id']]['r_rel'] 
                # r_rel of obj to self is r_rel of attachment to self + r_rel of obj to attachment
                if len(r_rel_att_self) < 3: # pad as needed
                    r_rel_att_self = np.hstack(r_rel_att_self,[0])
                if len(r_rel_att_obj) < 3: # pad as needed
                    r_rel_att_sub = np.hstack(r_rel_att_sub,[0])
                return(r_rel_att_self + r_rel_att_sub)
            # otherwise see if they are directly attached and return r_rel
            elif self.isattached(object):
                return self.attachments[object.id]['r_rel']
            else:
                raise Exception(f'Cannot determine how {self.id} and {object.id} are connected')
        else:
            raise Exception(f'{object} is not a Node or Edge')
              
                    
                



class Edge():
    '''A base class for objects that run between two Node-type objects.

    It has attached_to dictionaries for end A and B, with entries formatted as:
    id  : { 'ref' : object ) # simply a reference to the object it's attached to
    
    
    For grouped/multilevel edges, connections are stored at the highest level
    in the nodes they are attached to. But each edge object will store what it's
    attached to in its attached_to dict, even though that's repeated info between
    the levels.
    In general, edge methods will worry about their internal subcomponents, but won't
    "look up" to see if they are part of something bigger. Except getTopLevelEdge().
    '''

    def __init__(self, id):
        
        self.id = id  # id number or string, used as the key when attached to things
        
        self.attached_to = [None, None]  # object end [A, B] of this edge is attached to
        
        # End A and B locations
        self.rA = np.zeros(2)
        self.rB = np.zeros(2)
        
        # Some attributes related to super-edges used to group things
        self.part_of = None  # whether this object is part of an Edge group
        
        self.subcomponents = []  # chain of edges and nodes that make up this edge 
        # (e.g. sections of a mooring line, and connetors between them)
        
        self.subcons_A = [] # subcomponent for end A (can be multiple)
        self.subcons_B = [] # subcomponent for end B (can be multiple)
        
        whole = True  # false if there is a disconnect among the sub edges/nodes
    
        # installation dictionary [checks for installation status]
        self.inst = {'mobilized': False,
                     'installed': False,
                     'hookedUpA': False,
                     'hookedUpB': False}

    def isAttachedTo(self, query, end=None):
        '''Checks if this edge is attached to the Node object 'query'.
        It uses the node's isAttached method since that will also check
        if this edge in question is part of a higher-level edge that 
        might be what is stored in the attachments list.
        
        Parameters
        ----------
        query
            The Node we might be attached to.
        end
            Which end of the line being asked about, 'a' or 'b'.
        '''
        
        if isinstance(query, Node):  # call the node's method; it's easy
            return query.isAttached(self, end=end)
        else:
            raise Exception('Edges can only be attached to Nodes.')
    
    
    def attachTo(self, object, r_rel=[0,0], end=None):
        '''Attach an end of this edge to some node.
        
        Parameters
        ----------
        object
            The Node object to attach to.
        r_rel : list
            x and y coordinates of the attachment point [m].
        end
            Which end of the line is being attached, 'a' or 'b'.
        '''
        
        # Determine which end to attach
        i_end = endToIndex(end, estr='when attaching an edge to something.')
        
        # Error if this edge end is already attached to the object (or could pass?)
        #if object.isAttached(self, i_end):
        #    raise Exception(f"Edge {self.id} end {end} is already attached to Node {object.id}.")
        
        # Make sure this end isn't already attached to something
        if self.attached_to[i_end]:
            self.detachFrom(end=i_end)
        
        # Tell the Node in question to initiate the attachment 
        # (this will also call self._attach_to)
        object.attach(self, r_rel=r_rel, end=i_end)
    
    
    def _attach_to(self, object, end):
        '''Internal method to update the edge's attached_to registry when
        requested by a node. This doesn't do higher level attachments.
        
        Parameters
        ----------
        object
            The Node object to attach to.
        end
            Which end of the line is being attached, a=0, b=1.
        '''
        i_end = endToIndex(end) # 0 for end A, 1 for end B
        
        if not isinstance(object, Node):
            raise Exception('Edge objects can only be attached to Node objects.')
        
        # Could potentially check if it's already attached to something <<<
        
        # Add it to the attached_to registry
        self.attached_to[i_end] = object
        
        '''
        # Recursively attach any subcomponent at the end
        >>> do we really want to do this still? <<<
        if len(self.subcomponents) > 0:  # this edge has subcomponents
            for i_sub_end in self.end_inds[i_end]: # go through each subcomponent associated with this end
                subcon = self.subcomponents[i_sub_end]
            
                if isinstance(subcon, Node): # if the end subcomponent is a node
                    subcon._attach_to(object)
                    
                elif isinstance(subcon, Edge): # if it's an edge
                    subcon._attach_to(object, i_end)  # i_end will tell it whether to use end A or B
        ''' 
        
    
    
    def detachFrom(self, end):
        '''Detach the specified end of the edge from whatever it's attached to.
        
        Parameters
        ----------
        end
            Which end of the line is being dettached, 'a' or 'b'.
        '''
        # Determine which end to detach
        i_end = endToIndex(end, estr='when detaching an edge.')
        
        # Tell the Node the end is attached to to initiate detachment 
        # (this will then also call self._detach_from)
        self.attached_to[i_end].detach(self, end=i_end)
        
        
    def _detach_from(self, end):
        '''Internal method to update the edge's attached_to registry when
        requested by a node. In nested edges, it will recursively call any 
        lower objects.
        
        Parameters
        ----------
        end
            Which end of the line is being detached, a-false, b-true.
        '''
        i_end = endToIndex(end)
        # Delete the attachment(s) of the edge end (there should only be one)
        self.attached_to[i_end] = None
        
        # Go deeper >..... >>>
        if i_end == 0:
            end_subcons = self.subcons_A
        elif i_end == 1:
            end_subcons = self.subcons_B

        for subcon in end_subcons:
            if isinstance(subcon, Node):
                if subcon.attached_to:
                    subcon.attached_to.detach(subcon)
                # the above will then call subcon._detach_from()
                
            elif isinstance(subcon, Edge):
                if subcon.attached_to[i_end]:
                    subcon.detachFrom(end=i_end)
                # the above will eventually call subcon._detach_from(i_end)

        # could add a check that if it isn't the highest left edge and if it 
        # isn't called from another edge's _detach_from method, then error, 
        # or go up a level...
    
    
    def addSubcomponents(self, items, iA=[0], iB=[-1]):
        '''If items is a list: Adds a sequences of nodes and edges 
        (alternating) as subcomponents of this edge. It also connects 
        the sequence in the process, and saves it as a dict rather than list.??        
        If items is a list: Adds them, assuming they are already assembled.??
        iA, iB : index of the end subcomponent(s) - provided as lists
        '''
        
        # Attach the sequency of nodes and edges to each other
        assemble(items)
        
        # Store them as subcomponents of this edge
        self.subcomponents = items  # dict(enumerate(items))
        for item in items:
            if isinstance(item, list):
                for branch in item:
                    for subitem in branch:
                        subitem.part_of = self
            else:    
                item.part_of = self
        # Assign ends 
        if isinstance(items[0],list):
            self.subcons_A = [it[0] for it in items[0]] # subcomponent for end A (can be multiple)
        else:
            self.subcons_A = list([items[0]]) # subcomponents for end A (can be multiple)
        if isinstance(items[-1],list):
            self.subcons_B = [it[-1] for it in items[-1]] # subcomponent for end B (can be multiple)
        else:
            self.subcons_B = list([items[-1]]) # subcomponent for end B (can be multiple)
 
        
        '''
        # Make sure the subcomponents ends are connected appropriately
        # to whatever this Edge might be attached to
         >>> this seems like it shouldn't be done anymore! <<<
        for i in [0,1]:
            if self.attached_to[end_inds[i]]:
                if isinstance(self.attached_to[end_inds[i], Node):
                    self._attach_to(self.attached_to[end_inds[i]])
                else:  # it's an edge, so also tell it which end should be attached
                    self._attach_to(self.attached_to[end_inds[i]], end_ends[i])
        
        if self.attached_to[0]:
            for i in self.end_inds[0]:
                
                if isinstance(self.attached_to[end_inds[i], Node):
                    self._attach_to(self.attached_to[end_inds[i]])
                else:  # it's an edge, so also tell it which end should be attached
                    self._attach_to(self.attached_to[end_inds[i]], end_ends[i])
            #self._attach_to(self.attached_to[0], 0)
        if self.attached_to[1]:
            self._attach_to(self.attached_to[1], 1)
        '''
    """
    def getTopLevelObject(self, end):
        '''If this edge is part of a higher-level object, and the request
        corresponds to an end of that higher-level group, return the higher object,
        otherwise return this same object. Can be recursive. 
        A similar method exists for nodes.'''
        
        i_end = endToIndex(end)
        
        if self.part_of:  # if part of something higher
            supe = self.part_of  # shorthand for the super edge
            if supe.subcomponents[-i_end] == self:  # if we're at an end of supe
                return supe.getTopLevelEdge(i_end), i_end  # return supe, and which end
        
        else:
            return self, i_end  # if not part of something bigger, just return self
    """
    
    def findEnd(self, object):
        '''Checks if object is a subcomponent of self and which end it's at.'''
        
        if not object in self.subcomponents:
            obj_in = False
            for sub in self.subcomponents:
                if isinstance(sub,list):
                    for subsub in sub:
                        if object in subsub:
                            obj_in = True
                            break
            if not obj_in:                          
                raise Exception("This object is not a subcomponent of this edge!")
        
        if any([object is con for con in self.subcons_A]):
            end = 0
        elif any([object is con for con in self.subcons_B]):
            end = 1
        else:
            end = -1  # object isn't at the end of the higher level edge so do nothing
        
        return end
    
    
    def getSubcomponent(self, index):
        '''Returns the subcomponent of the edge corresponding to the provided
        index. An index with multiple entries can be used to refer to parallel
        subcomponents.
        
        Parameters
        ----------
        index: list
            The index of the subcomponent requested to be returned. Examples:
            [2]: return the third subcomponent in the series (assuming there
            are no parallel subcomponents).
            [2,1,0]: Return the first (or only) object along the second 
            parallel string at the third serial position.  Same as [2,1].
            [1,0,2]: Return the third object along the first paralle string
            at the first serial position.
        '''
    
        if np.isscalar(index):
            index = [index]  # put into a common list format if not already
        
        if len(index) == 2:  # assume length 2 is for a parallel branch with 
            index.append(0)  # with just one edge object, so add that 0 index.
        
        # Only one index means the simple case without a parallel string here
        if len(index) == 1:
            if isinstance(self.subcomponents[index[0]], list):
                raise Exception('There is a parallel string at the requested index.')
            
            object = self.subcomponents[index[0]]
        
        # Three indices means an object along a parallel string
        elif len(index) == 3:
            if not isinstance(self.subcomponents[index[0]], list):
                raise Exception('There is not a parallel string at the requested index.')
            if len(self.subcomponents[index[0]]) < index[1]+1:
                raise Exception('The number of parallel strings is less than the requested index.')
            if len(self.subcomponents[index[0]][index[1]]) < index[2]+1:
                raise Exception('The number of objects along the parallel string is less than the requested index.')
            
            object = self.subcomponents[index[0]][index[1]][index[2]]    
        
        else:  # other options are not yet supported
            raise Exception('Index must be length 1 or 3.')
        
        return object
    
    
    def setEndPosition(self, r, end):
        '''Set the position of an end of the edge. This method should only be
        called by a Node's setPosition method if the edge end is attached to
        the node.
        
        Parameters
        ----------
        r : list
            x and y coordinates to set the end at [m].
        end
            Which end of the edge is being positioned, 'a' or 'b'.
        '''
        
        if end in ['a', 'A', 0]:
            self.rA = np.array(r)
        elif end in ['b', 'B', 1]:
            self.rB = np.array(r)
        else:
            raise Exception('End A or B must be specified with either the letter, 0/1, or False/True.')


    def delete(self):
        '''Detach the point from anything it's attached to, then delete the 
        object (if such a thing is possible?).'''
        
        self.detachFrom(0)  # detach end A
        self.detachFrom(1)  # detach end B
        
        # next step would just be to separately remove any other references
        # to this object...



class Poly():
    '''A base class for objects that run between two OR MORE Node-type objects.

    It has attached_to dictionaries for ends 0-N, with entries formatted as:
    id  : { 'ref' : object ) # simply a reference to the object it's attached to
    
    With a Poly, end is no longer A/B (0/1) but instead 0:N where N is the 
    number of exposed ends of the poly.
    
    >>>>
    For grouped/multilevel edges, connections are stored at the highest level
    in the nodes they are attached to. But each edge object will store what it's
    attached to in its attached_to dict, even though that's repeated info between
    the levels.
    In general, edge methods will worry about their internal subcomponents, but won't
    "look up" to see if they are part of something bigger. Except getTopLevelEdge(). <<<< revise
    '''

    def __init__(self, id):
        
        self.id = id  # id number or string, used as the key when attached to things
        
        self.attached_to = [None, None]  # whether either end [A, B] of this object is bound to another object
        
        # End locations
        self.r = [[0,0]]
        
        # Some attributes related to super-polies used to group things
        self.part_of = None  # whether this object is part of a Poly group
        
        self.subcomponents = []  # collection of edges and nodes that make up this Poly
       
        self.sub_end_indices = [] # subcomponent index of each attachable end of this edge
        self.sub_end_ends = [] # if the subcomponent is an edge, which end corresponds to self Edge's end
        
        whole = True  # false if there are sub edges/nodes that aren't connected
    
    
    def isAttachedTo(self, query, end=None):
        '''Checks if this poly is attached to the Node object 'query'.
        It uses the node's isAttached method since that will also check
        if this poly in question is part of a higher-level poly that 
        might be what is stored in the attachments list.
        
        Parameters
        ----------
        query
            The Node we might be attached to.
        end
            Which end of the poly being asked about, 0:N.
        '''
        
        if isinstance(query, Node):  # call the node's method; it's easy
            return query.isAttached(self, end=end)
        else:
            raise Exception('Polies can only be attached to Nodes.')
    
    
    def attachTo(self, object, r_rel=[0,0], end=None):
        '''Attach an end of this poly to some node.
        
        Parameters
        ----------
        object
            The Node object to attach to.
        r_rel : list
            x and y coordinates of the attachment point [m].
        end
            Which end of the line is being attached, 0:N.
        '''
        
        # Determine which end to attach
        if end == None:
            raise Exception("Poly end must be given...")
        else:
            i_end = int(end) # <<>> endToIndex(end, estr='when attaching an edge to something.')
        
        # Make sure this end isn't already attached to something
        if self.attached_to[i_end]:
            self.detachFrom(end=i_end)
        
        # Tell the Node in question to initiate the attachment 
        # (this will also call self._attach_to)
        object.attach(self, r_rel=r_rel, end=i_end)
    
    
    def _attach_to(self, object, end):
        '''Internal method to update the poly's attached_to registry when
        requested by a node. In nested polies, expects to be called for the
        highest-level one first, then it will recursively call any lower ones.
        
        Parameters
        ----------
        object
            The Node object to attach to.
        end
            Which end of the Poly is being attached, 0:N.
        '''
        i_end = endToIndex(end) # <<< (all these are redundant for polies?) <<<
        
        if not isinstance(object, Node):
            raise Exception('Poly objects can only be attached to Node objects.')
        
        # Add it to the attached_to registry
        self.attached_to[i_end] = object
        
        # Recursively attach any subcomponent at the end
        if len(self.subcomponents) > 0:
            subcon = self.subcomponents[-i_end]
            
            if isinstance(subcon, Node): # if the end subcomponent is a node
                subcon._attach_to(object)
                
            elif isinstance(subcon, Edge): # if it's an edge
                subcon._attach_to(object, i_end)
    
    
    def detachFrom(self, end):
        '''Detach the specified end of the poly from whatever it's attached to.
        
        Parameters
        ----------
        end
            Which end of the line is being dettached, 'a' or 'b'.
        '''
        # Determine which end to detach
        #i_end = endToIndex(end, estr='when detaching a poly.')
        
        # Tell the Node the end is attached to to initiate detachment 
        # (this will then also call self._detach_from)
        self.attached_to[i_end].detach(self, end=i_end)
        
        
    def _detach_from(self, end):
        '''Internal method to update the poly's attached_to registry when
        requested by a node. In nested polies, expects to be called for the
        highest-level one first, then it will recursively call any lower ones.
        
        Parameters
        ----------
        end
            Which end of the line is being detached, a-false, b-true.
        '''
        i_end = endToIndex(end)
        # Delete the attachment(s) of the edge end (there should only be one)
        self.attached_to[i_end] = None
        
        # Recursively detach the ends of any sub-edges
        
        #>>> need to figure out which subcomponent would correspond to the requested end of the poly <<<
        
        if len(self.subcomponents) > 0:
            subcon = self.subcomponents[-i_end]  # index 0 for A, -1 for B
            
            if isinstance(subcon, Node):
                subcon._detach_from()
                
            elif isinstance(subcon, Edge):
                subcon._detach_from(i_end)
        
        '''
        if self.sub_edges:
            if end:
                self.sub_edges[-1]._detach_from(end)
            else:
                self.sub_edges[0]._detach_from(end)
        '''
        # could add a check that if it isn't the highest left edge and if it 
        # isn't called from another edge's _detach_from method, then error, 
        # or go up a level...
    
    
    def addSubcomponents(self, items):
        '''Adds a collection of nodes and edges as subcomponents of this Poly.
        These subcomponents should already be attached to each other.'''
        
        # >>> check if the items are already assembled? <<<
        
        # Store and register them as subcomponents of this Poly
        self.subcomponents = items
        for item in items:
            item.part_of = self
        
        
        # Make sure the subcomponents ends are connected appropriately
        # to whatever this Edge might be attached to
        for i in range(len(end_inds)):
            if self.attached_to[end_inds[i]]:
                if isinstance(self.attached_to[end_inds[i]], Node):
                    self._attach_to(self.attached_to[end_inds[i]])
                else:  # it's an edge, so also tell it which end should be attached
                    self._attach_to(self.attached_to[end_inds[i]], end_ends[i])
   
    """
    def getTopLevelEdge(self, end):
        '''If this edge is part of a higher-level edge group, and the request
        corresponds to an end of that higher-level group, return the higher edge,
        otherwise return this same object. Can be recursive. 
        A similar method exists for nodes.'''
        >>>
        i_end = endToIndex(end)
        
        if self.part_of:  # if part of something higher
            supe = self.part_of  # shorthand for the super edge
            if supe.subcomponents[-i_end] == self:  # if we're at an end of supe
                return supe.getTopLevelEdge(i_end), i_end  # return supe, and which end
            
            '''
            supe = self.part_of  # shorthand for the super edge
            if end:
                if supe.sub_edges[-1] == self:  # if we're at super end B
                    return supe.getTopLevelEdge(end)
            else:  # end A
                if supe.sub_edges[0] == self:  # if we're at super end A
                    return supe.getTopLevelEdge(end)
            '''
        else:
            return self, i_end  # if not part of something bigger, just return self
    
    
    def setEndPosition(self, r, end):
        '''Set the position of an end of the edge. This method should only be
        called by a Node's setPosition method if the edge end is attached to
        the node.
        
        Parameters
        ----------
        r : list
            x and y coordinates to set the end at [m].
        end
            Which end of the edge is being positioned, 'a' or 'b'.
        '''
        >>>
        if end in ['a', 'A', 0]:
            self.rA = np.array(r)
        elif end in ['b', 'B', 1]:
            self.rB = np.array(r)
        else:
            raise Exception('End A or B must be specified with either the letter, 0/1, or False/True.')


    def delete(self):
        '''Detach the point from anything it's attached to, then delete the 
        object (if such a thing is possible?).'''
        >>>>
        self.detachFrom(0)  # detach end A
        self.detachFrom(1)  # detach end B
        
        # next step would just be to separately remove any other references
        # to this object...
    """


# general functions

def areAttached(object1, object2):
    '''Check if two objects are attached to each other.
    If both are nodes, checks if one is in the attached list of the other.
    If both are edges, checks if both are attached to the same node.
    If one is an edge and one is a node, checks if the node is attached to 
    the edge.
    '''
    
    
    
def detach(object1, object2):
    '''Detaches two objects if they are attached in any way.
    Could optionally raise an error if they aren't already attached. '''
    
    
def attach(self, object1, object2, r_rel=[0,0], end=None):#end1=None, end2=None):
    '''Attached 1 to object2, as long as the object types are compatible.
    Ends can to be specified if either object is an edge type, otherwise
    if not specified then an available end will be connected.
    '''
    
    '''Attach something to this node.
    Attachments are noted with an entry in self.attachments
    consisting of a dictionary with ref, type, and end entries.
    
    Parameters
    ----------
    object
        The Node or Edge object being attached to this one.
    r_rel : list
        x and y coordinates of the attachment point [m].
    end
        If an Edge is being attached, which end of it, 'a' or 'b'.
    '''

    
# lower-level utility functions
def endToIndex(end, estr='', n=2):
    '''Converts an end specifier (a, b, A, B, 0, 1) to just 0, 1 for use
    in the attached_to indexing of an Edge-type object.'''
    
    if type(end) == str:
        if len(end) == 1:
            end = ord(end.lower())-97  # convert letter to integer (A=0, b=1, etc)
        else:
            raise Exception("When providing 'end' as a string, it must be a single letter.")

    if not type(end) == int:
        raise Exception('End must be provided as a character or integer.')

    if end < 0:
        raise Exception('End must be positive.')
    elif end > n-1:
        if n==2:
            raise Exception('End A/B must be specified (with a/b or 0/1) '+estr)
        else:
            raise Exception(f'The specified end value exceeds the limit of {n} values from 0/A' +estr)

    return end


def assemble(items):
    '''Strings together a sequence of nodes and edges'''
    
    # >>> TODO: adjust this so it can connect parallel elements .eg. for bridles <<<
    
    '''
    # If the provided items isn't a list, there's nothing to assemble so return
    if not isinstance(items[i], List):
        print('Not a list - returning!')
        return
    '''
    n = len(items)
    for i in range(n-1):
        if isinstance(items[i], list):
            for subitem in items[i]:  # go through each parallel subitem
                if isinstance(subitem, list):  # if it's a concatenation of multiple things
                
                    assemble(subitem) # make sure that any sublist is assembled
                    
                    # attach the end objects of the subitem to the nodes before and after
                    if i > 0 and isinstance(items[i-1], Node):  # attach to previous node
                        items[i-1].attach(subitem[0], end='a')
                    if i < n-1 and isinstance(items[i+1], Node):  # attach to next node
                        items[i+1].attach(subitem[-1], end='b')
                    # note: this requires the end objects to be edges
                
                elif isinstance(subitem, Edge): # if the subitem is just one edge
                    print("THIS CASE SHOULDN'T HAPPEN - the list should be nested more")
                    breakpoint()
                    if i > 0 and isinstance(items[i-1], Node):  # attach to previous node
                        items[i-1].attach(subitem, end='a')
                    if i < n-1 and isinstance(items[i+1], Node):  # attach to next node
                        items[i+1].attach(subitem, end='b')
                else:
                    raise Exception("Unsupported situation ... parallel subitems must be edges or concatenations")
                
        elif isinstance(items[i], Node) and isinstance(items[i+1], list):
            pass  # this node connects to a bridle or doubled section, 
            # so it will be hooked up in the next step
            
        elif isinstance(items[i], Node) and isinstance(items[i+1], Edge):
            items[i].attach(items[i+1], end='a')
        
        elif isinstance(items[i], Edge) and isinstance(items[i+1], Node):
            items[i+1].attach(items[i], end='b')
        
        else:
            raise Exception('sequences is not alternating between nodes and edges')
    # check if last item in items is a list (if length of items>1)
    # if it is a list, it won't have been attached/assembled previously, so 
    # attach and assemble now
    if n-1>0:        
        if isinstance(items[i+1], list):
            for subitem in items[i+1]:  # go through each parallel subitem
                if isinstance(subitem, list):  # if it's a concatenation of multiple things
                    assemble(subitem) # make sure that any sublist is assembled
                    
                    # attach the end objects of the subitem to the nodes before and after
                    if i > 0 and isinstance(items[i], Node):  # attach to previous node
                        items[i].attach(subitem[0], end='a')
                    if i < n-1 and isinstance(items[i+1], Node):  # attach to next node
                        items[i+1].attach(subitem[-1], end='b')
                    # note: this requires the end objects to be edges
                
                elif isinstance(subitem, Edge): # if the subitem is just one edge
                    print("THIS CASE SHOULDN'T HAPPEN - the list should be nested more")
                    breakpoint()
                    if i > 0 and isinstance(items[i], Node):  # attach to previous node
                        items[i].attach(subitem, end='a')
                    if i < n-1 and isinstance(items[i+1], Node):  # attach to next node
                        items[i+1].attach(subitem, end='b')
                else:
                    raise Exception("Unsupported situation ... parallel subitems must be edges or concatenations")

def rotationMatrix(x3,x2,x1):
    '''Calculates a rotation matrix based on order-z,y,x instrinsic (tait-bryan?) angles, meaning
    they are about the ROTATED axes. (rotation about z-axis would be (0,0,theta) )
    (Copied from MoorPy)
    
    Parameters
    ----------
    x3, x2, x1: floats
        The angles that the rotated axes are from the nonrotated axes. Normally roll,pitch,yaw respectively. [rad]

    Returns
    -------
    R : matrix
        The rotation matrix
    '''
    # initialize the sines and cosines
    s1 = np.sin(x1) 
    c1 = np.cos(x1)
    s2 = np.sin(x2) 
    c2 = np.cos(x2)
    s3 = np.sin(x3) 
    c3 = np.cos(x3)
    
    # create the rotation matrix
    R = np.array([[ c1*c2,  c1*s2*s3-c3*s1,  s1*s3+c1*c3*s2],
                  [ c2*s1,  c1*c3+s1*s2*s3,  c3*s1*s2-c1*s3],
                  [   -s2,           c2*s3,           c2*c3]])
    
    return R    
    

# test script
if __name__ == '__main__':
    
    
    # ----- base test of an edge with two nodes -----
    node1 = Node(id='node1')
    node2 = Node(id='node2')
    edge1 = Edge(id='edge1')
    
    node1.attach(edge1, r_rel=[10, 0], end='A')
    node2.attach(edge1, r_rel=[0, 10], end='B')
    
    node1.setPosition([0, 0])
    node2.setPosition([50, 50])
    
    
    # attaching the other end to the same object should cause an error:
    #edge1.attachTo(node2, end='a')
    
    # but this should work
    edge1.detachFrom(end='b')
    edge1.attachTo(node2, end='a')
    
    
    # ----- a test for super edges... -----
    
    e0 = Edge(id='e0')
    e1 = Edge(id='e1')
    e2 = Edge(id='e2')
    
    n0 = Node(id='n0')
    n1 = Node(id='n1')
    
    #assemble([e0,n0,e1,n1,e2])
    
    E = Edge(id='big edge')
    A = Node(id='Node A')
    B = Node(id='Node B')
    
    assemble([A, E, B])
    
    E.addSubcomponents([e0,n0,e1,n1,e2])    
    
    
    # ----- a test for bridles etc -----
    
    e0 = Edge(id='e0')
    n0 = Node(id='n0')
    e1 = Edge(id='e1')
    n1 = Node(id='n1')
    e2a1 = Edge(id='e2a')
    e2a2 = Node(id='e2a')
    e2a3 = Edge(id='e2a')
    e2b  = Edge(id='e2b')
    
    
    thelist = [e0, n0, e1, n1, [[e2a1, e2a2, e2a3], [e2b]]]
    
    E = Edge(id='big edge')
    
    E.addSubcomponents(thelist)
    
    s = E.getSubcomponent([4,0,2])
    
    # ----- try joining two nodes -----
    """
    A = Node(id='Node A')
    B = Node(id='Node B')
    A.join(B)
    
    
    # ----- tests connecting multi-level node and edge objects ----
    
    # --- Test 1 ---
    # platform and fairlead
    n1 = Node(id='n1')
    n2 = Node(id='n2')
    n1.attach(n2, r_rel=[20,0,-10])
    # mooring and contents
    e1 = Edge(id='e1')
    e1_e1 = Edge(id='e1_e1')
    e1_n2 = Node(id='e1_n2')
    e1_e3 = Edge(id='e1_e3')
    e1.addSubcomponents([e1_e1, e1_n2, e1_e3])
    # attach mooring to platfrom (by lower objects, then upper will be automatic)
    #n2.attach(e1_e1, end='A')
    n2.attach(e1_e3, end='B')
    
    # --- Test 2 ---
    # platform and fairlead
    n1 = Node(id='n1')
    n2 = Node(id='n2')
    n1.attach(n2, r_rel=[20,0,-10])
    # mooring and contents
    e1 = Edge(id='e1')
    e1_n1 = Node(id='e1_n1')
    e1_e2 = Edge(id='e1_e2')
    e1_n3 = Node(id='e1_n3')
    e1.addSubcomponents([e1_n1, e1_e2, e1_n3])
    # attach mooring to platfrom (by lower objects, then upper will be automatic)
    n2.attach(e1_n1)
    #n2.attach(e1_n3)
    #n2.join(e1_n1)
    #n2.join(e1_n3)
    
    # --- Test 3 ---
    # platform and fairlead
    n1 = Node(id='n1')
    # mooring and contents
    e1 = Edge(id='e1')
    e1_e1 = Edge(id='e1_e1')
    e1_n2 = Node(id='e1_n2')
    e1_e3 = Edge(id='e1_e3')
    e1.addSubcomponents([e1_e1, e1_n2, e1_e3])
    # attach mooring to platfrom (by lower objects, then upper will be automatic)
    #n1.attach(e1_e1, r_rel=[20,0,-10], end='A')
    n1.attach(e1_e3, r_rel=[20,0,-10], end='B')
    
    # --- Test 4 ---
    # platform and fairlead
    n1 = Node(id='n1')
    # mooring and contents
    e1 = Edge(id='e1')
    e1_n1 = Node(id='e1_n1')
    e1_e2 = Edge(id='e1_e2')
    e1_n3 = Node(id='e1_n3')
    e1.addSubcomponents([e1_n1, e1_e2, e1_n3])
    # attach mooring to platfrom (by lower objects, then upper will be automatic)
    #n1.attach(e1_n1, r_rel=[20,0,-10])
    n1.attach(e1_n3, r_rel=[20,0,-10])
    
    # --- Test 5 ---
    # platform and fairlead
    n1 = Node(id='n1')
    n2 = Node(id='n2')
    n1.attach(n2, r_rel=[20,0,-10])
    # mooring and contents
    e1_e1 = Edge(id='e1_e1')
    e1_n2 = Node(id='e1_n2')
    e1_n2.attach(e1_e1, end='B')
    # attach mooring to platfrom (by lower objects, then upper will be automatic)
    n2.attach(e1_n2)
    
    # --- Test 6 ---
    # platform and fairlead
    n1 = Node(id='n1')
    n2 = Node(id='n2')
    n1.attach(n2, r_rel=[20,0,-10])
    # mooring and contents
    e1_n1 = Node(id='e1_n1')
    e1_e2 = Edge(id='e1_e2')
    e1_n1.attach(e1_e2, end='A')
    # attach mooring to platfrom (by lower objects, then upper will be automatic)
    n2.attach(e1_e2, end='B')
    # --- done tests ---
    n1.setPosition(r=[0,0,0], theta=0)
    #print(n1.attachments)
    #print(e1.attached_to)
    """


         
    
    