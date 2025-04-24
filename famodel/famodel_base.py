import numpy as np

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
relationships are tracked with sub_edge and sub_node lists in the 
higher-level Edge object, and part_of entries in the lower-level node or edge
objects. 

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
        
        self.attachments = {}  # dictionary listing attached edges 
        # (key is the id, value is attachment object ref and other info)
        
        self.attached_to = None  # whether this object is bound to another object
        
        self.part_of = None  # whether this object is part of an Edge group
    
        # position/orientation variables
        self.r = np.zeros(2)  # position [m]
        self.theta = 0  # heading [rad]
        self.R = np.eye(2)  # rotation matrix
    
        # installation dictionary [checks for installation status]
        self.inst = {'mobilized': False,
                     'installed': False}

        # raft results dictionary
        self.raftResults = {}
        
    def isAttached(self, object, end=None):
        '''Check if something is attached to this node, even if it's part of
        a higher-level edge.
        
        Parameters
        ----------
        object
            The Node or Edge object being attached to this one.
        end
            If an Edge is being attached, which end of it, 'a' or 'b'.
        '''
        
        # find top-level edge end if applicable
        if isinstance(object, Node):
            object2, i_end = object.getTopLevelEdge()
        elif isinstance(object, Edge):
            i_end = endToIndex(end, estr='when checking if an edge is attached to a node.')
            object2, _ = object.getTopLevelEdge(i_end)
        else:
            raise Exception('Provided object is not an Edge or Node.')  
        
        # See if it's attached (might be a higher-level edge)
        # if object2.id in self.attachments:
        if object2.id in self.attachments:
        
            if isinstance(object2, Node):  # if it's a node, it's simple
                return True
                
            elif isinstance(object2, Edge):  # if it's an edge, end matters
                if endToIndex(self.attachments[object2.id]['end']) == i_end:
                    return True  # the end in question is attached
                else:
                    return False
        else:
            return False
        
    
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
        # find top-level edge end if applicable
        if isinstance(object, Node):
            object2, i_end = object.getTopLevelEdge()
        elif isinstance(object, Edge):
            i_end = endToIndex(end, estr='when attaching an edge to a node.')
            object2, _ = object.getTopLevelEdge(i_end)  
        else:
            raise Exception('Provided object is not an Edge or Node.')  
        
        
        # Make sure it's not already attached (note this doesn't distinguish end A/B)
        if object2.id in self.attachments:
            raise Exception(f"Object {object.id} is already attached to {self.id}")
            
            
        # Attach the object (might be a higher-level edge)
        if isinstance(object2, Node):
            self.attachments[object2.id] = dict(obj=object2, id=object2.id, 
                                            r_rel=np.array(r_rel), type='node')
            #object2._attach_to(self)  # tell it it's attached to this Node
            object2.attachments[self.id] = dict(obj=self,id=self.id,r_rel=np.array(r_rel),type='node')
            
        elif isinstance(object2, Edge):
            self.attachments[object2.id] = dict(obj=object2, id=object2.id, 
                                            r_rel=np.array(r_rel), type='edge', 
                                            end=['a', 'b'][i_end])
            object2._attach_to(self, i_end)  # tell it it's attached to this Node
            '''
            if end in ['a', 'A', 0]:
                new_entry['end'] = 'a'
                object._attach_to(self, 0)
            
            elif end in ['b', 'B', 1]:
                new_entry['end'] = 'b'
                object._attach_to(self, 1)
            
            else:
                raise Exception('End A or B must be specified when attaching an edge.')
            '''
        else:
            raise Exception('Unrecognized object type')
    
    
    def detach(self, object, end=None):
        '''Detach the specified object from this node.
        Note that this method doesn't search for highest-level edge 
        attachment because it should already be attached that way.
        
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
        
        # Remove it from the attachment registry
        del self.attachments[object.id]
        
        # Handle attachment type and the end if applicable, and record
        # the detachment in the subordinate object.
        if isinstance(object, Node):
            object._detach_from()  # tell the attached object to record things
            
        elif isinstance(object, Edge):
            
            i_end = endToIndex(end, estr='when detaching an edge from a node.')
            
            object._detach_from(i_end)
            '''
            if end in ['a', 'A', 0]:
                object._detach_from(0)
            
            elif end in ['b', 'B', 1]:
                object._detach_from(1)
            
            else:
                raise Exception('End A or B must be specified when detaching an edge.')
            '''
        else:
            raise Exception('Provided object is not an Edge or Node.')
        
        

    def _attach_to(self, object,sub=0):
        '''Internal method to update the Node's attached_to registry when
        requested by a node.
        
        Parameters
        ----------
        object
            The Node object to attach to.
        sub
            Boolean to mark if this node is subordinately attached to the object
        '''
        if not sub:
            if isinstance(object, Node):
                raise Exception('Node objects can only be attached subordinately to Node objects.')
        # Make sure it's not already attached to something else
        if self.attached_to:  # True if populated, False if empty dict
            # self.detach() # commented out because I don't think this would work for shared anchors
            # could optionally have a warning or error here
            print(f'Warning: node {self.id} is attached to 2 objects')
        
                
        # Add it to the attached_to registry
        self.attached_to = object

        
    def _detach_from(self):
        '''Internal method to update the Node's attached_to registry when
        requested by a node to detach.
        '''
        self.attached_to = None
        


    def getTopLevelEdge(self):
        '''If this node is part of a higher-level edge group, and the request
        corresponds to an end of that higher-level group, return the higher edge,
        otherwise return this same object. Can be recursive. 
        A similar method exists for edges.'''
        
        if self.part_of:  # if part of something higher
            supe = self.part_of  # shorthand for the super edge
            if supe.subcomponents[0] == self:  # if this node is at end A of supe
                return supe.getTopLevelEdge(0)  # return supe, and which end
            elif supe.subcomponents[-1] == self:  # if this node is at end B of supe
                return supe.getTopLevelEdge(1)
        else:
            return self, -1  # if not part of something bigger, just return self
    
    
    def setPosition(self, r, theta=0):
        '''Set the position of the node, as well as any attached objects.
        
        Parameters
        ----------
        r : list
            x and y coordinates to position the node at [m].
        theta, float (optional)
            The heading of the object [rad].
        '''
        
        # Store updated position and orientation
        self.r = np.array(r)
        self.theta = theta
        
        # Get rotation matrix...
        self.R = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
        
        
        # Update the position of any attach objects
        for att in self.attachments.values():
        
            # Compute the attachment's position
            r_att = self.r + np.matmul(self.R, att['r_rel'])
            
            if isinstance(att['obj'], Node):
                att['obj'].setPosition(r_att)
                
            elif isinstance(att['obj'], Edge):
                att['obj'].setEndPosition(r_att, att['end'])



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
        
        self.attached_to = [None, None]  # whether either end [A, B] of this object is bound to another object
        
        # End A and B locations
        self.rA = [0,0]
        self.rB = [0,0]
        
        # Some attributes related to super-edges used to group things
        self.part_of = None  # whether this object is part of an Edge group
        
        self.subcomponents = []  # chain of edges and nodes that make up this edge 
        # (e.g. sections of a mooring line, and connetors between them)
       
        whole = True  # false if there are sub edges/nodes that aren't connected
    
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
        
        '''
        # old version where query can be a reference to the object itself, or the object ID.
        # If query is a reference to an object, get its ID
        if isinstance(query, Node) or isinstance(query, Edge):
            id = query.id
        else:
            id = query
        
        if end == None:
             answer = bool(id == self.attached_to[0]['id'] or id == self.attached_to[1]['id'])
        else:
            if end in ['a', 'A', 0]:
                answer = bool(id == self.attached_to[0]['id'])
            elif end in ['b', 'B', 1]:
                answer = bool(id == self.attached_to[1]['id'])
            else:
                raise Exception("If an 'end' parameter is provided, it must be one of a,b,A,B,0,1,False,True.")
        '''
        # return answer 
        
    
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
        
        '''
        # Add it to the attachment registry
        new_entry = dict(ref=object, r_rel=np.array(r_rel))
        
        # Handle attachment type and the end if applicable
        if type(object) == Node:
            new_entry['type'] = 'Node'
            
        elif type(object) == Edge:
            new_entry['type'] = 'Edge'
            if not end:
                raise Exception('End A or B must be specified when attaching an edge.')
            if end.lower() in ['a', 'b']:
                new_entry['end'] = end.lower()
            else:
                raise Exception('End A or B must be specified when attaching an edge.')
        else:
            raise Exception('Provided object is not an Edge or Node.')
        
        # Since things have worked, add it to the list
        self.attachments[object['id']] = new_entry
        '''
    
    
    def _attach_to(self, object, end):
        '''Internal method to update the edge's attached_to registry when
        requested by a node. In nested edges, expects to be called for the
        highest-level one first, then it will recursively call any lower ones.
        
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
        
        # Add it to the attached_to registry
        self.attached_to[i_end] = object
        
        # Recursively attach any subcomponent at the end
        if len(self.subcomponents) > 0:
            subcon = self.subcomponents[-i_end]  # index 0 for A, -1 for B
            
            if isinstance(subcon, Node):
                # will be attaching a node to a node, tell _attach_to it is attaching the subcon as a subordinate node we don't throw an error
                #print(f'attaching {object.id} to {subcon.id}')
                subcon._attach_to(object,sub=1)
                # object._attach_to(subcon,sub=1)
                
            elif isinstance(subcon, Edge):
                #print(f'attaching {object.id} to {subcon.id}')
                subcon._attach_to(object, i_end)
                # object._attach_to(subcon)
            
        '''
        if end:
            self.sub_edges[-1]._attach_to(object, end)
        else:
            self.sub_edges[0]._attach_to(object, end)
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
        requested by a node. In nested edges, expects to be called for the
        highest-level on first, then it will recursively call any lower ones.
        
        Parameters
        ----------
        end
            Which end of the line is being detached, a-false, b-true.
        '''
        i_end = endToIndex(end)
        # Delete the attachment(s) of the edge end (there should only be one)
        self.attached_to[i_end] = None
        
        # Recursively detach the ends of any sub-edges
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
        '''Adds a sequences of nodes and edges (alternating) as subcomponents
        of this edge. It also connects the sequence in the process.'''
        
        
        assemble(items)
        self.subcomponents = items
        for item in items:
            item.part_of = self
        
        # Make sure the subcomponents ends are connected appropriately
        # to whatever this Edge might be attached to
        if self.attached_to[0]:
            self._attach_to(self.attached_to[0], 0)
        if self.attached_to[1]:
            self._attach_to(self.attached_to[1], 1)
        
    
    def getTopLevelEdge(self, end):
        '''If this edge is part of a higher-level edge group, and the request
        corresponds to an end of that higher-level group, return the higher edge,
        otherwise return this same object. Can be recursive. 
        A similar method exists for nodes.'''
        
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


# general functions

def are_attached(object1, object2):
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
    Ends can to be specified if either object is an end type, otherwise
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
def endToIndex(end, estr=''):
    '''Converts an end specifier (a, b, A, B, 0, 1) to just 0, 1 for use
    in the attached_to indexing of an Edge-type object.'''
            
    if end in ['a', 'A', 0]:
        return 0
    elif end in ['b', 'B', 1]:
        return 1
    else:
        raise Exception('End A/B must be specified (with a/b or 0/1) '+estr)


def assemble(items):
    '''Strings together a sequence of nodes and edges'''
    
    n = len(items)
    
    for i in range(n-1):
        if isinstance(items[i], Node) and isinstance(items[i+1], Edge):
            items[i].attach(items[i+1], end='a')
        
        elif isinstance(items[i], Edge) and isinstance(items[i+1], Node):
            items[i+1].attach(items[i], end='b')
        
        else:
            raise Exception('sequences is not alternating between nodes and edges')
    
    
    
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
    
    
    # ----- make a test for super edges... -----
    
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
    

