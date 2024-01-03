

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
'''

class Node():
	'''A base class for objects that have a single spatial reference point and
    can be interconnected by edge-type objects.
    
    Its attachments dictionary has the following format for each entry:
    
    id  : { 'ref'   :  ,  #
            'r_rel' :  ,  #
            'type'  :     # Node or Edge
    
    '''
	
	def __init__(self, id):
		
        self.id = id  # id number or string, used as the key when attached to things
        
        self.attachments = {}  # dictionary listing attached edges (key is edge ID)
        
        self.attached_to = {}  # whether this object is bound to another object
    
        # position/orientation variables
        self.r = np.zeros(2)  # position [m]
        self.theta = 0  # heading [rad]
        self.R = np.eye(2)  # rotation matrix
    
    
    def attach(self, object, r_rel=[0,0], end=None):
        '''Attach something to this node.
        
        Parameters
        ----------
        object
            The Node or Edge object being attached to this one.
        r_rel : list
            x and y coordinates of the attachment point [m].
        end
            If an Edge is being attached, which end of it, 'a' or 'b'.
        '''
        
        # Make sure it's not already attached
        if object['id'] in self.attachments:
            raise Exception(f"Object {id} is already attached to {self['id']})
        
        # Prepare its entry in the attachment registry
        new_entry = dict(ref=object, r_rel=np.array(r_rel))
        
        # Handle attachment type and the end if applicable, and record
        # the attachment in the subordinate object.
        if type(object) == Node:
            new_entry['type'] = 'node'
            object._attach_to(self)  # tell the attached object to record things
            
        elif type(object) == Edge:
            new_entry['type'] = 'edge'
            if not end:
                raise Exception("The 'end' argument must be provided when attaching an edge.")
            if end.lower() in ['a', 'b']:
                new_entry['end'] = end.lower()
                object._attach_to(self, new_entry['end']=='b')  # tell the attached object to record things
            else:
                raise Exception('End A or B must be specified when attaching an edge.')
        
        else:
            raise Exception('Provided object is not an Edge or Node.')
        
        # Since things have worked, add it to the list
        self.attachments[object.id] = new_entry
    
    
    def detach(self, object, end=None):
        '''Detach the specified object from this node.
        
        Parameters
        ----------
        object
            The object to detach (its ID should be in the node's attachments dict).
        end
            If an Edge is being detached, which end of it, 'a' or 'b'.
        '''
        
        # Make sure it's attached before trying to detach it
        if not object.id in self.attachments:
            raise Exception(f"Object {id} is not attached to {self['id']})
            # this exception could be optionally disabled
        
        # Remove it from the attachment registry
        del self.attachments[object.id]
        
        # Handle attachment type and the end if applicable, and record
        # the detachment in the subordinate object.
        if type(object) == Node:
            object._detach_from(self)  # tell the attached object to record things
            
        elif type(object) == Edge:
            if not end:
                raise Exception("The 'end' argument must be provided when detaching an edge.")
            if end.lower() in ['a', 'b']:
                object._detach_from(self, end.lower()=='b')  # tell the detached object to record things
            else:
                raise Exception('End A or B must be specified when detaching an edge.')
        
        else:
            raise Exception('Provided object is not an Edge or Node.')
        
        

    def _attach_to(self, object):
        '''Internal method to update the Node's attached_to registry when
        requested by a node.
        
        Parameters
        ----------
        object
            The Node object to attach to.
        '''
        
        if type(object) == Node:
            raise Exception('Node objects can only be attached subordinately to Node objects.')
        
        # Make sure it's not already attached to something else
        if self.attached_to:  # True if populated, False if empty dict
            self.detach()
            # could optionally have a warning or error here
                
        # Add it to the attached_to registry
        self.attached_to[object['id']] = dict(ref=object)
        
        
    def _detach_from(self):
        '''Internal method to update the Node's attached_to registry when
        requested by a node to detach.
        '''
        self.attached_to = {}


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
        for id, att in self.attachments.items():
        
            # Compute the attachment's position
            r_att = self.r + np.matmul(self.R, att['r_rel'])
            
            if att['type'] == 'node':
                att['ref'].setPosition(r_att)
                
            elif att['type'] == 'edge':
                att['ref'].setEndPosition(r_att, att['end'])



class Edge():
	'''A base class for objects that run between two Node-type objects.
    
    It has attached_to dictionaries for end A and B, with entries formatted as:
    id  : { 'ref' : object ) # simply a reference to the object it's attached to
    '''
	
	def __init__(self, id):
		
        self.id = id  # id number or string, used as the key when attached to things
        
        self.attached_to = [{},{}]  # whether either end [A, B] of this object is bound to another object
        # I need to decide whether these are dicts of multiple attachments, or just one attachment for each end.
        # For now let's assume they're dicts, but with only one entry.
    
    
    def isAttachedTo(self, query, end=None):
        '''Checks if the object 'query' is attached to this object.
        query can be a reference to the object itself, or the object ID.
        '''
        
        # If query is a reference to an object, get its ID
        if type(query) in [Node, Edge]:
            id = query.id
        else:
            id = query
        
        if end:
            if end in ['a', 'A', 0]:
                answer = bool(id in self.attached_to[0])
            elif end in ['b', 'B', 1]:
                answer = bool(id in self.attached_to[1])
            else:
                raise Exception("The 'end' paramter must be one of a,b,A,B,0,1,False,True.")
        else:
            answer = bool(id in self.attached_to[0] or id in self.attached_to[1])
        
        return answer
        
    
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
        if not end
            raise Exception('end must be specified (a or b).')
        if end.lower() == 'a':
            i_end = 0
        if end.lower() == 'b':
            i_end = 1
        else:
            raise Exception('End A or B must be specified when attaching an edge.')
        
        # Make sure it's not already attached
        if self.attached_to[i_end]:
            self.detach()
            
        if object['id'] in self.attachments:
            raise Exception(f"Object {id} is already attached to {self['id']})
        
        # Tell the Node in question to initiate the attachment 
        # (this will also call self._attach_to)
        object.attach(self, r_rel, end)
        
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
    
    
    def _attach_to(self, object, endB):
        '''Internal method to update the edge's attached_to registry when
        requested by a node.
        
        Parameters
        ----------
        object
            The Node object to attach to.
        endB
            Which end of the line is being attached, a-false, b-true.
        '''
        
        i_end = int(endB) # 0 for end A, 1 for end B
        
        if type(object) == Node:
            raise Exception('Edge objects can only be attached to Node objects.')
        
        # Make sure it's not already attached to something else
        if self.attached_to[i_end]:  # True if populated, False if empty dict
            self.detach()
            # could optionally have a warning or error here
                
        # Add it to the attached_to registry
        self.attached_to[i_end][object['id']] = dict(ref=object)
        
        
    def detachFrom(self, end=None):
        '''Detach the specified end of the edge from whatever it's attached to.
        
        Parameters
        ----------
        end
            Which end of the line is being dettached, 'a' or 'b'.
        '''
        
        # Determine which end to detach
        if not end
            raise Exception('end must be specified (a or b).')
        if end in ['a', 'A', 0]:
            i_end = 0
        elif end in ['b', 'B', 1]:
            i_end = 1
        else:
            raise Exception('End A or B must be specified when detaching an edge.')
        
        # Loop through the attachments of the selected end (there should only be one)
        # and tell the Node to initiate detachment (this will also call self._detach_from)
        for k, v in self.attached_to[i_end].items():
            v.detach(self)
        
        
    def _detach_from(self, endB):
        '''Internal method to update the edge's attached_to registry when
        requested by a node.
        
        Parameters
        ----------
        endB
            Which end of the line is being detached, a-false, b-true.
        '''
        # Delete the attachment(s) of the edge end (there should only be one)
        self.attached_to[i_end] = {}
    
    
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




# general functions

def detach(object1, object2):
    '''Detaches two objects if they are attached in any way.'''
    
    
def attach(object1, object2, r_rel=[0,0], end1=None, end2=None):
    '''Attached 1 to object2, as long as the object types are compatible.
    Ends can to be specified if either object is an end type, otherwise
    if not specified then an available end will be connected.
    '''
    
# lower-level utility functions
def endToIndex(end)
    '''Converts an end specifier (a, b, A, B, 0, 1) to just 0, 1 for use
    in the attached_to indexing of an Edge-type object.'''
            
    if end in ['a', 'A', 0]:
        return = 0
    elif end in ['b', 'B', 1]:
        return = 1
    else:
        raise Exception('End A or B must be specified with either the letter, 0/1, or False/True.')
