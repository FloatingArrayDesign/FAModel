# class for a mooring line in a floating array



class Mooring():
    '''
    Class for a floating array mooring line (anchored or shared).
    The idea is to have a common object that can be used for 2D
    layout optimization as well as (optionally) 3D design.
    Work in progress. Eventually will inherit from Edge.
    '''
    
    def __init__(self, subsystem=None, rA=[0,0,0], rB=[0,0,0]):
        '''
        Initialize an empty object for a mooring line.
        Eventually this will fully set one up from ontology inputs.
        '''
        
        self.subsystem = subsystem  # The MoorPy subsystem that corresponds to the mooring line
        
        # end points, to be set later
        self.rA = rA
        self.rB = rB
        
    
    def setEndPosition(self, r, end):
        '''Set the position of an end of the mooring.
        
        Parameters
        ----------
        r : list
            Cordinates to set the end at [m].
        end
            Which end of the edge is being positioned, 'a' or 'b'.
        '''
        
        if end in ['a', 'A', 0]:
            self.rA = np.array(r)
            
            if self.subsystem:
                self.subsystem.setEndPosition(self.rA, False)
            
        elif end in ['b', 'B', 1]:
            self.rB = np.array(r)
            
            if self.subsystem:
                self.subsystem.setEndPosition(self.rB, True)
                
        else:
            raise Exception('End A or B must be specified with either the letter, 0/1, or False/True.')
      
    
"""
class Platform():
    '''
    Class for a mooring floating platform.
    Eventually will inherit from Node.
    '''
    
    def __init__(self, r=[0,0], heading=0, mooring_headings=[60,180,300]):
        '''
        
        Parameters
        ----------
        r 
            x and y coordinates [m].
        theta, float (optional)
            The heading of the object [deg].
        mooring_headings (optional)
            relative headings of mooring lines [deg].
        '''
        
        self.r = np.array(r)
        
        self.theta = np.radians(heading)
        
        self.mooring_headings = np.radians(mooring_headings)
        
    
    def setPosition(self, r, heading=0):
        '''
        Set the position of the node, as well as any attached objects.
        
        Parameters
        ----------
        r : list
            x and y coordinates to position the node at [m].
        heading, float (optional)
            The heading of the object [deg].
        '''
        
        # Store updated position and orientation
        self.r = np.array(r)
        self.theta = np.radians(heading)
        
        
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
"""
 