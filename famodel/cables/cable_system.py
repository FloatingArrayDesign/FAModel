# class for the collection of power cables of an array

class CableSystem():
    def __init__(self):
        self.cableList = []
        
        
    # for calculations that involve both cables and seabed, 
    # should they be done in the cable classes (referencing/storing
    # the higher-level site/seabed data), or should it be done at
    # a higher level?
    
    
    
    
    def addCable(self, cableType, nSegs=40, pointA=0, pointB=0):
        '''Convenience function to add a Cable to a cable system.

        Parameters
        ----------
        lineType : string or dict
            string identifier of lineType for this line already added to the system, or dict specifying custom line type.
        nSegs : int, optional
            number of segments to split the line into. The default is 20.
        pointA int, optional
            Point number to attach end A of the line to.
        pointB int, optional
            Point number to attach end B of the line to.

        Returns
        -------
        None.
        '''
        
        if not isinstance(lineType, dict):                      # If lineType is not a dict, presumably it is a key for System.LineTypes.
            if lineType in self.lineTypes:                      # So make sure it matches up with a System.LineType
                lineType = self.lineTypes[lineType]             # in which case that entry will get passed to Line.init
            else:
                raise ValueError(f"The specified lineType name ({lineType}) does not correspond with any lineType stored in this MoorPy System")
        
        self.lineList.append( Line(self, len(self.lineList)+1, lUnstr, lineType, nSegs=nSegs, cb=cb) )
        
        if pointA > 0:
            if pointA <= len(self.pointList):
                self.pointList[pointA-1].attachLine(self.lineList[-1].number, 0)
            else:
                raise Exception(f"Provided pointA of {pointA} exceeds number of points.")
        if pointB > 0:
            if pointB <= len(self.pointList):
                self.pointList[pointB-1].attachLine(self.lineList[-1].number, 1)
            else:
                raise Exception(f"Provided pointB of {pointB} exceeds number of points.")
        