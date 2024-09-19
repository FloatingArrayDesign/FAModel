# Platform Class

The platform class contains properties and methods relating to a moored floating platform. The platform class inherits from the Node class, and stores information on connected Edge objects (mooring lines, cables).

Currently, RAFT information on the specific geometry and hydrostatic properties of the platform are not stored in the platform object. When a project is initialized, a RAFT model is created if directed to, and this information is stored in the RAFT model found in project.array.

## Platform Properties
 - **dd** : design dictionary of the platform, including 
 - **r** : [x,y] coordinates of the platform
 - **phi** : heading offset of the platform
 - **rFair** : fairlead radius [m]
 - **zFair** : fairlead depth [m]
 - **body** : moorpy body object associated with this platform
 - **mooring_headings** : headings of associated mooring lines
 - **n_mooring** : number of mooring lines
 - **endB** : dictionary with key as mooring object names and values of booleans for if the end B of that mooring object is attached to the platform
 - **rc** : [row,column] informating location in a uniform grid array
 - **envelopes** : dictionary of 2D motion envelopes, buffers, etc. Each entry is a dict with x,y or shape
- **loads** : dictionary of loads on platform
- **reliability** : dictionary of platform reliability information
- **cost** : dictionary of platform costs
- **failure_probability** : dictionary of platform failure probabilities


## Platform Methods
- **setPosition()** : Sets the position/orientation of the platform as well as the associated anchor points.
- **mooringSystem()** : Creates a moorpy system for the platform based on the mooring objects attached. This is different than the project moorpy system, which includes all platforms. This system only includes this platform and its associated moorings and anchors.
- **getWatchCircle()** : Computes watch circle of platform, as well as mooring and cable tension safety factors and cable sag safety factors based on rated thrust. The watch circle is stored in the envelopes dictionary, and 
- **getMoorings()** : Creates list of mooring objects connected to this platform
- **getCables()** : Creates list of cable objects connected to this platform
- **getAnchors()** : Creates list of anchor objects associated with this platform
- **getBufferZones()** : Calculates buffer zones around mooring lines and anchors, and saves buffer zones in envelopes dictionary.
