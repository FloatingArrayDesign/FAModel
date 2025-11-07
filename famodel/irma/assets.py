"Classes for vessels and port"
# Adapted from Rudy's original work and Ryan's updates

import yaml 
from copy import deepcopy

class Asset:
    '''
    Base class for vessel or port used in installation or maintenance.
    
    Attributes
    ----------
    name : str
        Name of the vessel.
    specs : dict
        Specifications of the vessel.
    state : dict
        Current state of the vessel including position, cargo, etc.
    '''
    
    def __init__(self, file = None, vesselDisc=None):
        """
        Initialize a Vessel object from a configuration file.
        
        Parameters
        ----------
        config_file : str
            Path to the vessel configuration file.
            
        Returns
        -------
        None
        """

        if vesselDisc is None and file is not None:
            with open(file) as f:
                vesselDisc = yaml.load(f, Loader=yaml.FullLoader)
        elif vesselDisc is not None and file is None:
            pass
        else:
            raise ValueError("Either vesselDisc or file must be provided.")
        
        # Set up general attributes from inputs
        self.name  = vesselDisc.get('name', "Unnamed Vessel")
        self.type  = vesselDisc.get('type', "Untyped Vessel")
        self.specs = vesselDisc['specs']
        self.state = {
            "spool_storage": self.specs['storage_specs']['max_spool_capacity'],
            "deck_storage": self.specs['storage_specs']['max_deck_space'],
            "cargo_mass": self.specs['storage_specs']['max_cargo'],
            "assigned_materials" : [],
            "log": []
        }
        
        # additional initialization should be done in the vessel or port subclass init
    
    
    def logState(self, time, new_state):
        """
        Log and update the vessel state.

        Parameters
        ----------
        time : float
            Current simulation time.
        new_state : dict
            New state information to update and log.

        Returns
        -------
        None
        """
        self.state.update(new_state)
        self.state["log"].append({"time": time, "state": new_state})


    def getState(self, t):
        """
        Retrieve vessel state at a specific time.

        Parameters
        ----------
        t : float
            Time at which to retrieve the vessel state.

        Returns
        -------
        state : dict
            The vessel state at time t, or None if no state exists before time t.
        """
        return next((log for log in reversed(self.state["log"]) if log["time"] <= t), None)
        
     




class Vessel(Asset):
    """
    Represents a vessel used in the installation process.
    
    Attributes
    ----------
    name : str
        Name of the vessel.
    specs : dict
        Specifications of the vessel.
    state : dict
        Current state of the vessel including position, cargo, etc.
    """

    def __init__(self, info):
        '''
        Initialize a Vessel object from a configuration file or dict.
        
        Parameters
        ----------
        config_file : str
            Path to the vessel configuration file.
        '''
        
        # Initialize the base class
        Asset.__init__(self, info)
        
        
        # Set up the action structures
        self.transit_to = Action("transit_to")
        self.transit_from = Action("transit_from")
        self.mobilize_material = Action("mobilize")
        self.install = Action("install")
        
    def get_mobilize_action(self, pkg):
        """
        Mobilize action for the vessel at port.

        This is a collection of code that looked duplicated between the vessel file and the install_helpers file.
        More than anything it helps give sense of how to calculate the mobilization time based on the vessel specs and the package of materials.

        Parameters
        ----------
        pkg : dict
            The package of materials to be mobilized.

        Returns
        -------
        Action
            Action for mobilizing the vessel.
        """

        # Old vessel mobilize action
        # mobilize_material.addItem("load spooling", duration=1)
        # mobilize_material.addItem("load line", duration=2, dependencies=[("self", "load spooling")])
        # mobilize_material.addItem("load anchor", duration=1)
        # mobilize_material.addItem("load gear", duration=2)
        # mobilize_material.addItem("seafasten", duration=3, dependencies=[
        #     ("self", "load spooling"), ("self", "load line"),
        #     ("self", "load anchor"), ("self", "load gear")
        # ])    
        
        # Old Mobilize function
        # mobilize_material.addItem("mobilize_vessel", duration=self.specs['vessel_specs']['mobilization_time'], dependencies=[])

        # Mobilize action from install_helpers
        winch_speed = self.specs['storage_specs']['winch_speed']*60  # m/hr
        anchor_loading_speed = self.specs['storage_specs']['anchor_loading_speed']
        
        self.mobilize_material.addItem("load spooling", duration=1, dependencies=[])
        self.mobilize_material.addItem("load line", duration=0, dependencies=[self.mobilize_material.items["load spooling"]]) # these need to be ActionItems in an Action object
        self.mobilize_material.addItem("load anchor", duration=0, dependencies=[])
        self.mobilize_material.addItem("load gear", duration=2, dependencies=[])
        self.mobilize_material.addItem("seafasten", duration=3, dependencies=[ # these need to be ActionItems in an Action object
            self.mobilize_material.items["load spooling"], self.mobilize_material.items["load line"],
            self.mobilize_material.items["load anchor"], self.mobilize_material.items["load gear"]
        ])

        for key, item in pkg.items():
            item['obj'].inst['mobilized'] = True
            if key.startswith("sec"):  # agnostic to line type
                self.mobilize_material.items["load line"].duration += item['length'] / winch_speed
                self.state['spool_storage'] -= item['length']
                
            elif key.startswith("anchor"):  # anchor
                if item['load'] > self.specs['storage_specs']['max_deck_load']:
                    raise ValueError(f"item {key} has a load higher than what the vessel can withstand.")
                
                self.mobilize_material.items["load anchor"].duration += anchor_loading_speed  # Assuming 1 anchor load = 1 * speed
                self.state['deck_storage'] -= item['space']
        
            self.state['cargo_mass'] -= item['mass'] # remaining capacity
            self.state['assigned_materials'].append(item['obj'])

        return self.mobilize_material

    def mob(self, time, **kwargs):
        """

        This function is not used yet. Example of what considering port location could look like. 
        
        Initialize the vessel and mobilize to port

        Parameters
        ----------
        time : float
            The current simulation time.
        location : str
            The target location for mobilization.

        Returns
        -------
        None
        """
        # Duration of the activity
        duration = self.specs['vessel_specs']['mobilization_time']
        
        # Vessel location at port
        portLocation = kwargs.get('port_r', [0, 0])
        self.r = portLocation
        self.state["location"] = self.r
        self.state["preceeds"] = "material_mob"

        # Get vessel latest activity
        log = self.getState(time)
        if not log:
            time=duration
        else:
            time = log["time"] + duration

        self.logState(time=time, new_state=self.state)

    def get_transit_to_action(self, distance2port):
        """
        Transit actions for the vessel to a destination from port.
        
        Parameters
        ----------
        distance2port : float
            The distance to the site from port.
            
        Returns
        -------
        transit_to : Action
            Action for transiting to the site from port.
        """

        self.transit_to.addItem("transit_to_site", duration=distance2port/self.specs['transport_specs']['transit_speed'], dependencies=[self.mobilize_material.items["seafasten"]]) # these need to be ActionItems in an Action object

        return self.transit_to


    def get_install_action(self, pkg):
        """
        Creates an action item for installing a materials package from a vessel.

        Parameters
        ----------
        pkg : dict
            The package of materials to be installed.

        Returns
        -------
        action : dict
            The action item for installing the materials.
        """
        
        # set up structure for filling in based on pkg
        self.install.addItem("position onsite", duration=0, dependencies=[])
        self.install.addItem("site survey", duration=0, dependencies=[self.install.items["position onsite"]])
        self.install.addItem("install anchor", duration=0, dependencies=[self.install.items["position onsite"], self.install.items["site survey"]])
        self.install.addItem("rerig deck", duration=0, dependencies=[self.install.items["position onsite"], self.install.items["install anchor"]])
        self.install.addItem("install line", duration=0, dependencies=[self.install.items["install anchor"], self.install.items["rerig deck"]])
        

        def installItem(key):
            '''
            NOT A PUBLIC FUNCTION
            This function installs an item and its dependencies.
            It checks if the item is already installed and if not, it installs its dependencies first.
            Then, it installs the item itself and updates the vessel state.

            Parameters
            ----------
            key : str
                The key of the item to be installed.

            Returns
            -------
            None
            '''
            item = pkg.get(key)
            for dep in item['dependencies']:
                if not pkg[dep]['obj'].inst['installed']:
                    installItem(dep)
            
            if key.startswith("anchor"):
                self.install.items["position onsite"].duration = 2  # from PPI (only once per anchor)           
                self.install.items["site survey"].duration = 2       # from PPI 
                if item['obj'].dd['design']['type']=='suction':
                    pile_fixed = self.specs["vessel_specs"]["pile_fixed_install_time"]
                    pile_depth = 0.005 * abs(item['obj'].r[-1])

                    self.install.items["install anchor"].duration = pile_fixed + pile_depth
                else:
                    # support for other anchor types
                    pass
                
                self.state['deck_storage'] += item.get('space', 0)

            elif key.startswith("sec"):
                if self.install.items["install line"].duration ==0:
                    # first line to install
                    self.install.items["rerig deck"].duration = self.specs['storage_specs'].get('rerig_deck', 0)
                winch_speed = self.specs['storage_specs']['winch_speed']*60  # m/hr
                line_fixed = self.specs["vessel_specs"]["line_fixed_install_time"]
                line_winch = item['length']/winch_speed
                self.install.items["install line"].duration += line_fixed + line_winch
                self.install.items["install line"].dependencies = [self.install.items["install anchor"], self.install.items["rerig deck"]]

                self.state['spool_storage'] += item.get('length', 0)

            item['obj'].inst['installed'] = True
            self.state['cargo_mass'] += item['mass']
            self.state['assigned_materials'].remove(item['obj'])

        for key in pkg.keys():
            installItem(key)
        
        return self.install
    
    def get_transit_from_action(self, distance2port, empty_factor=1.0):
        """
        Transit actions for the vessel from a destination to port.
        
        Parameters
        ----------
        distance2port : float
            The distance to the site from port.
        empty_factor : float, optional
            The factor to account for empty return trip.
            
        Returns
        -------
        transit_from : Action
            Action for transiting from the site to port.
        """

        self.transit_from.addItem("transit_from_site", duration= empty_factor * distance2port/self.specs['transport_specs']['transit_speed'], dependencies=[self.transit_to.items["transit_to_site"], self.install.items["install anchor"], self.install.items["install line"]])

        return self.transit_from


    def logState(self, time, new_state):
        """
        Log and update the vessel state.

        Parameters
        ----------
        time : float
            Current simulation time.
        new_state : dict
            New state information to update and log.

        Returns
        -------
        None
        """
        self.state.update(new_state)
        self.state["log"].append({"time": time, "state": new_state})


    def getState(self, t):
        """
        Retrieve vessel state at a specific time.

        Parameters
        ----------
        t : float
            Time at which to retrieve the vessel state.

        Returns
        -------
        state : dict
            The vessel state at time t, or None if no state exists before time t.
        """
        return next((log for log in reversed(self.state["log"]) if log["time"] <= t), None)


"Port base class"

__author__ = "Rudy Alkarem"

import yaml
from copy import deepcopy
import numpy as np

class Port:
    '''
    Represents a port for staging and logistics operations.
    
    Attributes
    ----------
    name : str
        Name of the port.
    capacity : dict
        Dictionary containing capacity parameters of the port.
    storage : dict
        Current storage state of the port.
    '''

    def __init__(self, file):
        '''
        Initialize a Port object from a configuration file.
        
        Parameters
        ----------
        config_file : str
            Path to the port configuration file.
        '''
        
        # Initialize the base class
        Asset.__init__(self, info)
        
        
        
        self.r = [portDisc['location']['lattitude'], portDisc['location']['longitude']]
        self.pkgs = {}

        # misc
        self.reel_refArea   = 13.5  # m^2
        self.reel_refCap    = 735   # m 
        self.chain_refArea  = 20.5  # m^2 
        self.chain_refLngth = 100   # m


    def staging(self, pkgs):
        """
        Perform staging and update port storage states.
        
        Parameters
        ----------
        pkgs : list
            List of packages to be staged at the port.
            
        Returns
        -------
        remaining_pkgs : list or None
            Packages that couldn't be staged due to capacity constraints,
            or None if all packages were staged successfully.
        """
        remainingPkgs = deepcopy(pkgs)
        # Get some information about polyester and chain lines 
        polyLineLength = 0
        chinLineLength = 0
        polyPkgs = {}
        chinPkgs = {}
        for pkgName, pkg in pkgs.items():
            if pkgName.startswith("sec"):
                if pkg["obj"]["type"]["material"]=="polyester":
                    polyLineLength += pkg["length"]
                    polyPkgs.append(pkg)
                if pkg["obj"]["type"]["material"]=="chain":
                    chinLineLength += pkg["length"]
                    chinPkgs.append(pkg)
        
        # TODO: can we generalize this beyond polyester and chain? Any number of lines and line types?
        # Store polyester lines
        # Compute number of reels required to roll all polyester lines in pkgs
        reelCount = np.ceil(polyLineLength/self.reel_refCap)
        reelAreaTot = reelCount * self.reel_refArea
        if self.state["yard_storage"] > reelAreaTot:
            self.state["yard_storage"] -= reelAreaTot
            for key in polyPkgs:
                remainingPkgs.pop(key, None)
            self.pkgs.update(polyPkgs)

        # Store chains
        # Compute area acquired by chains 
        area_per_unit_meter = self.chain_refArea/self.chain_refLngth  # for a pile of 1.5m tall [135mm chain nominal diameter]
        chinAreaTot = area_per_unit_meter * chinLineLength
        if self.state["yard_storage"] > chinAreaTot:
            self.state["yard_storage"] -= chinAreaTot
            for key in chinPkgs:
                remainingPkgs.pop(key, None)
            self.pkgs.update(chinPkgs)

        # remaining packages:
        for pkgName, pkg in pkgs.items():
            if pkgName.startswith("anchor"):
                if self.state["yard_storage"] > pkg["space"]: 
                    self.state["yard_storage"] -= pkg["space"]
                    remainingPkgs.pop(pkgName, None)
                    self.pkgs.append(pkg)
            if pkgName.startswith("conn"):
                'add logic to stage clump weights and buoys'
                pass
        
        if remainingPkgs=={}:
            remainingPkgs=None
        return remainingPkgs







