"Port base class"

__author__ = "Rudy Alkarem"

import yaml
from copy import deepcopy
import numpy as np

class Port:
    """
    Represents a port for staging and logistics operations.
    
    Attributes
    ----------
    name : str
        Name of the port.
    capacity : dict
        Dictionary containing capacity parameters of the port.
    storage : dict
        Current storage state of the port.
    """

    def __init__(self, file):
        """
        Initialize a Port object from a configuration file.
        
        Parameters
        ----------
        config_file : str
            Path to the port configuration file.
            
        Returns
        -------
        None
        """
        with open(file) as f:
            portDisc = yaml.load(f, Loader=yaml.FullLoader)
        
        self.name = portDisc['name']
        self.specs = portDisc['specs']
        self.state = {
            "yard_storage": self.specs['storage_specs']['yard_storage_capacity'],
            "wet_storage": self.specs['storage_specs']['wet_storage_capacity'], 
            "log": []
        }
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

    def logState(self, time, new_state):
        """
        Log and update the port state.

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
        Retrieve port state at a specific time.

        Parameters
        ----------
        t : float
            Time at which to retrieve the port state.

        Returns
        -------
        state : dict
            The port state at time t, or None if no state exists before time t.
        """
        return next((log for log in reversed(self.state["log"]) if log["time"] <= t), None)

