"Port base class"

__author__ = "Rudy Alkarem"

import yaml
from copy import deepcopy
import numpy as np

class Port:
    def __init__(self, file):
        with open(file) as f:
            portDisc = yaml.load(f, Loader=yaml.FullLoader)
        
        self.name = portDisc['name']
        self.specs = portDisc['specs']
        self.state = {
            "yard_storage": self.specs['storage_specs']['yard_storage_capacity'],
            "wet_storage": self.specs['storage_specs']['wet_storage_capacity'], 
            "log": []
        }
        
        self.r = [portDisc['lattitude'], portDisc['longitude']]
        self.pkgs = {}

        # misc
        self.reel_refArea   = 13.5  # m^2
        self.reel_refCap    = 735   # m 
        self.chain_refArea  = 20.5  # m^2 
        self.chain_refLngth = 100   # m

    def staging(self, pkgs):
        """Perform staging and update port storage states."""

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
        """Log and update port state."""
        self.state.update(new_state)
        self.state["log"].append({"time": time, "state": new_state})

    def getState(self, t):
        """Retrieve port state at time t."""
        return next((log for log in reversed(self.state["log"]) if log["time"] <= t), None)

