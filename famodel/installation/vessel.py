"Vessel base class"

__author__ = "Rudy Alkarem"

from .port import Port
from .action import Action
import yaml 
from copy import deepcopy

class Vessel:
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

    def __init__(self, file):
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
        with open(file) as f:
            vesselDisc = yaml.load(f, Loader=yaml.FullLoader)

        self.name  = vesselDisc['name']
        self.type  = vesselDisc['type']
        self.specs = vesselDisc['specs']
        self.state = {
            "spool_storage": self.specs['storage_specs']['max_spool_capacity'],
            "deck_storage": self.specs['storage_specs']['max_deck_space'],
            "cargo_mass": self.specs['storage_specs']['max_cargo'],
            "log": []
        }
        

    def mobilize(self):
        """
        Mobilize the vessel to a specified location.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        mobilize_material = Action("mobilize")

        mobilize_material.addItem("load spooling", duration=1)
        mobilize_material.addItem("load line", duration=2, dependencies=[("self", "load spooling")])
        mobilize_material.addItem("load anchor", duration=1)
        mobilize_material.addItem("load gear", duration=2)
        mobilize_material.addItem("seafasten", duration=3, dependencies=[
            ("self", "load spooling"), ("self", "load line"),
            ("self", "load anchor"), ("self", "load gear")
        ])        

    def transit(self):
        """
        Transit the vessel to a destination.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """
        transit_to   = Action("transit_to")
        transit_from = Action("transit_from")
        transit_to.addItem("transit_to_site", duration=18, dependencies=[("mobilize_material", "seafasten")])
        transit_from.addItem("transit_from_site", duration=18, dependencies=[("transit_to", "transit_to_site"), ("install", "finalizing")])


    def mob(self, time, **kwargs):
        """
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