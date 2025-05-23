"Installation manager for Mooring systems that could be imported from a Project Class in FAModel"

__author__ = "Rudy Alkarem"


import heapq
from .port import Port
from .vessel import Vessel


class InstallManager:
    """
    Manages the overall installation process for mooring systems.
    
    Attributes
    ----------
    vesselList : dict
        Dictionary of registered vessels.
    portList : dict
        Dictionary of registered ports.
    events : list
        Priority queue of scheduled events.
    now : float
        Current simulation time.
    totalDuration : float
        Total duration of the installation.
    logs : list
        List of logged events.
    allTasks : dict
        Dictionary tracking all actions by (agent_name, action_name, item_name).
    pkgs2Stage : list or None
        Packages that need staging later.
    """

    def __init__(self):
        """
        Initialize the installation manager.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.vesselList = {}
        self.portList = {}
        self.events = []
        self.now = 0
        self.totalDuration = 0
        self.logs = []
        self.allTasks = {}  # (agent_name, action_name, item_name) - I need this to keep track of all actions
        self.pkgs2Stage = None

    def registerVessel(self, file):
        """
        Register a vessel with the installation manager.
        
        Parameters
        ----------
        file : str
            Path to the vessel configuration file.
            
        Returns
        -------
        None
        """
        vessel = Vessel(file)
        self.vesselList[vessel.name] = vessel
    
    def registerPort(self, file):
        """
        Register a port with the installation manager.
        
        Parameters
        ----------
        file : str
            Path to the port configuration file.
            
        Returns
        -------
        None
        """
        port = Port(file)
        self.portList[port.name] = port

    def scheduleEvent(self, time, agent, action, params):
        """
        Schedule an event at a given time, checking weather conditions.
        
        Parameters
        ----------
        time : float
            Time at which the event should occur.
        agent : object
            Agent responsible for the action.
        action : str
            Name of the action to perform.
        params : dict
            Parameters for the action.
            
        Returns
        -------
        None
        """
        # check weather
        if self.weather_ok(time):
            heapq.heappush(self.events, (time, agent, action, params))
        else:
            delayTime = 0
            heapq.heappush(self.events, (time+delayTime, agent, action, params))

    def depSchedule(self, agent, action_name, item_name):
        """
        Schedule events based on dependencies of actions.
        
        Parameters
        ----------
        agent : object
            Agent responsible for the action.
        action_name : str
            Name of the action.
        item_name : str
            Name of the item.
            
        Returns
        -------
        None
        """
        key = (agent.name, action_name, item_name)
        item = self.allTasks(key)

        deps_met = True
        for dep in item.dependencies:
            dep_agent_name = agent.name if dep[0] == "self" else dep[0]
            dep_key = (dep_agent_name, dep[0] if dep[0] != "self" else action_name, dep[1])
            dep_item = self.allTasks.get(dep_key)
            if not dep_item or dep_item.status != "done":
                deps_met = False
                break

        if deps_met and item.status == "pending":
            duration = item.duration
            self.scheudleEvent(self.now + duration, agent, "completeItem", {
                "action_name": action_name,
                "item_name": item_name
            })

        
    def run(self):
        """
        Run the installation phase simulation.
        This method processes scheduled events in the priority queue, executing actions
        and updating the current time.
        It continues until all events have been processed.
        The method also logs each event as it is processed.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        None
        """

        print("Events:")
        print(self.events) # how is this determined? How do we loop through actions? What level is actions at? <-- each event has an action
        while self.events:
            t, agent, action, params = heapq.heappop(self.events) # each event is made of different agents with their own actions. 
            print(f"Event: {t}, {agent.name}, {action}, {params}")
            self.now = t
            func = getattr(agent, action) # get the function to call: agent.action (t, params)
            # TODO: what functions can be called here? How do we keep them consistent with the same output formatting?
                # all functions that are in the vessel class and port class. Do we want these all called? How does functions in install_helpers relate?
            completed, new_events = func(t, **params) # call funcion: agent.action(t, params)

            self.logs.append({"time": t, "agent": agent.name, "action": action})

            # TODO: what does this do? <-- schedules events that might be triggered by the current event? 
            for evt in new_events:
                self.scheduleEvent(*evt)
    
    def weather_ok(self, t):
        """
        Check if weather conditions are suitable for operations.
        
        Parameters
        ----------
        t : float
            Current time.
            
        Returns
        -------
        ok : bool
            True if weather conditions are suitable, False otherwise.
        """

        # TODO

        return True
    
    def createPkgs(self, project, stageMode):
        """
        Clusters components [or material Packages] from FAModel project class 
        that need to be mobilized/installed simultaneously.
        
        This function also updates the port staging based on stageMode.
        
        Parameters
        ----------
        project : Project
            The project class containing mooring systems.
        stageMode : int
            Determines what package goes into which port:
                stageMode=1: packages are stored in ports in a first-listed-first-served basis. Whenever the 
                    yard storage capacity of the first port is full, the remaining packages are stored in the second
                    port on the list and so on. If all ports are full and packages are still remaining, they are stored 
                    in self.pkgs2Stage for later staging. 
                stageMode=2: Not developed yet [but this could be used for more smart staging - for instance, 
                    suction piles can be stored separately in a different port and so on.]
            
        Returns
        -------
        None
        """
        # Stage Mode = 1
        def createPkg(moor):
            """
            NOT A PUBLIC FUNCTION
            Create a package of components from a mooring system.

            Parameters
            ----------
            moor : Mooring
                Mooring system to create a package from.

            Returns
            -------
            pkg : dict
                Package of components from the mooring system.
            """
            pkg = {}
            if moor.shared:
                # installation itemization for shared line
                for i, sec in enumerate(moor.dd['sections']):
                    pkg[f"sec_{i}"] = {
                        "obj": sec,
                        "mass": sec['type']['m']*sec['L']/1e3,  # mass [t]
                        "length": sec['L'],                     # length [m]
                        "dependencies": [],
                    }
                    if i>1:
                        pkg[f"sec_{i}"]['dependencies'].append(f"sec_{0}")
                        
                for i, conn in enumerate(moor.dd['connectors']):
                    if conn['m'] > 0:
                        pkg[f"conn_{i}"] = {
                            "obj": conn,
                            "mass": conn['m']/1e3, # mass [t]
                            # "load": ?,  # pressure [t/m^2]    # NOTE: I am not checking load because clump weights could be easily sized in a way to reduce load if needed.
                            "dependencies": [f"sec_{0}"], 
                        }
            else:
                for att in moor.attached_to:
                    if type(att).__name__ == "Anchor":
                        M = att.mass / 1e3        # mass [t]
                        D = att.dd['design']['D'] # diameter [m]
                        L = att.dd['design']['L'] # length [m]
                        area_xy = D * L
                        pkg = {
                            f"anchor_{att.id}": {
                                "obj": att,
                                "mass": M,
                                "load": M/area_xy,
                                "space": area_xy,
                                "dependencies": []
                            },
                        }
                        break

                for i, sec in enumerate(moor.dd['sections']):
                    pkg[f"sec_{i}"] = {
                        "obj": sec,
                        "mass": sec['type']['m']*sec['L']/1e3,
                        "length": sec['L'],
                        "dependencies": [f"anchor_{att.id}"],
                    }
                for i, conn in enumerate(moor.dd['connectors']):
                    if conn['m'] > 0:
                        pkg[f"conn_{i}"] = {
                            "obj": conn,
                            "mass": conn['m']/1e3,
                            "load": conn['m']/conn['v'],
                            "dependencies": [f"anchor_{att.id}"],
                        }         
            return pkg         
           
        pkgs = []

        for moor in project.mooringList.values():
            pkgs.append(createPkg(moor))

        # stage packages at ports
        
        portList = [port for port in self.portList.values()]
        for port in portList:
            remainingPkgs = port.staging(pkgs)
            if remainingPkgs==None:
                break
            else:
                pkgs = remainingPkgs

        self.pkgs2Stage = remainingPkgs  # If there are additional material packages that needs staging later, save here.