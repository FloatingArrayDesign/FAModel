"Installation manager for Mooring systems that could be imported from a Project Class in FAModel"

__author__ = "Rudy Alkarem"


import heapq
from fadesign.installation.port import Port
from fadesign.installation.vessel import Vessel


class InstallManager:

    def __init__(self):
        self.vesselList = {}
        self.portList = {}
        self.events = []
        self.now = 0
        self.totalDuration = 0
        self.logs = []
        self.allTasks = {}  # (agent_name, action_name, item_name) - I need this to keep track of all actions
        self.pkgs2Stage = None

    def registerVessel(self, file):
        """Register a vessel"""
        vessel = Vessel(file)
        self.vesselList[vessel.name] = vessel
    
    def registerPort(self, file):
        """Register a port"""
        port = Port(file)
        self.portList[port.name] = port

    def scheduleEvent(self, time, agent, action, params):
        """Schedule an event at a given time"""
        # check weather
        if self.weather_ok(time):
            heapq.heappush(self.events, (time, agent, action, params))
        else:
            delayTime = 0
            heapq.heappush(self.events, (time+delayTime, agent, action, params))

    def depSchedule(self, agent, action_name, item_name):
        """Schedule events based on dependencies of these actions"""
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
        """Run the installation phase"""
        while self.events:
            t, agent, action, params = heapq.heappop(self.events)
            self.now = t
            func = getattr(agent, action)
            completed, new_events = func(t, **params)
            self.logs.append({"time": t, "agent": agent.name, "action": action})
            for evt in new_events:
                self.scheduleEvent(*evt)
    
    def weather_ok(self, t):
        """Check for weather conditions"""
        return True
    
    def createPkgs(self, project, stageMode):
        """Clusters components [or material Packages] from FAModel project class that needs to bre mobilized/
        installed simultaneously.
        
        This function also updates the port staging based on stageMode.

        Params:
            project (Project): The project class
            stageMode (int): This determines what package goes into which port:
                stageMode=1: packages are stored in ports in a first-listed-first-served basis. Whenever the
                  yard storage capacity of the first port is full,
                  the remaining packages are stored in the second port on the list and so on. If all ports are
                  full and packages are yet remaining, they are stored in self.pkgs2Stage for later staging.
                stageMode=2: Not developed yet [but this could be used for more smart staging - for instance, 
                  suction piles can be stored separately in a different port and so on.]

        """

        # Stage Mode = 1
        def createPkg(moor):
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