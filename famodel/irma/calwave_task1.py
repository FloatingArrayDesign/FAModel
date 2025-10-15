# task1_calwave.py
# Build CalWave Task 1 (Anchor installation) from PDF routes and action table.

from action import Action
from task import Task
from calwave_irma import Scenario
# If you use Scenario/Project to add actions + hold vessels, import it:
from famodel.project import Project  

# ---- Constants from CalWave_IOM_Summary.pdf ----
# Transit times [h] (one-way)
SAN_DIEGO_NationalCity_to_TestSite = 4.5
JAG_NationalCity_to_TestSite       = 4.5
BEYSTER_PointLoma_to_TestSite      = 0.8

# At site support transit
JAG_TestSite                       = 7.5

# Mirror returns
SAN_DIEGO_TestSite_to_Home = SAN_DIEGO_NationalCity_to_TestSite
JAG_TestSite_to_Home       = JAG_NationalCity_to_TestSite
BEYSTER_TestSite_to_Home   = BEYSTER_PointLoma_to_TestSite

def build_task1_calwave(sc):
    """
    sc: scenario/project object that exposes:
        - sc.vessels['San_Diego'], sc.vessels['Jag'], sc.vessels['Beyster']
        - sc.addAction(action_type_name, name, objects=None) -> Action
    Returns:
        task (Task)
    """

    # --- Create Actions ---
    a_mob_sd   = sc.addAction('mobilize', 'mobilize_sandiego')
    a_mob_bys  = sc.addAction('mobilize', 'mobilize_beyster')

    a_load_cargo  = sc.addAction('load_cargo', 'load_cargo_task1',  objects=[])   # add anchors and moorings?

    a_tr_site_sd  = sc.addAction('transit_tug', 'transit_site_sandiego')
    a_tr_site_jag = sc.addAction('transit', 'transit_site_jag')
    a_tr_site_bys = sc.addAction('transit', 'transit_site_beyster')
    a_tr_at_site_jag = sc.addAction('at_site_support', 'at_site_jag')

    a_install_anchor  = sc.addAction('install_anchor',  'install_anchor_task1',  objects=[])
    a_install_mooring = sc.addAction('install_mooring', 'install_mooring_task1', objects=[])

    a_monitor = sc.addAction('monitor_installation', 'monitor_installation_task1')

    a_tr_home_sd  = sc.addAction('transit_tug', 'transit_home_sandiego')
    a_tr_home_jag = sc.addAction('transit', 'transit_home_jag')
    a_tr_home_bys = sc.addAction('transit', 'transit_home_beyster')

    # --- Assign assets and compute durations/costs where needed ---
    # Mobilize / Load: let evaluateAssets compute time/cost from capabilities
    a_mob_sd.evaluateAssets(  {'operator': sc.vessels['San_Diego']} )
    a_mob_bys.evaluateAssets( {'operator': sc.vessels['Beyster']} )
    a_load_cargo.evaluateAssets(  {'operator': sc.vessels['San_Diego']} )

    # Transit site: set duration from PDF table; still assign a vessel so costing works (if your calc uses day rate)
    a_tr_site_sd.duration  = SAN_DIEGO_NationalCity_to_TestSite
    a_tr_site_jag.duration = JAG_NationalCity_to_TestSite
    a_tr_site_bys.duration = BEYSTER_PointLoma_to_TestSite
    a_tr_at_site_jag.duration = JAG_TestSite
    # Optionally call evaluateAssets to pick up cost models:
    a_tr_site_sd.evaluateAssets(  {'carrier': sc.vessels['San_Diego']} )
    a_tr_site_jag.evaluateAssets( {'operator': sc.vessels['Jag']} )
    a_tr_site_bys.evaluateAssets( {'carrier': sc.vessels['Beyster']} )

    # Install: vessel + tug
    a_install_anchor.evaluateAssets(  {'operator': sc.vessels['San_Diego'], 'carrier': sc.vessels['Jag']} )
    a_install_mooring.evaluateAssets( {'operator': sc.vessels['San_Diego'], 'carrier': sc.vessels['Jag']} )
    a_tr_at_site_jag.evaluateAssets(  {'operator': sc.vessels['Jag']}) # Need to include this when the tug is included in install anchor and mooring

    # Monitor (Beyster)
    a_monitor.evaluateAssets( {'support': sc.vessels['Beyster']} )

    # Transit homeport: set duration from PDF and assign asset for costing
    a_tr_home_sd.duration  = SAN_DIEGO_TestSite_to_Home
    a_tr_home_jag.duration = JAG_TestSite_to_Home
    a_tr_home_bys.duration = BEYSTER_TestSite_to_Home
    # Optionally call evaluateAssets to pick up cost models:
    a_tr_home_sd.evaluateAssets(  {'carrier': sc.vessels['San_Diego']} )
    a_tr_home_jag.evaluateAssets( {'operator': sc.vessels['Jag']} )
    a_tr_home_bys.evaluateAssets( {'carrier': sc.vessels['Beyster']} )

    # --- Compose the action list for the Task ---
    actions = [
        a_mob_sd, a_mob_bys,
        a_load_cargo,
        a_tr_site_sd, 
        a_tr_site_jag, a_tr_site_bys,
        a_tr_at_site_jag,
        a_install_anchor, a_install_mooring,
        a_monitor,
        a_tr_home_sd, 
        a_tr_home_jag, a_tr_home_bys,
    ]

    # --- Define the sequencing (dependencies) ---
    # Keys are action names; values are lists of prerequisite action names.
    action_sequence = {
        # Mobilize and load
        a_mob_sd.name: [],
        a_mob_bys.name: [],
        a_load_cargo.name: [a_mob_sd.name],

        # Transit to site (can be parallel), but after loading
        #a_tr_site_sd.name:  [a_load_cargo.name],
        a_tr_site_jag.name: [a_mob_sd.name],   # tug stays on site for support
        a_tr_site_bys.name: [a_mob_bys.name],  # support can leave when ready

        # Installs: require everyone onsite
        a_install_anchor.name:  [a_tr_site_jag.name, a_tr_site_jag.name],
        a_install_mooring.name: [a_install_anchor.name],   # mooring after anchors

        # Monitoring: runs during/after install (simplify: after mooring)
        a_monitor.name: [a_install_mooring.name, a_tr_site_bys.name],

        # Transit homeport: everyone returns after monitoring
        #a_tr_home_sd.name:  [a_monitor.name],
        a_tr_home_jag.name: [a_monitor.name],
        a_tr_home_bys.name: [a_monitor.name],
    }

    # --- Build and return the Task (plots sequence with CPs) ---
    task = Task(actions=actions, action_sequence=action_sequence, name='CalWave_Task1_AnchorInstall')
    return task

if __name__ == '__main__':
    # Example skeleton: adjust to your actual Scenario/Project initializer
    project = Project(file='calwave_ontology.yaml', raft=False) # for Mac
    # create moorpy system of the array, include cables in the system
    project.getMoorPyArray(cables=1)
    sc = Scenario()
    task = build_task1_calwave(sc)
    # (The Task constructor will plot the sequence graph automatically)
    pass
