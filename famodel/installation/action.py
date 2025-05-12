"Action base class"

__author__ = "Rudy Alkarem"

class ActionItem: 
    def __init__(self, name, duration, cost=0, dependencies=None):
        self.name = name
        self.duration = duration
        self.cost = cost
        self.dependencies = dependencies if dependencies else []

    def add_dependency(self, action_name, item_name):
        """Add dependency to that action item"""
        self.dependencies.append((action_name, item_name))

class Action:
    def __init__(self, name):
        self.name = name
        self.items = {}

    def addItem(self, name, duration, cost=0, dependencies=None):
        """Add item to that action"""
        deps = []
        if dependencies:
            for dep_action, dep_item in dependencies:
                action = self.name if dep_action == "self" else dep_action
                deps.append((action, dep_item))

        item = ActionItem(name, duration, cost, deps)
        self.item[name] = item
    
    def get_dependencies(self, item_name):
        """Return dependencies of a given item name"""
        item = self.items[item_name]
        return item.dependencies
    
    def __repr__(self):
        out = f"Action({self.name}):\n"
        for name, item in self.items.items():
            out += f"  {name}: duration={item.duration}, cost={item.cost}, deps={item.dependencies}\n"
        return out    