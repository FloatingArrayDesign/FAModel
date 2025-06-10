import matplotlib.pyplot as plt
import networkx as nx

"Action base class"

__author__ = "Rudy Alkarem"

class ActionItem:
    """
    Represents an item within an action.

    Parameters
    ----------
    name : str
        The name of the item.
    duration : float
        The duration of the item in hours.
    cost : float, optional
        The cost of the item.
    dependencies : list, optional
        A list of dependencies for the item.

    Attributes
    ----------
    name : str
        The name of the item.
    duration : float
        The duration of the item in hours.
    cost : float
        The cost of the item.
    dependencies : list
        A list of dependencies for the item.
    """

    def __init__(self, name, duration, cost=0, dependencies=None):
        """
        Initializes an ActionItem object.

        Parameters
        ----------
        name : str
            The name of the item.
        duration : float
            The duration of the item in hours.
        cost : float, optional
            The cost of the item (default is 0).
        dependencies : list, optional
            A list of dependencies for the item (default is None).

        Returns
        -------
        None
        """
        self.name = name
        self.duration = duration
        self.cost = cost
        self.dependencies = dependencies if dependencies else []

    def add_dependency(self, action_name, item_name):
        """
        Adds a dependency to the item.

        Parameters
        ----------
        action_name : str
            The name of the action that the dependency belongs to.
        item_name : str
            The name of the item that the dependency belongs to.

        Returns
        -------
        None
        """
        self.dependencies.append((action_name, item_name))

class Action:
    """
    Represents an action in the installation process, such as transport or deployment.

    Parameters
    ----------
    name : str
        The name of the action.

    Attributes
    ----------
    name : str
        The name of the action.
    items : dict
        A dictionary of ActionItem objects associated with the action.
    """

    def __init__(self, name):
        """
        Initializes an Action object.

        Parameters
        ----------
        name : str
            The name of the action.

        Returns
        -------
        None
        """
        self.name = name
        self.items = {}

    def addItem(self, name, duration, cost=0, dependencies=None):
        """
        Adds an item to the action.

        Parameters
        ----------
        name : str
            The name of the item.
        duration : float
            The duration of the item in hours.
        cost : float, optional
            The cost of the item.
        dependencies : list, optional
            A list of dependencies for the item.

        Returns
        -------
        None
        """
        # Old code from when dependencies were tuples of (action, item). Now are just items themselves.
        # deps = []
        # if dependencies:
        #     for dep_action, dep_item in dependencies:
        #         action = self.name if dep_action == "self" else dep_action
        #         deps.append((action, dep_item))

        item = ActionItem(name, duration, cost, dependencies)
        self.items[name] = item

    def get_dependencies(self, item_name):
        """
        Returns the dependencies of a given item.

        Parameters
        ----------
        item_name : str
            The name of the item.

        Returns
        -------
        dependencies : list
            A list of dependencies for the item.
        """
        item = self.items[item_name]
        return item.dependencies

    def __repr__(self):
        """
        Returns a string representation of the action.

        Returns
        -------
        out : str
            A string representation of the action.
        """
        out = f"Action({self.name}):\n"
        for name, item in self.items.items():
            out += f"  {name}: duration={item.duration}, cost={item.cost}, deps={item.dependencies}\n"
        return out
    

    def visualize(self):
        """
        Visualizes the action items as a directed graph.

        Parameters
        ----------
        action : dict
            The action items to be visualized.

        Returns
        -------
        None
        """
        # Create the graph
        G = nx.DiGraph()
        for key in self.items.keys():
            item = self.items[key]
            for dep in item.dependencies:
                G.add_edge(dep, item, duration=item.duration)  # Store duration as edge attribute

        # Compute longest path & total duration
        longest_path = nx.dag_longest_path(G, weight='duration')
        longest_path_edges = list(zip(longest_path, longest_path[1:]))  # Convert path into edge pairs
        total_duration = sum(node.duration for node in longest_path)

        if len(longest_path)>=1:
            last_node = longest_path[-1]  # Identify last node of the longest path
            # Define layout
            pos = nx.shell_layout(G)        
            # Draw all nodes and edges (default gray)
            nx.draw(G, pos, with_labels=True, labels = {node: node.name for node in G.nodes()}, node_size=500, node_color='skyblue', font_size=10, font_weight='bold', font_color='black', edge_color='gray')

            # Highlight longest path in red
            nx.draw_networkx_edges(G, pos, edgelist=longest_path_edges, edge_color='red', width=2)

            # Annotate last node with total duration in red
            plt.text(pos[last_node][0], pos[last_node][1] - 0.1, f"{total_duration:.2f} hr", fontsize=12, color='red', fontweight='bold', ha='center') 

            plt.suptitle(f"Action: {self.name} - Total Duration: {total_duration:.2f} hr")         

        else:
            pass