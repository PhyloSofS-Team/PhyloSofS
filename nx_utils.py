import networkx as nx
from distutils.version import LooseVersion

# successors, predecessors and edges returns a generator/iterator after version 2.0 :
if LooseVersion(nx.__version__) < LooseVersion('2.0'):
    def successors(nx_graph, node):
        return nx_graph.successors(node)
    def predecessors(nx_graph, node):
        return nx_graph.predecessors(node)
    def get_edge_list(nx_graph):
        return nx_graph.edges()
else:
    def successors(nx_graph, node):
        return list(nx_graph.successors(node)) # Use list to get a list from the generator
    def predecessors(nx_graph, node):
        return list(nx_graph.predecessors(node))
    def get_edge_list(nx_graph):
        return [ e for e in nx_graph.edges() ]

# NetworkX deprecates max_flow in favor of maximum_flow_value in version 1.9 :
if LooseVersion(nx.__version__) < LooseVersion('1.9'):
    def maximum_flow_value(g, x, y):
        return nx.max_flow(g, x, y)
else:
    def maximum_flow_value(g, x, y):
        return nx.maximum_flow_value(g, x, y)
