# standard imports
import networkx as nx

# custom imports


def addEdges(edges):
    """
    Docstring for addEdges
    
    :param edges: Description
    """
    p = edges.shape[0]
    G = nx.Graph()
    for i in range(0, p):
        for j in range(0, p):
            cond = (j >= i) | (edges[i, j] == 0)
            if cond:
                continue
            G.add_edges_from([(j, i)])
    return G


def colorCoreGenes(edges):
    """
    Docstring for colorCoreGenes
    
    :param edges: Description
    """
    color_map = ["powderblue"]
    p = edges.shape[0]
    for node in range(1, p):
        nodek = edges[:, node]
        if nodek[0] == 1:
            color_map.append("tomato")
        else:
            color_map.append("lightgrey")
    return color_map


def mapNodeNames(edges, newlabs):
    """
    Docstring for mapNodeNames
    
    :param edges: Description
    :param newlabs: Description
    """
    p = edges.shape[0]
    if len(newlabs) != p:
        raise ValueError("size of edge matrix does not equal length of new node names")
    old = range(0, p)
    new = {}
    for _ in range(0, p):
        new[old[_]] = newlabs[_]
    return new
