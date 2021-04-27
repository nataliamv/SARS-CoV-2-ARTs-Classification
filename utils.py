### import libraries

import networkx as nx
from Bio import Phylo
import re 

from itertools import count
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
import matplotlib.colors as clr


#####################
### preprocessing ###
#####################

def get_tips(tree):
    '''creates a dictionary of tree
       tips (value) by name (key)'''
    labels = {}
    for clade in tree.get_terminals(order='preorder'):
        if clade.name: 
            if clade.name in labels:
                raise ValueError("Duplicate key: %s" % clade.name) #double check there are no duplicate names
            labels[clade.name] = clade
    return labels

def collapse_tips(tree):
    tips = get_tips(tree)
    for tip in tips:
        tree.collapse(tip)
    return tree

def add_names(tree):
    '''adds node name to name field in NEXUS tree
       and creates a dictionary of clades (value) by name (key)''' 
    names = {}
    for i, clade in enumerate(tree.find_clades()):
        if clade.name:
            clade.name = "%d_%s" % (i, clade.name)
        else:
            clade.name = str(i)  # if clade has no name 
        names[clade.name] = clade
    return names

def state_collapse(tree,k):
    '''merges a node with its parent only if they
       share the same amino acid state 
       adds 1 to the weight of the parent node
       removes nodes with state 'X'
       k: depth of branch ancestors'''
    n = 0
    for clade in tree.find_clades(): ### add weights
        clade.weight = 1
        n = n+1 
    tree_depths = tree.depths(unit_branch_lengths=True)
    for x in reversed(list(tree_depths.keys())):
        current = x ### change assign current max to name
        path = tree.get_path(current)
        if len(path) > k: ### stops at branch
            parent = path[-2]
            if 'X' in current.comment:
                tree.collapse(current)
            elif (current.comment == parent.comment): 
                if int(current.name) - int(parent.name) == 1:
                    tree.collapse(current)
                    parent.weight = parent.weight + current.weight 
                else:
                    sibling = parent.root.clades[0]
                    if current.comment == sibling.comment:
                        tree.collapse(current)
                        parent.weight = parent.weight + current.weight 
                    else:
                        pass
            else: 
                pass
            del tree_depths[current] 
    return tree

def get_root(D):
    return [n for n,d in D.in_degree() if d==0]

def phylo_to_digraph(tree): 
    '''convert tree to networkx digraph with amino acid state and weight as attributes'''
    tree_nx = Phylo.to_networkx(tree)
    node_weights = {clade.name:  clade.weight for clade in tree.find_clades()}
    node_comments = {clade.name:  clade.comment for clade in tree.find_clades()}
    clean_comments = {k:re.sub('[^a-zA-Z]+', '', v) for k, v in node_comments.items() if v!=None}
    mapping = {clade: clade.name for clade in tree.find_clades()}
    digraph = nx.relabel_nodes(tree_nx, mapping)
    nx.set_node_attributes(digraph, node_weights, 'weight')
    nx.set_node_attributes(digraph, clean_comments, 'comment')
    root = '0_NODE_0000000'
    T = nx.relabel_nodes(digraph, {root: '0'}, copy=False)
    D = nx.DiGraph(T) 
    arrows = list(D.edges())
    edgelist = [(x,y) for (x,y) in arrows if int(y)<int(x)]
    D.remove_edges_from(edgelist)
    return D


def low_freq_var(D,w):
    '''identify if low frequency variants appear within the reduced tree'''
    terminal_nodes = [n for n,d in D.out_degree() if d==0]
    if '3' in terminal_nodes:
        terminal_nodes.remove('3') ### hard coded special case
    weights = [n for n,k in D.nodes(data=True) if k['weight']<=w]
    low_freq = list(set(terminal_nodes) & set(weights))
    if len(low_freq) == 0: 
        return False, low_freq
    else:
        return True, low_freq

    
    
#####################
### clasification ###
#####################

def conserved(digraph):
    S = set()
    for node in digraph.nodes:
        if node != '0':
            state = digraph.nodes[node]['comment']
            S.add(state)
    if len(S) == 1:
        return True
    else:
        return False

def terminal_nodes(digraph):
    '''list of terminal nodes'''
    terminal_nodes = [x for x in digraph.nodes() if digraph.degree(x)==1]
    return terminal_nodes    
    
def states_path(digraph,u,v):
    '''generate list of states of path from node u to terminal node v'''
    amino_acid_states =[]
    for path in nx.all_simple_paths(digraph,u,v):
        for node in path:
            amino_acid_states.append(list(digraph.nodes[node]['comment'])[0])
    amino_acid_states = [state for state in amino_acid_states if state != '-'] 
    amino_acid_states = [state for state in amino_acid_states if state != 'X']
    return amino_acid_states
    
def all_states_paths(digraph,ancestors): ### de lo mas reciente hacia atras following coalescence
    '''returns all paths of amino acid states between a list
    of ancestor nodes and the terminal nodes of the digraph'''
    N = terminal_nodes(digraph)
    all_states_paths=[0] * len(N)
    i=0
    for u in ancestors:
        for v in N:
            path = states_path(digraph,u,v)
            if not path: 
                pass
            else:
                all_states_paths[i] = path 
                i=i+1
    return all_states_paths[:i]
    
def inter_homoplasy(D,p): 
    '''check if homoplasy occurs between lineages'''
    l1 = [a for path in p[0] for a in path]
    l2 = [a for path in p[1] for a in path]
    amino_acids = list(set(l1) & set(l2))
    if not list(set(l1) & set(l2)): 
        return False
    else:
        if (l1[0]==l2[0]) & (len(amino_acids)==1): 
            return False
        else:
            return True
    
def intra_homoplasy(D,p):
    '''check if homoplasy occurs within lineages'''
    p0 = [path for paths in p for path in paths]
    for path in p0:
        path = path[1:] #remove ancestor
        l1=['X']
        l2=[]
        i=0
        for x in path:
            if x != l1[i]:
                l1.append(x)
                i=i+1
        for y in l1:
            if y not in l2:
                l2.append(y)
            else:
                return True
                break
    return False

def stepwise(D,p):
    p0 = [path for paths in p for path in paths] #flatten 
    for path in p0:
        path = path[1:] #remove ancestor
        path = [x for i, x in enumerate(path) if i == 0 or x != path[i-1]]
        if (len(path)>=2) & (len(set(path)) == len(path)):
            return True
            break
    return False

def get_branches(D,k):
    '''k is the distance from the root'''
    b = list(D.successors('0'))
    if k == 0:
        return '0'
    elif k == 1:
        return b
    else:
        c=[]
        for x in b:
            c += list(D.successors(x))
        return c
    
    
#####################
### visualisation ###
#####################


def visualise_ART(D): ### add options save, attributes shown, make sure root is in right pos.
    nodes = D.nodes()
    # get unique groups
    groups = set(nx.get_node_attributes(D,'comment').values())
    c_mapping = dict(zip(sorted(groups),count()))
    colors = [c_mapping[D.nodes[n]['comment']] for n in nodes if n!='0']
    c = max(colors)+1
    colors.append(c)
    pos = graphviz_layout(D, prog='dot')  # positions for all nodes
    new_pos = {n:(-q,p) for n,(p,q) in pos.items()}
    pos = new_pos
    nx.draw_networkx_nodes(D, pos, nodelist=nodes, node_color=colors, node_size=400, cmap = plt.cm.get_cmap('rainbow'))
    nx.draw_networkx_edges(D, pos)
    #labels = {n: D.nodes[n]['weight'] for n in nodes if n!='0'}
    nx.draw_networkx_labels(D, pos)#, labels=labels)
    plt.gca().invert_yaxis()
    plt.axis("off")
    plt.show() 