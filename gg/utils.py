import networkx as nx
import numpy as np
import pandas as pd
from ase.neighborlist import NeighborList,natural_cutoffs
from ase.data import atomic_numbers, covalent_radii
from ase.io import read as read_atoms
from itertools import combinations
from numpy.linalg import norm
from ase import Atoms

#Function borrowed from surf graph
def node_symbol(atom):
    return "{}_{}".format(atom.symbol, atom.index)

#Function borrowed from surf graph
def relative_position(atoms, neighbor, offset):
    return atoms[neighbor].position + np.dot(offset, atoms.get_cell())

# Function to check if the nodes form a cycle
def is_cycle(G, nodes):
    start_node = next(iter(nodes))  # Get any node as starting point
    subgraph = G.subgraph(nodes)
    try:
        cycle = nx.find_cycle(subgraph, source=start_node)
        return True
    except nx.NetworkXNoCycle:
        return False

def are_points_collinear_with_tolerance(p1, p2, p3, tolerance=1e-7):
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    cross_product = np.cross(p2 - p1, p3 - p1)
    norm = np.linalg.norm(cross_product)
    
    return norm < tolerance
    
#Function to ad adsorbate to atoms object
def add_ads(atoms,ads,ref,offset):
    
    if isinstance(ads, str):
        _ads = read_atoms(ads)
    elif isinstance(ads, Atoms):
        _ads = ads.copy()
    else:
        print("Please provide proper adsorbate file")
    for atom in _ads:
        atom.position[0] = atom.position[0] + ref[0] + offset[0]    
        atom.position[1] = atom.position[1] + ref[1] + offset[1]
        atom.position[2] = atom.position[2] + ref[2] + offset[2]
    sub = atoms.copy()
    sub += _ads
    return sub

def check_contact(atoms,error=0.1,print_contact=False):
    nl=NeighborList(natural_cutoffs(atoms), self_interaction=False,  bothways=True)
    nl.update(atoms)
    close_contact = []
    for index,atom in enumerate(atoms):
        for neighbor, offset in zip(*nl.get_neighbors(index)):
            if sorted((index,neighbor)) not in close_contact:
                atom2 = atoms[neighbor]
                distance = np.linalg.norm(atom.position-relative_position(atoms, neighbor, offset))
                eqm_radii = (covalent_radii[atom.number] + covalent_radii[atom2.number])
                if distance < (1-error)*eqm_radii:
                    if print_contact:
                        print('Close Contact:', node_symbol(atom),node_symbol(atom2),round(distance,2))
                    close_contact.append(sorted((index,neighbor)))
    if close_contact:
        return True
    else:
        return False

def atoms_to_graph(atoms, nl, max_bond=0, max_bond_ratio = 0):
    if max_bond == 0 and max_bond_ratio == 0:
        print("Please Specify bond information")
        return

    G=nx.Graph()
    for index,atom in enumerate(atoms):
        if not G.has_node(node_symbol(atom)):
            G.add_node(node_symbol(atom),index=atom.index, symbol=atom.symbol)
        for neighbor, offset in zip(*nl.get_neighbors(index)):
            atom2 = atoms[neighbor]
            vector = atom.position-relative_position(atoms, neighbor, offset)
            distance = np.linalg.norm(vector)
            eqm_radii = (covalent_radii[atom.number] + covalent_radii[atom2.number])
            check = max(max_bond,eqm_radii*max_bond_ratio)
            if distance > check:
                    continue    
            if not G.has_node(node_symbol(atom2)):
                G.add_node(node_symbol(atom2),index=atom2.index, symbol=atom2.symbol)
            if not G.has_edge(node_symbol(atom),node_symbol(atom2)):
                G.add_edge(node_symbol(atom),node_symbol(atom2),weight = vector, start = index)
    return G

def get_normals(index,atoms,G):

    #Initially, get the right position vectors for the adsorption cluster
    ads_pos = np.zeros((len(index),3))
    initial = index[0]
    ads_pos[0] = np.array([0,0,0])
    for i,j in enumerate(index[1:]):
        atom = atoms[j]
        edge_data = G.get_edge_data(node_symbol(atoms[initial]), node_symbol(atom))
        vector = edge_data['weight']
        start = edge_data['start']
        if start == initial:
            ads_pos[i+1] = vector
        else:
            ads_pos[i+1] = -vector
            
    ads_pos_sum = np.sum(-ads_pos,axis=0)/len(index)
    
    #Construct a Matrix with vectors surrounding the adsorption cluster
    normals = []
    for i,j in enumerate(index):
        atom = atoms[j]
        for neighbor in G.neighbors(node_symbol(atom)):
                n_index = G.nodes[neighbor]['index']
                if n_index not in index:
                    edge_data = G.get_edge_data(node_symbol(atom), neighbor)
                    vector = edge_data['weight']
                    start = edge_data['start']
                    if start == j:
                        normal = vector
                    else:
                        normal = -vector
                    normal = normal + ads_pos_sum + ads_pos[i]
                    normals.append(normal/norm(normal))
    V = np.array(normals)

    #Find the direction that best represents the empty space around the adsorption cluster
    _, _, Vt = np.linalg.svd(V)
    Vt = Vt[-1]
    
    #The best normal could be +Vt or -Vt
    matrix1 = np.array([Vt, -Vt])
    dot_products = np.dot(matrix1, V.T)
    sums_per_vector = np.sum(dot_products, axis=1)
    max_index = np.argmax(sums_per_vector)
    vector_with_smallest_sum = matrix1[max_index]

    ref_pos = ads_pos_sum + atoms[initial].position
    return vector_with_smallest_sum, ref_pos

def generate_sites(atoms, ads, graph, index, coordination, ad_dist = 1.7, contact_error = 0.2):

    possible = list(combinations(index, coordination))
    valid = []

    for cycle in possible:
        if coordination == 1:
            valid.append(list(cycle))

        if coordination == 2:
            if graph.has_edge(node_symbol(atoms[cycle[0]]),node_symbol(atoms[cycle[1]])):
                valid.append(list(cycle))
                
        if coordination == 3:
            nodes = [node_symbol(atoms[i]) for i in cycle]
            pos = [atoms[i].position for i in cycle]
            if is_cycle(graph, nodes):
                if not are_points_collinear_with_tolerance(pos[0],pos[1],pos[2],tolerance=0.01):
                    valid.append(list(cycle))

    movie = []
    for cycle in valid:
        normal, ref_pos = get_normals(cycle,atoms,graph)
        offset = normal*ad_dist/norm(normal)
        ads_copy = ads.copy()
        ads_copy.rotate([0, 0, 1], normal, center=[0,0,0])
        new_atoms = add_ads(atoms,ads_copy,ref_pos,offset=offset)
        if check_contact(new_atoms,error = contact_error):
            print(f"Issue at {cycle}")    
        movie.append(new_atoms)
    return movie

class surface_sites:
    def __init__(self,max_coord, surf_atom_sym=None, max_bond_ratio=0, max_bond=0, contact_error=0.2, com=True):
        self.max_coord = max_coord
        if surf_atom_sym:
            self.surf_atom_sym = surf_atom_sym
        else:
            self.surf_atom_sym = list(set(self.atoms.symbols))
       
        self.max_bond_ratio = max_bond_ratio
        self.max_bond = max_bond
        self.contact_error = contact_error
        self.com = com
        

    def get_surface_sites(self, atoms, self_interaction=False,  bothways=True):
        for sym in atoms.symbols:
            if sym not in list(self.max_coord.keys()):
                print("Incomplete max_coord")
                return 
        nl=NeighborList(natural_cutoffs(atoms), self_interaction=self_interaction,  bothways=bothways)
        nl.update(atoms)
        G = atoms_to_graph(atoms, nl, max_bond_ratio = self.max_bond_ratio, max_bond = self.max_bond)
        sites = []
        for node in G.nodes():
            cord = len([edge for edge in G[node]])
            sum1 = sum([G[node][edge]['weight'] for edge in G[node]])
            index = G.nodes[node]['index']
            symbol = atoms[index].symbol
            diff_cord = self.max_coord[symbol] - cord
            sites.append({'ind': index, 'symbol': symbol,
                         'cord': cord,'diff_cord':diff_cord,'z_coord':atoms[index].position[2]})
            
        df = pd.DataFrame(sites)
        
        if self.com:
            df = df[df.z_coord > atoms.get_center_of_mass()[2]]
            
        df = df[df.diff_cord > 0].sort_values(by=['symbol','cord'])
        
        if isinstance(self.surf_atom_sym, str):
            df = df[df['symbol'] == self.surf_atom_sym]
        else:
            df = df[df['symbol'].isin(self.surf_atom_sym)]
        df = df.sort_values(by=['cord','z_coord'])
        return df, G

        