import os
import rdkit.Chem as Chem
from nxchem_utils import *
from scaff_decomp import *

def match_graph(molgraph, molgraph_list):
    for g_idx, G in enumerate(molgraph_list):
        if nx.is_isomorphic(G, molgraph, node_match=node_exact, edge_match=edge_exact):
            return g_idx

class Vocab(object):
    def __init__(self, filepath):
        self.molgraph_list = self._load_graph_list(filepath)

    def _load_graph_list(self, filepath):
        total_vocab_graphs=[]
        for filename in os.listdir(filepath):
            vocab_filepath = filepath + "/{}".format(filename)
            G = pickle.load(open(vocab_filepath, "rb"))
            total_vocab_graphs.append(G)
        return total_vocab_graphs

    def get_index(self, molgraph):
        return match_graph(molgraph, self.molgraph_list)
    
    def get_molgraph(self, idx):
        return self.molgraph_list[idx]


class MolScaffNode(object):

    def __init__(self, mol, clique=[], smiles=None):
        if mol and smiles: 
            self.smiles = smiles
            self.mol = mol

            self.tri_mol, self.graph = get_triangulated_graph(self.mol)
            self.clique = [x for x in clique] #copy
            self.neighbors = []
        else:
            self.smiles = get_fragments(mol , clique)
            self.smarts = get_smarts_fragments(mol , clique)

            self.mol = get_clique_mol(mol, clique) # use mol object directly
            # self.mol = data_to_mol(mol_to_data(self.mol))

            self.tri_mol, self.graph = get_triangulated_graph(self.mol)

            self.clique = [x for x in clique] #copy
            self.neighbors = []
        
    def add_neighbor(self, nei_node):
        self.neighbors.append(nei_node)

    def recover_G(self, original_graph):
        pass

class MolScaff(object):
    
    def __init__(self, smiles):
        self.smiles = smiles
        self.mol = get_mol(smiles)

        mol = Chem.MolFromSmiles(smiles)
        self.smiles3D = Chem.MolToSmiles(mol, isomericSmiles=True)
        self.smiles2D = Chem.MolToSmiles(mol)

        self.cliques, self.edges, self.triangulated_graph = scaff_decomp(mol)

        self.nodes = []
        self.nodes_dict = {}

        for i, clique in enumerate(self.cliques):
            if isinstance(clique, tuple):
                c = list(clique)
                m = MolScaffNode(mol, c)
                self.nodes_dict[i] = m
            self.nodes.append(m)

        for x,y in self.edges:
            self.nodes[x].add_neighbor(self.nodes[y])
            self.nodes[y].add_neighbor(self.nodes[x])
            
    def size(self):
        return len(self.nodes)
    
    def recover(self, mol_node_u, mol_node_v):
        u_clique = mol_node_u.clique
        v_clique = mol_node_v.clique

        atoms_subgraph = list(set(u_clique + v_clique))
        assemble_label = self.triangulated_graph.subgraph(atoms_subgraph)

        return assemble_label
    
    def nx_scaffold(self):
        graph = nx.Graph()
        saved_graph = nx.Graph()
        for u, v in self.edges:
            u_mol, v_mol = nx_to_mol(self.nodes_dict[u].graph), nx_to_mol(self.nodes_dict[v].graph)
            graph.add_node(u, smiles=Chem.MolToSmiles(u_mol))
            graph.add_node(v, smiles=Chem.MolToSmiles(v_mol))

            saved_graph.add_node(u, mol_node=self.nodes_dict[u])
            saved_graph.add_node(v, mol_node=self.nodes_dict[v])
            assemble_node = self.recover(self.triangulated_graph, self.nodes_dict[u], self.nodes_dict[v])
            saved_graph.add_edge(u, v, assemble=assemble_node)

    



