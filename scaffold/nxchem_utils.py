import networkx as nx
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
import itertools
import networkx.algorithms.isomorphism as iso
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor
import copy


KEY = "ghost"
def mol_to_nx(mol, skip_unattached=False):
    G = nx.Graph()

    for bond in mol.GetBonds():
        try: bond.GetBoolProp(KEY)
        except: bond.SetBoolProp(KEY, False)

    # if not mol.GetNumConformers():
    #     rdDepictor.Compute2DCoords(mol)
    # conf = mol.GetConformer()
    # Chem.WedgeMolBonds(mol,conf)

    for atom in mol.GetAtoms():
        if skip_unattached and not atom.GetNeighbors(): continue
        G.add_node(atom.GetIdx(),
                   symbol=atom.GetSymbol(),
                   formal_charge=atom.GetFormalCharge(),
                   chiral_tag=atom.GetChiralTag(),
                   hybridization=atom.GetHybridization(),
                   num_explicit_hs=atom.GetNumExplicitHs(),
                   is_aromatic=atom.GetIsAromatic(),
                   map_num=atom.GetAtomMapNum())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType(),
                   bond_dir=bond.GetBondDir(),
                   ghost=bond.GetBoolProp(KEY),
                   color='r' if bond.GetBoolProp(KEY) else 'b',
                   )
    return G

def nx_to_mol(G):
    mol = Chem.RWMol()
    # Chem.rdDepictor.Compute2DCoords(mol)
    atomic_nums = nx.get_node_attributes(G, 'symbol')
    chiral_tags = nx.get_node_attributes(G, 'chiral_tag')
    formal_charges = nx.get_node_attributes(G, 'formal_charge')
    node_is_aromatics = nx.get_node_attributes(G, 'is_aromatic')
    node_hybridizations = nx.get_node_attributes(G, 'hybridization')
    num_explicit_hss = nx.get_node_attributes(G, 'num_explicit_hs')
    map_nums = nx.get_node_attributes(G, 'map_num')
    node_to_idx = {}
    for node in G.nodes():
        # print(node, atomic_nums[node], num_explicit_hss[node], node_is_aromatics[node])
        a=Chem.Atom(atomic_nums[node])
        a.SetChiralTag(chiral_tags[node])
        a.SetFormalCharge(formal_charges[node])
        a.SetIsAromatic(node_is_aromatics[node])
        a.SetHybridization(node_hybridizations[node])
        a.SetNumExplicitHs(num_explicit_hss[node])
        a.SetAtomMapNum(map_nums[node])
        idx = mol.AddAtom(a)
        node_to_idx[node] = idx

    
    # mol.UpdatePropertyCache()
    bond_types = nx.get_edge_attributes(G, 'bond_type')
    ghost = nx.get_edge_attributes(G, 'ghost')
    bond_dir = nx.get_edge_attributes(G, 'bond_dir')
    for edge in G.edges():
        first, second = edge
        ifirst = node_to_idx[first]
        isecond = node_to_idx[second]
        bond_type = bond_types[first, second]

        # print(bond_type, first, second) 
        mol.AddBond(ifirst, isecond, bond_type)
        new_bond = mol.GetBondBetweenAtoms(ifirst, isecond)
        new_bond.SetBoolProp(KEY, ghost[first, second])
        # new_bond.SetBondDir(bond_dir[first, second])
        # if new_bond.GetBondDir() != Chem.BondDir.NONE: print(new_bond.GetBondDir())
        
    # Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    try: 
        Chem.SanitizeMol(mol)
        # Chem.rdmolops.AssignChiralTypesFromBondDirs(mol)
    except: pass
    return mol

import pickle

def mol_to_data(mol, filename="vocab.txt"):
    G = mol_to_nx(mol)
    pickle.dump(G, open(filename, "wb"))
    return filename

def data_to_mol(filename):
    G = pickle.load(open(filename, "rb"))
    mol = nx_to_mol(G)
    return mol

def save_mol_img(mol_ori, filename="mol1"):
    mol = Chem.Mol(mol_ori)
    label = "molAtomMapNumber"
    for atom in mol.GetAtoms():
        atom.SetProp(label, str(atom.GetIdx()))
    img = Draw.MolToImage(mol)
    img.save(f"{filename}.png")
# ON NETWORKX GRAPH

def node_equal_iso(node1, node2):
    return node1["symbol"] == node2["symbol"] and node1["formal_charge"] == node2["formal_charge"] \
        and node1["map_num"] == node2["map_num"] and node1["is_aromatic"] == node2["is_aromatic"] \
            and node1["num_explicit_hs"] == node2["num_explicit_hs"]

def node_equal_iso2(node1, node2): # honeycomb
    return node1["symbol"] == node2["symbol"] and node1["formal_charge"] == node2["formal_charge"] \

def node_exact(node1, node2):
    return node1["symbol"] == node2["symbol"] and node1["formal_charge"] == node2["formal_charge"] and node1["chiral_tag"] == node2["chiral_tag"] \
        and node1["is_aromatic"] == node2["is_aromatic"] and node1["num_explicit_hs"] == node2["num_explicit_hs"] and \
        node1["chiral_tag"] == node2["chiral_tag"] and node1["hybridization"] == node2["hybridization"]

def edge_exact(edge1, edge2):
    return edge1["bond_type"] == edge2["bond_type"] and \
        edge1["ghost"] == edge2["ghost"] and edge1["bond_dir"] == edge2["bond_dir"]

def ring_edge_equal_iso(edge1, edge2):
    return edge1["bond_type"] == edge2["bond_type"] and \
        edge1["ghost"] == edge2["ghost"]

def ring_edge_equal_iso2(edge1, edge2):
    return edge1["bond_type"] == edge2["bond_type"]


def copy_node_attr(G, idx):
    val = {
        "symbol": G.nodes[idx]["symbol"],
        "chiral_tag": G.nodes[idx]["chiral_tag"],
        "formal_charge": G.nodes[idx]["formal_charge"],
        "is_aromatic": G.nodes[idx]["is_aromatic"],
        "hybridization": G.nodes[idx]["hybridization"],
        "num_explicit_hs": G.nodes[idx]["num_explicit_hs"],
        "map_num": G.nodes[idx]["map_num"],
    }

    return val

def node_equal(a1, a2):
    return a1["symbol"] == a2["symbol"] and a1["formal_charge"] == a2["formal_charge"]


def ring_edge_equal(G1, G2, b1, b2, reverse=False):
    bond_prop = G1.get_edge_data(*b1)["ghost"] == G2.get_edge_data(*b2)["ghost"]
    # bond_prop = G1.get_edge_data(*b1) == G2.get_edge_data(*b2)
    if reverse: b2 = b2[::-1]

    return node_equal(G1.nodes[b1[0]], G2.nodes[b2[0]]) and node_equal(G1.nodes[b1[1]], G2.nodes[b2[1]]) and bond_prop

def draw_mol(cand_G, numb=0, attr=['symbol', 'bond_type', 'color'], folder="subgraph", label=""):
    plt.clf()
    symbol, bond_type, color = attr

    pos = nx.spring_layout(cand_G)
    nx.draw(cand_G, pos)
    node_labels = nx.get_node_attributes(cand_G, symbol)
    node_labels = {k : "       ({})".format(v) for k, v in node_labels.items()}
    nx.draw_networkx_labels(cand_G, pos, node_labels)
    edge_labels = nx.get_edge_attributes(cand_G, bond_type)
    nx.draw_networkx_edge_labels(cand_G, pos, edge_labels)
    colors_cand_G = nx.get_edge_attributes(cand_G, color).values()
    nx.draw(cand_G, pos, node_color="yellow", with_labels=True, edge_color=colors_cand_G)
    # plt.show()
    plt.savefig("{}/show{}_{}.png".format(folder, numb, label))
    plt.clf()

    return


def get_triangulated_graph(mol):
    if mol and mol.GetNumAtoms() == 3 and mol.GetNumBonds() <= 2:
        triangulated_mol = Chem.RWMol(Chem.MolFromSmiles(''))
        for atom in mol.GetAtoms():
            new_atom = copy_atom(atom)
            triangulated_mol.AddAtom(new_atom)

        subset_possible_bonds = list(itertools.combinations(mol.GetAtoms(), 2))
        subset = [(bond[0].GetIdx(), bond[1].GetIdx()) for bond in subset_possible_bonds]

        # G = nx.Graph()
        for bond in subset:
            a1, a2 = bond[0], bond[1]
            bond_obj = mol.GetBondBetweenAtoms(a1, a2)
            if bond_obj:
                triangulated_mol.AddBond(a1, a2, order=bond_obj.GetBondType())
                new_bond = triangulated_mol.GetBondBetweenAtoms(a1, a2)
                new_bond.SetBoolProp(KEY, False)
                # G.add_edge(a1, a2, order=bond_obj.GetBondTypeAsDouble())
            else:
                triangulated_mol.AddBond(a1, a2, order=Chem.BondType.SINGLE)
                new_bond = triangulated_mol.GetBondBetweenAtoms(a1, a2)
                new_bond.SetBoolProp(KEY, True)
                # G.add_edge(a1, a2, order=0.0)        

        mol = triangulated_mol.GetMol()
        mol.UpdatePropertyCache() # getNumImplicitHs() called without preceding call to calcImplicitValence()
        G = mol_to_nx(mol)

        return mol, G
    else:
        for bond in mol.GetBonds():
            bond.SetBoolProp(KEY, False)
        G = mol_to_nx(mol)
        return mol, G

def get_fragments(mol, atoms):
    try: return Chem.MolFragmentToSmiles(mol, atoms, kekuleSmiles=True)
    except: return Chem.MolFragmentToSmiles(mol, atoms)

def get_smarts_fragments(mol, atoms):
    return Chem.MolFragmentToSmarts(mol, atoms)

def get_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: 
        return None
    Chem.Kekulize(mol)
    return mol

def get_clique_mol(mol, clique):
    new_mol = Chem.RWMol()

    node_to_idx = {}

    for idx in clique:
        symbol=mol.GetAtomWithIdx(idx).GetSymbol()
        formal_charge=mol.GetAtomWithIdx(idx).GetFormalCharge()
        chiral_tag=mol.GetAtomWithIdx(idx).GetChiralTag()
        hybridization=mol.GetAtomWithIdx(idx).GetHybridization()
        num_explicit_hs=mol.GetAtomWithIdx(idx).GetNumExplicitHs()
        is_aromatic=mol.GetAtomWithIdx(idx).GetIsAromatic()

        a=Chem.Atom(symbol)
        a.SetChiralTag(chiral_tag)
        a.SetFormalCharge(formal_charge)
        a.SetIsAromatic(is_aromatic)
        a.SetHybridization(hybridization)
        a.SetNumExplicitHs(num_explicit_hs)
        out_idx = new_mol.AddAtom(a)
        node_to_idx[idx] = out_idx

    subset_possible_bonds = list(itertools.combinations(clique, 2))
    for sub in subset_possible_bonds:
        bond = mol.GetBondBetweenAtoms(*sub)
        if bond:
            ifirst = node_to_idx[sub[0]]
            isecond = node_to_idx[sub[1]]
            bond_type = bond.GetBondType()
            new_mol.AddBond(ifirst, isecond, bond_type)

    new_mol2 = new_mol.GetMol()
    new_mol2.UpdatePropertyCache()

    return new_mol2