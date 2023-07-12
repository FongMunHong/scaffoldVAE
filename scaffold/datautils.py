from torch.utils.data import Dataset
from mol_tree import MolTree
import numpy as np


class MoleculeDataset(Dataset):

    def __init__(self, data_file):
        with open(data_file) as f:
            self.data = [line.strip("\r\n ").split()[0] for line in f]

    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        smiles = self.data[idx]
        mol_tree = MolTree(smiles)
        mol_tree.recover() # recover tree from nodes?
        mol_tree.assemble() # enum_assemble for all nodes
        return mol_tree