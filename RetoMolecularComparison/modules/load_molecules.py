from rdkit import Chem
import os

def load_smi(path_to_smi):
    """
    Load a .smi file and return a list with the molecules.
    Args:
        - path_to_smi : path to the .smi file
    Returns:
        - molecules : list with the molecules
    """
    molecules = []
    supplier = Chem.SmilesMolSupplier(path_to_smi,titleLine=False)
    for mol in supplier:
        if mol is not None:
            molecules.append(mol)
    return molecules

class MolLoader:
    '''
    Class to load molecules in .dat format.

    Args:
        - target_QAOA_path : path to the target molecules to be used in the QAOA challenge
        - candidate_QAOA_path : path to the candidate molecules to be used in the QAOA challenge
        - target_DA_path : path to the target molecules to be used in the DA challenge
        - candidate_DA_path : path to the candidate molecules to be used in the DA challenge

    Returns:
        - target_mols_dict : dictionary with the target molecules (keys: molecule names, values: RDKit molecule objects)
        - candidate_molecules_dict : dictionary with the candidate molecules (keys: molecule names, values: RDKit molecule objects)

    Example:
        `mol_loader = MolLoader()` \\
        `target_mols_dict, candidate_molecules_dict = mol_loader.load_molecules_QAOA()`
    '''
    def __init__(self,target_QAOA_path='data/target_molecules/target_molecule_QAOA.dat',
                 candidate_QAOA_path='data/candidates/candidate_molecules_QAOA.dat',
                 target_DA_path='data/target_molecules/target_molecule_DA.dat',
                 candidate_DA_path='data/candidates/candidate_molecules_DA.dat') -> None:
        self.target_QAOA_path = target_QAOA_path
        self.candidate_QAOA_path = candidate_QAOA_path
        self.target_DA_path = target_DA_path
        self.candidate_DA_path = candidate_DA_path
        
        
    def load_molecules_QAOA(self):
        with open(self.target_QAOA_path, 'r') as file:
            target_molecules_lines = file.readlines()
        target_mols_dict = {}
        for line in target_molecules_lines:
            if line[0] == '#':
                continue
            line = line.split()
            target_mols_dict[line[0]] = Chem.MolFromSmiles(line[2])

        with open(self.candidate_QAOA_path, 'r') as file:
            candidate_molecules_lines = file.readlines()
        candidate_molecules_dict = {}
        for line in candidate_molecules_lines:
            if line[0] == '#':
                continue
            line = line.split()
            candidate_molecules_dict[line[0]] = Chem.MolFromSmiles(line[2])
        return target_mols_dict, candidate_molecules_dict
    
    def load_target_QAOA(self):
        with open(self.target_QAOA_path, 'r') as file:
            target_molecules_lines = file.readlines()
        for line in target_molecules_lines:
            if line[0] == '#':
                continue
            line = line.split()
            return line[0],Chem.MolFromSmiles(line[2])
    
    def load_molecules_DA(self):
        with open(self.target_DA_path, 'r') as file:
            target_molecules_lines = file.readlines()
        target_mols_dict = {}
        for line in target_molecules_lines:
            if line[0] == '#':
                continue
            line = line.split()
            target_mols_dict[line[0]] = Chem.MolFromSmiles(line[2])

        with open(self.candidate_DA_path, 'r') as file:
            candidate_molecules_lines = file.readlines()
        candidate_molecules_dict = {}
        for line in candidate_molecules_lines:
            if line[0] == '#':
                continue
            line = line.split()
            candidate_molecules_dict[line[0]] = Chem.MolFromSmiles(line[2])
        return target_mols_dict, candidate_molecules_dict
    
    def load_target_DA(self):
        with open(self.target_DA_path, 'r') as file:
            target_molecules_lines = file.readlines()
        for line in target_molecules_lines:
            if line[0] == '#':
                continue
            line = line.split()
            return line[0],Chem.MolFromSmiles(line[2])