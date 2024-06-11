import numpy as np
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, rdMolDescriptors
from rdkit.Avalon import pyAvalonTools

class Compute_FP:
    """
    This class contains a series of functions to obtain some kind of molecular descriptor or fingerprint as a numpy array.
    Rdkit 'molecules' or Pubchem 'cid' have to be passed into the functions, as appropriate.
    The function 'relate_fp_functions' organizes the different functions in a dictionary and runs the adequate on depending on input.
    """

    def compute_morgan_fp(self, mol, depth=2, nBits=2048):
        try:
            mor_fp = AllChem.GetMorganFingerprintAsBitVect(mol,depth,nBits)
        except:
            print('Something went wrong computing Morgan fingerprints')
            return None
        return np.array(mor_fp)

    def compute_maccskeys(self, mol):
        try:
            mkeys = MACCSkeys.GenMACCSKeys(mol)   
        except:
            print('Something went wrong computing MACCSKeys')
            return None
        return np.array(mkeys)

    def compute_atom_pair_fp(self, mol, nBits=2048):
        try:
            atom_pair_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits)
        except:
            print('Something went wrong computing Atom Pair fingerprints')
            return None
        return np.array(atom_pair_fp)
    
    def compute_atom_pair_fp_trials(self, mol, nBits=2048):
        atom_pair_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits)
        return np.array(atom_pair_fp)

    def compute_topological_torsion_fp(self, mol, nBits=2048):
        try:
            tt_fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol)
        except:
            print('Something went wrong computing Topological Torsion fingerprints')
            return None
        return np.array(tt_fp)

    def compute_avalon_fp(self, mol, nBits=2048):
        try:
            av_fp = pyAvalonTools.GetAvalonFP(mol, nBits)
        except:
            print('Something went wrong computing Avalon fingerprints')
            return None
        return np.array(av_fp)
    
    def compute_morgan_circular_fp(self, mol, depth=2, nBits=2048):
        try:
            mc_fp = AllChem.GetMorganFingerprintAsBitVect(mol, depth, nBits)
        except:
            print('Something went wrong computing Morgan circular fingerprints')
            return None
        return np.array(mc_fp)

    def compute_rdkit_fp(self, mol, maxPath=5, fpSize=2048):
        try:
            rdkit_fp = AllChem.RDKFingerprint(mol, maxPath, fpSize)
        except:
            print('Something went wrong computing RDKit fingerprints')
            return None
        return np.array(rdkit_fp)
    
    def compute_pubchem_fingerprints(self, mol):
        try:
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            comp = pcp.get_compounds(smiles, 'smiles')
            fp_bin = bin(int(comp[0].fingerprint, 16))[2:]   
        except:
            print('Something went wrong computing Pubchem fingerprints')
            return None
        return np.array(list(fp_bin)).astype('int')
    
    def relate_fp_functions(self, fp, mol):
        """
        This function takes the name of a fingerprint and an RDKit molecule object and returns the corresponding fingerprints.
        Input: string, RDKit molecule.
        Output: numpy array.
        """
        dic ={
            "Morgan2FP": self.compute_morgan_fp,
            "MACCSKeys": self.compute_maccskeys,
            "AtomPairFP": self.compute_atom_pair_fp,
            "TopTorFP": self.compute_topological_torsion_fp,
            "AvalonFP": self.compute_avalon_fp,
            "PubchemFP": self.compute_pubchem_fingerprints
        }
        return dic[fp](mol)
