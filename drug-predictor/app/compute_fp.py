import numpy as np
import pubchempy as pcp
from rdkit.Chem import AllChem, MACCSkeys, rdMolDescriptors
from rdkit.Avalon import pyAvalonTools

class Compute_FP:
    """
    This class contains a series of functions to obtain some kind of molecular descriptor or fingerprint as a numpy array.

    Rdkit 'molecules' or Pubchem 'cid' have to be passed into the functions, as appropriate.
    """

    def compute_connectivity_invariants(mol):
        try:
            con_inv_fp = rdMolDescriptors.GetConnectivityInvariants(mol)
        except:
            print('Something went wrong computing Feature Invariants')
            return None
        return np.array(con_inv_fp)

    def compute_feature_invariants(mol):
        try:
            inv_fp = rdMolDescriptors.GetFeatureInvariants(mol)
        except:
            print('Something went wrong computing Feature Invariants')
            return None
        return np.array(inv_fp)

    def compute_morgan_fp(mol, depth=2, nBits=2048):
        try:
            mor_fp = AllChem.GetMorganFingerprintAsBitVect(mol,depth,nBits)
        except:
            print('Something went wrong computing Morgan fingerprints')
            return None
        return np.array(mor_fp)

    def compute_maccskeys(mol):
        try:
            mkeys = MACCSkeys.GenMACCSKeys(mol)   
        except:
            print('Something went wrong computing MACCSKeys')
            return None
        return np.array(mkeys)

    def compute_atom_pair_fp(mol, nBits=2048):
        try:
            atom_pair_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits)
        except:
            print('Something went wrong computing Atom Pair fingerprints')
            return None
        return np.array(atom_pair_fp)

    def compute_topological_torsion_fp(mol, nBits=2048):
        try:
            tt_fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol)
        except:
            print('Something went wrong computing Topological Torsion fingerprints')
            return None
        return np.array(tt_fp)

    def compute_avalon_fp(mol, nBits=2048):
        try:
            av_fp = pyAvalonTools.GetAvalonFP(mol, nBits)
        except:
            print('Something went wrong computing Avalon fingerprints')
            return None
        return np.array(av_fp)
    
    def compute_morgan_circular_fp(mol, depth=2, nBits=2048):
        try:
            mc_fp = AllChem.GetMorganFingerprintAsBitVect(mol, depth, nBits)
        except:
            print('Something went wrong computing Morgan circular fingerprints')
            return None
        return np.array(mc_fp)

    def compute_rdkit_fp(mol, maxPath=5, fpSize=2048):
        try:
            rdkit_fp = AllChem.RDKFingerprint(mol, maxPath, fpSize)
        except:
            print('Something went wrong computing RDKit fingerprints')
            return None
        return np.array(rdkit_fp)

    def compute_pubchem_fingerprints(cid):
        try:
            comp = pcp.Compound.from_cid(int(cid))
            fp_bin = bin(int(comp.fingerprint, 16))[2:]   
        except:
            print('Something went wrong computing Pubchem fingerprints')
            return None
        return np.array(list(fp_bin)).astype('int')

    def compute_cactvs_fingerprints(cid):
        try:
            comp = pcp.Compound.from_cid(int(cid))
            cactvs_fp_bin = bin(int(comp.fingerprint, 16))[2:]
        except:
            print('Something went wrong computing Cactvs fingerprints')
            return None
        return np.array(list(cactvs_fp_bin)).astype('int')
    
    def relate_fp_functions(fp, mol):
        dic ={
            "FeatInvariants": Compute_FP.compute_feature_invariants,
            "ConnInvariants": Compute_FP.compute_connectivity_invariants,
            "Morgan2FP": Compute_FP.compute_morgan_fp,
            "MACCSKeys": Compute_FP.compute_maccskeys,
            "AtomPairFP": Compute_FP.compute_atom_pair_fp,
            "TopTorFP": Compute_FP.compute_topological_torsion_fp,
            "AvalonFP": Compute_FP.compute_avalon_fp,
            "PubchemFP": Compute_FP.compute_pubchem_fingerprints,
            "CactvsFP": Compute_FP.compute_cactvs_fingerprints
        }
        
        return dic[fp](mol)
