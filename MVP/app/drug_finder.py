import os
import numpy as np
import json
import pubchempy as pcp
from rdkit import Chem
from tensorflow.keras.models import load_model
from compute_fp import Compute_FP as CFP

query=input('Insert CID of molecule: ')

with open(os.path.join('..','res','config', 'class_codes.json'), 'r') as file:
    class_codes_dict = json.load(file)

CNN_model = load_model(os.path.join('..','compiled_models','checkpoints', '04-0.300.hdf5'))

   
def get_molecule(cid):
    """
    Function that obtains an RDKit molecule object from a Pubchem "CID" number.
    Input: Pubchem CID as an integer
    Output: RDKit molecule
    """
    comp_mol = pcp.Compound.from_cid(cid)
    smiles_mol = comp_mol.canonical_smiles
    mol = Chem.MolFromSmiles(smiles_mol)
    return mol

def get_fp(mol):
    """
    Function that obtains the fingerprints from an RDKit molecule and returns them as a Numpy array.
    Input: RDKit Molecule
    Output: Numpy array
    """
    mc_fp = CFP.compute_morgan_fp(mol)

    return mc_fp

def reshape_fp_array(fp_chain):
    """
    Function that reshapes a Numpy array to pass it into a 1D CNN model.
    Input: Numpy array with shape (n,)
    Output: Numpy array with shape (1, n, 1)
    """
    reshaped = np.array(fp_chain.reshape(1, fp_chain.shape[0], 1))
    return reshaped

def return_predictions(arr, query):
    """
    Function that takes an array and predicts its label using a CNN model.
    Input: Numpy array
    Output: Prints predictions of which class a molecule is assigned to and the probability with which it is assigned to that class.
    """
    y_pred = CNN_model.predict(arr)
    print(f"The compound with CID '{query}' is predicted as {class_codes_dict[str(np.argmax(y_pred))]}")
    print(f"with probability: {y_pred.max():.2f}")

def interface(query):
    """
    Checks the input and runs all the functions in order to obtain a prediction from a CID identifier.
    Input: user input
    Output: prints prediction
    """
    input_is_ok = False
    while not input_is_ok:
        try:
            mol = get_molecule(query)
            input_is_ok = True
            fp_chain = get_fp(mol)
            arr = reshape_fp_array(fp_chain)
            return_predictions(arr, query)
        except:
            query = input('Please, enter a valid CID: ')

interface(query)

