import streamlit as st
import os
import numpy as np
import json
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from tensorflow.keras.models import load_model
from compute_fp import Compute_FP as CFP

st.title('Drug predictor')

query=st.text_input("Insert molecule's CID or SMILES: ")
#query_smiles=[st.text_input('Insert SMILES'), 'SMILES']

with open(os.path.join('..','res','config', 'class_codes.json'), 'r') as file:
    class_codes_dict = json.load(file)

CNN_model = load_model(os.path.join('..','compiled_models','checkpoints', '13-0.330.hdf5'))


def get_molecule_from_cid(cid):
    """
    Function that obtains an RDKit molecule object from a Pubchem "CID" number.
    Input: Pubchem CID as an integer
    Output: RDKit molecule
    """
    comp_mol = pcp.Compound.from_cid(cid)
    smiles_mol = comp_mol.canonical_smiles
    mol = Chem.MolFromSmiles(smiles_mol)
    return mol

def get_molecule_from_smiles(smiles):
    """
    Function that obtains an RDKit molecule object from a SMILES.
    Input: SMILES as a string
    Output: RDKit molecule
    """
    mol = Chem.MolFromSmiles(smiles)
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
    st.markdown(f"##### The compound with {is_cid_or_or_smiles(query)} '{query}' is predicted as {class_codes_dict[str(np.argmax(y_pred))]} with probability: {y_pred.max():.2f}**")

def pred_from_cid(query):
    """
    Checks the input and runs all the functions in order to obtain a prediction from a CID identifier.
    Input: user input
    Output: prints prediction
    """
    input_is_ok = False
    while not input_is_ok:
        try:
            mol = get_molecule_from_cid(query)
            mol_image = Draw.MolToImage(mol)
            st.write("This is your compound's structure")
            st.image(mol_image)
            input_is_ok = True
            fp_chain = get_fp(mol)
            arr = reshape_fp_array(fp_chain)
            return_predictions(arr, query)
        except:
            st.write('Please, enter a valid CID')
            break

def pred_from_smiles(query):
    """
    Checks the input and runs all the functions in order to obtain a prediction from a CID identifier.
    Input: user input
    Output: prints prediction
    """
    input_is_ok = False
    while not input_is_ok:
        try:
            mol = get_molecule_from_smiles(query)
            mol_image = Draw.MolToImage(mol)
            st.write("This is your compound's structure")
            st.image(mol_image)
            input_is_ok = True
            fp_chain = get_fp(mol)
            arr = reshape_fp_array(fp_chain)
            return_predictions(arr, query)
        except:
            st.write('Please, enter a valid SMILES')
            break
            
def is_cid_or_or_smiles(user_input):
    try:
        pcp.Compound.from_cid(user_input)
        return 'CID'
    except:
        print('not a cid')
        pass
    try:
        Chem.MolFromSmiles(user_input)
        return 'SMILES'
    except:
        print('not a smiles')
        pass

query_type = is_cid_or_or_smiles(query)
if query_type == 'SMILES':
    pred_from_smiles(query)
if query_type == 'CID':
    pred_from_cid(query)



