from multiprocessing import Pool
import os
import pandas as pd
import pubchempy as pcp

def get_cid_from_inchi(inchi):
    """Function that obtains the CID of a molecule from its InChI
    Input: molecule's InChi
    Output: molecule's CID
    """
    try:
        comp = pcp.get_compounds(smiles, 'smiles')
    except:
        print('Something went wrong obtaining the CID')
        return None
    return comp[0].cid

drugs = pd.read_csv(os.path.join('Dataframes','Drugbank', 'drugbank_all_full_database.xml','recreation', 'drugbank_dataframe.csv'), low_memory=False)

sel_cols = ['H Bond Acceptor Count', 'H Bond Donor Count', 'Molecular Weight', 'logP', 'Rule of Five', 'SMILES', 'atc_code']
drugs_dataset = drugs[sel_cols]


if __name__ == '__main__':
    with Pool(8) as p:
        drugs_dataset['CID']= drugs['SMILES'].map(get_cid_from_inchi)

        drugs_dataset.to_csv('nuc_drugbank_dataset.csv', Index = False)

        
        


