from multiprocessing import Pool
import os
import pandas as pd
import pubchempy as pcp

drugs = pd.read_csv(os.path.join('Dataframes','Pubmed', 'pubchem_atc_1981'), low_memory=False)

if __name__ == '__main__':
    with Pool(8) as p:
        drugs['compound'] = drugs['cid'].map(pcp.Compound.from_cid)

    drugs.to_csv('pubchem_mapped_cores.csv', index=False)

        
        


