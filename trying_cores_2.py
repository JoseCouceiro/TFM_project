from multiprocessing import Pool
import requests as rq
import pandas as pd
import re

""" def f(x):
    return x*x

if __name__ == '__main__':
    with Pool(5) as p:
        print(p.map(f, [1, 2, 3]))
 """


def get_atc_code(response):
    pat = r"\"[A-Z]\d{2}[A-Z]{2}\d{2}\""
    atc_code_found = re.search(pat, response.text)
    if atc_code_found:
        return atc_code_found.group(0).strip('"')

def get_name(response):
    return response['Record']['RecordTitle']

df_keys = ['cid', 'name', 'atc_code']
parsing_errors = []



def get_atc(cid):
    df_dic = {}
    
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON'
    response_json = rq.get(url).json()
    response_text = rq.get(url)
    df_dic['cid'] = cid
    try:
        df_dic['atc_code'] = get_atc_code(response_text)
    except KeyError:
        print('There is some problem with the atc_code')
        parsing_errors.append(cid)
        df_dic['atc_code'] = None
    
    print(df_dic)
    return df_dic

def mount_df(df_dic, df):
    df.loc[len(df)] = df_dic
    return df

    
            

""" for cid in cids:
    get_data(cid)
print(df_atc)
print(df_no_atc)
df_no_atc.to_csv('foo') """

if __name__ == '__main__':
    with Pool(8) as p:
        cids = list(range(2001,2200))
        df_keys = ['cid', 'atc_code']
        df_atc = pd.DataFrame(columns = df_keys)
        df_no_atc = pd.DataFrame(columns = df_keys)
        
        add_dic = p.map(get_atc, cids)
        
        for dic in add_dic:
            if dic['atc_code']:
                print(dic)
                df_atc = mount_df(dic, df_atc)
            else:
                print(dic)
                df_no_atc = mount_df(dic, df_no_atc)
        print(df_atc)
        print(df_no_atc)
        print(parsing_errors)
        df_atc.to_csv('pubchem_atc_2001_3000')
        df_no_atc.to_csv('pubchem_no_atc_2001_3000')


