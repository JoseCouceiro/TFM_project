from multiprocessing import Pool
import requests as rq
import pandas as pd
import re


def get_atc_code(response):
    pat = r"\"[A-Z]\d{2}[A-Z]{2}\d{2}\""
    atc_code_found = re.search(pat, response.text)
    if atc_code_found:
        return atc_code_found.group(0).strip('"')
    else:
        return None

def get_name(response):
    return response.text.split('RecordTitle": ')[1].split('"')[1]
        
    """ except:
        print('There is some problem with the name')
        return None """
        
#df_keys = ['cid', 'name', 'atc_code']
parsing_errors = []

def get_atc(cid):
    df_dic = {}
    
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON'
    
    try:
        #response_json = rq.get(url).json()
        response_text = rq.get(url)
        print(response_text.text[0:200])
        #print(response_text.text)
        
        df_dic['cid'] = cid
        #df_dic['name'] = get_name(response_text)
        print(df_dic['name'])
        df_dic['atc_code'] = get_atc_code(response_text)

    except KeyError:
        print('There is some problem with the atc_code')
        parsing_errors.append(cid)
        df_dic['atc_code'] = None
    
    return df_dic

def mount_df(df_dic, df):
    df.loc[len(df)] = df_dic
    return df
     

if __name__ == '__main__':
    with Pool(8) as p:
        cids = 1983, 1983 #list(range(1983,1983))
        df_keys = ['cid', 'atc_code']
        df_atc = pd.DataFrame(columns = df_keys)
        df_no_atc = pd.DataFrame(columns = df_keys)
        
        add_dic = p.map(get_atc, cids)
        
        for dic in add_dic:
            if dic['atc_code']:
                df_atc = mount_df(dic, df_atc)
                df_atc.to_csv('pubchem_atc_1980_1985', index=False)
            else:
                df_no_atc = mount_df(dic, df_no_atc)
                df_no_atc.to_csv('pubchem_no_atc_1980_1985', index=False)

        print(df_atc)
        print(df_no_atc)
        print('parsing errors', parsing_errors)

        
        


