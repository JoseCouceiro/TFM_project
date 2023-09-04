from kedro.pipeline import Pipeline, node, pipeline

from .nodes import process_gitter, preprocess_pubchem, \
                   process_pubchem, preprocess_drugbank, \
                   process_drugbank, join_datasets

def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=process_gitter,
                inputs=["gitter",
                        'params:gitter',
                        'params:matc_gitter_conversion',
                        'params:matc_codes_explanation'],
                outputs="processed_gitter",
                name="process_gitter_node",
            ),
            node(
                func=preprocess_pubchem,
                inputs=["pubchem",
                        'params:pubchem'],
                outputs="preprocessed_pubchem",
                name="preprocess_pubchem_node",
            ),  
            node(
                func=process_pubchem,
                inputs=["preprocessed_pubchem",
                        'params:pubchem',
                        'params:matc_codes_explanation'],
                outputs="processed_pubchem",
                name="process_pubchem_node",
            ),
            node(
                func=preprocess_drugbank,
                inputs=["drugbank"],
                outputs='preprocessed_drugbank',
                name='preprocess_drugbank_node',
            ),
            node(
                func=process_drugbank,
                inputs=['preprocessed_drugbank',
                 'preprocessed_pubchem',
                 'params:drugbank',
                 'params:matc_codes_explanation'],
                outputs='processed_drugbank',
                name='process_drugbank_node'
            ), 
            node(
                func=join_datasets,
                inputs=['processed_pubchem',
                        'processed_drugbank',
                        'processed_gitter'],
                outputs='all_drugs_table',
                name='join_datasets_node'
            )            
        ]
    )
