"""
This is a boilerplate pipeline 'build_model'
generated using Kedro 0.18.10
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import determine_best_fp, prepare_data, \
                   obtain_model


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
    [ 
      node(
        func=determine_best_fp,
        inputs=["fingerprint_table"],
        outputs='fingerprints_accuracies',
        name="determine_best_fingerprints",
      ),
      node(
        func=prepare_data,
        inputs=['fingerprint_table', 'fingerprints_accuracies'],
        outputs=['split_col', 'model_data', 'selected_fingerprint'],
        name='prepare_data_node'      
      ),
      node(
        func=obtain_model,
        inputs=['split_col', 'model_data'],
        outputs=['def_model', 'history'],
        name='obtain_model_node'
      )
    ]
    )