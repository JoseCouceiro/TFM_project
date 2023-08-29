
from kedro.pipeline import Pipeline, node, pipeline
from .nodes import get_model_input

def create_pipeline(**kwargs) -> Pipeline:
  return pipeline(
    [
      node(
        func=get_model_input,
        inputs=["all_drugs_table", 'params:build_molecule_column', "params:table_features"],
        outputs='fingerprint_table',
        name="get_model_input_node",
      )
    ]
  )