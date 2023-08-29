"""
This is a boilerplate pipeline 'build_model'
generated using Kedro 0.18.10
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import determine_best_fp


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
    [
      node(
        func=determine_best_fp,
        inputs=["fingerprint_table"],
        outputs='fingerprints_accuracies',
        name="determine_best_fingerprints",
      )
    ]
  )
