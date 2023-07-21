from kedro.pipeline import Pipeline, node, pipeline

from .nodes import process_gitter, process_pubchem


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
                func=process_pubchem,
                inputs=["pubchem",
                        'params:pubchem',
                        'params:matc_codes_explanation'],
                outputs="processed_pubchem",
                name="process_pubchem_node",
            )
        ]
    )
