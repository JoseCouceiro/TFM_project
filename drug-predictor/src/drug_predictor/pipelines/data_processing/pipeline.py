from kedro.pipeline import Pipeline, node, pipeline

from .nodes import process_gitter


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
            )
        ]
    )
