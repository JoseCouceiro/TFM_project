from kedro.pipeline import Pipeline, node, pipeline
from .nodes import visualize_training, evaluate_model


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
    [
      node(
        func=evaluate_model,
        inputs=['def_model', 'model_data', 'split_col', 
                'validation_dataset', 'selected_fingerprint'],
        outputs=['train_predictions', 'train_classification_report',
                 'validation_predictions', 'validation_classification_report'],
        name='evaluation_node'
      ),
      node(
        func=visualize_training,
        inputs='history',
        outputs='training_fig',
        name='visualization_node'
      )
    ]
    )
