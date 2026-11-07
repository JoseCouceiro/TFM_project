# Refined and correct pipeline definition
from kedro.pipeline import Pipeline, node
from typing import Dict

# Only the main orchestration function needs to be imported here
from .nodes import train_model_on_partitions, train_multitask_model_on_partitions

def create_pipeline(**kwargs) -> Pipeline:
    return Pipeline(
        [
            node(
                func=train_model_on_partitions,
                inputs=[
                    "featurized_data",
                    "params:train_params",
                    "params:split_params",
                    # Directory behind `featurized_data`: lets the node read row
                    # counts and label columns from the parquet footers instead of
                    # materialising every shard (a full extra pass over ~7 GB).
                    "params:featurized_data_path",
                ],
                outputs=[
                    "lipinski_model",
                    "lipinski_training_history",
                    "lipinski_train_predictions",
                    "lipinski_train_report",
                    "lipinski_val_predictions",
                    "lipinski_val_report"
                ],
                name="train_model_on_partitions_node"
            ),
            node(
                func=train_multitask_model_on_partitions,
                inputs=[
                    "featurized_data",
                    "params:train_params_multitask",
                    "params:split_params_multitask",
                    "params:featurized_data_path",
                ],
                outputs=[
                    "lipinski_multitask_model",
                    "lipinski_multitask_training_history",
                    "lipinski_multitask_train_predictions",
                    "lipinski_multitask_train_report",
                    "lipinski_multitask_val_predictions",
                    "lipinski_multitask_val_report"
                ],
                name="train_multitask_model_on_partitions_node"
            )
        ]
    )