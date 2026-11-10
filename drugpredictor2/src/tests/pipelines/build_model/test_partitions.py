"""Tests for ``drugpredictor2.pipelines.build_model.partitions``.

Those helpers are deliberately free of Kedro/TensorFlow imports so the row
accounting can be tested without the ML stack. Importing the package still runs
``build_model/__init__.py``, which imports the Kedro pipeline, so that module is
stubbed when Kedro is not installed (same idea as the process_data tests).

Run with ``pytest src/tests/pipelines/build_model/test_partitions.py``.
"""

import sys
import types

import numpy as np
import pandas as pd
import pytest

_STUB = "drugpredictor2.pipelines.build_model.pipeline"
if _STUB not in sys.modules:
    try:  # a real install has the genuine module; stub only when Kedro is absent
        import kedro  # noqa: F401
    except ImportError:
        _module = types.ModuleType(_STUB)
        _module.create_pipeline = lambda **kwargs: None
        sys.modules[_STUB] = _module

from drugpredictor2.pipelines.build_model.partitions import (  # noqa: E402
    column_null_count,
    fingerprint_column_is_complete,
    labels_and_keep,
    parquet_row_count,
    partition_row_counts,
    partition_spans,
    predict_stream,
    predict_stream_multitask,
    resolve_partition_path,
    stack_batch,
    valid_rows,
)

PARAMS = {"X_column": "FP", "label": "RuleFive"}
MULTITASK_PARAMS = {
    "X_column": "FP",
    "label": "RuleFive",
    "labels": ["RuleFive", "QED"],
}


def _frame(rows, fp_nulls=(), qed_nulls=()):
    """Stand-in for a featurized shard: FP is a fixed-length list column."""
    return pd.DataFrame(
        {
            "RuleFive": [float(index % 2) for index in range(rows)],
            "QED": [np.nan if index in qed_nulls else 0.5 for index in range(rows)],
            "FP": [None if index in fp_nulls else [0.0, 1.0] for index in range(rows)],
        }
    )


def _write_shard(directory, name, rows, **kwargs):
    path = directory / f"{name}.parquet"
    _frame(rows, **kwargs).to_parquet(path)
    return path


def _loaders(directory, names):
    """Loader functions shaped like the ones ``PartitionedDataSet.load()`` gives."""
    return {
        name: (lambda path=directory / f"{name}.parquet": pd.read_parquet(path))
        for name in names
    }


def test_resolve_partition_path_accepts_both_spellings(tmp_path):
    _write_shard(tmp_path, "Balanced_000000001_000000002", 2)
    assert resolve_partition_path("Balanced_000000001_000000002", str(tmp_path)) is not None
    assert resolve_partition_path("Balanced_000000001_000000002.parquet", str(tmp_path)) is not None
    assert resolve_partition_path("Balanced_000000009_000000010", str(tmp_path)) is None
    assert resolve_partition_path("Balanced_000000001_000000002", None) is None


def test_footer_counts(tmp_path):
    path = _write_shard(tmp_path, "shard", 5, fp_nulls=(1,), qed_nulls=(3,))
    assert parquet_row_count(path) == 5
    assert column_null_count(path, "RuleFive") == 0
    assert column_null_count(path, "QED") == 1
    assert column_null_count(path, "not_a_column") is None
    assert fingerprint_column_is_complete(path, "FP", 5) is False
    assert fingerprint_column_is_complete(_write_shard(tmp_path, "clean", 5), "FP", 5) is True
    assert fingerprint_column_is_complete(path, "FP", 0) is False


def test_counts_return_none_for_non_parquet(tmp_path):
    csv = tmp_path / "shard.csv"
    csv.write_text("RuleFive\n0\n1\n", encoding="utf-8")
    assert parquet_row_count(csv) is None
    assert column_null_count(csv, "RuleFive") is None
    assert fingerprint_column_is_complete(csv, "FP", 2) is False


def test_spans_come_from_footers_without_loading(tmp_path):
    _write_shard(tmp_path, "shard_a", 3)
    _write_shard(tmp_path, "shard_b", 2)

    def _explode():
        raise AssertionError("partition was loaded - the footer path was not used")

    loaders = {"shard_a": _explode, "shard_b": _explode}
    assert partition_row_counts(loaders, str(tmp_path)) == {"shard_a": 3, "shard_b": 2}
    assert partition_spans(loaders, str(tmp_path)) == [("shard_a", 0, 3), ("shard_b", 3, 5)]


def test_spans_fall_back_to_the_loader(tmp_path):
    _write_shard(tmp_path, "shard_a", 3)
    loaders = {
        "shard_a": lambda: _frame(3),
        "not_on_disk": lambda: _frame(4),
    }
    assert partition_row_counts(loaders, str(tmp_path)) == {"shard_a": 3, "not_on_disk": 4}
    assert partition_spans(loaders, str(tmp_path)) == [
        ("shard_a", 0, 3),
        ("not_on_disk", 3, 7),
    ]
    assert partition_spans(loaders) == [("shard_a", 0, 3), ("not_on_disk", 3, 7)]


def test_labels_and_keep_skips_the_mask_when_the_footer_proves_it_is_all_true(tmp_path):
    _write_shard(tmp_path, "shard", 4)
    loaders = _loaders(tmp_path, ["shard"])
    labels, keep = labels_and_keep("shard", loaders["shard"], PARAMS, str(tmp_path))
    assert keep is None
    assert labels.tolist() == [0.0, 1.0, 0.0, 1.0]


def test_labels_and_keep_builds_the_mask_when_a_fingerprint_is_missing(tmp_path):
    _write_shard(tmp_path, "shard", 4, fp_nulls=(2,))
    loaders = _loaders(tmp_path, ["shard"])
    labels, keep = labels_and_keep("shard", loaders["shard"], PARAMS, str(tmp_path))
    assert keep is not None
    assert labels.tolist() == [0.0, 1.0, 0.0, 1.0]
    assert keep.tolist() == [True, True, False, True]


def test_valid_rows_ignores_auxiliary_nulls_for_the_single_label_split(tmp_path):
    _write_shard(tmp_path, "shard", 4, qed_nulls=(1,))
    loaders = _loaders(tmp_path, ["shard"])
    valid_idx, _ = valid_rows(loaders, PARAMS, str(tmp_path))
    assert valid_idx.tolist() == [0, 1, 2, 3]
    multitask_idx, _ = valid_rows(
        loaders, MULTITASK_PARAMS, str(tmp_path), extra_label_cols=["RuleFive", "QED"]
    )
    assert multitask_idx.tolist() == [0, 2, 3]


@pytest.mark.parametrize("params", [PARAMS, MULTITASK_PARAMS])
def test_valid_rows_matches_a_full_load(tmp_path, params):
    """The footer path and the load-everything path must agree exactly."""
    _write_shard(tmp_path, "shard_a", 3, qed_nulls=(0,))
    _write_shard(tmp_path, "shard_b", 2, fp_nulls=(1,))
    loaders = _loaders(tmp_path, ["shard_a", "shard_b"])

    fast_idx, fast_labels = valid_rows(loaders, params, str(tmp_path))
    slow_idx, slow_labels = valid_rows(loaders, params, None)

    assert np.array_equal(fast_idx, slow_idx)
    assert fast_labels.tolist() == slow_labels.tolist()


class _FakeModel:
    """Stands in for a Keras model: predicts a fixed function of each batch."""

    def __init__(self, predict):
        self._predict = predict
        self.batches = 0

    def predict_on_batch(self, X):
        self.batches += 1
        return self._predict(X)


def test_stack_batch_shape_and_values():
    frame = pd.DataFrame({"FP": [[0.0, 1.0], [2.0, 3.0], [4.0, 5.0]]})
    batch = stack_batch(frame, "FP")
    assert batch.shape == (3, 2, 1)
    assert batch.dtype == np.float32
    assert batch[:, :, 0].tolist() == [[0.0, 1.0], [2.0, 3.0], [4.0, 5.0]]


def test_predict_stream_stops_after_steps_and_keeps_order():
    batches = [
        (np.zeros((2, 2, 1), dtype=np.float32), np.array([0.0, 1.0], dtype=np.float32)),
        (np.ones((1, 2, 1), dtype=np.float32), np.array([1.0], dtype=np.float32)),
    ]
    model = _FakeModel(lambda X: X[:, :1, 0] + 0.25)

    y_true, y_pred = predict_stream(model, iter(batches * 5), steps=2)

    assert y_true.tolist() == [0.0, 1.0, 1.0]
    assert y_pred.tolist() == [0.25, 0.25, 1.25]
    assert model.batches == 2, "predict_stream must not walk past `steps` batches"


def test_predict_stream_handles_an_exhausted_generator():
    y_true, y_pred = predict_stream(_FakeModel(lambda X: X), iter([]), steps=3)
    assert y_true.shape == (0,)
    assert y_pred.shape == (0,)


def test_predict_stream_multitask_normalises_a_list_output():
    X = np.zeros((2, 2, 1), dtype=np.float32)
    batch = (X, {"RuleFive": np.array([0.0, 1.0]), "QED": np.array([0.5, 0.6])})
    model = _FakeModel(
        lambda _: [np.array([[0.1], [0.9]]), np.array([[0.5], [0.7]])]
    )

    y_true, y_pred = predict_stream_multitask(
        model, iter([batch]), steps=1, keys=["RuleFive", "QED"]
    )

    assert y_true["RuleFive"].tolist() == [0.0, 1.0]
    assert np.allclose(y_pred["RuleFive"], [0.1, 0.9])
    assert np.allclose(y_true["QED"], [0.5, 0.6])  # cast to float32 on the way in
    assert np.allclose(y_pred["QED"], [0.5, 0.7])


def test_predict_stream_multitask_accepts_a_dict_output():
    X = np.zeros((1, 2, 1), dtype=np.float32)
    batch = (X, {"RuleFive": np.array([1.0])})
    model = _FakeModel(lambda _: {"RuleFive": np.array([[0.8]])})

    y_true, y_pred = predict_stream_multitask(
        model, iter([batch]), steps=1, keys=["RuleFive"]
    )

    assert y_true["RuleFive"].tolist() == [1.0]
    assert y_pred["RuleFive"].tolist() == [0.8]
