"""Tests for ``build_model.targets`` (target scaling for the multi-task heads).

Run with ``pytest src/tests/pipelines/build_model/test_targets.py``.
"""

import numpy as np
import pandas as pd
import pytest

from drugpredictor2.pipelines.build_model.targets import (
    KINDS,
    TargetTransform,
    build_output_specs,
    fit_target_transform,
)


def _frame():
    """Two easy heads and one with the heavy tail that broke the RotB loss."""
    return pd.DataFrame(
        {
            "RuleFive": [0, 1, 0, 1, 0, 1],
            "Molecular Weight": [180.0, 320.0, 450.0, 700.0, 900.0, 2500.0],
            "NumRotatableBonds": [0.0, 1.0, 2.0, 5.0, 12.0, 240.0],
            "FractionCSP3": [0.1, 0.3, 0.5, 0.7, 0.9, 1.0],
        }
    )


def test_standard_transform_centres_and_scales():
    kinds = {"FractionCSP3": "standard"}
    transform = fit_target_transform(_frame(), kinds)
    scaled = transform.forward("FractionCSP3", _frame()["FractionCSP3"])

    assert abs(float(np.mean(scaled))) < 1e-5
    assert abs(float(np.std(scaled)) - 1.0) < 1e-5


def test_round_trip_returns_the_original_values_for_every_kind():
    kinds = {
        "FractionCSP3": "standard",
        "Molecular Weight": "log_standard",
        "NumRotatableBonds": "log1p_standard",
        "RuleFive": "none",
    }
    frame = _frame()
    transform = fit_target_transform(frame, kinds)

    for column, kind in kinds.items():
        original = frame[column].to_numpy(dtype=np.float32)
        restored = transform.inverse(column, transform.forward(column, original))
        assert np.allclose(restored, original, rtol=1e-4, atol=1e-4), (column, kind)


def test_log1p_standard_tames_the_heavy_tail():
    kinds = {"NumRotatableBonds": "log1p_standard"}
    transform = fit_target_transform(_frame(), kinds)
    scaled = transform.forward("NumRotatableBonds", _frame()["NumRotatableBonds"])

    # The 240-bond molecule must no longer dominate the squared error.
    assert scaled.max() < 3.0
    assert float(scaled[-1] ** 2) < 9.0


def test_binary_columns_pass_through_untouched():
    transform = fit_target_transform(_frame(), {"RuleFive": "none"})
    values = np.array([0, 1, 1, 0], dtype=np.float32)
    assert np.array_equal(transform.forward("RuleFive", values), values)
    assert np.array_equal(transform.inverse("RuleFive", values), values)


def test_statistics_are_fitted_on_the_requested_rows_only():
    frame = _frame()
    kinds = {"FractionCSP3": "standard"}

    everything = fit_target_transform(frame, kinds)
    first_three = fit_target_transform(frame, kinds, rows=[0, 1, 2])

    assert everything.stats["FractionCSP3"] != first_three.stats["FractionCSP3"]
    assert first_three.stats["FractionCSP3"][0] == pytest.approx(
        float(np.mean([0.1, 0.3, 0.5])), rel=1e-6
    )


def test_constant_column_does_not_divide_by_zero():
    frame = pd.DataFrame({"NumRotatableBonds": [3.0, 3.0, 3.0]})
    transform = fit_target_transform(frame, {"NumRotatableBonds": "log1p_standard"})

    assert transform.stats["NumRotatableBonds"][1] == 1.0
    scaled = transform.forward("NumRotatableBonds", frame["NumRotatableBonds"])
    assert np.all(scaled == 0.0)


def test_unknown_kind_and_missing_column_are_rejected():
    with pytest.raises(ValueError, match="unknown transform"):
        fit_target_transform(_frame(), {"QED": "zscore"})
    with pytest.raises(ValueError, match="not present"):
        fit_target_transform(_frame(), {"QED": "standard"})


def test_unlisted_columns_are_left_alone():
    transform = fit_target_transform(_frame(), {"FractionCSP3": "standard"})
    assert transform.kind_of("LogP") == "none"
    values = np.array([-1.5, 2.0], dtype=np.float32)
    assert np.array_equal(transform.forward("LogP", values), values)


def test_describe_lists_every_head():
    transform = fit_target_transform(
        _frame(), {"Molecular Weight": "log_standard", "RuleFive": "none"}
    )
    text = transform.describe()
    assert "Molecular Weight" in text and "log_standard" in text
    assert "RuleFive" in text


def test_all_documented_kinds_are_covered_by_the_constant():
    assert set(KINDS) == {"none", "standard", "log_standard", "log1p_standard"}
    assert TargetTransform().kind_of("anything") == "none"


LABEL_TYPES = {"RuleFive": "binary", "NumRotatableBonds": "regression", "QED": "regression"}


def test_output_specs_default_to_mse_and_bce():
    specs = build_output_specs(list(LABEL_TYPES), LABEL_TYPES)

    assert [spec["name"] for spec in specs] == list(LABEL_TYPES)
    assert specs[0]["loss"] == "binary_crossentropy"
    assert specs[1]["loss"] == "mse" and specs[2]["loss"] == "mse"
    assert all(spec["loss_weight"] == 1.0 for spec in specs)


def test_output_specs_apply_per_head_loss_weight_and_huber_delta():
    specs = build_output_specs(
        list(LABEL_TYPES),
        LABEL_TYPES,
        loss_weights={"NumRotatableBonds": 2.0},
        losses={"NumRotatableBonds": "huber"},
        huber_delta=3.0,
    )
    by_name = {spec["name"]: spec for spec in specs}

    assert by_name["NumRotatableBonds"]["loss"] == "huber"
    assert by_name["NumRotatableBonds"]["loss_weight"] == 2.0
    assert by_name["NumRotatableBonds"]["huber_delta"] == 3.0
    assert by_name["QED"]["loss"] == "mse"
    assert by_name["RuleFive"]["loss"] == "binary_crossentropy"


def test_output_specs_reject_unknown_loss_and_type():
    with pytest.raises(ValueError, match="unknown loss"):
        build_output_specs(["QED"], LABEL_TYPES, losses={"QED": "quantile"})
    with pytest.raises(ValueError, match="unknown label type"):
        build_output_specs(["QED"], {"QED": "ordinal"})
