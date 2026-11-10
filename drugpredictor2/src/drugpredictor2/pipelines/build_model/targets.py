"""Target scaling for the multi-task heads.

The raw targets have two pathologies that both showed up in the 07/11 run:

* **Mixed scales.** ``Molecular Weight`` reaches ~2 500 while ``FractionCSP3``
  lives inside [0, 1]. Summing one MSE per head therefore lets the widest-range
  head dominate: ``NumRotatableBonds`` alone owned 31 % of the composite
  ``val_loss`` that decided when training stopped.
* **Right tails.** ``NumRotatableBonds`` reaches 240; the worst 1 % of rows
  carried 47 % of that head's squared error, so a handful of polymer-sized
  molecules set the gradient.

:class:`TargetTransform` maps every head to roughly zero mean / unit variance
(optionally after a ``log`` or ``log1p`` compression) before the loss, and maps
predictions back to the original units for reporting. Binary heads are passed
through untouched. Statistics are fitted on the **training rows only**.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd

#: Supported per-head transforms.
KINDS = ("none", "standard", "log_standard", "log1p_standard")

_LOG_FLOOR = 1e-6


@dataclass
class TargetTransform:
    """Per-head forward/inverse scaling.

    Attributes:
        kinds: ``{column: kind}``, one of :data:`KINDS`.
        stats: ``{column: (mean, std)}`` in the *transformed* space, for the
            scaled kinds only.
    """

    kinds: Dict[str, str] = field(default_factory=dict)
    stats: Dict[str, Tuple[float, float]] = field(default_factory=dict)

    def kind_of(self, column: str) -> str:
        """Transform used for ``column`` (``"none"`` when unlisted)."""
        return self.kinds.get(column, "none")

    def forward(self, column: str, values) -> np.ndarray:
        """Raw target values -> model space."""
        kind = self.kind_of(column)
        array = np.asarray(values, dtype=np.float32)
        if kind == "none":
            return array
        if kind == "log_standard":
            array = np.log(np.maximum(array, _LOG_FLOOR))
        elif kind == "log1p_standard":
            array = np.log1p(np.maximum(array, 0.0))
        mean, std = self.stats.get(column, (0.0, 1.0))
        return ((array - mean) / std).astype(np.float32)

    def inverse(self, column: str, values) -> np.ndarray:
        """Model space -> raw target units."""
        kind = self.kind_of(column)
        array = np.asarray(values, dtype=np.float64)
        if kind == "none":
            return array.astype(np.float32)
        mean, std = self.stats.get(column, (0.0, 1.0))
        array = array * std + mean
        if kind == "log_standard":
            array = np.exp(array)
        elif kind == "log1p_standard":
            array = np.maximum(np.expm1(array), 0.0)
        return array.astype(np.float32)

    def describe(self) -> str:
        """One line per scaled head, for the run log."""
        lines = ["target scaling (fitted on the training rows):"]
        for column, kind in self.kinds.items():
            if kind == "none":
                lines.append(f"  {column:22s} {kind}")
            else:
                mean, std = self.stats[column]
                lines.append(f"  {column:22s} {kind:17s} mean={mean:+.4f} std={std:.4f}")
        return "\n".join(lines)


def build_output_specs(
    labels: Sequence[str],
    label_types: Mapping[str, str],
    loss_weights: Mapping[str, float] | None = None,
    losses: Mapping[str, str] | None = None,
    huber_delta: float = 2.0,
) -> list:
    """Per-head specs for :func:`build_multitask_model`.

    Kept here (not in ``nodes.py``) so the mapping from configuration to model
    heads can be tested without TensorFlow.

    Args:
        labels: Label columns, in output order.
        label_types: ``{column: "binary" | "regression"}``.
        loss_weights: Optional ``{column: weight}``; default 1.0.
        losses: Optional ``{column: "mse" | "mae" | "huber"}`` for regression
            heads. Huber is the bounded-gradient alternative to ``log1p`` for
            heavy-tailed targets: it keeps the fit in the original units (no
            inverse-transform bias) while stopping a handful of outliers from
            setting the gradient.
        huber_delta: Delta used by the Huber heads, in standardised units.

    Returns:
        ``[{"name", "type", "loss_weight", "loss", "huber_delta"}, ...]``.

    Raises:
        ValueError: an unknown label type or loss name.
    """
    weights = loss_weights or {}
    chosen = losses or {}
    specs = []
    for column in labels:
        kind = label_types[column]
        if kind not in ("binary", "regression"):
            raise ValueError(f"{column}: unknown label type {kind!r}")
        loss = chosen.get(column, "binary_crossentropy" if kind == "binary" else "mse")
        if loss not in ("binary_crossentropy", "mse", "mae", "huber"):
            raise ValueError(f"{column}: unknown loss {loss!r}")
        specs.append(
            {
                "name": column,
                "type": kind,
                "loss_weight": float(weights.get(column, 1.0)),
                "loss": loss,
                "huber_delta": float(huber_delta),
            }
        )
    return specs


def fit_target_transform(
    labels: pd.DataFrame,
    kinds: Mapping[str, str],
    rows: Sequence[int] | None = None,
) -> TargetTransform:
    """Fit the scaling statistics on the given rows (usually the training split).

    Args:
        labels: One column per head, in the same row order as the split indices.
        kinds: ``{column: kind}``; columns absent here are left unscaled.
        rows: Row positions to fit on. ``None`` uses every row.

    Returns:
        The fitted :class:`TargetTransform`.

    Raises:
        ValueError: an unknown transform kind, or a column missing from `labels`.
    """
    for column, kind in kinds.items():
        if kind not in KINDS:
            raise ValueError(f"{column}: unknown transform {kind!r}, expected one of {KINDS}")
        if column not in labels.columns:
            raise ValueError(f"{column}: not present in the label frame")

    subset = labels if rows is None else labels.iloc[list(rows)]
    transform = TargetTransform(kinds=dict(kinds))

    for column, kind in kinds.items():
        if kind == "none":
            continue
        values = subset[column].to_numpy(dtype=np.float64)
        if kind == "log_standard":
            values = np.log(np.maximum(values, _LOG_FLOOR))
        elif kind == "log1p_standard":
            values = np.log1p(np.maximum(values, 0.0))
        mean = float(np.mean(values))
        std = float(np.std(values))
        transform.stats[column] = (mean, std if std > 0 else 1.0)

    return transform
