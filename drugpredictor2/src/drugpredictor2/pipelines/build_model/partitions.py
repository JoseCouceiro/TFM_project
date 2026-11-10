"""Row accounting for the sharded model input, without Kedro or TensorFlow.

`featurized_data` is a ``PartitionedDataSet`` over ~60 parquet shards, each
holding ~10 000 molecules and a 3 004-float fingerprint (roughly 115 MB of
column data). Training only needs two things from a shard *before* the first
epoch:

* how many rows it holds, to map a global row index back to a shard;
* the label column, and whether the fingerprint / auxiliary columns are null.

Both used to be answered by materialising the whole shard, so the node walked
~7 GB of fingerprint arrays twice before a single epoch ran — the main reason
runs kept being OOM-killed. Both answers are available far more cheaply:

* row count  -> ``ParquetFile.metadata.num_rows``  (footer only)
* null count -> per-row-group column statistics    (footer only)

Every helper here degrades to calling the loader when a partition is not a
readable parquet file (CSV shards, in-memory dicts in tests, odd files), so the
behaviour is unchanged — just faster.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

LoaderMap = Dict[str, Callable[[], pd.DataFrame]]


def resolve_partition_path(name: str, directory: Optional[str]) -> Optional[Path]:
    """Locate the file that backs a partition key, or ``None`` if unknown.

    Partition keys carry the file suffix only when the catalog entry does not
    declare a ``filename_suffix``: with ``filename_suffix: ".parquet"`` Kedro
    strips it, so the key is ``Balanced_000000001_000010000``; without it the
    key keeps it, ``Balanced_000000001_000010000.parquet``. Both spellings are
    accepted here so the helpers work either way.
    """
    if not directory:
        return None
    base = Path(directory)
    for candidate in (name, f"{name}.parquet"):
        path = base / candidate
        if path.is_file():
            return path
    return None


def parquet_row_count(path) -> Optional[int]:
    """Rows in a parquet file, read from the footer (no column data)."""
    try:
        import pyarrow.parquet as pq

        return pq.ParquetFile(path).metadata.num_rows
    except Exception:
        return None


def _leaf_column_index(metadata, column: str) -> Optional[int]:
    """Index of ``column`` among the flattened leaf columns of the parquet schema.

    ``metadata.schema.names`` is not usable here: a ``LIST<double>`` written from
    a pandas column called ``FP`` shows up as ``item`` in that view (the element
    field), while the row-group column at that position is the one carrying the
    statistics. Leaf paths (``FP.list.item``) are matched by prefix instead.
    """
    try:
        row_group = metadata.row_group(0)
    except Exception:
        return None
    for index in range(row_group.num_columns):
        path = row_group.column(index).path_in_schema
        if path == column or path.startswith(f"{column}."):
            return index
    return None


def column_null_count(path, column: str) -> Optional[int]:
    """Nulls in a scalar ``column`` over the whole file, from footer statistics.

    ``None`` when pyarrow is unavailable, the file cannot be read, the column is
    missing, or any row group lacks statistics — the caller then falls back to
    reading the column.
    """
    try:
        import pyarrow.parquet as pq

        metadata = pq.ParquetFile(path).metadata
        index = _leaf_column_index(metadata, column)
        if index is None:
            return None
        total = 0
        for group in range(metadata.num_row_groups):
            stats = metadata.row_group(group).column(index).statistics
            if stats is None or stats.null_count is None:
                return None
            total += stats.null_count
        return total
    except Exception:
        return None


def fingerprint_column_is_complete(path, column: str, rows: int) -> bool:
    """True when the footer proves every row has a non-empty, non-null fingerprint.

    Uses leaf statistics only: the null count must be zero, and the number of
    stored values must divide evenly by the row count — a row whose list value is
    NULL stores no values at all, so it shows up as a non-divisible count (a
    ``NULL`` list is the one case plain null counts cannot see).
    """
    if rows <= 0:
        return False
    try:
        import pyarrow.parquet as pq

        metadata = pq.ParquetFile(path).metadata
        index = _leaf_column_index(metadata, column)
        if index is None:
            return False
        nulls = 0
        values = 0
        for group in range(metadata.num_row_groups):
            leaf = metadata.row_group(group).column(index)
            stats = leaf.statistics
            if stats is None or stats.null_count is None:
                return False
            nulls += stats.null_count
            values += leaf.num_values
        return nulls == 0 and values > 0 and values % rows == 0
    except Exception:
        return False


def read_label_columns(
    partition_loaders: LoaderMap,
    columns: Sequence[str],
    data_path: Optional[str] = None,
) -> pd.DataFrame:
    """Read only the label columns of every partition, in loader order.

    Used to fit the target scaling, so it deliberately never touches the
    fingerprint column: all nine label columns for ~584 k rows come to ~40 MB,
    against ~2 GB if the fingerprints were included.
    """
    wanted = list(columns)
    frames = []
    for name, load in partition_loaders.items():
        path = resolve_partition_path(name, data_path)
        frames.append(
            pd.read_parquet(path, columns=wanted) if path is not None else load()[wanted]
        )
    if not frames:
        return pd.DataFrame(columns=wanted)
    return pd.concat(frames, ignore_index=True)


def partition_row_counts(
    partition_loaders: LoaderMap, data_path: Optional[str] = None
) -> Dict[str, int]:
    """Row count of every partition, keyed by partition name, in loader order."""
    counts: Dict[str, int] = {}
    for name, load in partition_loaders.items():
        count = None
        path = resolve_partition_path(name, data_path)
        if path is not None:
            count = parquet_row_count(path)
        if count is None:
            count = len(load())
        counts[name] = count
    return counts


def partition_spans(
    partition_loaders: LoaderMap, data_path: Optional[str] = None
) -> List[Tuple[str, int, int]]:
    """``[(name, start, end), ...]`` global row span of every partition.

    Spans are half-open and cumulative in loader order, so a global row index
    can be resolved with ``start <= i < end``. Costs one footer read per shard
    instead of a full load.
    """
    spans, total = [], 0
    for name, count in partition_row_counts(partition_loaders, data_path).items():
        spans.append((name, total, total + count))
        total += count
    return spans


def labels_and_keep(
    name: str,
    load: Callable[[], pd.DataFrame],
    params: dict,
    data_path: Optional[str] = None,
    extra_label_cols: Sequence[str] = (),
) -> Tuple[pd.Series, Optional[pd.Series]]:
    """Primary labels for one partition, plus the rows that survive preprocessing.

    Returns ``(labels, keep)``. ``keep`` is ``None`` only on the cheap path,
    when the footer proves no row can be dropped — the caller can then treat it
    as "keep everything". The fallback always computes the mask, so the keep/drop
    rule is identical to ``load_and_preprocess_fingerprints``.
    """
    x_col = params["X_column"]
    y_col = params["label"]
    columns = list(dict.fromkeys([y_col, *extra_label_cols]))

    path = resolve_partition_path(name, data_path)
    if path is not None:
        rows = parquet_row_count(path) or 0
        if fingerprint_column_is_complete(path, x_col, rows):
            if all(column_null_count(path, col) == 0 for col in columns):
                frame = pd.read_parquet(path, columns=columns)
                return frame[y_col], None

    df = load()
    keep = df[x_col].notna()
    for col in extra_label_cols:
        keep = keep & df[col].notna()
    return df[y_col], keep


def valid_rows(
    partition_loaders: LoaderMap,
    params: dict,
    data_path: Optional[str] = None,
    extra_label_cols: Sequence[str] = (),
) -> Tuple[np.ndarray, pd.Series]:
    """``(valid_global_indices, valid_labels)`` over all partitions.

    A row is valid when its fingerprint — and every column in
    ``extra_label_cols`` — is present, matching the keep/drop rule applied by
    ``load_and_preprocess_fingerprints``. Indices are positions in the
    concatenation of all partitions, in loader order.
    """
    label_parts, keep_parts = [], []
    for name, load in partition_loaders.items():
        labels, keep = labels_and_keep(name, load, params, data_path, extra_label_cols)
        label_parts.append(labels.reset_index(drop=True))
        keep_parts.append(
            pd.Series(True, index=range(len(labels)))
            if keep is None
            else keep.reset_index(drop=True)
        )

    labels = (
        pd.concat(label_parts, ignore_index=True) if label_parts else pd.Series(dtype=float)
    )
    keep = pd.concat(keep_parts, ignore_index=True) if keep_parts else pd.Series(dtype=bool)

    valid_idx = np.flatnonzero(keep.to_numpy())
    valid_labels = labels.iloc[valid_idx].reset_index(drop=True)
    print(
        f"Scanned {len(partition_loaders)} partitions: {len(labels)} rows, "
        f"{len(valid_idx)} usable after keep/drop."
    )
    return valid_idx, valid_labels


def stack_batch(df: pd.DataFrame, x_col: str) -> np.ndarray:
    """``(n, n_features, 1)`` float32 batch from a frame of fingerprint arrays.

    Stacking only the rows of the batch at hand — instead of every row of the
    partition, which is what the generators used to do — keeps the largest
    transient allocation at a few MB instead of ~170 MB per partition.
    """
    X = np.stack(df[x_col].to_numpy()).astype(np.float32)
    return X.reshape((X.shape[0], X.shape[1], 1))


def predict_stream(model, gen, steps: int) -> Tuple[np.ndarray, np.ndarray]:
    """``(y_true, y_pred)`` for ``steps`` batches, without materialising the split.

    Concatenating every validation/test batch into one array (58 451 x 5 288
    float32 = 1.2 GB, plus the ``np.concatenate`` copy) is what pushed the
    process over its memory limit at the end of training. Predictions are tiny,
    so only those are kept; rows stay in generator order, which is what keeps
    ``y_true`` and ``y_pred`` aligned.
    """
    y_true, y_pred = [], []
    for _ in range(steps):
        try:
            X, y = next(gen)
        except StopIteration:
            break
        y_true.append(np.asarray(y, dtype=np.float32).ravel())
        y_pred.append(np.asarray(model.predict_on_batch(X)).ravel())

    if not y_true:
        return np.empty(0, dtype=np.float32), np.empty(0, dtype=np.float32)
    return np.concatenate(y_true), np.concatenate(y_pred)


def predict_stream_multitask(
    model, gen, steps: int, keys: Sequence[str]
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Multi-head version of :func:`predict_stream`.

    ``keys`` are the model's output names (the generator yields a dict keyed the
    same way). The model returns one array per head, in output order; both that
    shape and a dict are normalised to ``{key: values}``.
    """
    y_true = {key: [] for key in keys}
    y_pred = {key: [] for key in keys}

    for _ in range(steps):
        try:
            X, y = next(gen)
        except StopIteration:
            break
        predictions = model.predict_on_batch(X)
        if not isinstance(predictions, dict):
            predictions = {
                key: predictions[index] for index, key in enumerate(keys)
            }
        for key in keys:
            y_true[key].append(np.asarray(y[key], dtype=np.float32).ravel())
            y_pred[key].append(np.asarray(predictions[key]).ravel())

    empty = np.empty(0, dtype=np.float32)
    return (
        {key: (np.concatenate(v) if v else empty) for key, v in y_true.items()},
        {key: (np.concatenate(v) if v else empty) for key, v in y_pred.items()},
    )
