import math
import gc
import os
import shutil
import tempfile
import numpy as np
import pandas as pd
from typing import Dict, Callable, Iterator, Optional, Tuple, List

from sklearn import metrics
from sklearn.model_selection import train_test_split

from tensorflow import keras
from keras import layers, models, optimizers, regularizers
from keras_tuner import HyperParameters

from .partitions import (
    partition_spans,
    predict_stream,
    predict_stream_multitask,
    stack_batch,
    valid_rows,
)


def _rss_mb():
    """Resident set size of this process in MB, or None when unavailable.

    Used to log memory per epoch: the run that trained fine for 20 epochs and
    was then OOM-killed shows whether RSS is creeping (a leak / allocator
    fragmentation) or jumping (a single oversized allocation).
    """
    try:
        with open("/proc/self/status", encoding="utf-8") as status:
            for line in status:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024
    except OSError:
        pass
    return None


class MemoryLog(keras.callbacks.Callback):
    """Print RSS at the end of every epoch (no-op where /proc is missing)."""

    def on_epoch_end(self, epoch, logs=None):
        rss = _rss_mb()
        if rss is not None:
            print(f"  [mem] epoch {epoch + 1}: RSS {rss:.0f} MB")


# --------------------------
# Preprocessing
# --------------------------

def load_and_preprocess_fingerprints(df: pd.DataFrame, params: dict) -> pd.DataFrame:
    """
    Keep rows with non-null fingerprints; cast to float32 for Keras Conv1D.
    """
    x_col = params["X_column"]
    # drop rows where the fingerprint is missing
    df = df.dropna(subset=[x_col]).copy()
    # ensure numpy arrays of dtype float32 (from bool/ints)
    df[x_col] = df[x_col].apply(lambda fp: np.asarray(fp, dtype=np.float32))
    return df


# --------------------------
# Index discovery & splits
# --------------------------

def _discover_valid_rows(
    partition_loaders: Dict[str, Callable[[], pd.DataFrame]],
    params: dict,
    data_path: Optional[str] = None,
):
    """
    Return (valid_global_indices_sorted, valid_labels_series) where 'valid'
    means rows that survive the same preprocessing keep/drop rule.

    Delegates to ``partitions.valid_rows``, which reads each shard's parquet
    footer (row count, null counts) instead of materialising ~115 MB of
    fingerprints per shard just to look at one label column.
    """
    return valid_rows(partition_loaders, params, data_path)


def create_train_val_test_indices(
    partitioned_model_input: Dict[str, Callable[[], pd.DataFrame]],
    params: dict,
    data_path: Optional[str] = None,
) -> Tuple[List[int], List[int], List[int]]:
    """
    Stratified split on rows that will actually be used.
    Returns sorted lists of global indices.
    """
    val_size = params["val_size"]
    test_size = params["test_size"]

    valid_idx, valid_labels = _discover_valid_rows(partitioned_model_input, params, data_path)
    print(valid_labels.value_counts())

    train_val_idx, test_idx = train_test_split(
        valid_idx,
        test_size=test_size,
        stratify=valid_labels,
        random_state=25,
    )

    # labels for remaining (train+val)
    pos_in_valid = np.searchsorted(valid_idx, np.sort(train_val_idx))
    remain_labels = valid_labels.iloc[pos_in_valid]

    rel_val = val_size / (1.0 - test_size)
    train_idx, val_idx = train_test_split(
        train_val_idx,
        test_size=rel_val,
        stratify=remain_labels,
        random_state=42,
    )

    train_idx = sorted(train_idx.tolist())
    val_idx   = sorted(val_idx.tolist())
    test_idx  = sorted(test_idx.tolist())

    print(f"Dataset split (valid rows only): Train={len(train_idx)}, Val={len(val_idx)}, Test={len(test_idx)}")
    return train_idx, val_idx, test_idx


# --------------------------
# Generators (deterministic)
# --------------------------

def _partition_spans(
    partition_loaders: Dict[str, Callable[[], pd.DataFrame]],
    data_path: Optional[str] = None,
):
    """Return list of (name, start, end) global spans per partition, in insertion order.

    Row counts come from the parquet footers (see ``partitions.partition_spans``)
    — this used to load every shard in full just to call ``len()`` on it.
    """
    return partition_spans(partition_loaders, data_path)


def infinite_train_generator(
    partition_loaders: Dict[str, Callable[[], pd.DataFrame]],
    ordered_indices: List[int],
    params: Dict[str, str],
    batch_size: int,
    shuffle_each_epoch: bool = False,
    seed: Optional[int] = None,
    data_path: Optional[str] = None,
) -> Iterator[Tuple[np.ndarray, np.ndarray]]:
    """
    Infinite generator for training. Deterministic unless shuffling is enabled.
    If shuffle_each_epoch=True and seed is None, uses a fresh RNG each epoch.
    """
    x_col, y_col = params["X_column"], params["label"]
    spans = _partition_spans(partition_loaders, data_path)
    base_order = np.array(ordered_indices, dtype=int)

    # Create RNG once (don’t re-seed inside the loop)
    rng = np.random.default_rng(seed) if shuffle_each_epoch else None

    while True:
        order = base_order.copy()
        if rng is not None:
            rng.shuffle(order)

        for (name, start, end) in spans:
            mask = (order >= start) & (order < end)
            part_global = order[mask]
            if part_global.size == 0:
                continue

            local = np.sort(part_global - start)
            df = partition_loaders[name]()
            df = load_and_preprocess_fingerprints(df.iloc[local], params)
            if df.empty:
                continue

            # Stack one batch at a time: building X for the whole partition
            # (~170 MB) on top of the partition itself is what turned a single
            # epoch into a multi-GB spike.
            for i in range(0, len(df), batch_size):
                chunk = df.iloc[i:i + batch_size]
                yield stack_batch(chunk, x_col), chunk[y_col].to_numpy().astype(np.float32)

            # Release this partition's ~200 MB+ DataFrame before loading the
            # next one — lazy GC across a long fit() loop is what let memory
            # climb until the process was OOM-killed.
            del df
            gc.collect()


def repeating_eval_generator(
    partition_loaders: Dict[str, Callable[[], pd.DataFrame]],
    ordered_indices: List[int],
    params: Dict[str, str],
    batch_size: int,
    data_path: Optional[str] = None,
) -> Iterator[Tuple[np.ndarray, np.ndarray]]:
    """
    Infinite (repeating) generator for validation/test. Deterministic order each cycle.
    Safe to pass to Keras with `validation_steps`.
    """
    x_col, y_col = params["X_column"], params["label"]
    spans = _partition_spans(partition_loaders, data_path)
    order = np.array(ordered_indices, dtype=int)

    while True:
        for (name, start, end) in spans:
            mask = (order >= start) & (order < end)
            part_global = order[mask]
            if part_global.size == 0:
                continue

            local = np.sort(part_global - start)
            df = partition_loaders[name]()
            df = load_and_preprocess_fingerprints(df.iloc[local], params)
            if df.empty:
                continue

            for i in range(0, len(df), batch_size):
                chunk = df.iloc[i:i + batch_size]
                yield stack_batch(chunk, x_col), chunk[y_col].to_numpy().astype(np.float32)

            del df
            gc.collect()


# --------------------------
# Model
# --------------------------

def build_def_model_dynamic_Dense(hp, input_shape):
    model = models.Sequential([
        layers.Input(shape=input_shape, dtype="float32"),
        layers.Flatten(),
        layers.Dense(
            hp.Int('dense_1_units', 256, 512, step=64),
            activation=None,
            kernel_initializer='he_uniform'
        ),
        layers.BatchNormalization(),
        layers.LeakyReLU(alpha=0.1),
        layers.Dropout(hp.Float('dropout_1', 0.2, 0.4, step=0.1)),

        layers.Dense(
            hp.Int('dense_2_units', 128, 256, step=64),
            activation='leaky_relu',
            kernel_initializer='he_uniform',
            kernel_regularizer=regularizers.l2(1e-5)
        ),
        layers.Dropout(0.5),

        layers.Dense(1, activation='sigmoid')
    ])
    
    model.compile(
        optimizer=optimizers.Adam(
            learning_rate=hp.Choice('learning_rate', [1e-3, 5e-4, 1e-4])
        ),
        loss='binary_crossentropy',
        metrics=['accuracy']
    )
    return model

# This is the model used in the final experiments
def build_def_model_dynamic(hp: HyperParameters, input_shape: tuple) -> keras.Sequential:
    model = models.Sequential([
        layers.Input(shape=input_shape, dtype="float32"),

        # Single Conv1D layer
        layers.Conv1D(
            filters=hp.Int('conv_filters', 32, 128, step=16),
            kernel_size=hp.Choice('conv_kernel', [3,5]),
            activation='relu',
            padding='valid'
        ),

        layers.Flatten(),

        # Dense layers
        layers.Dense(
            units=hp.Int('dense_1_units', 256, 512, step=64),
            activation='relu',
            kernel_initializer='he_uniform'
        ),
        layers.Dropout(hp.Float('dropout_1', 0.3, 0.5, step=0.05)),

        layers.Dense(
            units=hp.Int('dense_2_units', 128, 256, step=32),
            activation='relu',
            kernel_initializer='he_uniform'
        ),
        layers.Dropout(hp.Float('dropout_2', 0.3, 0.5, step=0.05)),

        layers.Dense(1, activation='sigmoid')
    ])

    model.compile(
        optimizer=optimizers.Adam(
            learning_rate=hp.Choice('learning_rate', [1e-3, 5e-4])
        ),
        loss='binary_crossentropy',
        metrics=['accuracy']
    )

    return model


# --------------------------
# Orchestration
# --------------------------

def train_model_on_partitions(
    partitioned_model_input: Dict[str, Callable[[], pd.DataFrame]],
    train_params: dict,
    split_params: dict,
    data_path: Optional[str] = None,
):
    """
    Train + evaluate on partitioned data with correct splits and stable ordering.
    Expects split_params: {X_column, label, val_size, test_size}
            train_params: {epochs (int), batch_size (optional)}
    data_path is the directory backing `featurized_data`: row counts and label
    columns are then read from the parquet files directly instead of loading
    every shard into memory just to measure it.
    """
    batch_size = train_params['batch_size']

    # 1) Proper split on valid rows only (after keep/drop)
    train_idx, val_idx, test_idx = create_train_val_test_indices(
        partitioned_model_input, split_params, data_path
    )

    # 2) Generators
    train_gen = infinite_train_generator(
    partitioned_model_input, train_idx, split_params, batch_size, shuffle_each_epoch=True,
    data_path=data_path,
    )
    val_gen = repeating_eval_generator(
    partitioned_model_input, val_idx, split_params, batch_size, data_path=data_path
    )

    # 3) Steps
    steps_per_epoch  = math.ceil(len(train_idx) / batch_size)
    validation_steps = math.ceil(len(val_idx)   / batch_size)

    # 4) Build model
    sample_X, _ = next(train_gen)
    in_shape = sample_X.shape[1:]
    hp = HyperParameters()
    model = build_def_model_dynamic(hp, in_shape)

    # 5) Train. The best epoch is checkpointed to DISK instead of asking
    # EarlyStopping for restore_best_weights=True: that flag keeps a second full
    # copy of the weights alive for the whole run, at exactly the moment the
    # training buffers and validation outputs are also resident. Same fix the
    # multi-task run already uses.
    checkpoint_dir = tempfile.mkdtemp(prefix="lipinski_ckpt_")
    checkpoint_path = os.path.join(checkpoint_dir, "best.weights.h5")

    history = model.fit(
        train_gen,
        epochs=train_params['epochs'],
        steps_per_epoch=math.ceil(len(train_idx) / batch_size),
        validation_data=val_gen,
        validation_steps=math.ceil(len(val_idx) / batch_size),
        callbacks=[
            keras.callbacks.EarlyStopping(monitor="val_loss", patience=10, restore_best_weights=False),
            keras.callbacks.ModelCheckpoint(
                filepath=checkpoint_path, monitor="val_loss",
                save_best_only=True, save_weights_only=True, verbose=0,
            ),
            MemoryLog(),
        ],
        verbose=1,
    )

    # Load the best epoch's weights back once training is over and the
    # generators/buffers competing for memory are gone.
    if os.path.exists(checkpoint_path):
        model.load_weights(checkpoint_path)
    shutil.rmtree(checkpoint_dir, ignore_errors=True)


    # 6) Evaluation, predicted batch by batch. Stacking every val/test row
    # (58 451 x 5 288 float32 = 1.2 GB, plus the concatenate copy) is what got
    # the run OOM-killed at the end of training; only predictions are kept now,
    # still in generator order so y_true and y_pred stay aligned.
    def _evaluate(indices: List[int], steps: int):
        gen = repeating_eval_generator(
            partitioned_model_input, indices, split_params, batch_size, data_path=data_path
        )
        y_true, y_score = predict_stream(model, gen, steps)
        return y_true, (y_score > 0.5).astype(int)

    # Val
    yv, yv_pred = _evaluate(val_idx, validation_steps)
    n = min(len(yv), len(yv_pred))
    val_predictions = pd.DataFrame({"y_true": yv[:n], "y_pred": yv_pred[:n]})
    val_report = metrics.classification_report(yv[:n], yv_pred[:n])

    # Test
    test_steps = math.ceil(len(test_idx) / batch_size)
    yt, yt_pred = _evaluate(test_idx, test_steps)
    n2 = min(len(yt), len(yt_pred))
    test_predictions = pd.DataFrame({"y_true": yt[:n2], "y_pred": yt_pred[:n2]})
    test_report = metrics.classification_report(yt[:n2], yt_pred[:n2])

    return (
        model,
        history.history,
        test_predictions,
        test_report,
        val_predictions,
        val_report,
    )


# ============================================================
# Multi-task pretraining: RuleFive + QED + structural descriptors
# ============================================================
# Additive alongside the single-label pipeline above (nothing existing is
# changed) so both backbones can be compared head-to-head.

def _discover_valid_rows_multitask(
    partition_loaders: Dict[str, Callable[[], pd.DataFrame]],
    params: dict,
    data_path: Optional[str] = None,
):
    """Like _discover_valid_rows, but a row is only 'valid' if the fingerprint
    AND every auxiliary label column are present. Stratification still uses
    the primary binary label (params['label']) since stratifying jointly on
    several continuous targets isn't meaningful."""
    return valid_rows(
        partition_loaders, params, data_path, extra_label_cols=params["labels"]
    )


def create_train_val_test_indices_multitask(
    partitioned_model_input: Dict[str, Callable[[], pd.DataFrame]],
    params: dict,
    data_path: Optional[str] = None,
) -> Tuple[List[int], List[int], List[int]]:
    """Same stratified-split logic as create_train_val_test_indices, but the
    keep-mask additionally requires every auxiliary label column to be present."""
    val_size = params["val_size"]
    test_size = params["test_size"]

    valid_idx, valid_labels = _discover_valid_rows_multitask(
        partitioned_model_input, params, data_path
    )
    print(valid_labels.value_counts())

    train_val_idx, test_idx = train_test_split(
        valid_idx, test_size=test_size, stratify=valid_labels, random_state=25,
    )

    pos_in_valid = np.searchsorted(valid_idx, np.sort(train_val_idx))
    remain_labels = valid_labels.iloc[pos_in_valid]

    rel_val = val_size / (1.0 - test_size)
    train_idx, val_idx = train_test_split(
        train_val_idx, test_size=rel_val, stratify=remain_labels, random_state=42,
    )

    train_idx = sorted(train_idx.tolist())
    val_idx = sorted(val_idx.tolist())
    test_idx = sorted(test_idx.tolist())

    print(f"Dataset split (valid rows only): Train={len(train_idx)}, Val={len(val_idx)}, Test={len(test_idx)}")
    return train_idx, val_idx, test_idx


def _build_multitask_batch(df: pd.DataFrame, x_col: str, label_cols: List[str]):
    X = stack_batch(df, x_col)
    y = {col: df[col].to_numpy().astype(np.float32) for col in label_cols}
    return X, y


def infinite_train_generator_multitask(
    partition_loaders: Dict[str, Callable[[], pd.DataFrame]],
    ordered_indices: List[int],
    params: Dict[str, str],
    batch_size: int,
    shuffle_each_epoch: bool = False,
    seed: Optional[int] = None,
    data_path: Optional[str] = None,
):
    """Same traversal logic as infinite_train_generator, but yields a dict of
    label arrays (one per auxiliary target) instead of a single array."""
    x_col, label_cols = params["X_column"], params["labels"]
    spans = _partition_spans(partition_loaders, data_path)
    base_order = np.array(ordered_indices, dtype=int)
    rng = np.random.default_rng(seed) if shuffle_each_epoch else None

    while True:
        order = base_order.copy()
        if rng is not None:
            rng.shuffle(order)

        for (name, start, end) in spans:
            mask = (order >= start) & (order < end)
            part_global = order[mask]
            if part_global.size == 0:
                continue

            local = np.sort(part_global - start)
            df = partition_loaders[name]()
            df = load_and_preprocess_fingerprints(df.iloc[local], params)
            if df.empty:
                continue

            for i in range(0, len(df), batch_size):
                chunk = df.iloc[i:i + batch_size]
                yield _build_multitask_batch(chunk, x_col, label_cols)

            # Release this partition's ~200MB+ DataFrame before loading the
            # next one — relying on lazy GC across a long fit() loop is what
            # caused memory to climb until the process was OOM-killed.
            del df
            gc.collect()


def repeating_eval_generator_multitask(
    partition_loaders: Dict[str, Callable[[], pd.DataFrame]],
    ordered_indices: List[int],
    params: Dict[str, str],
    batch_size: int,
    data_path: Optional[str] = None,
):
    x_col, label_cols = params["X_column"], params["labels"]
    spans = _partition_spans(partition_loaders, data_path)
    order = np.array(ordered_indices, dtype=int)

    while True:
        for (name, start, end) in spans:
            mask = (order >= start) & (order < end)
            part_global = order[mask]
            if part_global.size == 0:
                continue

            local = np.sort(part_global - start)
            df = partition_loaders[name]()
            df = load_and_preprocess_fingerprints(df.iloc[local], params)
            if df.empty:
                continue

            for i in range(0, len(df), batch_size):
                chunk = df.iloc[i:i + batch_size]
                yield _build_multitask_batch(chunk, x_col, label_cols)

            del df
            gc.collect()


def build_multitask_model(hp: HyperParameters, input_shape: tuple, output_specs: List[Dict]) -> keras.Model:
    """Shared Conv1D backbone with one head per auxiliary target.

    output_specs: list of {"name": str, "type": "binary" | "regression", "loss_weight": float}.
    Binary heads use sigmoid + BCE; regression heads use a linear output + MSE.
    The backbone (Input -> Conv1D -> Flatten) has the same shape as
    build_def_model_dynamic's, so `_clone_conv1d_backbone` in
    multitask_model/nodes.py keeps working unchanged — it only looks for the
    first Conv1D/Flatten layers, regardless of what follows.
    """
    inputs = layers.Input(shape=input_shape, dtype="float32", name="fp_input")

    x = layers.Conv1D(
        filters=hp.Int('conv_filters', 32, 128, step=16),
        kernel_size=hp.Choice('conv_kernel', [3, 5]),
        activation='relu',
        padding='valid',
    )(inputs)
    x = layers.Flatten()(x)

    x = layers.Dense(
        units=hp.Int('dense_1_units', 256, 512, step=64),
        activation='relu',
        kernel_initializer='he_uniform',
    )(x)
    x = layers.Dropout(hp.Float('dropout_1', 0.3, 0.5, step=0.05))(x)

    shared = layers.Dense(
        units=hp.Int('dense_2_units', 128, 256, step=32),
        activation='relu',
        kernel_initializer='he_uniform',
    )(x)
    shared = layers.Dropout(hp.Float('dropout_2', 0.3, 0.5, step=0.05))(shared)

    outputs, losses, loss_weights, metrics_dict = {}, {}, {}, {}
    for spec in output_specs:
        name = spec["name"]
        head = layers.Dense(32, activation='relu', name=f"{name}_head_dense")(shared)
        if spec["type"] == "binary":
            outputs[name] = layers.Dense(1, activation='sigmoid', name=name)(head)
            losses[name] = 'binary_crossentropy'
            metrics_dict[name] = ['accuracy']
        else:
            outputs[name] = layers.Dense(1, activation='linear', name=name)(head)
            losses[name] = 'mse'
            metrics_dict[name] = ['mae']
        loss_weights[name] = spec.get("loss_weight", 1.0)

    model = models.Model(inputs=inputs, outputs=outputs)
    model.compile(
        optimizer=optimizers.Adam(learning_rate=hp.Choice('learning_rate', [1e-3, 5e-4])),
        loss=losses,
        loss_weights=loss_weights,
        metrics=metrics_dict,
    )
    return model


def train_multitask_model_on_partitions(
    partitioned_model_input: Dict[str, Callable[[], pd.DataFrame]],
    train_params: dict,
    split_params: dict,
    data_path: Optional[str] = None,
):
    """Train the Lipinski backbone with several auxiliary heads at once
    (RuleFive + QED + NumRotatableBonds + NumAromaticRings + FractionCSP3)
    instead of a single binary label, so the backbone learns a richer
    structural representation to transfer to the downstream drug/ATC heads.

    Expects split_params: {X_column, label, labels, label_types, loss_weights (optional),
                            val_size, test_size}
            train_params: {epochs, batch_size}
    data_path: directory backing `featurized_data` (see train_model_on_partitions).
    """
    batch_size = train_params['batch_size']
    label_cols = split_params["labels"]
    label_types = split_params["label_types"]
    loss_weights = split_params.get("loss_weights", {})
    output_specs = [
        {"name": c, "type": label_types[c], "loss_weight": loss_weights.get(c, 1.0)}
        for c in label_cols
    ]

    train_idx, val_idx, test_idx = create_train_val_test_indices_multitask(
        partitioned_model_input, split_params, data_path
    )

    train_gen = infinite_train_generator_multitask(
        partitioned_model_input, train_idx, split_params, batch_size, shuffle_each_epoch=True,
        data_path=data_path,
    )
    val_gen = repeating_eval_generator_multitask(
        partitioned_model_input, val_idx, split_params, batch_size, data_path=data_path
    )

    steps_per_epoch = math.ceil(len(train_idx) / batch_size)
    validation_steps = math.ceil(len(val_idx) / batch_size)

    sample_X, _ = next(train_gen)
    in_shape = sample_X.shape[1:]
    hp = HyperParameters()
    model = build_multitask_model(hp, in_shape, output_specs)

    # Checkpoint the best epoch to DISK instead of EarlyStopping's
    # restore_best_weights=True, which keeps an extra full in-memory copy of
    # all weights alive for the whole run. That extra copy — created right at
    # on_epoch_end, exactly when training buffers/validation results are also
    # all resident at once — is what pushed memory over the edge and got the
    # process OOM-killed right after a fully successful first epoch.
    checkpoint_dir = tempfile.mkdtemp(prefix="lipinski_multitask_ckpt_")
    checkpoint_path = os.path.join(checkpoint_dir, "best.weights.h5")

    history = model.fit(
        train_gen,
        epochs=train_params['epochs'],
        steps_per_epoch=steps_per_epoch,
        validation_data=val_gen,
        validation_steps=validation_steps,
        callbacks=[
            keras.callbacks.EarlyStopping(monitor="val_loss", patience=10, restore_best_weights=False),
            keras.callbacks.ModelCheckpoint(
                filepath=checkpoint_path, monitor="val_loss",
                save_best_only=True, save_weights_only=True, verbose=0,
            ),
            MemoryLog(),
        ],
        verbose=1,
    )

    # Load the best epoch's weights back in one shot, now that training is
    # done and the generators/buffers that were competing for memory are gone.
    if os.path.exists(checkpoint_path):
        model.load_weights(checkpoint_path)
    shutil.rmtree(checkpoint_dir, ignore_errors=True)

    def _report(y_dict, preds):
        report_lines, pred_cols = [], {}
        for c in label_cols:
            n = min(len(y_dict[c]), len(preds[c]))
            y_true = y_dict[c][:n]
            y_pred_raw = np.asarray(preds[c])[:n].ravel()
            pred_cols[f"{c}_true"] = y_true
            if label_types[c] == "binary":
                y_pred = (y_pred_raw > 0.5).astype(int)
                pred_cols[f"{c}_pred"] = y_pred
                report_lines.append(f"--- {c} (binary) ---")
                report_lines.append(metrics.classification_report(y_true, y_pred))
            else:
                pred_cols[f"{c}_pred"] = y_pred_raw
                mae = metrics.mean_absolute_error(y_true, y_pred_raw)
                r2 = metrics.r2_score(y_true, y_pred_raw)
                report_lines.append(f"--- {c} (regression) --- MAE={mae:.4f}  R2={r2:.4f}")
        return pd.DataFrame(pred_cols), "\n".join(report_lines)

    def _evaluate(indices: List[int], steps: int):
        gen = repeating_eval_generator_multitask(
            partitioned_model_input, indices, split_params, batch_size, data_path=data_path
        )
        y_true, y_pred = predict_stream_multitask(model, gen, steps, label_cols)
        return _report(y_true, y_pred)

    val_predictions, val_report = _evaluate(val_idx, validation_steps)

    # Use the small held-out test_idx (~10%) for the second report, NOT the
    # full training set — even streamed, evaluating all of train_idx would
    # walk ~80% of the corpus for a second time.
    test_steps = math.ceil(len(test_idx) / batch_size)
    train_predictions, train_report = _evaluate(test_idx, test_steps)

    return (
        model,
        history.history,
        train_predictions,
        train_report,
        val_predictions,
        val_report,
    )
