"""Tests for the partition duplicate-detection helpers in ``process_data.dedupe``.

Run with ``pytest src/tests/pipelines/process_data/test_dedupe.py``.
"""

from pathlib import Path

import pandas as pd
import pytest

from drugpredictor2.pipelines.process_data.dedupe import (
    KEY_COLUMN,
    apply_dedupe,
    count_duplicate_rows,
    find_identical_files,
    list_partitions,
    parse_partition_range,
    plan_dedupe,
)


def _write_csv(path: Path, smiles: list) -> None:
    """Write a minimal raw-style partition."""
    pd.DataFrame({"SMILES": smiles, "Molecular Weight": range(len(smiles))}).to_csv(path, index=False)


def _write_parquet(path: Path, smiles: list) -> None:
    """Write a minimal featurized-style partition (with a list column)."""
    pd.DataFrame(
        {"SMILES": smiles, "Molecular Weight": range(len(smiles)), "FP": [[0.0, 1.0] for _ in smiles]}
    ).to_parquet(path)


def test_parse_partition_range_accepts_the_shard_naming_scheme():
    assert parse_partition_range("Balanced_000000001_000010000") == (1, 10000)
    assert parse_partition_range("Balanced_000280001_000290000") == (280001, 290000)


def test_parse_partition_range_rejects_other_names():
    assert parse_partition_range("Balanced") is None
    assert parse_partition_range("notes_000000001") is None
    assert parse_partition_range(".tmp_dedupe_123") is None


def test_list_partitions_sorts_by_numeric_range(tmp_path):
    for name in ("Balanced_000000011_000000020", "Balanced_000000001_000000010", "notes"):
        _write_csv(tmp_path / f"{name}.csv", ["C"])
    (tmp_path / "ignored.parquet").write_bytes(b"")

    names = [path.name for path in list_partitions(tmp_path, ".csv")]

    assert names == [
        "Balanced_000000001_000000010.csv",
        "Balanced_000000011_000000020.csv",
        "notes.csv",
    ]


def test_find_identical_files_groups_byte_identical_copies(tmp_path):
    _write_csv(tmp_path / "Balanced_000000001_000000002.csv", ["CCO", "CCN"])
    _write_csv(tmp_path / "Balanced_000000011_000000012.csv", ["CCO", "CCN"])
    _write_csv(tmp_path / "Balanced_000000021_000000022.csv", ["CCC", "CCCl"])

    groups = find_identical_files(list_partitions(tmp_path, ".csv"))

    assert len(groups) == 1
@pytest.fixture()
def duplicated_csv_dir(tmp_path):
    """Two one-off partitions plus a byte-identical copy and cross-file dupes."""
    _write_csv(tmp_path / "Balanced_000000001_000000003.csv", ["CCO", "CCN", "CCC"])
    _write_csv(tmp_path / "Balanced_000000011_000000013.csv", ["CCO", "CCN", "CCC"])  # identical copy
    _write_csv(tmp_path / "Balanced_000000021_000000025.csv", ["CCO", "CCCl", "CCO", "CCBr"])
    return tmp_path


def test_plan_dedupe_reports_redundant_files_and_repeated_rows(duplicated_csv_dir):
    plan = plan_dedupe(duplicated_csv_dir, ".csv")

    assert plan.files_scanned == 3
    assert plan.redundant_files == ["Balanced_000000011_000000013.csv"]
    # Position 0 repeats a molecule from the first partition, position 2 repeats
    # position 0 inside the same file.
    assert plan.rows_to_drop == {"Balanced_000000021_000000025.csv": [0, 2]}
    assert plan.rows_before == 7
    assert plan.rows_after == 5
    assert plan.duplicate_rows == 2
    assert plan.repeated_keys == 1
    assert not plan.is_clean


def test_apply_dedupe_dry_run_leaves_the_directory_untouched(duplicated_csv_dir):
    before = {path.name: path.read_bytes() for path in duplicated_csv_dir.iterdir()}

    plan = apply_dedupe(duplicated_csv_dir, ".csv", dry_run=True)

    assert plan.duplicate_rows == 2
    assert {path.name: path.read_bytes() for path in duplicated_csv_dir.iterdir()} == before


def test_apply_dedupe_cleans_csv_partitions(duplicated_csv_dir):
    target = duplicated_csv_dir / "Balanced_000000021_000000025.csv"
    original_lines = target.read_text().splitlines()

    apply_dedupe(duplicated_csv_dir, ".csv", dry_run=False)

    assert sorted(path.name for path in duplicated_csv_dir.iterdir()) == [
        "Balanced_000000001_000000003.csv",
        "Balanced_000000021_000000025.csv",
    ]
    # Kept rows keep their original bytes; dropped rows are the duplicated ones.
    assert target.read_text().splitlines() == [original_lines[0], original_lines[2], original_lines[4]]
    assert count_duplicate_rows(duplicated_csv_dir, ".csv") == {
        "files": 2,
        "rows": 5,
        "distinct": 5,
        "duplicate_rows": 0,
        "repeated_molecules": 0,
    }


def test_apply_dedupe_cleans_parquet_partitions(tmp_path):
    _write_parquet(tmp_path / "Balanced_000000001_000000003.parquet", ["CCO", "CCN", "CCC"])
    _write_parquet(tmp_path / "Balanced_000000011_000000013.parquet", ["CCO", "CCN", "CCC"])
    _write_parquet(tmp_path / "Balanced_000000021_000000025.parquet", ["CCO", "CCCl", "CCO", "CCBr"])

    apply_dedupe(tmp_path, ".parquet", dry_run=False)

    partitions = list_partitions(tmp_path, ".parquet")
    assert [path.name for path in partitions] == [
        "Balanced_000000001_000000003.parquet",
        "Balanced_000000021_000000025.parquet",
    ]
    # The fingerprint column must survive the partial rewrite.
    frame = pd.read_parquet(partitions[1])
    assert frame[KEY_COLUMN].tolist() == ["CCCl", "CCBr"]
    assert [len(fingerprint) for fingerprint in frame["FP"]] == [2, 2]
    # No temporary files left behind by the atomic writes.
    assert list(tmp_path.glob(".tmp_dedupe_*")) == []


def test_apply_dedupe_is_a_no_op_on_a_clean_directory(tmp_path):
    _write_csv(tmp_path / "Balanced_000000001_000000002.csv", ["CCO", "CCN"])

    plan = apply_dedupe(tmp_path, ".csv", dry_run=False)

    assert plan.is_clean
    assert count_duplicate_rows(tmp_path, ".csv")["duplicate_rows"] == 0
