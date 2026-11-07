"""Tests for the partition renumbering helpers in ``process_data.renumber``.

Run with ``pytest src/tests/pipelines/process_data/test_renumber.py``.
"""

from pathlib import Path

import pandas as pd
import pytest

from drugpredictor2.pipelines.process_data.renumber import (
    TEMP_PREFIX,
    apply_renumbering,
    count_rows,
    list_partitions,
    parse_partition_name,
    plan_renumbering,
)


def _write_csv(directory: Path, name: str, rows: int) -> Path:
    """Write a partition holding `rows` distinct molecules."""
    path = directory / name
    pd.DataFrame({"SMILES": [f"C{'C' * index}" for index in range(rows)]}).to_csv(path, index=False)
    return path


def _write_parquet(directory: Path, name: str, rows: int) -> Path:
    path = directory / name
    pd.DataFrame(
        {"SMILES": [f"C{'C' * index}" for index in range(rows)], "FP": [[0.0, 1.0]] * rows}
    ).to_parquet(path)
    return path


def test_parse_partition_name():
    assert parse_partition_name("Balanced_000080001_000090000.csv") == ("Balanced", 80001, 90000)
    assert parse_partition_name("Balanced_000080001_000090000") == ("Balanced", 80001, 90000)
    assert parse_partition_name("notes.csv") is None
    assert parse_partition_name(".tmp_renumber_0000_Balanced_000000001_000000002.csv") is None


def test_list_partitions_sorts_by_range_and_counts_rows(tmp_path):
    _write_csv(tmp_path, "Balanced_000000011_000000020.csv", 2)
    _write_csv(tmp_path, "Balanced_000000001_000000010.csv", 3)
    (tmp_path / "notes.csv").write_text("SMILES\nCCO\n", encoding="utf-8")

    partitions = list_partitions(tmp_path, ".csv")

    assert [path.name for path in partitions] == [
        "Balanced_000000001_000000010.csv",
        "Balanced_000000011_000000020.csv",
    ]
    assert [count_rows(path) for path in partitions] == [3, 2]


def test_plan_renumbering_uses_the_real_row_counts(tmp_path):
    _write_csv(tmp_path, "Balanced_000000001_000010000.csv", 4)
    _write_csv(tmp_path, "Balanced_000010001_000020000.csv", 3)
    _write_csv(tmp_path, "Balanced_000020001_000030000.csv", 5)

    plan = plan_renumbering(tmp_path, ".csv")

    assert plan.rows == 12
    assert plan.renames == [
        ("Balanced_000000001_000010000.csv", "Balanced_000000001_000000004.csv"),
        ("Balanced_000010001_000020000.csv", "Balanced_000000005_000000007.csv"),
        ("Balanced_000020001_000030000.csv", "Balanced_000000008_000000012.csv"),
    ]
    assert len(plan.changes) == 3


def test_plan_renumbering_leaves_already_correct_names_alone(tmp_path):
    _write_csv(tmp_path, "Balanced_000000001_000000004.csv", 4)
    _write_csv(tmp_path, "Balanced_000000005_000000008.csv", 4)

    plan = plan_renumbering(tmp_path, ".csv")

    assert plan.is_clean
    assert plan.renames == [
        ("Balanced_000000001_000000004.csv", "Balanced_000000001_000000004.csv"),
        ("Balanced_000000005_000000008.csv", "Balanced_000000005_000000008.csv"),
    ]


def test_apply_renumbering_dry_run_touches_nothing(tmp_path):
    _write_csv(tmp_path, "Balanced_000000001_000010000.csv", 2)
    before = sorted(path.name for path in tmp_path.iterdir())

    plan = apply_renumbering(tmp_path, ".csv", dry_run=True)

    assert len(plan.changes) == 1
    assert sorted(path.name for path in tmp_path.iterdir()) == before


def test_apply_renumbering_renames_and_keeps_the_content(tmp_path):
    first = _write_csv(tmp_path, "Balanced_000000001_000010000.csv", 2)
    second = _write_csv(tmp_path, "Balanced_000010001_000020000.csv", 3)
    first_bytes = first.read_bytes()
    second_bytes = second.read_bytes()

    apply_renumbering(tmp_path, ".csv", dry_run=False)

    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "Balanced_000000001_000000002.csv",
        "Balanced_000000003_000000005.csv",
    ]
    assert (tmp_path / "Balanced_000000001_000000002.csv").read_bytes() == first_bytes
    assert (tmp_path / "Balanced_000000003_000000005.csv").read_bytes() == second_bytes
    # A second pass has nothing left to do.
    assert plan_renumbering(tmp_path, ".csv").is_clean


def test_apply_renumbering_handles_parquet(tmp_path):
    _write_parquet(tmp_path, "Balanced_000000001_000010000.parquet", 2)
    _write_parquet(tmp_path, "Balanced_000010001_000020000.parquet", 1)
    _write_parquet(tmp_path, "Balanced_000020001_000030000.parquet", 4)

    apply_renumbering(tmp_path, ".parquet", dry_run=False)

    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "Balanced_000000001_000000002.parquet",
        "Balanced_000000003_000000003.parquet",
        "Balanced_000000004_000000007.parquet",
    ]
    frame = pd.read_parquet(tmp_path / "Balanced_000000004_000000007.parquet")
    assert frame["SMILES"].tolist() == ["C", "CC", "CCC", "CCCC"]
    assert [len(fingerprint) for fingerprint in frame["FP"]] == [2, 2, 2, 2]


def test_plan_renumbering_rejects_leftover_temporary_files(tmp_path):
    _write_csv(tmp_path, "Balanced_000000001_000010000.csv", 2)
    (tmp_path / f"{TEMP_PREFIX}0000_Balanced_000000001_000010000.csv").write_text("x", encoding="utf-8")

    with pytest.raises(RuntimeError, match="interrupted run"):
        plan_renumbering(tmp_path, ".csv")


def test_plan_renumbering_rejects_mixed_prefixes(tmp_path):
    _write_csv(tmp_path, "Balanced_000000001_000000002.csv", 2)
    _write_csv(tmp_path, "Extra_000000003_000000004.csv", 2)

    with pytest.raises(ValueError, match="mixes partition prefixes"):
        plan_renumbering(tmp_path, ".csv")


def test_plan_renumbering_rejects_empty_partitions(tmp_path):
    (tmp_path / "Balanced_000000001_000010000.csv").write_text("SMILES\n", encoding="utf-8")

    with pytest.raises(ValueError, match="holds no rows"):
        plan_renumbering(tmp_path, ".csv")


def test_plan_renumbering_requires_at_least_one_partition(tmp_path):
    with pytest.raises(FileNotFoundError, match="no .csv partitions"):
        plan_renumbering(tmp_path, ".csv")
