"""Duplicate detection and cleanup for dataset partitions.

The raw partitions in ``data/01_raw/csv_files`` and the featurized partitions in
``data/05_model_input/featurized_data`` are both sharded with the naming scheme
``<prefix>_<start:09d>_<end:09d>`` (1-based row ranges), one molecule per row,
identified by its ``SMILES`` string.

The sharding loop that produced them used an in-memory ``last_index`` carried in
the notebook kernel, so re-running it wrote the same molecules again under fresh
file names. Two independent kinds of duplication are therefore possible:

1. **Wholly duplicated partitions** - two files with byte-identical content.
   These are redundant files, so the later copy is deleted outright.
2. **Molecules repeated across distinct partitions** - the same ``SMILES``
   present in more than one file. Only the offending *rows* are dropped, so the
   partitions stay intact for Kedro's ``PartitionedDataSet``.

Both cleanups keep the **first** occurrence in deterministic partition order
(sorted by the numeric range in the file name) and drop the later ones. Because
the rule only depends on the file names and row order, cleaning the raw CSVs and
cleaning the featurized parquet files with this module produces a consistent
pair of datasets.
"""

from __future__ import annotations

import hashlib
import os
import re
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

#: Column that identifies a molecule.
KEY_COLUMN = "SMILES"

#: ``<prefix>_<start:09d>_<end:09d>``, e.g. ``Balanced_000000001_000010000``.
PARTITION_RE = re.compile(r"^(?P<prefix>[A-Za-z][A-Za-z0-9]*)_(?P<start>\d{9})_(?P<end>\d{9})$")


@dataclass
class DedupePlan:
    """What a dedupe run on a single directory would do (or did).

    Attributes:
        directory: The scanned directory.
        suffix: File suffix that was scanned (``.csv`` or ``.parquet``).
        files_scanned: Number of partition files found.
        redundant_files: Byte-identical copies that get deleted (later copies).
        rows_to_drop: ``{file name: [row positions]}`` for files that need a
            partial rewrite. Row positions are zero-based positions in the file
            as it currently exists on disk.
        rows_before: Total row count across all scanned files.
        rows_after: Total row count once the plan is applied.
        duplicate_rows: ``rows_before - rows_after``.
        repeated_keys: Number of distinct molecules that appear more than once.
    """

    directory: str
    suffix: str
    files_scanned: int
    redundant_files: List[str] = field(default_factory=list)
    rows_to_drop: Dict[str, List[int]] = field(default_factory=dict)
    rows_before: int = 0
    rows_after: int = 0
    duplicate_rows: int = 0
    repeated_keys: int = 0

    @property
    def is_clean(self) -> bool:
        """True when there is nothing to remove."""
        return not self.redundant_files and not self.rows_to_drop

    def to_dict(self) -> dict:
        """JSON-serialisable summary (row position lists are summarised)."""
        return {
            "directory": self.directory,
            "suffix": self.suffix,
            "files_scanned": self.files_scanned,
            "redundant_files": self.redundant_files,
            "files_rewritten": {name: len(rows) for name, rows in self.rows_to_drop.items()},
            "rows_before": self.rows_before,
            "rows_after": self.rows_after,
            "duplicate_rows": self.duplicate_rows,
            "repeated_molecules": self.repeated_keys,
        }


def parse_partition_range(stem: str) -> Optional[Tuple[int, int]]:
    """Parse the ``(start, end)`` row range out of a partition file stem.

    Args:
        stem: File name without suffix, e.g. ``Balanced_000000001_000010000``.

    Returns:
        The ``(start, end)`` pair, or ``None`` when the name does not follow the
        partition naming scheme.
    """
    match = PARTITION_RE.match(stem)
    if match is None:
        return None
    return int(match.group("start")), int(match.group("end"))


def list_partitions(directory: str | Path, suffix: str) -> List[Path]:
    """Return the partition files in a directory, in deterministic order.

    Ordering is by the numeric range in the file name (falling back to the plain
    name for files that do not match), so "first occurrence wins" is stable
    across runs and across the raw/featurized copies of the dataset.

    Args:
        directory: Directory to scan.
        suffix: File suffix to select, e.g. ``".csv"``.

    Returns:
        Sorted list of matching paths.
    """
    directory = Path(directory)
    if not directory.is_dir():
        return []

    def sort_key(path: Path) -> Tuple[int, int, str]:
        parsed = parse_partition_range(path.stem)
        if parsed is None:
            return (10**12, 10**12, path.name)
        return (parsed[0], parsed[1], path.name)

    return sorted(
        (p for p in directory.iterdir() if p.is_file() and p.suffix == suffix),
        key=sort_key,
    )


def file_digest(path: str | Path, chunk_size: int = 1 << 20) -> str:
    """Streaming SHA-256 of a file, so large partitions never load into memory."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def find_identical_files(paths: Sequence[Path]) -> List[List[Path]]:
    """Group files whose content is byte-for-byte identical.

    Args:
        paths: Files to compare.

    Returns:
        One list per duplicated group, ordered by first appearance. Only groups
        with two or more members are returned, so a clean directory yields ``[]``.
    """
    by_digest: Dict[str, List[Path]] = {}
    for path in paths:
        by_digest.setdefault(file_digest(path), []).append(path)
    return [group for group in by_digest.values() if len(group) > 1]


def read_key_column(path: str | Path, key: str = KEY_COLUMN) -> pd.Series:
    """Read only the molecule-identifying column of a partition.

    Reading a single column keeps the peak memory of a full scan low (the
    featurized partitions carry a 5,288-element fingerprint column).

    Args:
        path: Partition file (``.csv`` or ``.parquet``).
        key: Column to read.

    Returns:
        The column as a ``Series``.
    """
    path = Path(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path, columns=[key])[key]
    return pd.read_csv(path, usecols=[key])[key]


def plan_dedupe(directory: str | Path, suffix: str, key: str = KEY_COLUMN) -> DedupePlan:
    """Work out which files to delete and which rows to drop. Read-only.

    Args:
        directory: Directory holding the partitions.
        suffix: File suffix to scan, ``".csv"`` or ``".parquet"``.
        key: Molecule-identifying column.

    Returns:
        A :class:`DedupePlan` describing the cleanup. Applying it removes
        ``redundant_files`` completely and drops the listed row positions.
    """
    directory = Path(directory)
    plan = DedupePlan(directory=str(directory), suffix=suffix, files_scanned=0)

    partitions = list_partitions(directory, suffix)
    plan.files_scanned = len(partitions)

    # 1) Whole files that are byte-identical copies: keep the first, delete the rest.
    redundant = set()
    for group in find_identical_files(partitions):
        redundant.update(group[1:])
    plan.redundant_files = sorted(p.name for p in redundant)

    # 2) Molecules repeated across the surviving partitions: keep first, drop later.
    #    Keys are added to `seen` row by row so a molecule repeated *within* one
    #    partition is caught as well.
    seen: set = set()
    repeated: set = set()
    for path in partitions:
        if path in redundant:
            continue
        keys = read_key_column(path, key)
        plan.rows_before += len(keys)

        positions = []
        for position, value in enumerate(keys):
            if value in seen:
                positions.append(position)
                repeated.add(value)
            else:
                seen.add(value)

        if positions:
            plan.rows_to_drop[path.name] = positions

        plan.rows_after += len(keys) - len(positions)

    plan.duplicate_rows = plan.rows_before - plan.rows_after
    plan.repeated_keys = len(repeated)
    return plan


def _atomic_replace(target: Path, writer) -> None:
    """Write via a sibling temp file and ``os.replace`` so a crash cannot truncate.

    Args:
        target: Final file path.
        writer: Callable receiving a temp path and persisting the data there.
    """
    handle, tmp_name = tempfile.mkstemp(dir=str(target.parent), suffix=target.suffix, prefix=".tmp_dedupe_")
    os.close(handle)
    tmp_path = Path(tmp_name)
    try:
        writer(tmp_path)
        os.replace(tmp_path, target)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def _rewrite_partition(path: Path, drop_positions: Sequence[int], key: str = KEY_COLUMN) -> None:
    """Drop the given row positions from one partition, in place and atomically.

    CSV partitions are rewritten line-by-line so that every kept row keeps its
    original bytes (a pandas round-trip re-formats floats and churns the whole
    file). The line-based rewrite is verified against a column read before
    falling back to a full pandas rewrite, which keeps it safe if a field ever
    contains an embedded newline.

    Args:
        path: Partition to rewrite.
        drop_positions: Zero-based row positions to remove.
        key: Molecule-identifying column, used to verify the rewrite.
    """
    if path.suffix != ".parquet":
        expected = (
            pd.read_csv(path, usecols=[key])[key].drop(index=list(drop_positions)).tolist()
        )

        def write_lines(tmp_path: Path) -> None:
            with open(path, encoding="utf-8", newline="") as handle:
                lines = handle.readlines()
            # Line 0 is the header, so row position `p` lives on line `p + 1`.
            drop_lines = {position + 1 for position in drop_positions}
            with open(tmp_path, "w", encoding="utf-8", newline="") as handle:
                handle.writelines(line for index, line in enumerate(lines) if index not in drop_lines)
            if pd.read_csv(tmp_path, usecols=[key])[key].tolist() != expected:
                raise ValueError(f"line-based rewrite of {path.name} did not reproduce the expected rows")

        try:
            _atomic_replace(path, write_lines)
            return
        except ValueError:
            pass  # fall through to the generic pandas rewrite

    frame = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
    cleaned = frame.drop(index=list(drop_positions))

    if path.suffix == ".parquet":
        _atomic_replace(path, lambda tmp: cleaned.to_parquet(tmp))
    else:
        _atomic_replace(path, lambda tmp: cleaned.to_csv(tmp, index=False))


def apply_dedupe(
    directory: str | Path,
    suffix: str,
    key: str = KEY_COLUMN,
    dry_run: bool = True,
) -> DedupePlan:
    """Clean a partition directory: delete duplicate files, drop duplicate rows.

    Args:
        directory: Directory holding the partitions.
        suffix: File suffix to scan, ``".csv"`` or ``".parquet"``.
        key: Molecule-identifying column.
        dry_run: When True (the default) nothing on disk is touched.

    Returns:
        The :class:`DedupePlan` that was applied (or would be applied).
    """
    directory = Path(directory)
    plan = plan_dedupe(directory, suffix, key)
    if dry_run or plan.is_clean:
        return plan

    for name in plan.redundant_files:
        (directory / name).unlink()

    for name, positions in plan.rows_to_drop.items():
        _rewrite_partition(directory / name, positions, key)

    return plan


def count_duplicate_rows(directory: str | Path, suffix: str, key: str = KEY_COLUMN) -> dict:
    """Report how many molecules repeat in a directory. Read-only verification.

    Args:
        directory: Directory holding the partitions.
        suffix: File suffix to scan.
        key: Molecule-identifying column.

    Returns:
        ``{"files", "rows", "distinct", "duplicate_rows", "repeated_molecules"}``.
    """
    frames = [read_key_column(path, key) for path in list_partitions(directory, suffix)]
    if not frames:
        return {"files": 0, "rows": 0, "distinct": 0, "duplicate_rows": 0, "repeated_molecules": 0}

    all_keys = pd.concat(frames, ignore_index=True)
    counts = all_keys.value_counts()
    return {
        "files": len(frames),
        "rows": int(len(all_keys)),
        "distinct": int(all_keys.nunique()),
        "duplicate_rows": int(all_keys.duplicated().sum()),
        "repeated_molecules": int((counts > 1).sum()),
    }
