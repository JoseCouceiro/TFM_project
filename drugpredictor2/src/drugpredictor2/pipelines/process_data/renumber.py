"""Renumber dataset partitions so every file name matches its contents.

Partitions are named ``<prefix>_<start:09d>_<end:09d>`` where the range is meant
to be the 1-based span of the rows stored in that file. Nothing enforced it:

* the notebook's old write loop advanced a fixed ``chunk_size`` per file, so the
  last chunk (and every chunk after a short one) claimed more rows than it held;
* removing duplicated rows (see :mod:`drugpredictor2.pipelines.process_data.dedupe`)
  shortened files without touching their names.

A name such as ``Balanced_000080001_000090000.csv`` can therefore hold 4,303
rows. This module recomputes the ranges from the real row counts, walking the
partitions in name order, and renames them so the names describe the data::

    Balanced_000000001_000010000.csv  (10,000 rows) -> unchanged
    Balanced_000010001_000020000.csv  ( 4,303 rows) -> Balanced_000010001_000014303.csv
    Balanced_000020001_000030000.csv  ( 5,697 rows) -> Balanced_000014304_000020000.csv

Rows are never moved between files and no data is rewritten: only the names
change. Running the plan twice is a no-op.

The module depends on nothing but pandas (plus pyarrow, optionally, to read
parquet row counts from the footer), so ``rename_partitions.py`` can load it even
when Kedro is not installed in the maintenance interpreter.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd

#: Column that identifies a molecule.
KEY_COLUMN = "SMILES"

#: ``<prefix>_<start:09d>_<end:09d>``, e.g. ``Balanced_000000001_000010000``.
PARTITION_RE = re.compile(r"^(?P<prefix>[A-Za-z][A-Za-z0-9]*)_(?P<start>\d{9})_(?P<end>\d{9})$")

#: Prefix used for the intermediate names of the two-phase rename.
TEMP_PREFIX = ".tmp_renumber_"


def parse_partition_name(name: str) -> Optional[Tuple[str, int, int]]:
    """Split ``<prefix>_<start>_<end>`` into its three parts.

    Args:
        name: File name, with or without its suffix.

    Returns:
        ``(prefix, start, end)``, or ``None`` when the name is not a partition.
    """
    match = PARTITION_RE.match(Path(name).stem)
    if not match:
        return None
    return match.group("prefix"), int(match.group("start")), int(match.group("end"))


def list_partitions(directory: str | Path, suffix: str) -> List[Path]:
    """Return the partition files in `directory`, ordered by their row range.

    Args:
        directory: Directory to scan.
        suffix: File suffix, e.g. ``".csv"``.

    Returns:
        Partition paths sorted by ``(start, end)``.
    """
    partitions = [
        path
        for path in Path(directory).iterdir()
        if path.is_file() and path.suffix == suffix and parse_partition_name(path.name)
    ]
    return sorted(partitions, key=lambda path: parse_partition_name(path.name)[1:])


def count_rows(path: str | Path) -> int:
    """Number of rows stored in a partition, read from a single column.

    Args:
        path: Partition file (``.csv`` or ``.parquet``).

    Returns:
        The row count.
    """
    path = Path(path)
    if path.suffix == ".parquet":
        try:
            import pyarrow.parquet as pq  # optional fast path: footer metadata only

            return int(pq.ParquetFile(path).metadata.num_rows)
        except ImportError:  # pragma: no cover - pyarrow ships with pandas[parquet]
            return int(len(pd.read_parquet(path, columns=[KEY_COLUMN])))
    return int(len(pd.read_csv(path, usecols=[KEY_COLUMN])))


@dataclass
class RenumberPlan:
    """The renames a renumbering run would perform (or performed).

    Attributes:
        directory: The scanned directory.
        suffix: File suffix that was scanned.
        renames: ``[(current name, correct name)]`` in partition order.
        rows: Total row count across the partitions.
    """

    directory: str
    suffix: str
    renames: List[Tuple[str, str]] = field(default_factory=list)
    rows: int = 0

    @property
    def changes(self) -> List[Tuple[str, str]]:
        """Only the entries whose name actually has to change."""
        return [(old, new) for old, new in self.renames if old != new]

    @property
    def is_clean(self) -> bool:
        """True when every name already matches its contents."""
        return not self.changes

    def to_dict(self) -> dict:
        """JSON-serialisable summary."""
        return {
            "directory": self.directory,
            "suffix": self.suffix,
            "files": len(self.renames),
            "rows": self.rows,
            "renamed": len(self.changes),
            "renames": [{"from": old, "to": new} for old, new in self.changes],
        }


def plan_renumbering(directory: str | Path, suffix: str) -> RenumberPlan:
    """Work out the correct name for every partition in `directory`.

    Args:
        directory: Directory holding the partitions.
        suffix: File suffix to scan, ``".csv"`` or ``".parquet"``.

    Returns:
        The :class:`RenumberPlan`, in partition order.

    Raises:
        FileNotFoundError: No partition matches `suffix` in `directory`.
        RuntimeError: A previous run left temporary files behind.
        ValueError: Mixed prefixes, or an empty partition.
    """
    directory = Path(directory)
    partitions = list_partitions(directory, suffix)
    if not partitions:
        raise FileNotFoundError(f"no {suffix} partitions found in {directory}")

    leftovers = sorted(path.name for path in directory.iterdir() if path.name.startswith(TEMP_PREFIX))
    if leftovers:
        raise RuntimeError(
            f"{len(leftovers)} temporary file(s) left behind by an interrupted run, e.g. "
            f"{leftovers[0]!r}. Rename them back before running this again."
        )

    plan = RenumberPlan(str(directory), suffix)
    prefixes = set()
    running = 0
    for path in partitions:
        prefix, _start, _end = parse_partition_name(path.name)
        prefixes.add(prefix)

        rows = count_rows(path)
        if rows <= 0:
            raise ValueError(f"{path.name} holds no rows; delete it before renumbering")

        # Rows already sit in partition order, so the correct range for this file
        # is simply [running + 1, running + rows].
        new_name = f"{prefix}_{running + 1:09d}_{running + rows:09d}{suffix}"
        plan.renames.append((path.name, new_name))
        running += rows

    if len(prefixes) > 1:
        raise ValueError(
            f"{directory} mixes partition prefixes {sorted(prefixes)}; renumber them separately"
        )

    plan.rows = running
    return plan


def apply_renumbering(directory: str | Path, suffix: str, dry_run: bool = True) -> RenumberPlan:
    """Rename the partitions in `directory` to their correct ranges.

    The rename happens in two phases (everything to a temporary name, then to the
    final name) so a new name can never collide with a name still in place.

    Args:
        directory: Directory holding the partitions.
        suffix: File suffix to scan, ``".csv"`` or ``".parquet"``.
        dry_run: When True (the default) nothing on disk is touched.

    Returns:
        The :class:`RenumberPlan` that was applied (or would be applied).
    """
    directory = Path(directory)
    plan = plan_renumbering(directory, suffix)
    if dry_run or plan.is_clean:
        return plan

    staged: List[Tuple[Path, str]] = []
    for index, (old, new) in enumerate(plan.renames):
        if old == new:
            continue
        temporary = directory / f"{TEMP_PREFIX}{index:04d}_{old}"
        (directory / old).rename(temporary)
        staged.append((temporary, new))

    for temporary, new in staged:
        temporary.rename(directory / new)

    return plan
