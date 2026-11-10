"""Remove duplicated molecules from the dataset partitions.

Cleans both layers of the dataset with the same deterministic rule (keep the
first occurrence in partition-range order):

* ``data/01_raw/csv_files`` - the raw sharded CSVs.
* ``data/05_model_input/featurized_data`` - the featurized parquet shards.

Two kinds of duplication are removed:

* whole partitions that are byte-identical copies (deleted outright),
* molecules repeated across distinct partitions (only the repeated rows are
  dropped, so Kedro's ``PartitionedDataSet`` keeps working).

Nothing is written unless ``--apply`` is passed, so the default run is a
report-only dry run.

Usage:
    python clean_dataset.py                # dry run, print what would change
    python clean_dataset.py --apply        # clean the data in place
    python clean_dataset.py --apply --skip-raw
"""

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent
DEDUPE_MODULE_PATH = PROJECT_ROOT / "src" / "drugpredictor2" / "pipelines" / "process_data" / "dedupe.py"

try:  # normal case: the project package is installed (`pip install -e src`)
    from drugpredictor2.pipelines.process_data.dedupe import (  # noqa: E402
        KEY_COLUMN,
        apply_dedupe,
        count_duplicate_rows,
        list_partitions,
        read_key_column,
    )
except ImportError:
    # Kedro (and therefore the package) may not be installed in the interpreter
    # running this maintenance script; it only needs pandas, so load the module
    # straight from the source tree instead of importing the package.
    sys.path.insert(0, str(PROJECT_ROOT / "src"))
    _spec = importlib.util.spec_from_file_location("drugpredictor2_dedupe", DEDUPE_MODULE_PATH)
    _dedupe = importlib.util.module_from_spec(_spec)
    sys.modules[_spec.name] = _dedupe  # dataclasses resolves annotations via sys.modules
    _spec.loader.exec_module(_dedupe)
    KEY_COLUMN = _dedupe.KEY_COLUMN
    apply_dedupe = _dedupe.apply_dedupe
    count_duplicate_rows = _dedupe.count_duplicate_rows
    list_partitions = _dedupe.list_partitions
    read_key_column = _dedupe.read_key_column

RAW_DIR = PROJECT_ROOT / "data" / "01_raw" / "csv_files"
FEATURIZED_DIR = PROJECT_ROOT / "data" / "05_model_input" / "featurized_data"
DEFAULT_REPORT = PROJECT_ROOT / "data" / "08_reporting" / "dedupe_report.json"


def _summarise(label: str, plan) -> dict:
    """Print the plan for one directory and return its JSON summary."""
    print(f"\n=== {label} ===")
    print(f"directory            : {plan.directory}")
    print(f"partitions scanned   : {plan.files_scanned}")
    print(f"redundant file copies: {len(plan.redundant_files)}")
    for name in plan.redundant_files:
        print(f"    - delete {name}")
    print(f"molecules repeated   : {plan.repeated_keys}")
    print(f"rows before          : {plan.rows_before}")
    print(f"rows after           : {plan.rows_after}")
    print(f"rows to drop         : {plan.duplicate_rows}")
    if plan.rows_to_drop:
        for name, positions in plan.rows_to_drop.items():
            print(f"    - {name}: drop {len(positions)} row(s)")
    return plan.to_dict()


def _check_consistency(raw_dir: Path, featurized_dir: Path) -> dict:
    """Confirm the raw and featurized layers describe the same molecules."""
    raw = pd.concat(
        [read_key_column(path, KEY_COLUMN) for path in list_partitions(raw_dir, ".csv")],
        ignore_index=True,
    )
    featurized = pd.concat(
        [read_key_column(path, KEY_COLUMN) for path in list_partitions(featurized_dir, ".parquet")],
        ignore_index=True,
    )
    result = {
        "raw_rows": int(len(raw)),
        "featurized_rows": int(len(featurized)),
        "raw_distinct": int(raw.nunique()),
        "featurized_distinct": int(featurized.nunique()),
        "same_molecule_sets": bool(set(raw.unique()) == set(featurized.unique())),
        "same_row_counts": bool(len(raw) == len(featurized)),
    }
    print("\n=== raw vs featurized consistency ===")
    for name, value in result.items():
        print(f"{name:<22}: {value}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--apply", action="store_true", help="write the changes (default is a dry run)")
    parser.add_argument("--raw-dir", type=Path, default=RAW_DIR, help=f"raw CSV partitions (default: {RAW_DIR})")
    parser.add_argument(
        "--featurized-dir",
        type=Path,
        default=FEATURIZED_DIR,
        help=f"featurized parquet partitions (default: {FEATURIZED_DIR})",
    )
    parser.add_argument("--skip-raw", action="store_true", help="do not touch the raw CSV partitions")
    parser.add_argument(
        "--skip-featurized", action="store_true", help="do not touch the featurized parquet partitions"
    )
    parser.add_argument(
        "--report", type=Path, default=DEFAULT_REPORT, help=f"JSON report path (default: {DEFAULT_REPORT})"
    )
    args = parser.parse_args()

    dry_run = not args.apply
    if dry_run:
        print("DRY RUN - nothing will be written. Re-run with --apply to clean the data.")
    else:
        print("APPLYING CHANGES")

    report = {"dry_run": dry_run, "directories": {}}

    if not args.skip_raw:
        report["directories"]["raw"] = _summarise(
            "RAW CSV PARTITIONS", apply_dedupe(args.raw_dir, ".csv", dry_run=dry_run)
        )
    if not args.skip_featurized:
        report["directories"]["featurized"] = _summarise(
            "FEATURIZED PARTITIONS", apply_dedupe(args.featurized_dir, ".parquet", dry_run=dry_run)
        )

    if args.apply and not (args.skip_raw or args.skip_featurized):
        print("\n=== verification after cleanup ===")
        report["verification"] = {}
        for label, directory, suffix in (
            ("raw", args.raw_dir, ".csv"),
            ("featurized", args.featurized_dir, ".parquet"),
        ):
            stats = count_duplicate_rows(directory, suffix)
            report["verification"][label] = stats
            print(f"{label + ' duplicates remaining':<34}: {stats['duplicate_rows']}")
        report["consistency"] = _check_consistency(args.raw_dir, args.featurized_dir)

    args.report.parent.mkdir(parents=True, exist_ok=True)
    with open(args.report, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    print(f"\nReport written to {args.report}")

    leftover = sum(entry.get("duplicate_rows", 0) for entry in report.get("verification", {}).values())
    return 1 if leftover else 0


if __name__ == "__main__":
    raise SystemExit(main())
