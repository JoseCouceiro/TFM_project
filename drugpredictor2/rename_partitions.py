"""Rename the dataset partitions so each file name spans exactly its own rows.

Partitions are named ``<prefix>_<start:09d>_<end:09d>``, e.g.
``Balanced_000080001_000090000.csv``. The range is supposed to describe the rows
stored in the file, but nothing kept it honest: the notebook's old write loop
advanced a fixed chunk size per file, and removing duplicated rows shortened
files without renaming them. ``Balanced_000080001_000090000.csv`` therefore only
holds 4,303 molecules.

This script walks the partitions in name order, re-derives every range from the
real row counts, and renames the files:

    Balanced_000080001_000090000.csv  (4,303 rows) -> Balanced_000080001_000084303.csv

Both layers are renumbered with the same rule, so the raw CSVs and the featurized
parquet shards keep matching names (``process_new_partitions`` copies the CSV stem
to the parquet file). No row is moved between files and no content is rewritten.

Nothing on disk is touched unless ``--apply`` is passed.

Usage:
    python rename_partitions.py                # dry run, print the renames
    python rename_partitions.py --apply        # rename in place
    python rename_partitions.py --apply --skip-raw
"""

import argparse
import importlib.util
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
RENUMBER_MODULE_PATH = (
    PROJECT_ROOT / "src" / "drugpredictor2" / "pipelines" / "process_data" / "renumber.py"
)

try:  # normal case: the project package is installed (`pip install -e src`)
    from drugpredictor2.pipelines.process_data.renumber import (  # noqa: E402
        apply_renumbering,
        plan_renumbering,
    )
except ImportError:
    # Kedro (and therefore the package) may not be installed in the interpreter
    # running this maintenance script; it only needs pandas, so load the module
    # straight from the source tree instead of importing the package.
    _spec = importlib.util.spec_from_file_location("drugpredictor2_renumber", RENUMBER_MODULE_PATH)
    _renumber = importlib.util.module_from_spec(_spec)
    sys.modules[_spec.name] = _renumber  # dataclasses resolves annotations via sys.modules
    _spec.loader.exec_module(_renumber)
    apply_renumbering = _renumber.apply_renumbering
    plan_renumbering = _renumber.plan_renumbering

RAW_DIR = PROJECT_ROOT / "data" / "01_raw" / "csv_files"
FEATURIZED_DIR = PROJECT_ROOT / "data" / "05_model_input" / "featurized_data"
DEFAULT_REPORT = PROJECT_ROOT / "data" / "08_reporting" / "rename_report.json"

LAYERS = (
    ("raw", "RAW CSV PARTITIONS", RAW_DIR, ".csv"),
    ("featurized", "FEATURIZED PARTITIONS", FEATURIZED_DIR, ".parquet"),
)


def _summarise(label: str, plan) -> dict:
    """Print the plan for one directory and return its JSON summary."""
    print(f"\n=== {label} ===")
    print(f"directory : {plan.directory}")
    print(f"partitions: {len(plan.renames)}")
    print(f"rows      : {plan.rows}")
    print(f"renames   : {len(plan.changes)}")
    for old, new in plan.changes:
        print(f"    {old} -> {new}")
    if plan.is_clean:
        print("    (every name already matches its contents)")
    return plan.to_dict()


def _check_layers(raw_dir: Path, featurized_dir: Path) -> dict:
    """Confirm both layers ended up with the same partition names."""
    raw_stems = {Path(new).stem for _old, new in plan_renumbering(raw_dir, ".csv").renames}
    featurized_stems = {
        Path(new).stem for _old, new in plan_renumbering(featurized_dir, ".parquet").renames
    }

    result = {
        "raw_partitions": len(raw_stems),
        "featurized_partitions": len(featurized_stems),
        "same_names": bool(raw_stems == featurized_stems),
        "missing_in_featurized": sorted(raw_stems - featurized_stems)[:5],
        "missing_in_raw": sorted(featurized_stems - raw_stems)[:5],
    }
    print("\n=== raw vs featurized names ===")
    for name, value in result.items():
        print(f"{name:<22}: {value}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--apply", action="store_true", help="rename the files (default is a dry run)")
    parser.add_argument("--raw-dir", type=Path, default=RAW_DIR, help=f"raw CSV partitions (default: {RAW_DIR})")
    parser.add_argument(
        "--featurized-dir",
        type=Path,
        default=FEATURIZED_DIR,
        help=f"featurized parquet partitions (default: {FEATURIZED_DIR})",
    )
    parser.add_argument("--skip-raw", action="store_true", help="do not rename the raw CSV partitions")
    parser.add_argument(
        "--skip-featurized", action="store_true", help="do not rename the featurized parquet partitions"
    )
    parser.add_argument(
        "--report", type=Path, default=DEFAULT_REPORT, help=f"JSON report path (default: {DEFAULT_REPORT})"
    )
    args = parser.parse_args()

    dry_run = not args.apply
    print("DRY RUN - nothing will be renamed. Re-run with --apply to rename the files." if dry_run else "RENAMING")
    if args.skip_raw or args.skip_featurized:
        print(
            "WARNING: only one layer is being renumbered. `process_new_partitions` matches the "
            "featurized shards to the raw ones by file stem, so the layers must share names - "
            "run both unless you know why."
        )

    report = {"dry_run": dry_run, "directories": {}}
    for key, label, directory, suffix in LAYERS:
        if (key == "raw" and args.skip_raw) or (key == "featurized" and args.skip_featurized):
            continue
        report["directories"][key] = _summarise(label, apply_renumbering(directory, suffix, dry_run=dry_run))

    if not (args.skip_raw or args.skip_featurized):
        report["cross_layer"] = _check_layers(args.raw_dir, args.featurized_dir)

    args.report.parent.mkdir(parents=True, exist_ok=True)
    with open(args.report, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    print(f"\nReport written to {args.report}")

    inconsistent = report.get("cross_layer", {}).get("same_names", True)
    return 0 if inconsistent else 1


if __name__ == "__main__":
    raise SystemExit(main())
