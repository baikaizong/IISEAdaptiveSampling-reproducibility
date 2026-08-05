"""Prepare, run, and verify the published solar-monitoring case study."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = (
    REPOSITORY_ROOT
    / "AdaptiveSamplingFigures"
    / "ExternalData"
    / "solar_frames"
    / "data.mat"
)
DEFAULT_CACHE = REPOSITORY_ROOT / "cache" / "case-study" / "norm_data.bin"
DEFAULT_RESULTS = REPOSITORY_ROOT / "results" / "case-study"
PUBLISHED_RESULTS = (
    REPOSITORY_ROOT
    / "AdaptiveSamplingFigures"
    / "ExternalData"
    / "solar_statistics"
)
METHOD_FILES = (
    "mM-KST_ChartingStatistics-5.txt",
    "SASAM_ChartingStatistics-5.txt",
    "TOPR_ChartingStatistics-5.txt",
    "TSSPR_ChartingStatistics-5.txt",
    "NAS_ChartingStatistics-5.txt",
)
EXPECTED_SHAPE = (232, 292, 300)
IC_SIZE = 100
TIME_WINDOW = 299
EXPECTED_CACHE_BYTES = 232 * 292 * TIME_WINDOW * np.dtype("<f8").itemsize


def prepare_cache(source: Path, output: Path, force: bool = False) -> Path:
    """Normalize the solar frames and write Fortran-order float64 data."""

    source = source.resolve()
    output = output.resolve()
    if output.exists() and not force:
        if output.stat().st_size != EXPECTED_CACHE_BYTES:
            raise ValueError(
                f"Existing cache has {output.stat().st_size} bytes; "
                f"expected {EXPECTED_CACHE_BYTES}. Use --force to replace it."
            )
        print(f"Using existing cache: {output}")
        return output

    if not source.is_file():
        raise FileNotFoundError(f"Solar frame data not found: {source}")

    variables = loadmat(source, variable_names=["data"])
    if "data" not in variables:
        raise KeyError(f"MAT file has no 'data' variable: {source}")

    raw = np.asarray(variables["data"], dtype=np.float64)
    if raw.shape != EXPECTED_SHAPE:
        raise ValueError(f"Expected data shape {EXPECTED_SHAPE}, found {raw.shape}")

    differences = np.diff(raw, axis=2)
    baseline = differences[:, :, :IC_SIZE]
    means = baseline.mean(axis=2, keepdims=True)
    standard_deviations = baseline.std(axis=2, ddof=1, keepdims=True)
    if np.any(standard_deviations == 0):
        raise ValueError("At least one in-control pixel has zero standard deviation.")

    differences -= means
    differences /= standard_deviations
    np.round(differences, decimals=6, out=differences)
    if not np.isfinite(differences).all():
        raise ValueError("Normalization produced a non-finite value.")

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    try:
        with temporary.open("wb") as stream:
            for frame in range(TIME_WINDOW):
                values = differences[:, :, frame].ravel(order="F").astype("<f8")
                values.tofile(stream)
        if temporary.stat().st_size != EXPECTED_CACHE_BYTES:
            raise OSError(
                f"Prepared cache has {temporary.stat().st_size} bytes; "
                f"expected {EXPECTED_CACHE_BYTES}."
            )
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)

    print(f"Prepared normalized cache: {output}")
    return output


def verify_results(results: Path, reference: Path) -> bool:
    """Compare five generated series with the canonical plotting inputs."""

    ok = True
    for filename in METHOD_FILES:
        candidate_path = results / filename
        reference_path = reference / filename
        if not candidate_path.is_file() or not reference_path.is_file():
            print(f"MISSING {filename}")
            ok = False
            continue

        candidate = np.loadtxt(candidate_path, dtype=np.float64)
        expected = np.loadtxt(reference_path, dtype=np.float64)
        if candidate.shape != expected.shape:
            print(
                f"FAIL    {filename}: shape {candidate.shape}, "
                f"expected {expected.shape}"
            )
            ok = False
            continue

        max_difference = float(np.max(np.abs(candidate - expected)))
        equal = np.array_equal(candidate, expected)
        print(
            f"{'PASS' if equal else 'DIFF'}    {filename}: "
            f"n={candidate.size}, max_abs_diff={max_difference:.6g}"
        )
        ok = ok and equal
    return ok


def run_case_study(
    executable: Path,
    source: Path,
    cache: Path,
    results: Path,
    seed: int | None,
    force_prepare: bool,
) -> None:
    cache = prepare_cache(source, cache, force=force_prepare)
    executable = executable.resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"Case-study executable not found: {executable}")

    results = results.resolve()
    results.mkdir(parents=True, exist_ok=True)
    command = [str(executable), str(cache), str(results)]
    if seed is not None:
        command.append(str(seed))
    print("Running:", subprocess.list2cmdline(command))
    subprocess.run(command, cwd=REPOSITORY_ROOT, check=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare", help="prepare normalized binary data")
    prepare.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    prepare.add_argument("--output", type=Path, default=DEFAULT_CACHE)
    prepare.add_argument("--force", action="store_true")

    run = subparsers.add_parser("run", help="prepare data and run case_study")
    run.add_argument(
        "--executable",
        type=Path,
        default=REPOSITORY_ROOT / "build" / "windows-debug" / "bin" / "case_study.exe",
    )
    run.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    run.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    run.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    run.add_argument("--seed", type=int, default=20260804)
    run.add_argument("--force-prepare", action="store_true")

    verify = subparsers.add_parser(
        "verify", help="compare generated results with the plotting inputs"
    )
    verify.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    verify.add_argument("--reference", type=Path, default=PUBLISHED_RESULTS)
    return parser


def main() -> int:
    arguments = build_parser().parse_args()
    if arguments.command == "prepare":
        prepare_cache(arguments.source, arguments.output, force=arguments.force)
        return 0
    if arguments.command == "run":
        run_case_study(
            arguments.executable,
            arguments.source,
            arguments.cache,
            arguments.results,
            arguments.seed,
            arguments.force_prepare,
        )
        return 0
    return 0 if verify_results(arguments.results, arguments.reference) else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, KeyError, OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2) from error
