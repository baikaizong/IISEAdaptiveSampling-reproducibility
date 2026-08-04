"""Validate PNG resolution and EPS structure for the final figure set."""

from __future__ import annotations

from pathlib import Path
import argparse
import sys

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / "FinalVision"
MAX_EPS_BYTES = 5_000_000
EXPECTED = {
    "ABS",
    "Bspline_anomaly",
    "ConAnomaly",
    "Flowchart",
    "NAS",
    "NAS_log",
    "NAS_log_points",
    "Q",
    "SASAM",
    "SASAM_log",
    "SASAM_log_points",
    "TRAS",
    "TRAS_log",
    "TRAS_log_points",
    "TSSRP",
    "TSSRP_log",
    "TSSRP_log_points",
    "exploration",
    "frame",
    "mM_KSTD",
    "mM_KSTD_log",
    "mM_KSTD_log_points",
    "shape_figure_1",
    "shape_figure_2",
    "sparsecase",
}
KNOWN_MISSING_EPS: set[str] = set()


def main() -> int:
    """Validate every expected figure pair and return a shell-friendly status."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--strict", action="store_true", help="Treat documented missing sources as failures.")
    args = parser.parse_args()
    failures: list[str] = []
    warnings: list[str] = []
    for stem in sorted(EXPECTED):
        png = FINAL / f"{stem}.png"
        eps = FINAL / f"{stem}.eps"
        if not png.exists():
            failures.append(f"missing PNG: {png.name}")
        else:
            with Image.open(png) as image:
                width, height = image.size
            if max(width, height) < 2000:
                warnings.append(f"PNG below 2000 px on long edge: {png.name} ({width}x{height})")
        if not eps.exists():
            message = f"missing EPS: {eps.name} (source data unavailable)"
            if stem in KNOWN_MISSING_EPS and not args.strict:
                warnings.append(message)
            else:
                failures.append(message)
        elif not eps.read_bytes().startswith(b"%!PS-Adobe"):
            failures.append(f"invalid EPS header: {eps.name}")
        elif eps.stat().st_size > MAX_EPS_BYTES:
            failures.append(
                f"EPS above 5 MB: {eps.name} ({eps.stat().st_size / 1_000_000:.2f} MB)"
            )

    for message in warnings:
        print(f"[warning] {message}")
    for message in failures:
        print(f"[missing] {message}")
    print(f"Checked {len(EXPECTED)} figure stems: {len(failures)} failure(s), {len(warnings)} warning(s).")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
