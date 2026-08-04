"""Regenerate every figure whose complete source data is stored in this project."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / "FinalVision"

ALWAYS_REPRODUCIBLE = [
    "solar_charts.py",
    "bspline_anomaly.py",
    "convolution_anomaly.py",
    "sparse_case.py",
    "monitoring_area.py",
]

EXTERNAL = ROOT / "ExternalData"

CONDITIONAL = {
    "solar_frames.py": [EXTERNAL / "solar_frames" / "data.mat"],
    "exploration.py": [
        EXTERNAL / "exploration" / "Mmdis_K-2.0_P-0.95.txt",
        EXTERNAL / "exploration" / "Mmdis_K-5.0_P-0.95.txt",
        EXTERNAL / "exploration" / "Mmdis_K-10.0_P-0.95.txt",
        EXTERNAL / "exploration" / "Mmdis_K-20.0_P-0.95.txt",
        EXTERNAL / "exploration" / "Exploration_Random.txt",
    ],
    "shape_figures.py": [
        EXTERNAL / "shape" / "mM_KSTD_ellipse.txt",
        EXTERNAL / "shape" / "SASAM_ellipse.txt",
        EXTERNAL / "shape" / "mM_KSTD_crescent.txt",
        EXTERNAL / "shape" / "SASAM_crescent.txt",
    ],
    "q_boxplot.py": [
        EXTERNAL / "q_boxplot" / f"OC_K-10.0_P-0.90_Q-{q}.txt"
        for q in range(5, 105, 5)
    ],
}


def input_exists(path: Path) -> bool:
    """Large text inputs may be kept raw or losslessly gzip-compressed."""
    return path.exists() or Path(f"{path}.gz").exists()


def run_script(name: str) -> None:
    """Run one renderer with a deterministic headless Matplotlib environment."""
    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    env.setdefault("MPLCONFIGDIR", "/tmp/adaptive-sampling-mplconfig")
    env.setdefault("XDG_CACHE_HOME", "/tmp/adaptive-sampling-xdg-cache")
    print(f"[run] {name}")
    subprocess.run([sys.executable, ROOT / "src" / name], check=True, env=env)


def render_flowchart_png() -> None:
    """Rasterize the vector flowchart at 600 dpi when Ghostscript is available."""
    eps = FINAL / "Flowchart.eps"
    png = FINAL / "Flowchart.png"
    gs = shutil.which("gs")
    if not eps.exists():
        print("[skip] Flowchart.png: missing FinalVision/Flowchart.eps")
        return
    if not gs:
        print("[skip] Flowchart.png: Ghostscript is not installed")
        return
    subprocess.run(
        [
            gs,
            "-dSAFER",
            "-dBATCH",
            "-dNOPAUSE",
            "-dEPSCrop",
            "-q",
            "-sDEVICE=pngalpha",
            "-r600",
            f"-sOutputFile={png}",
            str(eps),
        ],
        check=True,
    )
    print("[render] Flowchart.png")


def restore_source_assets() -> None:
    """Restore immutable source artwork before derived artifacts are rebuilt."""
    for name in ("Flowchart.eps", "Q.png"):
        source = ROOT / "assets" / name
        if source.exists():
            shutil.copy2(source, FINAL / name)


def main() -> None:
    FINAL.mkdir(exist_ok=True)
    restore_source_assets()
    for script in ALWAYS_REPRODUCIBLE:
        run_script(script)

    for script, required in CONDITIONAL.items():
        missing = [path for path in required if not input_exists(path)]
        if missing:
            print(f"[skip] {script}: missing {', '.join(str(p.relative_to(ROOT)) for p in missing)}")
            continue
        run_script(script)

    render_flowchart_png()
    run_script("optimize_eps.py")
    if not input_exists(CONDITIONAL["solar_frames.py"][0]):
        print("[note] frame cannot be regenerated until ExternalData/solar_frames/data.mat is supplied.")
    if any(not input_exists(path) for path in CONDITIONAL["q_boxplot.py"]):
        print("[note] Q.eps cannot be regenerated until ExternalData/q_boxplot is supplied.")


if __name__ == "__main__":
    main()
