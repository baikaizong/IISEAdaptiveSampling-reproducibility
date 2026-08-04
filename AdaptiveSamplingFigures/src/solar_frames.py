"""Render the four solar frames and their temporal differences."""

from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_FILE = PROJECT_ROOT / "ExternalData" / "solar_frames" / "data.mat"
OUTPUT = PROJECT_ROOT / "FinalVision" / "frame"
TIMESTAMPS = (220, 221, 222, 223)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "axes.unicode_minus": False,
})


def load_frames() -> np.ndarray:
    """Load and validate the 232 x 292 x 300 MATLAB image stack.

    A transposed ``300 x 232 x 292`` representation is accepted for
    compatibility with MATLAB export conventions and normalized in memory.
    """
    if not DATA_FILE.exists():
        raise FileNotFoundError(f"Missing solar frame data: {DATA_FILE}")
    mat_data = loadmat(DATA_FILE)
    if "data" not in mat_data:
        raise KeyError(f"Expected variable 'data' in {DATA_FILE}")
    data = np.asarray(mat_data["data"], dtype=float)
    if data.shape == (300, 232, 292):
        data = np.transpose(data, (1, 2, 0))
    if data.shape != (232, 292, 300):
        raise ValueError(f"Unexpected solar data shape: {data.shape}")
    return data


def main() -> None:
    """Render four consecutive frames and their first temporal differences."""
    data = load_frames()
    differences = data[:, :, 1:] - data[:, :, :-1]
    original_frames = [data[:, :, timestamp] for timestamp in TIMESTAMPS]
    difference_frames = [differences[:, :, timestamp - 1] for timestamp in TIMESTAMPS]

    original_norm = mcolors.Normalize(
        vmin=np.percentile(original_frames, 1),
        vmax=np.percentile(original_frames, 99),
    )
    difference_norm = mcolors.Normalize(vmin=0, vmax=100)
    difference_cmap = mcolors.LinearSegmentedColormap.from_list(
        "difference_hot",
        [
            (0.00, 0.00, 0.00),
            (0.80, 0.00, 0.00),
            (1.00, 0.40, 0.00),
            (1.00, 1.00, 0.00),
            (1.00, 1.00, 1.00),
        ],
    )

    fig, axes = plt.subplots(2, 4, figsize=(15, 6), facecolor="white")
    for column, timestamp in enumerate(TIMESTAMPS):
        axes[0, column].imshow(original_frames[column], cmap="hot", norm=original_norm)
        axes[0, column].axis("off")
        axes[0, column].text(
            0.5, -0.05, rf"original frame at $t={timestamp}$",
            transform=axes[0, column].transAxes,
            fontsize=16, ha="center", va="top",
        )

        axes[1, column].imshow(
            difference_frames[column], cmap=difference_cmap, norm=difference_norm
        )
        axes[1, column].axis("off")
        axes[1, column].text(
            0.5, -0.05, rf"difference frame at $t={timestamp}$",
            transform=axes[1, column].transAxes,
            fontsize=16, ha="center", va="top",
        )

    fig.subplots_adjust(
        left=0.02, right=0.92, top=0.95, bottom=0.10, wspace=0.02, hspace=0.25
    )
    OUTPUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUTPUT.with_suffix(".eps"), format="eps", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {OUTPUT.with_suffix('.png')} and {OUTPUT.with_suffix('.eps')}")


if __name__ == "__main__":
    main()
