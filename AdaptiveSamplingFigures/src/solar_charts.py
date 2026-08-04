"""Render solar monitoring statistics for five detection methods.

The script produces three publication variants for each method: the original
scale, a logarithmic scale with point markers, and a compact logarithmic scale.
All variants share the same display window and highlighted flare intervals so
method-to-method comparisons remain aligned in time.
"""

import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# ============================================================
# Repository-relative paths; no workstation-specific configuration is required.
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "ExternalData" / "solar_statistics"
OUTPUT_DIR = PROJECT_ROOT / "FinalVision"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================
# 2. Input files
# ============================================================

FILES = {
    "mM-KSTD": "mM-KST_ChartingStatistics-5.txt",
    "SASAM": "SASAM_ChartingStatistics-5.txt",
    "TRAS": "TOPR_ChartingStatistics-5.txt",
    "TSSRP": "TSSPR_ChartingStatistics-5.txt",
    "NAS": "NAS_ChartingStatistics-5.txt",
}


# ============================================================
# Method-specific control limits used by the monitoring overlays.
# ============================================================

CONTROL_LIMITS = {
    "mM-KSTD": 25,
    "SASAM": 25,
    "TRAS": 25,
    "TSSRP": 25,
    "NAS": 25,
}


# ============================================================
# 4. Display settings
# ============================================================

FULL_SERIES_FIRST_FRAME = 1

DISPLAY_START = 180
DISPLAY_END = 300

FLARE_INTERVALS = [
    (191, 194),
    (216, 237),
    (257, 258),
]

USE_LOG_FOR = {
    "mM-KSTD": False,
    "SASAM": False,
    "TRAS": False,
    "TSSRP": True,
    "NAS": False,
}

YLIMS = {
    "mM-KSTD": None,
    "SASAM": None,
    "TRAS": None,
    "TSSRP": None,
    "NAS": None,
}


# ============================================================
# 5. Style
# ============================================================

plt.rcParams.update({
    "font.family": "Times New Roman",
    "mathtext.fontset": "stix",
    "axes.linewidth": 0.9,
    "axes.edgecolor": "#333333",
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 10,
    "figure.dpi": 300,
    "savefig.dpi": 600,
})


# ============================================================
# 6. Functions
# ============================================================

def load_series(file_path):
    """Load one charting-statistics file as a flat floating-point series."""
    data = np.loadtxt(file_path)
    return np.asarray(data, dtype=float).ravel()


def get_display_data(y):
    """Crop a complete time series to the manuscript display-frame interval."""
    n = len(y)

    x_all = np.arange(
        FULL_SERIES_FIRST_FRAME,
        FULL_SERIES_FIRST_FRAME + n
    )

    mask = (x_all >= DISPLAY_START) & (x_all <= DISPLAY_END)

    return x_all[mask], y[mask]


def add_flare_regions(ax):
    """Shade the documented solar-flare intervals behind the data curve."""
    for start, end in FLARE_INTERVALS:
        ax.axvspan(
            start,
            end,
            color="#B58BEA",
            alpha=0.25,
            linewidth=0,
            zorder=0
        )


def style_axis(ax):
    """Apply the common low-ink axis treatment used by all solar charts."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.grid(
        axis="y",
        color="#D9DEE8",
        linewidth=0.6,
        alpha=0.85
    )

    ax.grid(axis="x", visible=False)

    ax.tick_params(
        direction="out",
        length=3.5,
        width=0.8
    )

    ax.set_facecolor("white")


def auto_ylim(y, control_limit, use_log=False):
    """Derive readable limits while retaining both data and control threshold."""
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]

    if use_log:
        lower = max(min(np.min(y) * 0.8, control_limit * 0.6), 1e-4)
        upper = max(np.max(y) * 1.25, control_limit * 2.0)
        return lower, upper

    lower = 0
    upper = max(np.max(y) * 1.15, control_limit * 1.8)
    return lower, upper


def plot_single_chart(method, y, save_path):
    """Render and export one original-scale monitoring chart."""
    x, y_display = get_display_data(y)

    fig, ax = plt.subplots(
        figsize=(7.2, 4.3),
        facecolor="white"
    )

    add_flare_regions(ax)

    ax.plot(
        x,
        y_display,
        color="black",
        linewidth=1.15,
        zorder=3
    )

    ax.axhline(
        CONTROL_LIMITS[method],
        color="#D62728",
        linewidth=1.15,
        linestyle=(0, (5, 4)),
        zorder=2
    )

    if USE_LOG_FOR[method]:
        ax.set_yscale("log")
        ax.set_ylabel("Charting statistic (log scale)")
    else:
        ax.set_ylabel("Charting statistic")

    ylim = YLIMS[method]
    if ylim is None:
        ylim = auto_ylim(
            y_display,
            CONTROL_LIMITS[method],
            use_log=USE_LOG_FOR[method]
        )

    ax.set_ylim(ylim)
    ax.set_xlim(DISPLAY_START, DISPLAY_END)

    ax.set_xlabel("Time (frame index)")

    ax.set_title(
        method,
        fontweight="bold",
        color="#003C8F" if method == "mM-KSTD" else "#222222",
        pad=8
    )

    style_axis(ax)

    legend_elements = [
        Line2D([0], [0], color="black", lw=1.2, label="Charting statistic"),
        Line2D([0], [0], color="#D62728", lw=1.2,
               linestyle=(0, (5, 4)), label="Control limit"),
        Patch(facecolor="#B58BEA", alpha=0.25,
              edgecolor="none", label="Flare interval"),
    ]

    ax.legend(
        handles=legend_elements,
        frameon=False,
        loc="upper right"
    )

    plt.tight_layout()

    fig.savefig(
        save_path,
        dpi=600,
        bbox_inches="tight",
        facecolor="white"
    )
    fig.savefig(
        Path(save_path).with_suffix(".eps"),
        format="eps",
        bbox_inches="tight",
        facecolor="white"
    )

    plt.close(fig)


# ============================================================
# 7. Main
# ============================================================

series = {}

for method, filename in FILES.items():
    path = os.path.join(DATA_DIR, filename)

    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    series[method] = load_series(path)

    print(
        f"{method}: length={len(series[method])}, "
        f"min={np.min(series[method]):.4f}, "
        f"max={np.max(series[method]):.4f}, "
        f"max_frame={FULL_SERIES_FIRST_FRAME + np.argmax(series[method])}"
    )


for method, y in series.items():
    save_name = method.replace("-", "_").replace(" ", "_") + ".png"
    save_path = os.path.join(OUTPUT_DIR, save_name)

    plot_single_chart(
        method=method,
        y=y,
        save_path=save_path
    )

    print(f"Saved: {save_path}")


print("\nAll corrected single charts have been generated.")


import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# ============================================================
# 1. Paths
# ============================================================

DATA_DIR = PROJECT_ROOT / "ExternalData" / "solar_statistics"
OUTPUT_DIR = PROJECT_ROOT / "FinalVision"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================
# 2. Input files
# ============================================================

FILES = {
    "mM-KSTD": "mM-KST_ChartingStatistics-5.txt",
    "SASAM": "SASAM_ChartingStatistics-5.txt",
    "TRAS": "TOPR_ChartingStatistics-5.txt",
    "TSSRP": "TSSPR_ChartingStatistics-5.txt",
    "NAS": "NAS_ChartingStatistics-5.txt",
}


# ============================================================
# 3. Control limits
# Fill in original-scale control limits here
# ============================================================

CONTROL_LIMITS = {
    "mM-KSTD": 23,
    "SASAM": 30,
    "TRAS": 25,
    "TSSRP": 80000,
    "NAS": 2.05,
}


# ============================================================
# 4. Display settings
# ============================================================

FULL_SERIES_FIRST_FRAME = 1

DISPLAY_START = 180
DISPLAY_END = 300

FLARE_INTERVALS = [
    (191, 194),
    (216, 237),
    (257, 258),
]

LOG_TYPE = "ln"   # "ln" or "log10"

YLIMS = {
    # Frozen from the original auto-ranging calculation so the Nature-style
    # redesign cannot change any method's vertical comparison range.
    "mM-KSTD": (1.5679919107668387, 7.06910567646332),
    "SASAM": (2.285592311348565, 7.236653027380592),
    "TRAS": (2.517231259836763, 6.358633975651549),
    "TSSRP": (5.891268032077202, 25.822795587957014),
    "NAS": (0.4449947647261698, 1.5287147370419463),
}


# ============================================================
# 5. Style
# ============================================================

plt.rcParams.update({
    "font.family": "Times New Roman",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "axes.linewidth": 0.9,
    "axes.edgecolor": "#333333",
    "axes.labelsize": 15,
    "axes.titlesize": 17,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 12,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

INK = "#252525"
IN_CONTROL = "#3F3F3F"
ALERT = "#D55E00"
LIMIT = "#9E2A2B"
FLARE = "#D9D2E9"


# ============================================================
# 6. Functions
# ============================================================

def load_series(file_path):
    data = np.loadtxt(file_path)
    return np.asarray(data, dtype=float).ravel()


def log_transform(y):
    y = np.asarray(y, dtype=float)

    if np.any(y <= 0):
        raise ValueError("Log transform requires all charting statistics to be positive.")

    if LOG_TYPE == "log10":
        return np.log10(y)
    if LOG_TYPE == "ln":
        return np.log(y)

    raise ValueError("LOG_TYPE must be either 'ln' or 'log10'.")


def log_control_limit(control_limit):
    if control_limit <= 0:
        raise ValueError("Control limit must be positive.")

    if LOG_TYPE == "log10":
        return np.log10(control_limit)

    return np.log(control_limit)


def get_display_data(y):
    n = len(y)

    x_all = np.arange(
        FULL_SERIES_FIRST_FRAME,
        FULL_SERIES_FIRST_FRAME + n
    )

    mask = (x_all >= DISPLAY_START) & (x_all <= DISPLAY_END)

    return x_all[mask], y[mask]


def add_flare_regions(ax):
    for start, end in FLARE_INTERVALS:
        ax.axvspan(
            start,
            end,
            color=FLARE,
            alpha=0.42,
            linewidth=0,
            zorder=0
        )


def style_axis(ax):
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#333333")
        spine.set_linewidth(0.9)
    ax.grid(
        axis="y",
        color="#D9DEE8",
        linewidth=0.6,
        alpha=0.85
    )

    ax.grid(axis="x", visible=False)

    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        top=True,
        right=True,
        bottom=True,
        left=True,
        labeltop=False,
        labelright=False,
        labelbottom=True,
        labelleft=True,
        length=4.0,
        width=0.8,
        colors="#333333",
    )

    ax.set_facecolor("white")
def auto_ylim(y_log, control_limit_log):

    y_log = np.asarray(y_log, dtype=float)
    y_log = y_log[np.isfinite(y_log)]

    ymin = min(np.min(y_log), control_limit_log)
    ymax = max(np.max(y_log), control_limit_log)

    yrange = ymax - ymin

    if yrange <= 0:
        yrange = 1.0

    lower = ymin - 0.05 * yrange
    upper = ymax + 0.4 * yrange

    return lower, upper


def plot_single_chart(method, y_raw, save_path):
    x, y_display_raw = get_display_data(y_raw)

    y_display_log = log_transform(y_display_raw)

    cl_raw = CONTROL_LIMITS[method]
    cl_log = log_control_limit(cl_raw)

    exceed = y_display_raw > cl_raw

    fig, ax = plt.subplots(
        figsize=(7.2, 4.3),
        facecolor="white"
    )

    add_flare_regions(ax)

    # Connecting line
    ax.plot(
        x,
        y_display_log,
        color=IN_CONTROL,
        linewidth=0.85,
        alpha=0.92,
        zorder=1
    )

    # In-control points
    ax.scatter(
        x[~exceed],
        y_display_log[~exceed],
        s=15,
        marker="o",
        facecolor=IN_CONTROL,
        edgecolor="white",
        linewidth=0.4,
        zorder=3
    )

    # Out-of-control points
    ax.scatter(
        x[exceed],
        y_display_log[exceed],
        s=23,
        marker="o",
        facecolor=ALERT,
        edgecolor="white",
        linewidth=0.45,
        zorder=5
    )

    # Control limit
    ax.axhline(
        cl_log,
        color=LIMIT,
        linewidth=1.0,
        linestyle=(0, (4, 3)),
        zorder=2
    )

    if LOG_TYPE == "log10":
        ax.set_ylabel(r"$\log_{10}(\mathrm{charting\ statistic})$")
    else:
        ax.set_ylabel(r"$\ln(\mathrm{charting\ statistic})$")

    ylim = YLIMS[method]

    if ylim is None:
        ylim = auto_ylim(y_display_log, cl_log)

    ax.set_ylim(ylim)
    ax.set_xlim(DISPLAY_START, DISPLAY_END)

    ax.set_xticks(np.arange(DISPLAY_START, DISPLAY_END + 1, 20))
    ax.set_xlabel("Time (frame index)")

    style_axis(ax)

    legend_elements = [
        Line2D(
            [0], [0],
            color=IN_CONTROL,
            lw=0.85,
            marker="o",
            markersize=5,
            markerfacecolor=IN_CONTROL,
            markeredgecolor="white",
            label="In-control statistic"
        ),
        Line2D(
            [0], [0],
            color="none",
            marker="o",
            markersize=6,
            markerfacecolor=ALERT,
            markeredgecolor="white",
            label="Out-of-control statistic"
        ),
        Line2D(
            [0], [0],
            color=LIMIT,
            lw=1.0,
            linestyle=(0, (4, 3)),
            label="Control limit"
        ),
        Patch(
            facecolor=FLARE,
            alpha=0.42,
            edgecolor="none",
            label="Flare interval"
        ),
    ]

    ax.legend(
        handles=legend_elements,
        loc="upper right",
        ncol=1,
        frameon=False,
        borderaxespad=0.55,
        handlelength=2.4,
        handletextpad=0.7,
        fontsize=12,
    )

    fig.tight_layout()

    fig.savefig(
        save_path,
        dpi=600,
        bbox_inches="tight",
        facecolor="white"
    )
    fig.savefig(
        Path(save_path).with_suffix(".eps"),
        format="eps",
        bbox_inches="tight",
        facecolor="white"
    )

    plt.close(fig)


# ============================================================
# 7. Main
# ============================================================

series = {}

for method, filename in FILES.items():
    path = os.path.join(DATA_DIR, filename)

    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    series[method] = load_series(path)

    print(
        f"{method}: length={len(series[method])}, "
        f"min={np.min(series[method]):.4f}, "
        f"max={np.max(series[method]):.4f}, "
        f"max_frame={FULL_SERIES_FIRST_FRAME + np.argmax(series[method])}"
    )


for method, y in series.items():
    save_name = method.replace("-", "_").replace(" ", "_") + "_log_points.png"
    save_path = os.path.join(OUTPUT_DIR, save_name)

    plot_single_chart(
        method=method,
        y_raw=y,
        save_path=save_path
    )

    print(f"Saved: {save_path}")


print("\nAll log-transformed charts with point markers have been generated.")


# The Nature-style context above applies only to the five *_log_points figures.
# Restore the retained compact log-figure typography before rendering *_log.
plt.rcParams.update({
    "font.family": "Times New Roman",
    "mathtext.fontset": "stix",
    "axes.linewidth": 0.9,
    "axes.edgecolor": "#333333",
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 10,
})


def plot_log_chart(method, y_raw, save_path):
    """Render the compact log-scale variant retained in FinalVision."""
    x, y_display_raw = get_display_data(y_raw)
    y_display_log = log_transform(y_display_raw)
    cl_log = log_control_limit(CONTROL_LIMITS[method])
    exceed = y_display_raw > CONTROL_LIMITS[method]

    fig, ax = plt.subplots(figsize=(7.2, 4.3), facecolor="white")
    add_flare_regions(ax)
    ax.plot(x, y_display_log, color="black", linewidth=0.9, alpha=0.75, zorder=1)
    ax.scatter(
        x[~exceed], y_display_log[~exceed], s=24, color="black",
        edgecolor="white", linewidth=0.45, zorder=3,
    )
    ax.scatter(
        x[exceed], y_display_log[exceed], s=34, color="#D62728",
        edgecolor="white", linewidth=0.55, zorder=4,
    )
    ax.axhline(cl_log, color="#D62728", linewidth=1.2, linestyle=(0, (5, 4)))
    ax.set_xlim(DISPLAY_START, DISPLAY_END)
    ax.set_ylim(YLIMS[method] or auto_ylim(y_display_log, cl_log))
    ax.set_xlabel("Time (frame index)")
    ax.set_ylabel(r"$\log(\mathrm{charting\ statistic})$")
    ax.set_title(
        method, fontweight="bold",
        color="#003C8F" if method == "mM-KSTD" else "#222222", pad=8,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#D9DEE8", linewidth=0.6, alpha=0.85)
    ax.grid(axis="x", visible=False)
    ax.legend(
        handles=[
            Line2D([0], [0], color="black", lw=0.9, marker="o", markersize=5,
                   markerfacecolor="black", markeredgecolor="white", label="Statistic"),
            Line2D([0], [0], color="none", marker="o", markersize=6,
                   markerfacecolor="#D62728", markeredgecolor="white", label="Exceeds limit"),
            Line2D([0], [0], color="#D62728", lw=1.2,
                   linestyle=(0, (5, 4)), label="Control limit"),
            Patch(facecolor="#B58BEA", alpha=0.25, edgecolor="none", label="Flare interval"),
        ],
        loc="upper right", frameon=False,
    )
    fig.tight_layout()
    fig.savefig(save_path, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(Path(save_path).with_suffix(".eps"), format="eps",
                bbox_inches="tight", facecolor="white")
    plt.close(fig)


for method, y in series.items():
    save_name = method.replace("-", "_").replace(" ", "_") + "_log.png"
    save_path = os.path.join(OUTPUT_DIR, save_name)
    plot_log_chart(method, y, save_path)
    print(f"Saved: {save_path}")
