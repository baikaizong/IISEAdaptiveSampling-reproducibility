"""Plot out-of-control run-length distributions by sampling budget.

The script reads ``OC_K-*_P-*_Q-{q}.txt`` files, selects configured anomaly
amplitude/radius groups, and draws grouped boxplots across sampling budget
``q``. A line through group means complements the quartile and whisker summary.
Serif typography and STIX math provide a compact publication-ready rendering.
"""

import os
from pathlib import Path
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colormaps
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple

# ============================================================
# Publication typography and drawing defaults
# ============================================================
plt.rcParams.update({
    # Serif text and STIX mathematical glyphs approximate LaTeX typography.
    "font.family":           "serif",
    "font.serif":            ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset":      "stix",
    "axes.titlesize":        18,
    "axes.labelsize":        16,
    "xtick.labelsize":       14,
    "ytick.labelsize":       14,
    "legend.fontsize":       14,
    "legend.title_fontsize": 14,
    # Antialias text and geometry; embed editable TrueType outlines where used.
    "text.antialiased":      True,
    "lines.antialiased":     True,
    "patch.antialiased":     True,
    "pdf.fonttype":          42,
    "ps.fonttype":           42,
    # High-contrast axes and consistent tick hierarchy.
    "axes.edgecolor":        "#000000",
    "axes.linewidth":        1.3,
    "xtick.color":           "#000000",
    "ytick.color":           "#000000",
    "xtick.major.width":     1.2,
    "ytick.major.width":     1.2,
    "text.color":            "#000000",
    # Transparent output allows placement on the manuscript page background.
    "figure.facecolor":      "none",
    "axes.facecolor":        "none",
    "savefig.facecolor":     "none",
    "savefig.transparent":   True,
})

# ============================================================
# User-adjustable input and group configuration
# ============================================================

# Resolve data from the repository, independent of the working directory.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "ExternalData" / "q_boxplot"

FILE_PATTERN = "OC_K-*_P-*_Q-*.txt"  # One file per sampling budget.

# Each entry selects one anomaly amplitude/radius combination. Additional
# groups are supported; colors extend deterministically from ``tab10``.
GROUPS = [
    # {"Amp": 0.5, "Radius": 0.02, "label": r"$\delta=0.5,\ r=0.02$"},
    # {"Amp": 0.5, "Radius": 0.05, "label": r"$\delta=0.5,\ r=0.05$"},
    {"Amp": 0.5, "Radius": 0.1,  "label": r"$\delta=0.5,\ r=0.10$"},
    # {"Amp": 1.0, "Radius": 0.02, "label": r"$\delta=1.0,\ r=0.02$"},
    {"Amp": 1.0, "Radius": 0.05,  "label": r"$\delta=1.0,\ r=0.05$"},
    # {"Amp": 1.0, "Radius": 0.1, "label": r"$\delta=1.0,\ r=0.10$"},
    {"Amp": 2.0, "Radius": 0.02,  "label": r"$\delta=2.0,\ r=0.02$"},
    # {"Amp": 2.0, "Radius": 0.05, "label": r"$\delta=2.0,\ r=0.05$"},
    # {"Amp": 2.0, "Radius": 0.1, "label": r"$\delta=2.0,\ r=0.10$"},
]

Q_VALUES = list(range(5, 105, 5))  # 5, 10, 15, ..., 100

# Restrained high-contrast base palette, extended only when required.
_BASE_COLORS = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52",
    "#8172B2", "#937860", "#DA8BC3", "#8C8C8C",
]
_tab10 = colormaps["tab10"]
BOX_COLORS = _BASE_COLORS + [_tab10(k % 10) for k in range(len(_BASE_COLORS), max(len(GROUPS), len(_BASE_COLORS)))]


# ============================================================
# Input parsing and statistical preparation
# ============================================================

def compute_offsets(n_groups: int, total_width: float = 0.75):
    """Return a box width and centered offsets for grouped boxplots.

    ``total_width`` is the fraction of one categorical x interval assigned to
    the complete group. Equal offsets keep all boxes centered on the sampling
    budget they represent.
    """
    bar_w   = total_width / n_groups
    half    = total_width / 2 - bar_w / 2
    offsets = [-half + i * bar_w for i in range(n_groups)]
    return bar_w, offsets


def parse_q_from_filename(fname: str):
    """Extract the numeric sampling budget following ``_Q-`` in a filename."""
    m = re.search(r"_Q-([\d.]+)", os.path.basename(fname))
    return float(m.group(1)) if m else None


def load_all_files(data_dir: str, file_pattern: str) -> dict:
    """Load all matching tables and return a mapping ``{q: DataFrame}``."""
    files = glob.glob(os.path.join(data_dir, file_pattern))
    data  = {}
    for f in sorted(files):
        q = parse_q_from_filename(f)
        if q is None:
            continue
        try:
            df = pd.read_csv(f, sep=r"\s+")
            df.columns = df.columns.str.strip()
            data[q] = df
        except Exception as exc:
            print(f"[warning] Could not read {f}: {exc}")
    return data


def extract_group_stats(all_data: dict, shift_val: float,
                        radius_val: float, q_values: list) -> pd.DataFrame:
    """Select one anomaly group and build boxplot-ready statistics.

    Whiskers follow Tukey's 1.5-IQR rule but are clipped to the recorded
    minimum and maximum. The returned mean is plotted as the overlaid trend.
    """
    records = []
    for q in q_values:
        df = all_data.get(q)
        if df is None or df.empty:
            continue
        # Match the canonical public column names used by every Q input table.
        sub = df[(df["Amp"] == shift_val) & (df["Radius"] == radius_val)]
        if sub.empty:
            continue
        row = sub.iloc[0]
        iqr = row["OC_Q3"] - row["OC_Q1"]
        records.append({
            "q":       q,
            "ARL":     row["OC_Mean"],
            "Std":     row["OC_Std"],
            "Q1":      row["OC_Q1"],
            "Med":     row["OC_Med"],
            "Q3":      row["OC_Q3"],
            "Max":     row["OC_Max"],
            "Min":     row["OC_Min"],
            "CalcMin": max(row["OC_Min"], row["OC_Q1"] - 1.5 * iqr),
            "CalcMax": min(row["OC_Max"], row["OC_Q3"] + 1.5 * iqr),
        })
    return pd.DataFrame(records)


# ============================================================
# Grouped boxplot renderer
# ============================================================

def draw_panel(ax, stats_list_per_group: list, q_values_present: list,
               title: str = "", xlabel: str = "Sampling Budget", ylabel: str = "ARL"):
    """Draw grouped distribution summaries and mean trajectories on one axis.

    Box position is categorical rather than numerical so equal visual spacing
    is preserved across the uniformly sampled ``q`` values. The mapping back
    to the numeric budget is carried by explicit tick labels.
    """

    n_groups        = len(stats_list_per_group)
    bar_w, offsets  = compute_offsets(n_groups)
    all_q           = sorted(q_values_present)
    q_to_idx        = {q: i for i, q in enumerate(all_q)}

    for grp_i, df in enumerate(stats_list_per_group):
        if df is None or df.empty:
            continue
        color  = BOX_COLORS[grp_i % len(BOX_COLORS)]
        offset = offsets[grp_i]

        bxp_stats, positions, line_x, line_y = [], [], [], []

        for _, row in df.iterrows():
            q = row["q"]
            if q not in q_to_idx:
                continue
            pos = q_to_idx[q] + offset
            positions.append(pos)
            line_x.append(pos)
            line_y.append(row["ARL"])
            bxp_stats.append({
                "med":    row["Med"],
                "q1":     row["Q1"],
                "q3":     row["Q3"],
                "whislo": row["CalcMin"],
                "whishi": row["CalcMax"],
                "mean":   row["ARL"],
            })

        if bxp_stats:
            box = ax.bxp(bxp_stats, positions=positions,
                         widths=bar_w * 0.85,
                         patch_artist=True, showfliers=False,
                         manage_ticks=False)
            for patch in box["boxes"]:
                patch.set(facecolor=color, alpha=0.75,
                          edgecolor="#000000", linewidth=1.4)
            for el in ["whiskers", "caps", "medians"]:
                c = "white" if el == "medians" else "#000000"
                plt.setp(box[el], color=c,
                         linewidth=2.0 if el == "medians" else 1.3)

        if line_x:
            ax.plot(line_x, line_y,
                    color=color, linewidth=2.2,
                    marker="o", markersize=5.5,
                    markeredgecolor="white", markeredgewidth=0.8,
                    label="_nolegend_")

    # Horizontal tick labels remain legible at the selected wide aspect ratio.
    ax.set_xticks(range(len(all_q)))
    ax.set_xticklabels([f"{int(q)}" for q in all_q],
                       rotation=0, ha="center")
    # Half-category side padding prevents the outer boxes from touching axes.
    ax.set_xlim(-0.65, len(all_q) - 0.35)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=10)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.8, alpha=0.6, color="#555555")
    ax.set_axisbelow(True)

    # A complete thin frame is retained for journal-style panel definition.
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.2)
        spine.set_color("#000000")

    # Keep ticks on the data-bearing left and bottom axes only.
    ax.tick_params(axis="both", which="major",
                   direction="in",
                   left=True,  bottom=True,
                   right=False, top=False,
                   length=5, width=1.2, color="#000000")


# ============================================================
# Main rendering workflow
# ============================================================

def main():
    """Load all Q-budget inputs, render the figure, and export PNG/EPS pairs."""
    print(f"[read] Loading Q-budget data from {DATA_DIR}")
    all_data = load_all_files(DATA_DIR, FILE_PATTERN)

    if not all_data:
        raise FileNotFoundError(f"Missing Q-boxplot inputs in {DATA_DIR}")

    found_qs = sorted(all_data.keys())
    print(
        f"[read] Found {len(found_qs)} files; "
        f"q range: {found_qs[0]} to {found_qs[-1]}"
    )
    q_list = [q for q in Q_VALUES if q in all_data]

    # Prepare one boxplot-statistics table per configured anomaly group.
    group_stats = []
    for g in GROUPS:
        stats = extract_group_stats(all_data, g["Amp"], g["Radius"], q_list)
        group_stats.append(stats)
        print(f"[read] Group {g['label']}: {len(stats)} records")

    # Scale width with category count while preserving a compact fixed height.
    fig_w = max(12, len(q_list) * 0.55)
    fig, ax = plt.subplots(figsize=(fig_w, 5.5))

    draw_panel(
        ax,
        stats_list_per_group=group_stats,
        q_values_present=q_list,
        #title="OC Performance vs. Sampling Budget",
        xlabel="Sampling Budget $q$",
        ylabel="ARL",
    )

    # Each legend entry combines the distribution box and mean-line encodings.
    legend_entries, labels = [], []
    for i, g in enumerate(GROUPS):
        color = BOX_COLORS[i % len(BOX_COLORS)]
        patch = mpatches.Patch(facecolor=color, alpha=0.6,
                               edgecolor="black", linewidth=1.0)
        line  = Line2D([0], [0], color=color, linewidth=2.0,
                       marker="o", markersize=5,
                       markeredgecolor="white")
        legend_entries.append((patch, line))
        labels.append(g["label"])

    ax.legend(
        legend_entries,
        labels,
        handler_map={tuple: HandlerTuple(ndivide=None)},
        framealpha=0.85,
        loc="upper right",
    )

    plt.tight_layout()
    out_path = PROJECT_ROOT / "FinalVision" / "Q"
    plt.savefig(out_path.with_suffix('.png'), dpi=600, bbox_inches="tight", transparent=True)
    plt.savefig(out_path.with_suffix('.eps'), format='eps', bbox_inches="tight", transparent=True)
    print(
        f"[saved] {out_path.with_suffix('.png')} and "
        f"{out_path.with_suffix('.eps')}"
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
