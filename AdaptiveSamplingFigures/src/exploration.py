"""Plot exploration-distance summaries across neighborhood sizes.

Each input contains repeated distance measurements for 100 candidate columns.
The figure uses the first 40 columns, applies the square-root transform used in
the manuscript, and plots the median with an interquartile band. Three insets
expose local differences that would be difficult to read at the full scale.
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import warnings
import gzip
from pathlib import Path
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from scipy.ndimage import gaussian_filter1d

# ==========================================
# Repository paths and input configuration
# ==========================================
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "ExternalData" / "exploration"
OUTPUT_DIR = PROJECT_ROOT / "FinalVision"


def resolve_data_file(name):
    """Resolve an input name, preferring its lossless gzip representation."""
    raw = DATA_DIR / name
    compressed = DATA_DIR / f"{name}.gz"
    return compressed if compressed.exists() else raw


file_paths = [
    resolve_data_file("Mmdis_K-2.0_P-0.95.txt"),
    resolve_data_file("Mmdis_K-5.0_P-0.95.txt"),
    resolve_data_file("Mmdis_K-10.0_P-0.95.txt"),
    resolve_data_file("Mmdis_K-20.0_P-0.95.txt"),
    resolve_data_file("Exploration_Random.txt"),
]

group_names = ['$K=2$', '$K=5$', '$K=10$', '$K=20$', 'Random']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

FILE_COLS = 100
KEEP_COLS = 40
ROWS = 10000
SMOOTH_SIGMA = 1.0


def smooth(arr, sigma=SMOOTH_SIGMA):
    """Apply light one-dimensional Gaussian smoothing along window size."""
    if sigma <= 0:
        return arr
    return gaussian_filter1d(arr, sigma=sigma)


# ==========================================
# Input parsing and robust summary statistics
# ==========================================
def get_stats_from_file(file_path):
    """Load one experiment matrix and return smoothed median and quartiles.

    Infinite values are treated as missing after reshaping. If the raw element
    count differs from the documented ``ROWS * FILE_COLS`` shape, complete rows
    are retained and a warning reports the inferred row count.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Missing required exploration data: {file_path}")
    else:
        print(f"[read] {file_path}")
        try:
            opener = gzip.open if str(file_path).endswith('.gz') else open
            with opener(file_path, 'rt') as stream:
                raw_data = np.fromstring(stream.read(), sep=' ')
            if raw_data.size != ROWS * FILE_COLS:
                actual_rows = raw_data.size // FILE_COLS
                print(f"[warning] Input size mismatch; using {actual_rows} complete rows")
                raw_matrix = raw_data[:actual_rows*FILE_COLS].reshape(actual_rows, FILE_COLS)
            else:
                raw_matrix = raw_data.reshape(ROWS, FILE_COLS)
        except Exception as e:
            print(f"[error] Failed to read {file_path}: {e}")
            return None, None, None

    raw_matrix = raw_matrix[:, :KEEP_COLS]
    raw_matrix[np.isinf(raw_matrix)] = np.nan

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', r'invalid value encountered in sqrt')
        raw_matrix = np.sqrt(raw_matrix)

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
        medians = np.nanmedian(raw_matrix, axis=0)
        q1s = np.nanpercentile(raw_matrix, 25, axis=0)
        q3s = np.nanpercentile(raw_matrix, 75, axis=0)

    medians = smooth(medians)
    q1s     = smooth(q1s)
    q3s     = smooth(q3s)

    return medians, q1s, q3s


# ==========================================
# Load all experimental groups in the display order defined above.
# ==========================================
all_stats = []
for fp in file_paths:
    med, q1, q3 = get_stats_from_file(fp)
    if med is not None:
        all_stats.append({'med': med, 'q1': q1, 'q3': q3})

if not all_stats:
    raise ValueError("No valid exploration data were loaded")

# ==========================================
# Main figure styling and summary curves
# ==========================================
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
plt.rcParams['font.size'] = 17
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.top']       = False
plt.rcParams['ytick.right']     = False
plt.rcParams['xtick.bottom']    = True
plt.rcParams['ytick.left']      = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['mathtext.fontset'] = 'stix'

# A wide aspect ratio supports three insets without excessive vertical whitespace.
fig, ax = plt.subplots(figsize=(14, 6.0), dpi=300)

for spine in ax.spines.values():
    spine.set_visible(True)
ax.tick_params(axis='both', which='both', direction='in',
               left=True, bottom=True, top=False, right=False,
               labelsize=13)

x_axis = np.arange(1, KEEP_COLS + 1)

for i, stats in enumerate(all_stats):
    label = group_names[i]
    color = colors[i % len(colors)]
    ax.plot(x_axis, stats['med'],
            label=label, color=color, linewidth=1.5, linestyle='-', zorder=10 + i)
    ax.fill_between(x_axis, stats['q1'], stats['q3'],
                    color=color, alpha=0.15, edgecolor=None, zorder=1)

# ==========================================
# Inset panels for local comparisons
# ==========================================
def add_inset_plot(ax, all_stats, x_axis, x_range, insert_pos, mark_locs):
    """Add a data-driven zoom panel and mark its source region.

    The vertical limits are computed from the pooled quartile envelopes inside
    ``x_range`` with 10% padding. This keeps every method visible without using
    independently chosen manual ranges.
    """
    zoom_start, zoom_end = x_range

    y_vals = []
    for stats in all_stats:
        mask = (x_axis >= zoom_start) & (x_axis <= zoom_end)
        if np.any(mask):
            y_vals.extend(stats['q3'][mask])
            y_vals.extend(stats['q1'][mask])

    if y_vals:
        y_min, y_max = np.min(y_vals), np.max(y_vals)
        pad = (y_max - y_min) * 0.1
        zoom_ylim = (y_min - pad, y_max + pad)
    else:
        zoom_ylim = (0, 1)

    axins = ax.inset_axes(insert_pos)
    axins.patch.set_facecolor('none')
    for spine in axins.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1.0)
    axins.tick_params(axis='both', which='both', direction='in',
                      left=True, bottom=True, top=False, right=False)

    for i, stats in enumerate(all_stats):
        color = colors[i % len(colors)]
        axins.plot(x_axis, stats['med'], color=color, linewidth=1.5)
        axins.fill_between(x_axis, stats['q1'], stats['q3'], color=color, alpha=0.15, edgecolor=None)

    axins.set_xlim(zoom_start, zoom_end)
    axins.set_ylim(zoom_ylim)
    axins.grid(True, linestyle=':', alpha=0.5, linewidth=0.5)
    # Slightly smaller inset labels prevent crowding at publication width.
    axins.tick_params(axis='both', which='major', labelsize=11)

    mark_inset(ax, axins, loc1=mark_locs[0], loc2=mark_locs[1],
               fc="none", ec="0.5", linestyle='--', linewidth=0.8)


# Axes-relative inset rectangles are tuned for the wide, compact main panel.
add_inset_plot(ax, all_stats, x_axis,
               x_range=(3, 7),
               insert_pos=[0.10, 0.50, 0.20, 0.46],
               mark_locs=(2, 3))

add_inset_plot(ax, all_stats, x_axis,
               x_range=(9, 13),
               insert_pos=[0.40, 0.27, 0.20, 0.46],
               mark_locs=(2, 4))

add_inset_plot(ax, all_stats, x_axis,
               x_range=(22, 25),
               insert_pos=[0.70, 0.20, 0.20, 0.46],
               mark_locs=(3, 4))
# Final labels, grid hierarchy, legend, and export
ax.set_xlabel('window size', fontsize=18)
ax.set_ylabel(r'$\sqrt{d_{\mathrm{void}}(k)}$', fontsize=18)
ax.set_xlim(1, KEEP_COLS)

ax.grid(which='major', linestyle='--', linewidth=0.5, alpha=0.6)
ax.minorticks_on()
ax.grid(which='minor', linestyle=':', linewidth=0.3, alpha=0.3)

legend = ax.legend(loc='upper right', fontsize=15, frameon=True, fancybox=False,
                   edgecolor='black', shadow=False)
legend.get_frame().set_facecolor('none')
legend.get_frame().set_linewidth(1.0)

# Compact padding preserves label clearance while minimizing unused margins.
plt.tight_layout(pad=0.8)

fig.patch.set_facecolor('none')
ax.patch.set_facecolor('none')

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color('black')
    spine.set_linewidth(1.0)

save_name = OUTPUT_DIR / 'exploration'
plt.savefig(save_name.with_suffix('.png'), dpi=600, bbox_inches='tight')
plt.savefig(save_name.with_suffix('.eps'), format='eps', bbox_inches='tight')

print(f"Saved: {save_name}.png and {save_name}.eps")
plt.close(fig)
