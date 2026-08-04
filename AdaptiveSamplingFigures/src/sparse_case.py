"""Illustrate exploration-to-exploitation transitions in a sparse case.

For six ``(xi, alpha)`` parameter combinations, the left panel shows the
objective over distance and the right panel tracks the maximizing distance as
``x`` changes. Shaded regimes and LaTeX-style threshold annotations identify
the transitions between global exploration and local exploitation.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.ticker as ticker

# Serif typography and STIX math provide a LaTeX-like publication appearance.
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.unicode_minus'] = False

# Numerical domain and model constants.
hs = 0.02
rho = 0.0005
d_fine = np.linspace(0, 0.1, 5000)
x_scan = np.linspace(0.1, 5.0, 500)
x_examples = [1.2, 1.5, 1.6, 1.7, 1.8, 2.0]

def get_f(d, x, xi, alpha):
    """Evaluate the sparse-case objective over candidate distance ``d``.

    The compact-support first term rewards nearby exploitation when
    ``d <= hs``. The logarithmic term supplies the global exploration trend.
    Vectorized masking avoids evaluating the compact term outside its support.
    """
    term1 = np.zeros_like(d)
    mask = d <= hs
    if np.any(mask):
        u_sq = (d[mask] / hs) ** 2
        core = (1 - u_sq) ** 2
        numerator = alpha * core * (x**2)
        denominator = core + xi
        term1[mask] = numerator / denominator
    term2 = np.log(d**2 + rho)
    return term1 + term2

# Paired-panel renderer used for each parameter combination.
def draw_pair(pair_axes, xi, alpha, label):
    """Draw objective curves and the corresponding optimal-distance path."""
    ax1, ax2 = pair_axes
    
    # Left panel: objective profiles and their discrete maximizing locations.
    for x_val in x_examples:
        y = get_f(d_fine, x_val, xi, alpha)
        idx = np.argmax(y)
        ax1.plot(d_fine, y, label=r'$x=%.1f$' % x_val)
        ax1.scatter(d_fine[idx], y[idx], s=20, zorder=5)
    
    ax1.set_xlabel(r'$d$', fontsize=10)
    ax1.set_ylabel(r'$\eta(d)$', fontsize=10)
    ax1.legend(fontsize=7, ncol=2, loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    # Right panel: maximizing distance as the state variable x increases.
    peak_locations = []
    for x_val in x_scan:
        y = get_f(d_fine, x_val, xi, alpha)
        peak_locations.append(d_fine[np.argmax(y)])
    peak_locations = np.array(peak_locations)

    idxL = np.where(peak_locations < 0.0999)[0]
    x_L = x_scan[idxL[0]] if len(idxL) > 0 else 0
    idxU = np.where(peak_locations < 0.0005)[0]
    x_U = x_scan[idxU[0]] if len(idxU) > 0 else 0

    ax2.plot(x_scan, peak_locations, color='black', linewidth=1.5, zorder=10)
    ax2.axvline(x_L, color='red', linestyle='--', alpha=0.6)
    ax2.axvline(x_U, color='red', linestyle='--', alpha=0.6)
    ax2.axvspan(x_scan[0], x_L, color='gray', alpha=0.1)
    ax2.axvspan(x_L, x_U, color='blue', alpha=0.05)
    ax2.axvspan(x_U, x_scan[-1], color='green', alpha=0.1)

    # Center vertical labels inside the three transition regimes.
    y_center = 0.05
    mid_x1 = (x_scan[0] + x_L) / 2
    mid_x2 = (x_L + x_U) / 2
    mid_x3 = (x_U + x_scan[-1]) / 2

    ax2.text(mid_x1, y_center, 'Prevailed by\nglobal exploration', 
             rotation=90, ha='center', va='center', fontsize=8, fontweight='bold')
    ax2.text(mid_x2, y_center, 'Exploitation\nin the neighborhood', 
             rotation=90, ha='center', va='center', fontsize=8, fontweight='bold')
    ax2.text(mid_x3, y_center, 'Continuous exploitation\nat the same location', 
             rotation=90, ha='center', va='center', fontsize=8, fontweight='bold')

    # Offset threshold labels from the transition line to avoid overlap.
    ax2.annotate(r'$x_{\mathrm{L}}$', xy=(x_L, 0.06), xytext=(x_L-0.4, 0.07),
                 arrowprops=dict(arrowstyle='->', color='red', lw=0.7), fontsize=10)
    ax2.annotate(r'$x_{\mathrm{U}}$', xy=(x_U, 0.01), xytext=(x_U+0.2, 0.03),
                 arrowprops=dict(arrowstyle='->', color='red', lw=0.7), fontsize=10)

    ax2.set_xlabel(r'$x$', fontsize=10)
    ax2.set_ylabel(r'$d_{\mathrm{opt}}$', fontsize=10)
    ax2.set_ylim(-0.01, 0.11)
    ax2.grid(True, alpha=0.2)

    # Place one shared caption below the paired panels using the panel label.
    shared_title = r'%s $\xi=%.2f, \alpha=%.1f$' % (label, xi, alpha)
    ax1.set_title(shared_title, x=1.1, y=-0.25, fontsize=10, fontweight='bold')

# Arrange six paired examples as a dense 3 x 4 panel grid.
param_sets = [
    (0.20, 1.0), (0.20, 2.0),
    (0.10, 1.0), (0.10, 2.0), 
    (0.05, 1.0), (0.05, 2.0)
]
labels = ['(a)', '(b)', '(c)', '(d)','(e)', '(f)']

fig, axes = plt.subplots(3, 4, figsize=(20, 10))
# Preserve clear axis labels while keeping related panel pairs visually close.
plt.subplots_adjust(hspace=0.35, wspace=0.25)

for i, (xi_val, alpha_val) in enumerate(param_sets):
    row = i // 2
    col_start = (i % 2) * 2
    draw_pair(axes[row, col_start : col_start + 2], xi_val, alpha_val, labels[i])

# Export the complete multi-panel figure in raster and vector formats.
save_name = Path(__file__).resolve().parents[1] / 'FinalVision' / 'sparsecase'
plt.savefig(save_name.with_suffix('.png'), dpi=600, bbox_inches='tight')
plt.savefig(save_name.with_suffix('.eps'), format='eps', dpi=600, bbox_inches='tight')

plt.close(fig)
