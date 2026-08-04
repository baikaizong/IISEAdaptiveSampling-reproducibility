"""Illustrate an anomaly with and without a smooth B-spline background.

The figure decomposes the observed field into a localized high-resolution
B-spline anomaly, independent Gaussian noise, and a low-resolution smooth
background. A fixed random seed makes the illustrative fields reproducible.
All heatmaps share one symmetric color scale so color has the same numerical
meaning in every panel.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
from matplotlib.patches import Circle
from scipy.interpolate import BSpline

# ---------------------------------------------------------------------------
# Tensor-product B-spline basis construction
# ---------------------------------------------------------------------------
def get_1d_bspline_matrix(n_control, degree, n_pixel=100):
    """Evaluate an open-uniform one-dimensional B-spline basis.

    Parameters
    ----------
    n_control : int
        Number of control coefficients and therefore basis functions.
    degree : int
        Polynomial degree of each B-spline basis function.
    n_pixel : int
        Number of evenly spaced evaluation points on the unit interval.

    Returns
    -------
    numpy.ndarray
        Matrix with shape ``(n_pixel, n_control)``. Small negative values from
        floating-point evaluation are clipped to zero.
    """
    n_internal = n_control - degree
    internal_knots = np.linspace(0, 1, n_internal + 1)
    knots = np.concatenate(([0] * degree, internal_knots, [1] * degree))
    x = np.linspace(0, 1, n_pixel)
    basis_mat = np.zeros((n_pixel, n_control))
    for i in range(n_control):
        c = np.zeros(n_control)
        c[i] = 1.0
        spl = BSpline(knots, c, degree)
        basis_mat[:, i] = spl(x)
    return np.maximum(basis_mat, 0)

def create_2d_basis_matrix(n_ctrl_x, n_ctrl_y, degree, n_pixel=100):
    """Create flattened tensor-product basis images on a square grid.

    Each output row is one outer product of a y-direction and x-direction
    basis vector. The row-major flattening convention is retained when the
    coefficient-weighted field is reshaped below.
    """
    Bx = get_1d_bspline_matrix(n_ctrl_x, degree, n_pixel)
    By = get_1d_bspline_matrix(n_ctrl_y, degree, n_pixel)
    n_basis = n_ctrl_x * n_ctrl_y
    n_sq = n_pixel * n_pixel
    basis_matrix = np.zeros((n_basis, n_sq))
    row_idx = 0
    for i in range(n_ctrl_y):
        for j in range(n_ctrl_x):
            basis_img = np.outer(By[:, i], Bx[:, j])
            basis_matrix[row_idx, :] = basis_img.flatten()
            row_idx += 1
    return basis_matrix

# ---------------------------------------------------------------------------
# Reproducible synthetic fields
# ---------------------------------------------------------------------------
np.random.seed(42)
N_PIXEL = 100

# A coarse quadratic basis models the slowly varying background, while a
# finer cubic basis localizes a single anomaly without pixel-level edges.
B1_mat = create_2d_basis_matrix(5, 5, degree=2, n_pixel=N_PIXEL)
B2_mat = create_2d_basis_matrix(20, 20, degree=3, n_pixel=N_PIXEL)

anom_idx = 250
coef_anom = np.zeros(400)
coef_anom[anom_idx] = 5.0
anomaly_img = (coef_anom @ B2_mat).reshape(N_PIXEL, N_PIXEL)

# Draw one background coefficient vector and one independent noise field.
coef_bg = np.random.normal(loc=0, scale=3.0, size=25)
bg_img = (coef_bg @ B1_mat).reshape(N_PIXEL, N_PIXEL)
noise_img = np.random.normal(loc=0, scale=1.0, size=(N_PIXEL, N_PIXEL))

# The five panels encode the additive decomposition explicitly:
# 1. Anomaly (pure anomaly signal)
# 2. Noise
# 3. Anomaly without background = anomaly + noise
# 4. Smoothing background = bg_img
# 5. Anomaly with background = anomaly + noise + bg

img1 = anomaly_img                        # Anomaly
img2 = noise_img                          # Noise
img3 = anomaly_img + noise_img            # Anomaly without background
img4 = bg_img                             # Smoothing background
img5 = anomaly_img + noise_img + bg_img   # Anomaly with background

# Use one zero-centered range across all panels to preserve visual comparability.
all_vals = np.concatenate([img1.flatten(), img2.flatten(), img3.flatten(),
                            img4.flatten(), img5.flatten()])
abs_max = np.max(np.abs(all_vals))
vmin_global, vmax_global = -abs_max, abs_max

# Convert the anomaly maximum from array coordinates to the unit-square axes.
y_idx, x_idx = np.unravel_index(np.argmax(anomaly_img), anomaly_img.shape)
anom_center_x = (x_idx + 0.5) / N_PIXEL
anom_center_y = (y_idx + 0.5) / N_PIXEL

# ---------------------------------------------------------------------------
# Axis styling helper
# ---------------------------------------------------------------------------
def format_ax(ax, title=None):
    """Apply consistent unit-square ticks, typography, and light panel borders."""
    if title:
        ax.set_title(title, fontsize=11, pad=8)
    ticks = [0, 0.5, 1.0]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(['0', '0.5', '1'])
    ax.set_yticklabels(['0', '0.5', '1'])
    ax.tick_params(axis='both', direction='out', length=4, width=0.8,
                   colors='black', labelsize=8, pad=2)
    for spine in ax.spines.values():
        spine.set_edgecolor('#B0B0B0')
        spine.set_linewidth(0.8)

# ---------------------------------------------------------------------------
# Figure assembly and export
# ---------------------------------------------------------------------------
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'mathtext.fontset': 'stix'
})

fig = plt.figure(figsize=(16, 4))

# Alternate equal-width image panels with narrow arithmetic-symbol columns.
ratios = [1, 0.12, 1, 0.12, 1, 0.12, 1, 0.12, 1]
gs = gridspec.GridSpec(1, 9, width_ratios=ratios, wspace=0.02,
                       left=0.05, right=0.90, top=0.88, bottom=0.12)

extent = [0, 1, 0, 1]
cmap = 'RdBu_r'

images = [img1, img2, img3, img4, img5]
titles = ['Anomaly', 'Noise', 'Anomaly without\nbackground',
          'Smoothing\nbackground', 'Anomaly with\nbackground']
symbols = ['+', '=', '+', '=']  # between images
sym_cols = [1, 3, 5, 7]         # gridspec column indices for symbols
img_cols = [0, 2, 4, 6, 8]      # gridspec column indices for images

im_last = None
for k, (img, title, col) in enumerate(zip(images, titles, img_cols)):
    ax = fig.add_subplot(gs[0, col])
    im = ax.imshow(img, extent=extent, cmap=cmap,
                   vmin=vmin_global, vmax=vmax_global, origin='lower')
    format_ax(ax, title=title)
    im_last = im

    # Mark every panel containing the anomaly, including the pure signal panel.
    if k in [0, 2, 4]:
        circle = Circle((anom_center_x, anom_center_y), radius=0.08,
                        edgecolor='#FF0000', facecolor='none',
                        linewidth=1.2, linestyle='-')
        ax.add_patch(circle)

# Render arithmetic symbols in their own axes so image geometry stays square.
for sym, col in zip(symbols, sym_cols):
    ax_sym = fig.add_subplot(gs[0, col])
    ax_sym.axis('off')
    ax_sym.text(0.5, 0.5, sym, fontsize=28, fontweight='bold',
                ha='center', va='center', transform=ax_sym.transAxes)

# A single colorbar documents the shared amplitude encoding.
cax = fig.add_axes([0.915, 0.12, 0.012, 0.76])
cb = plt.colorbar(im_last, cax=cax)
cb.set_label('Amplitude', fontsize=10)
cb.ax.tick_params(labelsize=9)
cb.outline.set_edgecolor('#B0B0B0')
cb.outline.set_linewidth(0.8)

output = Path(__file__).resolve().parents[1] / 'FinalVision' / 'Bspline_anomaly'
plt.savefig(output.with_suffix('.eps'),
            format='eps', dpi=300, bbox_inches='tight')
plt.savefig(output.with_suffix('.png'),
            format='png', dpi=600, bbox_inches='tight')
print(f"Saved: {output.with_suffix('.png')} and {output.with_suffix('.eps')}")
plt.close(fig)
