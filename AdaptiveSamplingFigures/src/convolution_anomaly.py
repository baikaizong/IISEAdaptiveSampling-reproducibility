"""Visualize the convolution-based anomaly observation model.

The single-row equation shows a smooth Gaussian background plus a localized
impulse anomaly, convolution with a normalized Gaussian kernel, addition of
measurement noise, and the resulting observed field. A shared robust color
normalization makes amplitudes comparable across all field panels.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Circle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy.signal import convolve2d

# ============================================================
# Central typography controls for every label and annotation in the figure.
# ============================================================
GLOBAL_FONT_SIZE = 28
SMALL_FONT_SIZE  = GLOBAL_FONT_SIZE
TITLE_FONT_SIZE  = GLOBAL_FONT_SIZE
SYMBOL_FONT_SIZE = GLOBAL_FONT_SIZE
LABEL_FONT_SIZE  = GLOBAL_FONT_SIZE
TICK_FONT_SIZE   = GLOBAL_FONT_SIZE
KERNEL_TEXT_SIZE = GLOBAL_FONT_SIZE

# ============================================================
# Model and grid parameters. The fixed seed preserves anomaly location and noise.
# ============================================================
np.random.seed(42)

num_x         = 20
windowsize    = 6
sigma_b       = 5.0
anom_strength = 100.0
sigma_obs     = 1.0
sigma_kernel  = 1.0

N_large = num_x + windowsize
N_small = num_x
N_win   = 1 + windowsize

# ============================================================
# Synthetic field construction
# ============================================================

# Gaussian background defined on the padded grid required by valid convolution.
bg_large = np.random.normal(0, sigma_b, (N_large, N_large))

# Place one impulse anomaly far enough from the boundary to survive valid cropping.
anom_large = np.zeros((N_large, N_large))
margin = N_win // 2 + 1
ry = np.random.randint(margin, N_large - margin)
rx = np.random.randint(margin, N_large - margin)
anom_large[ry, rx] = anom_strength

# Construct a normalized discrete Gaussian kernel; normalization preserves mass.
cy_k, cx_k = (N_win - 1) / 2, (N_win - 1) / 2
ys, xs = np.mgrid[0:N_win, 0:N_win]

kernel_raw = np.exp(
    -((xs - cx_k) ** 2 + (ys - cy_k) ** 2)
    / (2 * sigma_kernel ** 2)
)
kernel = kernel_raw / kernel_raw.sum()

# ``valid`` convolution reduces the padded field to the target 20 x 20 grid.
combined_large = bg_large + anom_large
img_conv = convolve2d(combined_large, kernel, mode='valid')

# Independent observation noise is added after the convolution step.
noise_obs = np.random.normal(0, sigma_obs, (N_small, N_small))

# Final observed field in the monitoring model.
img_final = img_conv + noise_obs

# ============================================================
# Robust shared color scale based on pooled 1st and 99th percentiles.
# ============================================================
all_data = np.concatenate([
    bg_large.flatten(),
    anom_large.flatten(),
    img_conv.flatten(),
    noise_obs.flatten(),
    img_final.flatten()
])

p_low  = np.percentile(all_data, 1)
p_high = np.percentile(all_data, 99)

abs_bound = max(abs(p_low), abs(p_high))
vmin_g, vmax_g = -abs_bound, abs_bound

# ============================================================
# Map array indices to the normalized [0, 1] coordinates used by ``imshow``.
# ============================================================
anom_cx_norm = (rx + 0.5) / N_large
anom_cy_norm = (ry + 0.5) / N_large

conv_rx = rx - (N_win // 2)
conv_ry = ry - (N_win // 2)

conv_cx_norm = np.clip((conv_rx + 0.5) / N_small, 0.05, 0.95)
conv_cy_norm = np.clip((conv_ry + 0.5) / N_small, 0.05, 0.95)

# ============================================================
# Publication typography and mathematical text rendering.
# ============================================================
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif', 'Georgia', 'serif'],
    'mathtext.fontset': 'stix',
    'axes.unicode_minus': False,
    'font.size': GLOBAL_FONT_SIZE,
    'axes.titlesize': GLOBAL_FONT_SIZE,
    'axes.labelsize': GLOBAL_FONT_SIZE,
    'xtick.labelsize': GLOBAL_FONT_SIZE,
    'ytick.labelsize': GLOBAL_FONT_SIZE,
    'legend.fontsize': GLOBAL_FONT_SIZE,
})

cmap = 'RdBu_r'
extent = [0, 1, 0, 1]

# ============================================================
# Reusable panel helpers
# ============================================================

def format_ax(ax, title=None):
    """Apply consistent titles, normalized ticks, and unobtrusive borders."""
    if title:
        ax.set_title(title, fontsize=TITLE_FONT_SIZE, pad=12)

    ax.set_xticks([0, 1.0])
    ax.set_yticks([0, 1.0])

    ax.set_xticklabels(['0', '1'], fontsize=TICK_FONT_SIZE)
    ax.set_yticklabels(['0', '1'], fontsize=TICK_FONT_SIZE)

    ax.tick_params(
        axis='both',
        direction='out',
        length=4,
        width=0.9,
        colors='black',
        pad=4
    )

    for sp in ax.spines.values():
        sp.set_edgecolor('#B0B0B0')
        sp.set_linewidth(0.8)


def add_size_label(ax, text):
    """Place a centered grid-size or parameter annotation below a panel."""
    ax.text(
        0.5, -0.22, text,
        transform=ax.transAxes,
        fontsize=LABEL_FONT_SIZE,
        ha='center',
        va='top',
        color='#444444',
        style='italic'
    )


def sym_ax(ax, sym):
    """Render an arithmetic operator in a dedicated layout column."""
    ax.axis('off')
    ax.text(
        0.5, 0.5, sym,
        fontsize=SYMBOL_FONT_SIZE,
        fontweight='bold',
        ha='center',
        va='center',
        transform=ax.transAxes
    )


def draw_kernel_panel(ax, kernel, title=None):
    """Render the Gaussian kernel as a labeled matrix.

    Cell luminance encodes weight magnitude, while printed values expose the
    exact normalized coefficients. Text switches from dark to white when the
    cell exceeds 55% of the maximum to maintain contrast.
    """
    n = kernel.shape[0]

    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_aspect('equal')
    ax.axis('off')

    if title:
        ax.set_title(title, fontsize=TITLE_FONT_SIZE, pad=12)

    norm = Normalize(vmin=0, vmax=kernel.max())
    cmap_k = plt.get_cmap('Blues')

    for row in range(n):
        for col in range(n):
            val = kernel[row, col]
            color = cmap_k(norm(val))

            rect = plt.Rectangle(
                (col - 0.5, row - 0.5),
                1.0,
                1.0,
                facecolor=color,
                edgecolor='white',
                linewidth=0.8
            )
            ax.add_patch(rect)

            txt_color = 'white' if val > kernel.max() * 0.55 else '#333333'

            ax.text(
                col, row, f'{val:.3f}',
                ha='center',
                va='center',
                fontsize=KERNEL_TEXT_SIZE * 0.45,
                color=txt_color,
                fontfamily='serif'
            )

    ax.text(
        n / 2 - 0.5,
        -1.25,
        f'$(1+{windowsize})\\times(1+{windowsize})$,  $\\sigma_k={sigma_kernel}$',
        ha='center',
        va='top',
        fontsize=LABEL_FONT_SIZE,
        style='italic',
        color='#444444',
        transform=ax.transData
    )

# ============================================================
# One-row equation layout. Narrow columns contain operators and parentheses.
# ============================================================
col_ratios = [
    0.08, 1.0, 0.15, 1.0, 0.08, 0.15,
    0.85, 0.15, 1.0, 0.15, 1.0, 0.15, 1.0
]

fig = plt.figure(figsize=(30, 5.6))

gs = gridspec.GridSpec(
    1, 13,
    width_ratios=col_ratios,
    left=0.02,
    right=0.915,
    top=0.86,
    bottom=0.25,
    wspace=0.04
)

# Opening parenthesis for the pre-convolution sum.
ax_lp = fig.add_subplot(gs[0, 0])
sym_ax(ax_lp, '(')

# Background field.
ax1 = fig.add_subplot(gs[0, 1])
ax1.imshow(
    bg_large,
    extent=extent,
    cmap=cmap,
    vmin=vmin_g,
    vmax=vmax_g,
    origin='lower'
)
format_ax(ax1, title='Background noise')
add_size_label(
    ax1,
    f'$({num_x}+{windowsize})\\times({num_x}+{windowsize})$,  $\\sigma_b={sigma_b}$'
)

# Addition operator between background and anomaly.
ax_p1 = fig.add_subplot(gs[0, 2])
sym_ax(ax_p1, '+')

# Sparse anomaly field with a red location marker.
ax2 = fig.add_subplot(gs[0, 3])
ax2.imshow(
    anom_large,
    extent=extent,
    cmap=cmap,
    vmin=vmin_g,
    vmax=vmax_g,
    origin='lower'
)
format_ax(ax2, title='Anomaly')
add_size_label(
    ax2,
    f'$({num_x}+{windowsize})\\times({num_x}+{windowsize})$,  strength $={anom_strength}$'
)
ax2.add_patch(
    Circle(
        (anom_cx_norm, anom_cy_norm),
        radius=0.07,
        edgecolor='#EE0000',
        facecolor='none',
        linewidth=1.8
    )
)

# Closing parenthesis for the pre-convolution sum.
ax_rp = fig.add_subplot(gs[0, 4])
sym_ax(ax_rp, ')')

# Convolution operator.
ax_cs = fig.add_subplot(gs[0, 5])
sym_ax(ax_cs, '$\\circledast$')

# Explicit Gaussian-kernel matrix.
ax3 = fig.add_subplot(gs[0, 6])
draw_kernel_panel(ax3, kernel, title='Gaussian kernel')

# Equality after convolution.
ax_eq1 = fig.add_subplot(gs[0, 7])
sym_ax(ax_eq1, '=')

# Convolved signal with the transformed anomaly location marked.
ax4 = fig.add_subplot(gs[0, 8])
ax4.imshow(
    img_conv,
    extent=extent,
    cmap=cmap,
    vmin=vmin_g,
    vmax=vmax_g,
    origin='lower'
)
format_ax(ax4, title='Convolved signal')
add_size_label(ax4, f'${num_x}\\times{num_x}$')
ax4.add_patch(
    Circle(
        (conv_cx_norm, conv_cy_norm),
        radius=0.12,
        edgecolor='#EE0000',
        facecolor='none',
        linewidth=1.8
    )
)

# Addition operator before observation noise.
ax_p2 = fig.add_subplot(gs[0, 9])
sym_ax(ax_p2, '+')

# Independent observation-noise field.
ax5 = fig.add_subplot(gs[0, 10])
ax5.imshow(
    noise_obs,
    extent=extent,
    cmap=cmap,
    vmin=vmin_g,
    vmax=vmax_g,
    origin='lower'
)
format_ax(ax5, title='Noise')
add_size_label(
    ax5,
    f'${num_x}\\times{num_x}$,  $\\sigma_{{\\mathrm{{obs}}}}={sigma_obs}$'
)

# Equality before the final observation.
ax_eq2 = fig.add_subplot(gs[0, 11])
sym_ax(ax_eq2, '=')

# Final observed field.
ax6 = fig.add_subplot(gs[0, 12])
ax6.imshow(
    img_final,
    extent=extent,
    cmap=cmap,
    vmin=vmin_g,
    vmax=vmax_g,
    origin='lower'
)
format_ax(ax6, title='Final observation')
add_size_label(ax6, f'${num_x}\\times{num_x}$')
ax6.add_patch(
    Circle(
        (conv_cx_norm, conv_cy_norm),
        radius=0.12,
        edgecolor='#EE0000',
        facecolor='none',
        linewidth=1.8
    )
)

# ============================================================
# Shared colorbar for every diverging field panel.
# ============================================================
cax = fig.add_axes([0.926, 0.25, 0.010, 0.61])

sm = ScalarMappable(
    cmap=cmap,
    norm=Normalize(vmin=vmin_g, vmax=vmax_g)
)
sm.set_array([])

cb = plt.colorbar(sm, cax=cax)
cb.set_label('Amplitude', fontsize=GLOBAL_FONT_SIZE)
cb.ax.tick_params(labelsize=GLOBAL_FONT_SIZE)

cb.outline.set_edgecolor('#B0B0B0')
cb.outline.set_linewidth(0.8)

# ============================================================
# Export a high-resolution preview and a publication-oriented EPS artifact.
# ============================================================
output = Path(__file__).resolve().parents[1] / 'FinalVision' / 'ConAnomaly'
plt.savefig(
    output.with_suffix('.png'),
    format='png',
    dpi=600,
    bbox_inches='tight'
)

plt.savefig(
    output.with_suffix('.eps'),
    format='eps',
    dpi=300,
    bbox_inches='tight'
)

print(f"Saved: {output.with_suffix('.png')} and {output.with_suffix('.eps')}")
plt.close(fig)
