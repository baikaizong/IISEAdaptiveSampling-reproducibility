"""Render monitoring-area fields before and after an anomaly emerges.

The two-panel figure combines a reproducible spatial noise field, a smooth
near-elliptical anomaly region, regular sampling locations, and publication
annotations. Both panels use one percentile-based normalization so differences
reflect the data rather than independent color rescaling.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.patheffects as path_effects
from matplotlib import rcParams
from matplotlib.colors import Normalize

# ── Times New Roman everywhere ─────────────────────────────────────────────
rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "axes.labelsize": 15,
    "xtick.labelsize": 12.5,
    "ytick.labelsize": 12.5,
})

# ── Color and style settings ──────────────────────────────────────────────
DOMAIN_EDGE = "#1f4e79"

# Sampling nodes: high contrast against blue background.
GRID_FACE = "#ffb000"
GRID_EDGE = "#111111"

ANOM_E = "#0b2c4d"
BBOX_C = "#4a4a8a"
OBS_C = "#0b2c4d"

TEXT_C = "#111111"
ANN_C = "#111111"
ANN_FS = 13.5

# ── Intensity controls ────────────────────────────────────────────────────
NOISE_LEVEL = 1.0
ANOMALY_MAIN = 3.0
ANOMALY_WAVE = 0.50
ANOMALY_LOCAL = 0.20

# Adaptive colorbar range.
CBAR_LOW_Q = 2.0
CBAR_HIGH_Q = 98.0

rng = np.random.default_rng(2026)


# ── Helper: point-in-polygon ───────────────────────────────────────────────
def pip(px, py, poly_x, poly_y):
    """Classify points with the vectorized ray-crossing polygon test.

    Boundary coordinates ``poly_x`` and ``poly_y`` describe a closed polygon;
    ``px`` and ``py`` may contain many query points. The parity accumulator
    toggles whenever a horizontal ray crosses an edge.
    """
    n = len(poly_x)
    inside = np.zeros(len(px), dtype=bool)
    j = n - 1

    for i in range(n):
        xi, yi = poly_x[i], poly_y[i]
        xj, yj = poly_x[j], poly_y[j]

        cond = ((yi > py) != (yj > py)) & (
            px < (xj - xi) * (py - yi) / (yj - yi + 1e-12) + xi
        )
        inside ^= cond
        j = i

    return inside


# ── Spatial grid for scalar fields ────────────────────────────────────────
n_field = 560
x = np.linspace(0, 1, n_field)
y = np.linspace(0, 1, n_field)
X, Y = np.meshgrid(x, y)

# ── Smooth anomaly region A: near-elliptical shape ─────────────────────────
cx, cy = 0.50, 0.50
theta = np.linspace(0, 2 * np.pi, 600)

# Ellipse axes.
a = 0.145
b = 0.115

# Low-frequency, low-amplitude deformation only.
smooth_deform = (
    1.0
    + 0.030 * np.cos(2 * theta + 0.35)
    + 0.018 * np.sin(3 * theta - 0.20)
)

ax_pts = cx + a * smooth_deform * np.cos(theta)
ay_pts = cy + b * smooth_deform * np.sin(theta)

inside_field = pip(
    X.ravel(),
    Y.ravel(),
    ax_pts,
    ay_pts
).reshape(X.shape)

# ── Fully independent noise over the whole monitoring area ────────────────
global_noise = rng.normal(0, 1, size=(n_field, n_field))
global_noise = (global_noise - global_noise.mean()) / global_noise.std()

background_field = 0.30 + NOISE_LEVEL * global_noise

# ── Smooth anomaly signal for t = tau + 1 ─────────────────────────────────
smooth_bump = np.exp(
    -(
        ((X - cx) ** 2) / (2 * 0.090 ** 2)
        + ((Y - cy) ** 2) / (2 * 0.072 ** 2)
    )
)

smooth_wave = (
    0.12 * np.sin(7.0 * X + 5.0 * Y)
    + 0.09 * np.cos(6.5 * X - 7.5 * Y)
    + 0.06 * np.sin(10.0 * np.sqrt((X - cx) ** 2 + (Y - cy) ** 2))
)

local_anomaly_variation = rng.normal(0, 1, size=(n_field, n_field))
for _ in range(20):
    local_anomaly_variation = (
        local_anomaly_variation
        + np.roll(local_anomaly_variation, 1, axis=0)
        + np.roll(local_anomaly_variation, -1, axis=0)
        + np.roll(local_anomaly_variation, 1, axis=1)
        + np.roll(local_anomaly_variation, -1, axis=1)
    ) / 5.0

local_anomaly_variation = (
    local_anomaly_variation - local_anomaly_variation.mean()
) / local_anomaly_variation.std()

anomaly_signal = (
    ANOMALY_MAIN * smooth_bump
    + ANOMALY_WAVE * smooth_wave
    + ANOMALY_LOCAL * local_anomaly_variation
)

# ── Two fields: before and after anomaly emerges ──────────────────────────
Z_tau = background_field.copy()

Z_tau_plus_1 = background_field.copy()
Z_tau_plus_1[inside_field] += anomaly_signal[inside_field]

# ── Adaptive shared normalization for both subplots ───────────────────────
Z_all = np.concatenate([
    Z_tau[np.isfinite(Z_tau)].ravel(),
    Z_tau_plus_1[np.isfinite(Z_tau_plus_1)].ravel()
])

vmin = np.percentile(Z_all, CBAR_LOW_Q)
vmax = np.percentile(Z_all, CBAR_HIGH_Q)

# Avoid degenerate color scale.
if np.isclose(vmin, vmax):
    vmin = np.min(Z_all)
    vmax = np.max(Z_all)

# Small padding; range remains adaptive.
pad = 0.015 * (vmax - vmin)
vmin_adapt = vmin - pad
vmax_adapt = vmax + pad

norm = Normalize(vmin=vmin_adapt, vmax=vmax_adapt)
cmap = "Blues"

# ── Regular grid sampling locations R ─────────────────────────────────────
grid_vals = np.linspace(0.05, 0.95, 19)
gx, gy = np.meshgrid(grid_vals, grid_vals)
rx = gx.ravel()
ry = gy.ravel()

inside_mask = pip(rx, ry, ax_pts, ay_pts)
in_rx = rx[inside_mask]
in_ry = ry[inside_mask]

# ── Figure layout: two subplots in one row ────────────────────────────────
fig, axes = plt.subplots(
    1, 2,
    figsize=(13.2, 6),
    dpi=180,
    constrained_layout=True
)

fig.patch.set_alpha(0)


def add_sampling_legend(ax, h_grid):
    legend = ax.legend(
        handles=[h_grid],
        loc="upper right",
        fontsize=11.5,
        framealpha=1.0,
        edgecolor="#111111",
        facecolor="white",
        frameon=True,
        borderpad=0.75,
        labelspacing=0.5,
        handlelength=1.7
    )

    for text in legend.get_texts():
        text.set_color(TEXT_C)
        text.set_fontweight("semibold")

    legend.get_frame().set_linewidth(1.0)
    legend.get_frame().set_alpha(1.0)
    legend.set_zorder(1000)

    return legend


def draw_base_panel(ax, Z, panel_title):
    ax.set_facecolor("none")

    im = ax.imshow(
        Z,
        extent=(0, 1, 0, 1),
        origin="lower",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        alpha=0.96,
        zorder=0
    )

    ax.grid(False)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xticks(np.arange(0, 1.01, 0.1))
    ax.set_yticks(np.arange(0, 1.01, 0.1))

    ax.set_xlabel("$r_1$", fontsize=16, color=TEXT_C)
    ax.set_ylabel("$r_2$", fontsize=16, color=TEXT_C)
    ax.tick_params(colors=TEXT_C, labelsize=13)

    for spine in ax.spines.values():
        spine.set_edgecolor(TEXT_C)
        spine.set_linewidth(1.2)

    domain_edge = plt.Rectangle(
        (0, 0), 1, 1,
        facecolor="none",
        edgecolor=DOMAIN_EDGE,
        linewidth=2.0,
        zorder=4
    )
    ax.add_patch(domain_edge)

    h_grid = ax.scatter(
        rx, ry,
        s=19,
        facecolor=GRID_FACE,
        edgecolor=GRID_EDGE,
        linewidth=0.55,
        alpha=1.0,
        zorder=8,
        label="sampling candidate nodes"
    )

    title = ax.set_title(
        panel_title,
        fontsize=17,
        color=TEXT_C,
        pad=10,
        fontweight="semibold"
    )
    title.set_path_effects([
        path_effects.Stroke(linewidth=3.0, foreground="white"),
        path_effects.Normal()
    ])

    add_sampling_legend(ax, h_grid)

    return im, h_grid


# ── Left panel: t = tau ───────────────────────────────────────────────────
im0, h_grid0 = draw_base_panel(
    axes[0],
    Z_tau,
    r"$t=\tau$"
)

# ── Right panel: t = tau + 1 ──────────────────────────────────────────────
im1, h_grid1 = draw_base_panel(
    axes[1],
    Z_tau_plus_1,
    r"$t=\tau+1$"
)

ax = axes[1]

# ── Light and smooth anomaly boundary ─────────────────────────────────────
h_anom = plt.Polygon(
    list(zip(ax_pts, ay_pts)),
    facecolor="none",
    edgecolor=ANOM_E,
    linewidth=1.65,
    alpha=0.96,
    zorder=18
)
ax.add_patch(h_anom)

anom_text = ax.text(
    cx, cy + 0.005,
    "Anomaly region $\\mathcal{A}$",
    fontsize=13.5,
    color="#111111",
    ha="center",
    va="center",
    fontweight="bold",
    zorder=30
)
anom_text.set_path_effects([
    path_effects.Stroke(linewidth=3.8, foreground="white"),
    path_effects.Normal()
])

# Highlight sampling nodes inside anomaly region.
ax.scatter(
    in_rx,
    in_ry,
    s=33,
    color="white",
    edgecolors=OBS_C,
    linewidths=1.9,
    zorder=20
)

# ── Annotation style helper ───────────────────────────────────────────────
ann_kw = dict(
    fontsize=ANN_FS,
    color=ANN_C,
    bbox=dict(
        boxstyle="round,pad=0.34",
        fc="white",
        ec="#111111",
        lw=0.9,
        alpha=1.0
    ),
    zorder=35
)

arrow_kw = dict(
    arrowstyle="->",
    color="#111111",
    lw=1.35,
    shrinkA=2,
    shrinkB=2
)

# ── (i) Spatial concentration: arrow points to anomaly boundary ───────────
margin = 0.012
xmin_b = ax_pts.min() - margin
xmax_b = ax_pts.max() + margin
ymin_b = ay_pts.min() - margin
ymax_b = ay_pts.max() + margin

bbox_rect = plt.Rectangle(
    (xmin_b, ymin_b),
    xmax_b - xmin_b,
    ymax_b - ymin_b,
    facecolor="none",
    edgecolor=BBOX_C,
    linewidth=1.65,
    linestyle=":",
    zorder=19
)
ax.add_patch(bbox_rect)

edge_idx = np.argmin(ax_pts + ay_pts)
edge_target = (ax_pts[edge_idx], ay_pts[edge_idx])

ann1 = ax.annotate(
    "(i) Spatial concentration",
    xy=edge_target,
    xytext=(0.07, 0.22),
    ha="left",
    arrowprops=dict(
        **arrow_kw,
        connectionstyle="arc3,rad=-0.22"
    ),
    **ann_kw
)

# ── (ii) Local similarity: arrow points to anomaly interior ────────────────
interior_target = (cx + 0.030, cy + 0.020)

ann2 = ax.annotate(
    "(ii) Local similarity",
    xy=interior_target,
    xytext=(0.68, 0.82),
    ha="center",
    arrowprops=dict(
        **arrow_kw,
        connectionstyle="arc3,rad=-0.18"
    ),
    **ann_kw
)

for ann in [ann1, ann2]:
    ann.set_path_effects([
        path_effects.Stroke(linewidth=2.6, foreground="white"),
        path_effects.Normal()
    ])

# Re-raise legends so both remain visually on top.
for ax_i in axes:
    leg = ax_i.get_legend()
    if leg is not None:
        leg.set_zorder(1000)

# ── Shared adaptive colorbar ──────────────────────────────────────────────
cbar = fig.colorbar(
    im1,
    ax=axes,
    fraction=0.035,
    pad=0.025,
    extend="both"
)

# ── Coarse colorbar ticks with 0.5 spacing ────────────────────────────────
tick_step = 0.5

tick_min = np.ceil(vmin_adapt / tick_step) * tick_step
tick_max = np.floor(vmax_adapt / tick_step) * tick_step

ticks = np.arange(tick_min, tick_max + 0.5 * tick_step, tick_step)

# If no 0.5-spaced tick falls inside the adaptive range, show endpoints.
if len(ticks) == 0:
    ticks = np.array([vmin_adapt, vmax_adapt])

cbar.set_ticks(ticks)
cbar.set_ticklabels([f"{t:.1f}" for t in ticks])

cbar.set_label("Observation intensity", fontsize=13.5, color=TEXT_C)
cbar.ax.tick_params(labelsize=11.5, colors=TEXT_C)
cbar.outline.set_edgecolor(TEXT_C)
cbar.outline.set_linewidth(1.0)

plt.savefig(
    Path(__file__).resolve().parents[1] / "FinalVision" / "ABS.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True
)
plt.savefig(
    Path(__file__).resolve().parents[1] / "FinalVision" / "ABS.eps",
    format="eps",
    dpi=300,
    bbox_inches="tight",
    transparent=True
)

plt.close(fig)

print(f"Noise level used: {NOISE_LEVEL:.3f}")
print(f"Anomaly main intensity used: {ANOMALY_MAIN:.3f}")
print(f"Adaptive shared colorbar range: [{norm.vmin:.3f}, {norm.vmax:.3f}]")
print(f"Colorbar ticks: {ticks}")
