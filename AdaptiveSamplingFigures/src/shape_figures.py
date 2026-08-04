"""Compare shape-dependent anomaly fields and detection performance.

Two publication figures pair geometric definitions with representative noisy
anomaly fields and out-of-control average run-length summaries. Ellipse panels
vary ellipticity ``kappa = a / b``; crescent panels vary normalized center
separation ``zeta = R_d / R``. Shared color limits preserve comparability
between the heatmaps, and fixed random seeds reproduce the illustrative noise.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib import ticker
import pandas as pd
import numpy as np
import os
from pathlib import Path

# ============================================================
# Publication-wide typography and stroke defaults
# ============================================================
plt.rcParams.update({
    'font.family':        'serif',
    'font.serif':         ['Liberation Serif', 'DejaVu Serif'],
    'mathtext.fontset':   'stix',
    'axes.unicode_minus': False,
    'font.size':          17,       # Default text size in points.
    'axes.linewidth':     0.8,      # Thin journal-style axis frames.
    'lines.linewidth':    1.0,      # Default curve width in points.
    'figure.dpi':         300,      # Interactive DPI; savefig controls export DPI.
})

# ============================================================
# Parametric anomaly-shape generator
# ============================================================
class AnomalyVisualizer:
    """Rasterize binary ellipse and crescent anomalies on a unit-square grid."""

    def __init__(self, grid_size=100, domain_size=1.0):
        self.N = grid_size  # Spatial resolution along each square-grid axis.
        x = np.linspace(0, domain_size, self.N)
        self.X, self.Y = np.meshgrid(x, x)

    def _d(self, v1, v2):
        """Return elementwise absolute displacement on the Cartesian grid."""
        return np.abs(v1 - v2)

    def generate_ellipse(self, area, ellipticity, center, value=1.0):
        """Generate a filled ellipse with fixed area and axis ratio.

        ``ellipticity`` is the semi-axis ratio ``a / b``. Solving the area
        equation ``pi * a * b = area`` keeps anomaly mass comparable while the
        shape parameter changes.
        """
        cx, cy = center
        sb = np.sqrt(area / (np.pi * ellipticity))
        sa = sb * ellipticity
        mask = (self._d(self.X, cx)**2 / sa**2 +
                self._d(self.Y, cy)**2 / sb**2) <= 1.0
        g = np.zeros_like(self.X); g[mask] = value
        return g

    def generate_crescent(self, area, shift_ratio, center, value=1.0):
        """Generate a crescent as the difference of two equal-radius disks.

        ``shift_ratio`` is the center separation divided by disk radius. The
        outer radius is derived so the retained crescent area remains fixed.
        Values approaching two reduce to a full disk in this parameterization.
        """
        cx, cy = center
        ratio = np.clip(shift_ratio, 0.05, 2.0)
        g = np.zeros_like(self.X)
        if ratio > 1.99:
            Ro = np.sqrt(area / np.pi)
            g[(self._d(self.X, cx)**2 + self._d(self.Y, cy)**2) <= Ro**2] = value
        else:
            t  = min(ratio * 0.5, 1.0)
            ai = 2.0 * np.arccos(t) - t * np.sqrt(4.0 - ratio**2)
            Ro = np.sqrt(area / (np.pi - ai))
            cxs = cx + ratio * Ro
            dx  = self._d(self.X, cx); dy = self._d(self.Y, cy)
            dxi = self._d(self.X, cxs)
            mask = (dx**2 + dy**2 <= Ro**2) & (dxi**2 + dy**2 > Ro**2)
            g[mask] = value
        return g


# ============================================================
# Statistics input filtering
# ============================================================
def load_filter_data(file_path, target_amp, target_area, shape_col='Ellipse'):
    """Load one statistics table and select a target anomaly scenario.

    Column names are normalized only when a case-insensitive match exists.
    Tukey whiskers are derived from the interquartile range after filtering by
    amplitude and area. Unexpected or missing columns return an empty frame so
    the calling panel can display an explicit no-data state.
    """

    # Fixed statistical columns plus the selected dynamic shape parameter.
    USE_COLS = ['Amp', 'Area', shape_col,
                'OC_Mean', 'OC_Std', 'OC_Q1', 'OC_Med', 'OC_Q3', 'OC_Max', 'OC_Min']

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Missing required shape statistics: {file_path}")
    else:
        try:
            # Read the header first so legacy case differences can be normalized.
            df_raw = pd.read_csv(file_path, sep=r'\s+')
            # A legacy SASAM crescent file uses an uppercase shape-column name.
            if shape_col not in df_raw.columns:
                matched = next(
                    (column for column in df_raw.columns if column.casefold() == shape_col.casefold()),
                    None,
                )
                if matched is not None:
                    df_raw.rename(columns={matched: shape_col}, inplace=True)
        except Exception as exc:
            print(f"[error] Could not read {file_path}: {exc}")
            return pd.DataFrame()

        missing = [c for c in USE_COLS if c not in df_raw.columns]
        if missing:
            print(f"[error] Missing required columns in {file_path}: {missing}")
            return pd.DataFrame()

        # Retain only fields consumed by the final figure pipeline.
        df = df_raw[USE_COLS].copy()

    if 'Amp' not in df.columns or 'Area' not in df.columns:
        return pd.DataFrame()

    mask = (np.isclose(df['Amp'],  target_amp,  atol=1e-5) &
            np.isclose(df['Area'], target_area, atol=1e-5))
    df_f = df[mask].copy()
    if not df_f.empty:
        df_f.sort_values(shape_col, inplace=True)
        iqr = df_f['OC_Q3'] - df_f['OC_Q1']
        df_f['CalcMin'] = (df_f['OC_Q1'] - 1.5 * iqr).clip(lower=0.0)
        df_f['CalcMax'] =  df_f['OC_Q3'] + 1.5 * iqr
    return df_f


# ============================================================
# Reproducible heatmap fields
# ============================================================
np.random.seed(42)
viz = AnomalyVisualizer(grid_size=100)  # 100 x 100 publication raster grid.
target_area_img, anomaly_val, noise_std = 0.03, 2.0, 1.0
img_center = (0.5, 0.5)  # Center every anomaly to isolate shape effects.

ellipse_params  = [0.1, 0.2, 0.5, 0.8, 1.0]
crescent_params = [0.4, 0.8, 1.2, 1.6, 2.0]
KEPT_COLS = [0, 2, 4]  # Show low, middle, and high shape-parameter examples.

ellipse_noisy, crescent_noisy = [], []
for p in ellipse_params:
    c = viz.generate_ellipse(target_area_img, p, img_center, anomaly_val)
    ellipse_noisy.append(c + np.random.normal(0, noise_std, c.shape))
for p in crescent_params:
    c = viz.generate_crescent(target_area_img, p, img_center, anomaly_val)
    crescent_noisy.append(c + np.random.normal(0, noise_std, c.shape))

ellipse_sel  = [ellipse_noisy[i]  for i in KEPT_COLS]
crescent_sel = [crescent_noisy[i] for i in KEPT_COLS]
ellipse_params_sel  = [ellipse_params[i]  for i in KEPT_COLS]
crescent_params_sel = [crescent_params[i] for i in KEPT_COLS]

cmap_name = 'RdBu_r'
n_img = len(KEPT_COLS)  # Three representative heatmaps per row.

# ============================================================
# Detection-statistics inputs and method colors
# ============================================================
box_colors = ['#1f77b4', '#d62728']  # Blue: mM-KSTD; red: SASAM.
SCENARIOS  = [(2.0, 0.03), (0.5, 0.03)]

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SHAPE_DATA = PROJECT_ROOT / 'ExternalData' / 'shape'
ef_file1 = SHAPE_DATA / 'mM_KSTD_ellipse.txt'
ef_file2 = SHAPE_DATA / 'SASAM_ellipse.txt'

gh_file1 = SHAPE_DATA / 'mM_KSTD_crescent.txt'
gh_file2 = SHAPE_DATA / 'SASAM_crescent.txt'

# ============================================================
# Geometric-definition panels
# ============================================================
def draw_ellipticity_geo(ax):
    """Draw the semi-axes that define ellipse ellipticity ``kappa = a / b``."""
    ax.set_aspect('equal')
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    ax.axis('off')

    a, b = 3.0, 1.8

    ax.add_patch(patches.Ellipse(
        (0, 0), 2*a, 2*b,
        edgecolor='black', facecolor='#EE7942',
        linewidth=1.0, alpha=0.9))

    ax.plot([-a, a], [0, 0], color='gray', linestyle='--', linewidth=0.8)
    ax.plot([0, 0], [-b, b], color='gray', linestyle='--', linewidth=0.8)
    ax.plot(0, 0, 'ko', markersize=3)

    ax.text(0.15, -0.38, '$O$', fontsize=17)

    ax.annotate('', xy=(0, 0.12), xytext=(a, 0.12),
        arrowprops=dict(arrowstyle='<->', color='black', lw=0.8))
    ax.text(a/2, 0.35, '$a$', ha='center', fontsize=18)

    ax.annotate('', xy=(-0.12, 0), xytext=(-0.12, b),
        arrowprops=dict(arrowstyle='<->', color='black', lw=0.8))
    ax.text(-0.35, b/2, '$b$', va='center', ha='right', fontsize=18)

    ax.text(0, -3.3, r'$\kappa = \dfrac{a}{b}$', ha='center', fontsize=18)


def draw_eccentricity_geo(ax):
    """Draw the equal-radius construction defining crescent eccentricity.

    The visible crescent is the set difference of two disks separated by
    ``R_d``. Dashed arcs expose the occluded disk boundaries needed to interpret
    the geometric definition.
    """
    ax.set_aspect('equal')
    ax.set_xlim(-3.5, 3.5)
    ax.set_ylim(-3.5, 3.5)
    ax.axis('off')

    R, Rd = 2.0, 1.0

    xi = Rd / 2
    theta_rad = np.arctan2(np.sqrt(R**2 - xi**2), xi)
    theta_deg = np.degrees(theta_rad)
    ao = np.linspace(theta_rad, 2*np.pi - theta_rad, 150)
    ai = np.linspace(np.pi + theta_rad, np.pi - theta_rad, 150)

    ax.fill(np.concatenate([R*np.cos(ao),  Rd + R*np.cos(ai)]),
            np.concatenate([R*np.sin(ao),       R*np.sin(ai)]),
            color='#EE7942', alpha=0.9)

    ax.add_patch(patches.Arc((0,  0), 2*R, 2*R, angle=0,
        theta1=theta_deg,     theta2=360-theta_deg,
        color='black', linewidth=1.0))
    ax.add_patch(patches.Arc((0,  0), 2*R, 2*R, angle=0,
        theta1=-theta_deg,    theta2=theta_deg,
        color='gray',  linewidth=0.8, linestyle='--'))
    ax.add_patch(patches.Arc((Rd, 0), 2*R, 2*R, angle=0,
        theta1=180-theta_deg, theta2=180+theta_deg,
        color='black', linewidth=1.0))
    ax.add_patch(patches.Arc((Rd, 0), 2*R, 2*R, angle=0,
        theta1=180+theta_deg, theta2=180-theta_deg+360,
        color='gray',  linewidth=0.8, linestyle='--'))

    ax.plot(0,  0, 'ko', markersize=3)
    ax.text(-0.28, -0.38, '$O_1$', fontsize=17)
    ax.plot(Rd, 0, 'ko', markersize=3)
    ax.text(Rd+0.1, -0.38, '$O_2$', fontsize=17)

    ax.annotate('', xy=(0, 0), xytext=(Rd, 0),
        arrowprops=dict(arrowstyle='<->', color='black', lw=1.0))
    ax.text(Rd/2, 0.18, r'$R_d$', ha='center', va='bottom', fontsize=18)

    th1 = 110 * np.pi / 180
    rx1, ry1 = R*np.cos(th1), R*np.sin(th1)
    ax.plot([0, rx1], [0, ry1], 'k-', lw=0.8)
    ax.annotate('', xy=(0, 0), xytext=(rx1, ry1),
        arrowprops=dict(arrowstyle='<->', color='black', lw=0.8))
    ax.text(rx1/2 - 0.28, ry1/2, '$R$', fontsize=18, ha='right')

    th2 = -60 * np.pi / 180
    rx2, ry2 = Rd + R*np.cos(th2), R*np.sin(th2)
    ax.plot([Rd, rx2], [0, ry2], 'k-', lw=0.8)
    ax.annotate('', xy=(Rd, 0), xytext=(rx2, ry2),
        arrowprops=dict(arrowstyle='<->', color='black', lw=0.8))
    ax.text(Rd + (rx2-Rd)/2 + 0.28, ry2/2, '$R$', fontsize=18, ha='left')

    ax.text(0, -3.0, r'$\zeta = \dfrac{R_d}{R}$', ha='center', fontsize=18)


# ============================================================
# Heatmap-row rendering
# ============================================================
def draw_img_row(fig, gs_specs, data_list, params, param_sym, vmin, vmax,
                 float_fmt='.2f'):
    """Render representative anomaly fields with shared limits and math labels."""
    axes_img = []
    last_im  = None
    for i in range(n_img):
        ax = fig.add_subplot(gs_specs[i])
        im = ax.imshow(data_list[i], origin='lower', extent=[0, 1, 0, 1],
                       cmap=cmap_name, vmin=vmin, vmax=vmax)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.6)
        param_str = format(params[i], float_fmt)
        ax.set_xlabel(rf'${param_sym} = {param_str}$', fontsize=15)
        axes_img.append(ax)
        last_im = im
    return axes_img, last_im


def add_img_colorbar(fig, axes_img, last_im):
    """Place a colorbar aligned exactly with the final heatmap row.

    Positions are measured after the canvas layout is resolved; this prevents
    the colorbar height from drifting when figure margins or aspect ratios vary.
    """
    fig.canvas.draw()

    bb_last  = axes_img[-1].get_position()
    bb_first = axes_img[0].get_position()

    cb_left   = bb_last.x1 + 0.010
    cb_bottom = bb_last.y0
    cb_height = bb_last.y1 - bb_last.y0
    cb_width  = 0.013

    cax = fig.add_axes([cb_left, cb_bottom, cb_width, cb_height])
    cb  = fig.colorbar(last_im, cax=cax)
    cb.ax.tick_params(labelsize=12)
    cb.locator = ticker.MaxNLocator(nbins=5, symmetric=True)
    cb.update_ticks()

    x_center   = (bb_first.x0 + bb_last.x1) / 2
    y_subtitle = bb_last.y0 - 0.030

    return x_center, y_subtitle, bb_last.y0


# ============================================================
# Detection-performance boxplots
# ============================================================
bar_width = 0.35
offsets   = [-bar_width / 2, bar_width / 2]

def draw_one_boxplot(ax, file1, file2, t_amp, t_area, label, xlabel,
                     float_fmt='.2f', shape_col='Ellipse'):
    """Draw one method-comparison boxplot panel.

    Canonical public columns are ``Amp``, ``Area``, the selected ``shape_col``,
    and the ``OC_*`` summary statistics. ``OC_Mean`` supplies the overlaid mean
    trajectory; ``OC_Q1``, ``OC_Med``, and ``OC_Q3`` define each box; and
    ``CalcMin``/``CalcMax`` are Tukey whiskers derived from the quartiles.
    Unused legacy columns are intentionally ignored.
    """
    df1 = load_filter_data(file1, t_amp, t_area, shape_col=shape_col)
    df2 = load_filter_data(file2, t_amp, t_area, shape_col=shape_col)

    vals1 = df1[shape_col].values if not df1.empty else np.array([])
    vals2 = df2[shape_col].values if not df2.empty else np.array([])
    all_x = np.unique(np.concatenate([vals1, vals2]))

    if len(all_x) == 0:
        ax.text(0.5, 0.5, "No Data", ha='center', va='center', fontsize=14)
        ax.set_title(rf'{label} $\delta={t_amp}$',
                     fontsize=15, y=-0.18)
        return

    val_to_idx = {v: i for i, v in enumerate(all_x)}

    for i, df in enumerate([df1, df2]):
        if df.empty: continue
        color  = box_colors[i]
        offset = offsets[i]
        stats_list, positions, line_x, line_y = [], [], [], []

        for _, row in df.iterrows():
            val = row[shape_col]
            if val in val_to_idx:
                pos = val_to_idx[val] + offset
                positions.append(pos); line_x.append(pos); line_y.append(row['OC_Mean'])
                stats_list.append({
                    'med':     row['OC_Med'],
                    'q1':      row['OC_Q1'],
                    'q3':      row['OC_Q3'],
                    'whislo':  row['CalcMin'],
                    'whishi':  row['CalcMax'],
                })

        if stats_list:
            box = ax.bxp(stats_list, positions=positions,
                         widths=bar_width * 0.8,
                         patch_artist=True, showfliers=False, manage_ticks=False)
            for patch in box['boxes']:
                patch.set(facecolor=color, alpha=0.6,
                          edgecolor='black', linewidth=1.0)
            for el in ['whiskers', 'caps', 'medians']:
                c = 'white' if el == 'medians' else 'black'
                plt.setp(box[el], color=c,
                         linewidth=1.5 if el == 'medians' else 1.0)

        ax.plot(line_x, line_y, color=color,
                linewidth=2.0, marker='o', markersize=5,
                markeredgecolor='white', label='_nolegend_')

    ax.set_xlim(-0.5, len(all_x) - 0.5)

    step = 2
    sparse_idx = np.arange(0, len(all_x), step)
    ax.set_xticks(sparse_idx)
    ax.set_xticklabels([format(all_x[i], float_fmt) for i in sparse_idx],
                       rotation=0, fontsize=13)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)

    ax.tick_params(axis='both', which='major',
                   direction='in', length=5, width=1,
                   bottom=True, top=False, left=True, right=False,
                   labelbottom=True, labelleft=True)
    ax.tick_params(axis='y', labelsize=13)
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel('OC ARL values', fontsize=15)
    ax.set_title(rf'{label} $\delta={t_amp}$, $S_\delta={t_area}$',
                 fontsize=15, y=-0.22)
    ax.grid(True, axis='y', linestyle='--', alpha=0.5)
    ax.legend(handles=[Patch(facecolor=box_colors[0], edgecolor='black', label='mM-KSTD'),
                        Patch(facecolor=box_colors[1], edgecolor='black', label='SASAM')],
              loc='upper right', fontsize=12,
              frameon=True, edgecolor='black', framealpha=0.9)


# ============================================================
# Shared layout constants for both multi-panel figures
# ============================================================
geo_w = 1.45
img_w = 1.0
FIG_W = 20
FIG_H = 16

# ============================================================
# Align row captions only after Matplotlib has resolved the axes geometry.
# ============================================================
def place_row_titles(fig, ax_geo, axes_img, title_a, title_b, img_bottom_y):
    """Center panel captions beneath the geometry and heatmap groups."""
    bb_geo   = ax_geo.get_position()
    bb_first = axes_img[0].get_position()
    bb_last  = axes_img[-1].get_position()

    y_title  = img_bottom_y - 0.048
    x_a = (bb_geo.x0 + bb_geo.x1) / 2
    x_b = (bb_first.x0 + bb_last.x1) / 2

    fig.text(x_a, y_title, title_a, ha='center', va='top', fontsize=17)
    fig.text(x_b, y_title, title_b, ha='center', va='top', fontsize=16, fontweight='normal')


# ============================================================
# Main figure 1: ellipse geometry, fields, and performance.
# Expected boxplot columns: Amp, Area, Ellipse,
#   OC_Mean, OC_Std, OC_Q1, OC_Med, OC_Q3, OC_Max, OC_Min
# ============================================================
vmin_shared = min(np.min(ellipse_sel), np.min(crescent_sel))
vmax_shared = max(np.max(ellipse_sel), np.max(crescent_sel))

fig1 = plt.figure(figsize=(FIG_W, FIG_H))

gs1_top = gridspec.GridSpec(
    1, 4, figure=fig1,
    left=0.02, right=0.95,
    top=0.97,  bottom=0.15,
    wspace=0.08,
    width_ratios=[geo_w, img_w, img_w, img_w])

gs1_bot = gridspec.GridSpec(
    1, 2, figure=fig1,
    left=0.05, right=0.97,
    top=0.35,  bottom=0.07,
    wspace=0.18)

ax_geo1 = fig1.add_subplot(gs1_top[0, 0])
draw_ellipticity_geo(ax_geo1)

img_specs1  = [gs1_top[0, 1], gs1_top[0, 2], gs1_top[0, 3]]
axes_img1, last_im1 = draw_img_row(fig1, img_specs1, ellipse_sel,
                                    ellipse_params_sel, r'\kappa', vmin_shared, vmax_shared,
                                    float_fmt='.2f')

x_sub1, y_sub1, img_bot1 = add_img_colorbar(fig1, axes_img1, last_im1)

place_row_titles(fig1, ax_geo1, axes_img1,
                 "(a) Ellipticity geometry", "(b) Ellipse Anomalies", img_bot1)

ax_ef0 = fig1.add_subplot(gs1_bot[0, 0])
ax_ef1 = fig1.add_subplot(gs1_bot[0, 1])
# ``shape_col='Ellipse'`` selects Amp / Area / Ellipse / OC_* fields.
draw_one_boxplot(ax_ef0, ef_file1, ef_file2,
                 SCENARIOS[0][0], SCENARIOS[0][1], '(c)', r'Ellipticity $\kappa$',
                 float_fmt='.2f', shape_col='Ellipse')
draw_one_boxplot(ax_ef1, ef_file1, ef_file2,
                 SCENARIOS[1][0], SCENARIOS[1][1], '(d)', r'Ellipticity $\kappa$',
                 float_fmt='.2f', shape_col='Ellipse')

fig1.savefig(PROJECT_ROOT / 'FinalVision' / 'shape_figure_1.png', dpi=600, bbox_inches='tight')
fig1.savefig(PROJECT_ROOT / 'FinalVision' / 'shape_figure_1.eps', format='eps', bbox_inches='tight')
print("Main figure 1 saved.")


# ============================================================
# Main figure 2: crescent geometry, fields, and performance.
# Expected boxplot columns: Amp, Area, Crescent,
#   OC_Mean, OC_Std, OC_Q1, OC_Med, OC_Q3, OC_Max, OC_Min
# ============================================================
fig2 = plt.figure(figsize=(FIG_W, FIG_H))

gs2_top = gridspec.GridSpec(
    1, 4, figure=fig2,
    left=0.05, right=0.95,
    top=0.97,  bottom=0.15,
    wspace=0.08,
    width_ratios=[geo_w, img_w, img_w, img_w])

gs2_bot = gridspec.GridSpec(
    1, 2, figure=fig2,
    left=0.08, right=0.97,
    top=0.35,  bottom=0.07,
    wspace=0.18)

ax_geo2 = fig2.add_subplot(gs2_top[0, 0])
draw_eccentricity_geo(ax_geo2)

img_specs2  = [gs2_top[0, 1], gs2_top[0, 2], gs2_top[0, 3]]
axes_img2, last_im2 = draw_img_row(fig2, img_specs2, crescent_sel,
                                    crescent_params_sel, r'\zeta', vmin_shared, vmax_shared,
                                    float_fmt='.1f')

x_sub2, y_sub2, img_bot2 = add_img_colorbar(fig2, axes_img2, last_im2)

place_row_titles(fig2, ax_geo2, axes_img2,
                 "(a) Eccentricity geometry", "(b) Crescent Anomalies", img_bot2)

ax_gh0 = fig2.add_subplot(gs2_bot[0, 0])
ax_gh1 = fig2.add_subplot(gs2_bot[0, 1])
# ``shape_col='Crescent'`` selects Amp / Area / Crescent / OC_* fields.
draw_one_boxplot(ax_gh0, gh_file1, gh_file2,
                 SCENARIOS[0][0], SCENARIOS[0][1], '(c)', r'Eccentricity $\zeta$',
                 float_fmt='.1f', shape_col='Crescent')
draw_one_boxplot(ax_gh1, gh_file1, gh_file2,
                 SCENARIOS[1][0], SCENARIOS[1][1], '(d)', r'Eccentricity $\zeta$',
                 float_fmt='.1f', shape_col='Crescent')

fig2.savefig(PROJECT_ROOT / 'FinalVision' / 'shape_figure_2.png', dpi=600, bbox_inches='tight')
fig2.savefig(PROJECT_ROOT / 'FinalVision' / 'shape_figure_2.eps', format='eps', bbox_inches='tight')
print("Main figure 2 saved.")

plt.close(fig1)
plt.close(fig2)
