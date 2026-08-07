import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
from matplotlib.patches import Circle
from scipy.interpolate import BSpline

# ================= 0. 统一出版级绘图参数 =================
GLOBAL_FONT_SIZE = 32
TITLE_FONT_SIZE = GLOBAL_FONT_SIZE
SYMBOL_FONT_SIZE = GLOBAL_FONT_SIZE
TICK_FONT_SIZE = GLOBAL_FONT_SIZE
FIGURE_SIZE = (30, 5.6)

# ================= 1. B样条基底生成 =================
def get_1d_bspline_matrix(n_control, degree, n_pixel=100):
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

# ================= 2. 数据生成 =================
np.random.seed(42)
N_PIXEL = 100
B1_mat = create_2d_basis_matrix(5, 5, degree=2, n_pixel=N_PIXEL)
B2_mat = create_2d_basis_matrix(20, 20, degree=3, n_pixel=N_PIXEL)

anom_idx = 250
coef_anom = np.zeros(400)
coef_anom[anom_idx] = 5.0
anomaly_img = (coef_anom @ B2_mat).reshape(N_PIXEL, N_PIXEL)

# Generate one background and one noise
coef_bg = np.random.normal(loc=0, scale=3.0, size=25)
bg_img = (coef_bg @ B1_mat).reshape(N_PIXEL, N_PIXEL)
noise_img = np.random.normal(loc=0, scale=1.0, size=(N_PIXEL, N_PIXEL))

# The five images:
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

# Global color scale across all images
all_vals = np.concatenate([img1.flatten(), img2.flatten(), img3.flatten(),
                            img4.flatten(), img5.flatten()])
abs_max = np.max(np.abs(all_vals))
vmin_global, vmax_global = -abs_max, abs_max

# Anomaly center for circle annotation
y_idx, x_idx = np.unravel_index(np.argmax(anomaly_img), anomaly_img.shape)
anom_center_x = (x_idx + 0.5) / N_PIXEL
anom_center_y = (y_idx + 0.5) / N_PIXEL

# ================= 3. 辅助函数 =================
def format_ax(ax, title=None):
    ax.set_aspect('equal', adjustable='box')

    if title:
        title_artist = ax.set_title(
            title,
            fontsize=TITLE_FONT_SIZE,
            pad=8,
        )
        title_artist.set_linespacing(0.9)

    ticks = [0, 1.0]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(['0', '1'], fontsize=TICK_FONT_SIZE)
    ax.set_yticklabels(['0', '1'], fontsize=TICK_FONT_SIZE)
    ax.tick_params(
        axis='both',
        direction='out',
        length=4,
        width=0.9,
        colors='black',
        pad=4,
    )
    for spine in ax.spines.values():
        spine.set_edgecolor('#B0B0B0')
        spine.set_linewidth(0.8)

# ================= 4. 绘图 =================
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif', 'Georgia', 'serif'],
    'mathtext.fontset': 'stix',
    'axes.unicode_minus': False,
    'font.size': GLOBAL_FONT_SIZE,
    'axes.titlesize': TITLE_FONT_SIZE,
    'axes.labelsize': GLOBAL_FONT_SIZE,
    'xtick.labelsize': TICK_FONT_SIZE,
    'ytick.labelsize': TICK_FONT_SIZE,
})

fig = plt.figure(figsize=FIGURE_SIZE)

# Layout: img(1) + (1) img(2) =(1) img(3) +(1) img(4) =(1) img(5)
# Column width ratios: [image, symbol, image, symbol, image, symbol, image, symbol, image]
ratios = [1, 0.15, 1, 0.15, 1, 0.15, 1, 0.15, 1]
gs = gridspec.GridSpec(
    1,
    9,
    width_ratios=ratios,
    wspace=0.04,
    left=0.02,
    right=0.915,
    top=0.79,
    bottom=0.11,
)

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

    # Add red circle on images that contain the anomaly (img3 and img5)
    if k in [0, 2, 4]:
        circle = Circle((anom_center_x, anom_center_y), radius=0.08,
                        edgecolor='#EE0000', facecolor='none',
                        linewidth=1.8, linestyle='-')
        ax.add_patch(circle)

# Symbols
for sym, col in zip(symbols, sym_cols):
    ax_sym = fig.add_subplot(gs[0, col])
    ax_sym.axis('off')
    ax_sym.text(0.5, 0.5, sym, fontsize=SYMBOL_FONT_SIZE, fontweight='bold',
                ha='center', va='center', transform=ax_sym.transAxes)

# Colorbar
cax = fig.add_axes([0.926, 0.11, 0.010, 0.68])
cb = plt.colorbar(im_last, cax=cax)
cb.set_label('Amplitude', fontsize=GLOBAL_FONT_SIZE)
cb.ax.tick_params(labelsize=GLOBAL_FONT_SIZE)
cb.outline.set_edgecolor('#B0B0B0')
cb.outline.set_linewidth(0.8)

output = Path(__file__).resolve().parents[1] / 'FinalVision' / 'Bspline_anomaly'
plt.savefig(output.with_suffix('.eps'),
            format='eps', dpi=300)
plt.savefig(output.with_suffix('.png'),
            format='png', dpi=600)
print("Saved!")
plt.close(fig)
