import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Circle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy.signal import convolve2d

# ============================================================
# 全局字体设置：只修改这一处即可改变所有文字大小
# ============================================================
GLOBAL_FONT_SIZE = 32
TITLE_FONT_SIZE  = GLOBAL_FONT_SIZE
SYMBOL_FONT_SIZE = GLOBAL_FONT_SIZE
TICK_FONT_SIZE   = GLOBAL_FONT_SIZE
# 高斯核单元格内三位小数的最大无重叠字号。
KERNEL_TEXT_SIZE = 12.5
FIGURE_SIZE = (30, 5.2)

# ============================================================
# 参数设置
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
# 数据生成
# ============================================================

# 背景高斯噪声
bg_large = np.random.normal(0, sigma_b, (N_large, N_large))

# 单点异常场
anom_large = np.zeros((N_large, N_large))
margin = N_win // 2 + 1
ry = np.random.randint(margin, N_large - margin)
rx = np.random.randint(margin, N_large - margin)
anom_large[ry, rx] = anom_strength

# 高斯卷积核
cy_k, cx_k = (N_win - 1) / 2, (N_win - 1) / 2
ys, xs = np.mgrid[0:N_win, 0:N_win]

kernel_raw = np.exp(
    -((xs - cx_k) ** 2 + (ys - cy_k) ** 2)
    / (2 * sigma_kernel ** 2)
)
kernel = kernel_raw / kernel_raw.sum()

# 卷积信号
combined_large = bg_large + anom_large
img_conv = convolve2d(combined_large, kernel, mode='valid')

# 观测噪声
noise_obs = np.random.normal(0, sigma_obs, (N_small, N_small))

# 最终观测
img_final = img_conv + noise_obs

# ============================================================
# 统一色标范围
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
# 异常位置归一化
# ============================================================
anom_cx_norm = (rx + 0.5) / N_large
anom_cy_norm = (ry + 0.5) / N_large

conv_rx = rx - (N_win // 2)
conv_ry = ry - (N_win // 2)

conv_cx_norm = np.clip((conv_rx + 0.5) / N_small, 0.05, 0.95)
conv_cy_norm = np.clip((conv_ry + 0.5) / N_small, 0.05, 0.95)

# ============================================================
# 绘图风格设置
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
# 辅助函数
# ============================================================

def format_ax(ax, title=None):
    """统一设置子图标题、坐标轴刻度和边框样式。"""
    ax.set_aspect('equal', adjustable='box')

    if title:
        title_artist = ax.set_title(title, fontsize=TITLE_FONT_SIZE, pad=8)
        title_artist.set_linespacing(0.9)

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


def sym_ax(ax, sym):
    """绘制运算符号，例如 +, =, ⊛, 括号。"""
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
    """绘制高斯卷积核矩阵。"""
    n = kernel.shape[0]

    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_aspect('equal')
    ax.axis('off')

    if title:
        title_artist = ax.set_title(title, fontsize=TITLE_FONT_SIZE, pad=8)
        title_artist.set_linespacing(0.9)

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
                fontsize=KERNEL_TEXT_SIZE,
                color=txt_color,
                fontfamily='serif'
            )

# ============================================================
# 一行布局
# ============================================================
col_ratios = [
    0.08, 1.0, 0.15, 1.0, 0.08, 0.15,
    1.0, 0.15, 1.0, 0.15, 1.0, 0.15, 1.0
]

fig = plt.figure(figsize=FIGURE_SIZE)

gs = gridspec.GridSpec(
    1, 13,
    width_ratios=col_ratios,
    left=0.02,
    right=0.915,
    top=0.82,
    bottom=0.10,
    wspace=0.04
)

# 左括号
ax_lp = fig.add_subplot(gs[0, 0])
sym_ax(ax_lp, '(')

# 背景噪声
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

# 加号
ax_p1 = fig.add_subplot(gs[0, 2])
sym_ax(ax_p1, '+')

# 异常场
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
ax2.add_patch(
    Circle(
        (anom_cx_norm, anom_cy_norm),
        radius=0.07,
        edgecolor='#EE0000',
        facecolor='none',
        linewidth=1.8
    )
)

# 右括号
ax_rp = fig.add_subplot(gs[0, 4])
sym_ax(ax_rp, ')')

# 卷积符号
ax_cs = fig.add_subplot(gs[0, 5])
sym_ax(ax_cs, '$\\circledast$')

# 高斯核
ax3 = fig.add_subplot(gs[0, 6])
draw_kernel_panel(ax3, kernel, title='Gaussian kernel')

# 等号
ax_eq1 = fig.add_subplot(gs[0, 7])
sym_ax(ax_eq1, '=')

# 卷积信号
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
ax4.add_patch(
    Circle(
        (conv_cx_norm, conv_cy_norm),
        radius=0.12,
        edgecolor='#EE0000',
        facecolor='none',
        linewidth=1.8
    )
)

# 加号
ax_p2 = fig.add_subplot(gs[0, 9])
sym_ax(ax_p2, '+')

# 观测噪声
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

# 等号
ax_eq2 = fig.add_subplot(gs[0, 11])
sym_ax(ax_eq2, '=')

# 最终观测
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
# 共用色条
# ============================================================
cax = fig.add_axes([0.926, 0.10, 0.010, 0.72])

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
# 保存图片
# ============================================================
output = Path(__file__).resolve().parents[1] / 'FinalVision' / 'ConAnomaly'
plt.savefig(
    output.with_suffix('.png'),
    format='png',
    dpi=600
)

plt.savefig(
    output.with_suffix('.eps'),
    format='eps',
    dpi=300
)

print("Saved!")
plt.close(fig)
