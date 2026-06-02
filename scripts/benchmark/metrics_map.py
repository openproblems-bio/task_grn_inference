"""
Generate a figure showing the mapping from metric code names to display names,
grouped by category. Saves to local output and docs images folder.
"""

import sys
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path

matplotlib.rcParams["font.family"] = "Arial"

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.helper import load_env, surrogate_names

env = load_env()
TASK_GRN_INFERENCE_DIR = env["TASK_GRN_INFERENCE_DIR"]
sys.path.insert(0, TASK_GRN_INFERENCE_DIR)
from task_grn_inference.src.utils.config import METRICS

CATEGORY_TO_METRICS = {
    "Regression":            ["r_precision", "r_recall"],
    "WS distance":           ["ws_precision", "ws_recall"],
    "TF recovery":           ["t_rec_precision", "t_rec_recall"],
    "TF binding":            ["tfb_f1"],
    "SEM":                   ["sem"],
    "Gene sets recovery":    ["gs_f1"],
    "Virtual cell":          ["vc"],
    "Replicate consistency": ["replicate_consistency"],
}

rows = []
for category, metrics in CATEGORY_TO_METRICS.items():
    for i, m in enumerate(metrics):
        if m not in METRICS:
            continue
        rows.append((category if i == 0 else "", m, surrogate_names.get(m, m)))

# ── Layout ─────────────────────────────────────────────────────────────────────
FIG_W    = 4.0
FIG_H    = 4.0
N        = len(rows)           # 11
HEADER_H = FIG_H * 0.10        # 10% for header
ROW_H    = (FIG_H - HEADER_H) / N

COL_FRACS  = [0.30, 0.36, 0.34]
COL_LABELS = ["Category", "Raw name", "Display name"]

HDR_BG  = "#1a3a4a";  HDR_FG  = "#ffffff"
CAT_BG  = "#d6e8f0";  CAT_FG  = "#1a3a4a"
ODD_BG  = "#f0f5f8";  EVEN_BG = "#ffffff"
CODE_FG = "#1a6080";  DISP_FG = "#222222"
GRID    = "#b0c8d8"

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, FIG_W)
ax.set_ylim(0, FIG_H)
ax.axis("off")

# Column left edges in data (inch) coords
col_x = [sum(COL_FRACS[:i]) * FIG_W for i in range(len(COL_FRACS))]

def draw_rect(x, y, w, h, color):
    ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=color, edgecolor="none", zorder=2))

# ── Header (sits at the TOP: y from FIG_H-HEADER_H to FIG_H) ──────────────────
hy = FIG_H - HEADER_H
for ci, (label, frac) in enumerate(zip(COL_LABELS, COL_FRACS)):
    cw = frac * FIG_W
    draw_rect(col_x[ci], hy, cw, HEADER_H, HDR_BG)
    ax.text(col_x[ci] + cw / 2, hy + HEADER_H / 2, label,
            ha="center", va="center", fontsize=7.5,
            color=HDR_FG, fontweight="bold", zorder=3)

# ── Data rows (top row just below header, going downward) ─────────────────────
for ri, (cat, code, display) in enumerate(rows):
    # row top = hy - ri*ROW_H, row bottom = hy - (ri+1)*ROW_H
    y_bottom = hy - (ri + 1) * ROW_H
    bg = CAT_BG if cat else (ODD_BG if ri % 2 == 1 else EVEN_BG)

    for ci, (text, frac) in enumerate(zip([cat, code, display], COL_FRACS)):
        cw = frac * FIG_W
        draw_rect(col_x[ci], y_bottom, cw, ROW_H, bg)
        if text:
            if ci == 0:
                fg, fw = CAT_FG, "bold"
            elif ci == 1:
                fg, fw = CODE_FG, "normal"
            else:
                fg, fw = DISP_FG, "normal"
            ax.text(col_x[ci] + 0.08, y_bottom + ROW_H / 2, text,
                    ha="left", va="center", fontsize=7,
                    color=fg, fontweight=fw, zorder=3)

    # Horizontal separator
    lw = 0.7 if cat and ri > 0 else 0.3
    ax.plot([0, FIG_W], [y_bottom + ROW_H, y_bottom + ROW_H],
            color=GRID, linewidth=lw, zorder=4)
    if ri == N - 1:
        ax.plot([0, FIG_W], [y_bottom, y_bottom],
                color=GRID, linewidth=0.3, zorder=4)

# ── Vertical dividers ──────────────────────────────────────────────────────────
for x in col_x[1:]:
    ax.plot([x, x], [0, FIG_H], color=GRID, linewidth=0.5, zorder=4)

# ── Outer border ───────────────────────────────────────────────────────────────
ax.add_patch(plt.Rectangle((0, 0), FIG_W, FIG_H,
                             facecolor="none", edgecolor="#5a8aa0",
                             linewidth=1.0, zorder=5))

plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

# ── Save ───────────────────────────────────────────────────────────────────────
out_local = Path(env["genernbi_supp_DIR"]) / "output" / "figs" / "metrics_map.png"
out_docs  = Path(TASK_GRN_INFERENCE_DIR) / "docs" / "source" / "images" / "metrics_map.png"

out_local.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out_local, dpi=200, bbox_inches="tight", facecolor="white")
print(f"Saved: {out_local}")
fig.savefig(out_docs, dpi=200, bbox_inches="tight", facecolor="white")
print(f"Saved: {out_docs}")
plt.close()


