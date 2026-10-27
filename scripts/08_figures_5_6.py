# Figures 5-6
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Ellipse
from matplotlib import gridspec
import os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "07_figures_1_to_4.py"), encoding="utf-8").read().split("A = pd.read_csv")[0])   # reuse style + helpers
S = pd.read_csv(R + "S3_vKO_summary.csv")
fig = plt.figure(figsize=(175 * MM, 140 * MM))
gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[1, 1.1], height_ratios=[1, 0.62], wspace=0.9, hspace=0.75)
v = pd.read_csv(R + "S3_vKO_LTB_lymphocytes.csv"); v = v[v.gene != "LTB"].sort_values("z_KO")
top = pd.concat([v.head(12), v.tail(6)]).iloc[::-1]
axis_genes = {"MAP3K14", "CCL19", "CCR7", "GPR183", "CD44", "RELB", "NFKB2", "IL7R", "CXCL13", "CCL21"}
ax = fig.add_subplot(gs[0, 0]); letter(ax, "A", x=-0.42)
cols = [("#D55E00" if g in axis_genes else (C_CTRL if s < 0 else GREY)) for g, s in zip(top.gene, top.z_KO)]
ax.barh(range(len(top)), top.z_KO, color=cols, height=0.7)
ax.set_yticks(range(len(top))); ax.set_yticklabels(top.gene, style="italic", fontsize=6.5); ax.axvline(0, c="k", lw=0.6)
ax.set_xlabel("Predicted change after virtual $\\it{LTB}$\nknockout (z vs 1,000 random knockouts)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#D55E00", label="Lymphotoxin / lymphoid-homing gene"), Patch(color=C_CTRL, label="Other, decreased"), Patch(color=GREY, label="Increased")],
          frameon=False, loc="upper left", bbox_to_anchor=(0.40, 1.0), fontsize=6, handlelength=1.0, handletextpad=0.4)
ax = fig.add_subplot(gs[0, 1]); letter(ax, "B", x=-0.75)
lab = {"LTB_lymphocytes": "LTB | lymphocytes", "LTBR_stroma": "LTBR | fibroblast/mesangial", "RELB_immune": "RELB | immune",
       "CD44_immune": "CD44 | immune", "LTBR_PT_discovery": "LTBR | proximal tubule*", "MAP3K14_PT_discovery": "MAP3K14 | proximal tubule*"}
S["lab"] = S.label.map(lab); S = S.iloc[::-1]
ax.barh(range(len(S)), S.reversal_rho_KO, color=[C_DKD if p < 0.05 else ("#E69F00" if p < 0.1 else GREY) for p in S.empirical_P_KO_reverses], height=0.6)
for i, (r, p) in enumerate(zip(S.reversal_rho_KO, S.empirical_P_KO_reverses)):
    ax.text(min(r, -0.03) - 0.02, i, f"P = {p:.3f}", ha="right", va="center", fontsize=6)
ax.legend(handles=[Patch(color=C_DKD, label="P < 0.05"), Patch(color="#E69F00", label="P < 0.10"), Patch(color=GREY, label="P ≥ 0.10")],
          frameon=False, loc="lower left", fontsize=6, handlelength=1.0, handletextpad=0.4)
ax.set_yticks(range(len(S))); ax.set_yticklabels(S.lab, fontsize=6.5); ax.axvline(0, c="k", lw=0.6); ax.set_xlim(-0.9, 0.1); ax.set_ylim(-1.3, len(S) - 0.5)

ax.set_xlabel("Spearman ρ, knockout shift vs\nprogression log$_2$FC (negative = reversal)")
# C stability of the virtual LTB knockout
ST = pd.read_csv(R + "S6_vKO_LTB_stability.csv")
ax = fig.add_subplot(gs[1, :]); letter(ax, "C", x=-0.08)
lab_map = {"baseline (1,000 genes, 3 PCs, 600 cells, seed 1)": "Baseline", "seed 2": "Seed 2", "seed 3": "Seed 3", "500 genes": "500 genes",
           "1,500 genes": "1,500 genes", "2 PCs": "2 PCs", "5 PCs": "5 PCs", "1,000 cells per network": "1,000 cells", "expression-matched null": "Matched null"}
x = np.arange(len(ST))
ax.bar(x, ST.rho, color=[C_DKD if p_ < 0.05 else GREY for p_ in ST.P], width=0.6)
for i, (r_, p_) in enumerate(zip(ST.rho, ST.P)):
    ax.text(i, min(r_, 0) - 0.008, f"P = {p_:.3f}", ha="center", va="top", fontsize=5.8)
ax.axhline(0, c="k", lw=0.6); ax.set_ylim(-0.135, 0.02)
ax.set_xticks(x); ax.set_xticklabels([lab_map.get(v, v) for v in ST.setting], rotation=25, ha="right")
ax.set_ylabel("Reversal ρ, virtual\nLTB knockout")
ax.legend(handles=[Patch(color=C_DKD, label="Empirical P < 0.05"), Patch(color=GREY, label="P ≥ 0.05")], frameon=False, fontsize=6, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2)
save(fig, "Figure5")

# Figure 6 model
fig = plt.figure(figsize=(120 * MM, 100 * MM)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 120); ax.set_ylim(0, 100); ax.axis("off")
def b(x, y, w, h, t, fc, fs=7):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=2", fc=fc, ec="#404040", lw=0.6)); ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=fs, linespacing=1.3)
def ar(x1, y1, x2, y2, **k):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", lw=0.8, color="#404040", shrinkA=0, shrinkB=0, mutation_scale=8, **k))
# compartments
ax.add_patch(FancyBboxPatch((2, 30), 50, 66, boxstyle="round,pad=0,rounding_size=3", fc="none", ec="#D55E00", lw=0.8, ls="--"))
ax.add_patch(FancyBboxPatch((68, 30), 50, 66, boxstyle="round,pad=0,rounding_size=3", fc="none", ec="#0072B2", lw=0.8, ls="--"))
ax.text(27, 92, "Lymphocytes (T / B)", ha="center", fontsize=7.5, fontweight="bold", color="#D55E00")
ax.text(93, 92, "Fibroblast / mesangial cells", ha="center", fontsize=7.5, fontweight="bold", color="#0072B2")
b(8, 72, 38, 14, "↑ LTB", "#FBE3D6")
b(74, 72, 38, 14, "LTβR", "#DDEBF7")
b(74, 53, 38, 12, "NIK (MAP3K14)", "#DDEBF7")
b(74, 35, 38, 12, "RELB / NFKB2\n→ CCL19 / CCL21", "#DDEBF7")
b(8, 35, 38, 14, "CCR7$^{+}$ lymphocyte recruitment\n(possible lymphoid organisation)", "#FBE3D6")
b(14, 3, 92, 16, "Inflammation, matrix deposition and loss of tubular\nmetabolic programmes → DKD progression", "#EDEDED")
ar(46, 79, 74, 79); ax.text(60, 81, "LTα$_1$β$_2$", ha="center", fontsize=6.5)
ar(93, 72, 93, 65); ar(93, 53, 93, 47)
ar(74, 41, 46, 41); ax.text(60, 43, "chemokines", ha="center", fontsize=6.5, style="italic")
ar(27, 49, 27, 72); ax.text(25, 60.5, "proposed\nfeedback", ha="right", va="center", fontsize=6.5, style="italic")
ar(27, 35, 27, 19); ar(93, 35, 93, 19)
save(fig, "Figure6")
print(S[["lab","reversal_rho_KO","empirical_P_KO_reverses"]])
