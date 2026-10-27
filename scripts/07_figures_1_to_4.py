# Figures 1-4 from the result tables
import os, sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib import gridspec
from PIL import Image

BASE = os.path.dirname(os.path.abspath(__file__)) + "/../"
R = BASE + "results/"
F = BASE + "figures/"
MR = R + "S4_MR_candidates_kidney_disease.csv"
OUT = BASE + "paper_figures/"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7, "axes.linewidth": 0.8,
                     "xtick.direction": "out", "ytick.direction": "out", "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
C_CTRL, C_DKD, GREY = "#0072B2", "#D55E00", "#9E9E9E"   # Okabe-Ito
MM = 1 / 25.4

def letter(ax, s, x=-0.18, y=1.06):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")

def save(fig, name):
    fig.savefig(OUT + name + ".pdf", bbox_inches="tight")
    fig.savefig(OUT + name + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT + name + "_tmp.png", dpi=600, bbox_inches="tight")
    im = Image.open(OUT + name + "_tmp.png").convert("RGB")
    im.save(OUT + name + ".tiff", compression="tiff_lzw", dpi=(600, 600)); os.remove(OUT + name + "_tmp.png")
    plt.close(fig)

def img(ax, path, crop=None):
    im = Image.open(path).convert("RGB")
    if crop: w, h = im.size; im = im.crop((int(crop[0] * w), int(crop[1] * h), int(crop[2] * w), int(crop[3] * h)))
    ax.imshow(im); ax.axis("off")

def box(ax, x, y, w, h, txt, fc="#F2F2F2", fs=6.3):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.015", fc=fc, ec="#404040", lw=0.6))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs, wrap=True)

def arrow(ax, x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", lw=0.7, color="#404040"))

def pstar(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."

A = pd.read_csv(R + "ALL_genes_integrated.csv", index_col=0)
cand = pd.read_csv(R + "progression_genes_tiered.csv", index_col=0)

# ============================================================ Figure 1 (bulk progression programme; study flowchart moved to Figure S1)
def fletter(fig, ax, s, dx=0.0):
    bb = ax.get_position(); fig.text(bb.x0 + dx, bb.y1 + 0.015, s, fontsize=10, fontweight="bold", va="bottom", ha="left")
fig = plt.figure(figsize=(175 * MM, 170 * MM))
gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1.35, 0.85], width_ratios=[1.1, 1], hspace=0.42, wspace=0.62,
                       left=0.09, right=0.92, top=0.95, bottom=0.07)
# A volcano
d = pd.read_csv(R + "DE_GSE142025_ADVvsEARLY.csv", index_col=0)
ax = fig.add_subplot(gs[0, 0])
sig = d.FDR < 0.05; up = sig & (d.logFC > 0); dn = sig & (d.logFC < 0)
ax.scatter(d.logFC[~sig], -np.log10(d.P[~sig]), s=2, c="#C8C8C8", lw=0, rasterized=True, label="FDR ≥ 0.05")
ax.scatter(d.logFC[up], -np.log10(d.P[up]), s=2.5, c=C_DKD, lw=0, rasterized=True, label=f"Up (n = {up.sum():,})")
ax.scatter(d.logFC[dn], -np.log10(d.P[dn]), s=2.5, c=C_CTRL, lw=0, rasterized=True, label=f"Down (n = {dn.sum():,})")
ax.axhline(-np.log10(d.P[sig].max()), c="#606060", lw=0.5, ls="--")
lab = {"IL7R": (3.95, 7.7), "CCL21": (3.95, 6.3), "CCL19": (4.05, 2.9), "COL1A1": (3.4, 1.0), "RELB": (0.15, 8.35), "CD44": (2.2, 8.7), "EGF": (-3.4, 6.3)}
for g, (tx, ty) in lab.items():
    x0, y0 = d.at[g, "logFC"], -np.log10(d.at[g, "P"])
    ax.scatter([x0], [y0], s=9, facecolors="none", edgecolors="k", lw=0.6, zorder=4)
    ax.annotate(g, (x0, y0), xytext=(tx, ty), textcoords="data", fontsize=6.5, style="italic", ha="center", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.45, color="k", shrinkA=0, shrinkB=2), zorder=5)
ax.set_xlim(-4.6, 4.9); ax.set_ylim(0, 9.3)
ax.set_xlabel("log$_2$ fold change, advanced vs early DKD"); ax.set_ylabel("$-$log$_{10}$ P")
lg = ax.legend(frameon=False, loc="upper left", fontsize=6, markerscale=2.5, handletextpad=0.1, borderaxespad=0.1)
fletter(fig, ax, "A", dx=-0.075)
# B heatmap of representative Tier-1 genes
t1 = cand[cand.tier.str.startswith("Tier1")].copy()
t1["score"] = A.loc[t1.index, "ADVvsEARLY_P"].rank() + A.loc[t1.index, "glom_meta_FDR"].rank() + A.loc[t1.index, "tub_meta_FDR"].rank()
top = list(t1.sort_values("score").index[:22]) + [g for g in ["CCL19", "CCL21", "IL7R", "TNFRSF1B", "EGF"] if g in t1.index]
top = list(dict.fromkeys(top))[:26]
M = A.loc[top, ["glom_meta_g", "tub_meta_g", "ADVvsEARLY_g"]].astype(float)
M = M.loc[M.mean(axis=1).sort_values(ascending=False).index]
ax = fig.add_subplot(gs[0, 1])
im = ax.imshow(M.clip(-3, 3).values, cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto")
ax.set_yticks(range(len(M))); ax.set_yticklabels(M.index, fontsize=6.3, style="italic")
ax.set_xticks(range(3)); ax.set_xticklabels(["Glom.\nDKD vs\ncontrol", "Tub.-int.\nDKD vs\ncontrol", "Biopsy\nadvanced\nvs early"], fontsize=6.3)
ax.tick_params(length=0); [s_.set_visible(False) for s_ in ax.spines.values()]
ax.set_xticks(np.arange(-0.5, 3, 1), minor=True); ax.set_yticks(np.arange(-0.5, len(M), 1), minor=True)
ax.grid(which="minor", color="white", lw=0.4); ax.tick_params(which="minor", length=0)
bb = ax.get_position(); cax = fig.add_axes([bb.x1 + 0.015, bb.y0 + bb.height * 0.3, 0.012, bb.height * 0.4])
cb = fig.colorbar(im, cax=cax, ticks=[-3, -2, -1, 0, 1, 2, 3]); cb.set_label("Hedges' g", fontsize=7); cb.ax.tick_params(labelsize=6.3, length=2); cb.outline.set_linewidth(0.5)
fletter(fig, ax, "B", dx=-0.11)
# C enrichment of progression-up genes
e = pd.read_csv(R + "ENRICH_progression_up.csv")
pick = ["Immune System", "Extracellular matrix organization", "Integrin cell surface interactions", "Cytokine-cytokine receptor interaction",
        "Neutrophil extracellular trap formation", "ECM-receptor interaction", "Complement and coagulation cascades", "Phagosome",
        "Signaling by Interleukins", "Degradation of the extracellular matrix"]
ee = e[e.term.isin(pick)].drop_duplicates("term").sort_values("P_adj").iloc[::-1]
gsC = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, :], width_ratios=[0.62, 1], wspace=0.0)
ax = fig.add_subplot(gsC[0, 1])
src_col = {"REAC": C_DKD, "KEGG": "#E69F00"}
ax.barh(range(len(ee)), -np.log10(ee.P_adj), color=[src_col.get(s_, GREY) for s_ in ee.source], height=0.65)
ax.set_yticks(range(len(ee))); ax.set_yticklabels(list(ee.term), fontsize=6.5)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=c_, label={"REAC": "Reactome", "KEGG": "KEGG"}[s_]) for s_, c_ in src_col.items()], frameon=False, loc="lower right", fontsize=6.3)
ax.set_xticks([0, 5, 10, 15, 20])
ax.set_xlabel("$-$log$_{10}$ adjusted P (progression-up genes)")
fig.text(0.015, ax.get_position().y1 + 0.015, "C", fontsize=10, fontweight="bold", va="bottom")
fig.savefig(OUT + "Figure1.pdf"); fig.savefig(OUT + "Figure1.png", dpi=300)
fig.savefig(OUT + "Figure1_tmp.png", dpi=600); Image.open(OUT + "Figure1_tmp.png").convert("RGB").save(OUT + "Figure1.tiff", compression="tiff_lzw", dpi=(600, 600)); os.remove(OUT + "Figure1_tmp.png"); plt.close(fig)

# ============================================================ Figure 2 (discovery snRNA)
fig = plt.figure(figsize=(175 * MM, 175 * MM))
gs = gridspec.GridSpec(3, 4, figure=fig, height_ratios=[1.45, 1, 1], hspace=0.6, wspace=0.75)
ax = fig.add_subplot(gs[0, :3]); img(ax, F + "umap_hires/UMAP_discovery_snRNA.png"); ax.set_anchor("W"); letter(ax, "A", x=-0.02, y=0.98)
comp = pd.read_csv(R + "S2_composition_per_donor.csv", index_col=0)
ax = fig.add_subplot(gs[0, 3]); letter(ax, "B", x=-0.45)
for i, (g, col) in enumerate([("CTRL", C_CTRL), ("DKD", C_DKD)]):
    v = comp.loc[comp.condition == g, "IMM"]
    ax.boxplot(v, positions=[i], widths=0.55, showfliers=False, medianprops=dict(color="k", lw=0.8), boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    ax.scatter(np.full(len(v), i) + np.linspace(-0.12, 0.12, len(v)), v, c=col, s=14, marker="o" if g == "CTRL" else "s", zorder=3)
cs = pd.read_csv(R + "S2_composition_DKDvsCTRL.csv").set_index("cell_type")
ax.set_xticks([0, 1]); ax.set_xticklabels(["Control", "DKD"]); ax.set_ylabel("Immune nuclei (% per donor)")
ax.set_title(f"P = {cs.at['IMM','P_MWU']:.3f} (Mann–Whitney)", fontsize=7)
# C immune subtypes, per donor (% of all nuclei)
from scipy import stats as _st
dd = pd.read_csv(R + "S2_immune_subtypes_pct_of_all_nuclei_per_donor.csv")
ax = fig.add_subplot(gs[1, :2]); letter(ax, "C", x=-0.14)
order = ["T cell", "B cell", "Plasma", "NK", "Macrophage", "Monocyte"]
rng_ = np.random.default_rng(3)
for i, c in enumerate(order):
    for off, g, col, mk in [(-0.18, "CTRL", C_CTRL, "o"), (0.18, "DKD", C_DKD, "s")]:
        v = dd.loc[dd.condition == g, c].values
        ax.scatter(i + off + rng_.uniform(-0.06, 0.06, len(v)), v, s=12, c=col, marker=mk, lw=0, zorder=3)
        ax.hlines(np.median(v), i + off - 0.13, i + off + 0.13, color="k", lw=0.9, zorder=4)
    p_ = _st.mannwhitneyu(dd.loc[dd.condition == "DKD", c], dd.loc[dd.condition == "CTRL", c]).pvalue
    ax.text(i, 1.62, f"{p_:.3f}", ha="center", fontsize=6)
ax.text(-0.5, 1.62, "P =", ha="right", fontsize=6)
ax.set_yscale("function", functions=(lambda v: np.sqrt(np.clip(v, 0, None)), lambda v: v ** 2))
ax.set_yticks([0, 0.05, 0.25, 0.5, 1.0, 1.5]); ax.set_yticklabels(["0", "0.05", "0.25", "0.5", "1.0", "1.5"]); ax.set_ylim(0, 1.75)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=30, ha="right"); ax.set_xlim(-1.15, len(order) - 0.4)
ax.set_ylabel("% of all nuclei per donor")
from matplotlib.lines import Line2D as _L2
ax.legend(handles=[_L2([], [], ls="", marker="o", ms=4, mfc=C_CTRL, mec=C_CTRL, label="Control (6 donors)"),
                   _L2([], [], ls="", marker="s", ms=4, mfc=C_DKD, mec=C_DKD, label="DKD (5 donors)")], frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.9), fontsize=6.5)
# D signature by cell type
sg = pd.read_csv(R + "S2_signature_by_celltype.csv")
ax = fig.add_subplot(gs[1, 2:]); letter(ax, "D", x=-0.12)
up = sg[sg.signature == "Prog_up_T1"].set_index("ct"); dn = sg[sg.signature == "Prog_down_T1"].set_index("ct")
cts = ["IMM", "FIB", "PEC", "PT", "EC", "POD", "tL", "TAL", "DCT", "PC"]
cts = [c for c in cts if c in up.index]
x = np.arange(len(cts))
ax.bar(x - 0.2, up.loc[cts, "diff"], 0.4, color=C_DKD, label="Progression-up score")
ax.bar(x + 0.2, dn.loc[cts, "diff"], 0.4, color=C_CTRL, label="Progression-down score")
for i, c in enumerate(cts):
    for off, t in [(-0.2, up), (0.2, dn)]:
        p = t.at[c, "P_t"]; v = t.at[c, "diff"]
        if p < 0.05: ax.text(i + off, v + (0.008 if v > 0 else -0.02), pstar(p), ha="center", fontsize=7)
ax.axhline(0, c="k", lw=0.6); ax.set_xticks(x); ax.set_xticklabels(cts, rotation=0); ax.set_ylabel("DKD − control\n(donor-mean module score)")
ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=6, handlelength=1.2, columnspacing=0.8, borderaxespad=0.2)
# E TF activities
tf = pd.read_csv(R + "S2_TF_DKDvsCTRL.csv")
gsb = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[2, :], wspace=0.95, width_ratios=[1, 1])
ax = fig.add_subplot(gsb[0, 0]); letter(ax, "E", x=-0.4)
rows = [("IMM", "RELB"), ("IMM", "NFKB2"), ("IMM", "IRF1"), ("IMM", "NFKB1"), ("IMM", "RELA"), ("PT", "PPARA"), ("EC", "PPARA"), ("PEC", "PPARA")]
vals = [tf[(tf.ct == c) & (tf.TF == t)].iloc[0] for c, t in rows]
ax.barh(range(len(rows))[::-1], [v.t for v in vals], color=[C_DKD if v.t > 0 else C_CTRL for v in vals], height=0.6)
ax.set_yticks(range(len(rows))[::-1]); ax.set_yticklabels([f"{t} ({c})" for c, t in rows])
for i, v in enumerate(vals): ax.text(v.t + (0.2 if v.t > 0 else -0.2), len(rows) - 1 - i, f"P = {v.P:.3f}", va="center", ha="left" if v.t > 0 else "right", fontsize=6)
ax.axvline(0, c="k", lw=0.6); ax.set_xlabel("TF activity, t statistic (DKD vs control)"); ax.set_xlim(-8, 10)
# F gained L-R fibroblast -> immune
g = pd.read_csv(R + "S2_liana_gained_in_DKD.csv")
g = g[(g.source == "FIB") & (g.target == "IMM")].sort_values("delta_lr_means", ascending=False).head(8)
ax = fig.add_subplot(gsb[0, 1]); letter(ax, "F", x=-0.6)
ax.barh(range(len(g))[::-1], g.delta_lr_means, color="#009E73", height=0.6)
ax.set_yticks(range(len(g))[::-1]); ax.set_yticklabels([f"{l} → {r}" for l, r in zip(g.ligand_complex, g.receptor_complex)], fontsize=6.5, style="italic")
ax.set_xlabel("Δ mean ligand–receptor expression\n(fibroblast → immune, DKD − control)")
save(fig, "Figure2")

# ============================================================ Figure 3 (lymphotoxin axis in bulk tissue + genetic triangulation)
from matplotlib.lines import Line2D
fig = plt.figure(figsize=(175 * MM, 175 * MM))
ax = fig.add_axes([0.085, 0.56, 0.205, 0.40])
genes = ["LTB", "LTBR", "MAP3K14", "RELB", "NFKB2", "CCL19", "CCL21", "CCR7", "IL7R", "CD3E", "MS4A1", "JCHAIN"]
cols = [("glom_meta_g", "glom_meta_FDR", "Glom.\nDKD vs\ncontrol"), ("tub_meta_g", "tub_meta_FDR", "Tub.-int.\nDKD vs\ncontrol"), ("ADVvsEARLY_g", "ADVvsEARLY_FDR", "Advanced\nvs early\nDKD")]
def dsize(q): return 8 + 62 * min(-np.log10(max(q, 1e-10)), 6) / 6
for j, (gc, fc, lab) in enumerate(cols):
    for i, gname in enumerate(genes):
        if gname not in A.index: continue
        v = A.at[gname, gc]; q = A.at[gname, fc]
        if pd.isna(v): ax.scatter(j, i, marker="x", s=12, c=GREY, lw=0.6); continue
        ax.scatter(j, i, s=dsize(q), c=[plt.cm.RdBu_r((np.clip(v, -3, 3) + 3) / 6)], edgecolors="k" if q < 0.05 else "none", linewidths=0.6)
ax.set_xticks(range(3)); ax.set_xticklabels([c[2] for c in cols], fontsize=6)
ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, style="italic", fontsize=6.8)
ax.set_ylim(len(genes) - 0.5, -0.5); ax.set_xlim(-0.6, 2.6)
fig.text(0.01, 0.965, "A", fontsize=10, fontweight="bold", va="bottom")
cax = fig.add_axes([0.315, 0.80, 0.010, 0.13])
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-3, 3)); cb = fig.colorbar(sm, cax=cax, ticks=[-3, 0, 3])
cb.ax.tick_params(labelsize=6, length=2); cb.outline.set_linewidth(0.5); cb.set_label("Hedges' g", fontsize=6.3, labelpad=1)
hs = [Line2D([], [], ls="", marker="o", ms=np.sqrt(dsize(10 ** -k)), mfc="#BDBDBD", mec="k" if k > 1.3 else "none", mew=0.6, label=l)
      for k, l in [(1, "0.1"), (2, "0.01"), (6, "≤10$^{-6}$")]] + [Line2D([], [], ls="", marker="x", ms=4, mec=GREY, label="not tested")]
fig.legend(handles=hs, title="FDR", loc="upper left", bbox_to_anchor=(0.305, 0.76), frameon=False, fontsize=6, title_fontsize=6, handletextpad=0.3, labelspacing=0.6, alignment="left")
# B: MHC proteins - DN effect tracks T1D effect
S5 = pd.read_csv(R + "S5_MR_MHC_vs_T1D.tsv", sep="\t").dropna(subset=["FG_DN", "FG_T1D"])
ax = fig.add_axes([0.53, 0.56, 0.42, 0.40])
from scipy import stats as _st
rho = _st.spearmanr(S5.FG_T1D, S5.FG_DN); sl = np.polyfit(S5.FG_T1D, S5.FG_DN, 1)
xx = np.linspace(S5.FG_T1D.min() - 0.2, S5.FG_T1D.max() + 0.2, 10); ax.plot(xx, np.polyval(sl, xx), c="#606060", lw=0.7, ls="--")
sig = S5.FG_DN_P < 0.05
ax.scatter(S5.FG_T1D[~sig], S5.FG_DN[~sig], s=14, c=GREY, lw=0, zorder=3)
ax.scatter(S5.FG_T1D[sig], S5.FG_DN[sig], s=16, c="#E69F00", lw=0, zorder=3)
l = S5[S5.Gene == "LTB"].iloc[0]; ax.scatter([l.FG_T1D], [l.FG_DN], s=30, c=C_DKD, edgecolors="k", lw=0.6, zorder=4)
offs = {"LTB": (8, 6), "AGER": (6, -8), "CFB": (6, -6), "HCG22": (6, -2), "BTN2A1": (-8, 7), "TNXB": (6, -5), "ATP6V1G2": (6, -4)}
for _, r in S5.iterrows():
    if r.Gene in offs:
        ax.annotate(r.Gene, (r.FG_T1D, r.FG_DN), xytext=offs[r.Gene], textcoords="offset points", fontsize=6.2, style="italic",
                    fontweight="bold" if r.Gene == "LTB" else "normal", ha="left" if offs[r.Gene][0] > 0 else "right")
ax.axhline(0, c="k", lw=0.4); ax.axvline(0, c="k", lw=0.4)
ax.set_xlabel("MR estimate for type 1 diabetes\n(log OR per SD protein, FinnGen R12)"); ax.set_ylabel("MR estimate for diabetic nephropathy\n(log OR per SD protein, FinnGen R12)")
ax.text(0.03, 0.97, f"{len(S5)} MHC-region Olink proteins\nSpearman ρ = {rho.correlation:.2f}, P = {rho.pvalue:.1e}".replace("e-0", "×10$^{-").replace("e-", "×10$^{-") + ("}$" if "e-" in f"{rho.pvalue:.1e}" else ""),
        transform=ax.transAxes, va="top", fontsize=6.3)
ax.legend(handles=[Line2D([], [], ls="", marker="o", ms=4, mfc="#E69F00", mec="#E69F00", label="DN P < 0.05"), Line2D([], [], ls="", marker="o", ms=4, mfc=GREY, mec=GREY, label="DN P ≥ 0.05"),
                   Line2D([], [], ls="", marker="o", ms=5, mfc=C_DKD, mec="k", label="LTB")], loc="lower right", frameon=False, fontsize=6)
fig.text(0.44, 0.965, "B", fontsize=10, fontweight="bold", va="bottom")
# C: forest of key genetic estimates
RR = pd.read_csv(R + "S5_MR_robustness_all.tsv", sep="\t")
OL = {"FG_DN": "Diabetic nephropathy", "FG_T2D": "Type 2 diabetes", "RAPID3": "Rapid eGFR decline", "CKDI25": "Incident CKD"}
spec = [("LTB", "rs2395470", "LTB, lead pQTL (rs2395470, MHC)", ["FG_DN", "FG_T2D"]),
        ("LTB", "rs3093553", "LTB, gene-proximal pQTL (rs3093553, MHC)", ["FG_DN", "FG_T2D"]),
        ("LTBR", "rs10849449", "LTBR, pQTL (rs10849449)", ["FG_DN", "RAPID3"]),
        ("CCL19", "rs10972202", "CCL19, pQTL (rs10972202)", ["FG_DN", "RAPID3"]),
        ("CCL21", "rs10814138", "CCL21, pQTL (rs10814138)", ["FG_DN", "RAPID3"]),
        ("RELB", "rs8111826", "RELB, blood cis-eQTL (rs8111826)", ["FG_DN", "RAPID3", "CKDI25"])]
rows = []
for gname, rs, lab, outs in spec:
    for o in outs:
        x = RR[(RR.Gene == gname) & (RR.rsid == rs) & (RR.outcome == o)]
        if len(x): x = x.iloc[0]; rows.append((lab, OL[o], x.wald, x.wald_se, x.wald_p))
ax = fig.add_axes([0.40, 0.10, 0.30, 0.38])
ypos, yy_ = [], 0; heads = []; prev = None
for lab, ol, b, se, p in rows:
    if lab != prev: heads.append((lab, yy_)); yy_ += 1; prev = lab
    ypos.append(yy_); yy_ += 1
ymax = yy_; ypos = [ymax - 1 - v for v in ypos]; heads = [(l, ymax - 1 - v) for l, v in heads]
ylab = []
for yy, (lab, ol, b, se, p) in zip(ypos, rows):
    c = C_DKD if p < 0.05 / len(rows) else ("#E69F00" if p < 0.05 else GREY)
    ax.errorbar(np.exp(b), yy, xerr=[[np.exp(b) - np.exp(b - 1.96 * se)], [np.exp(b + 1.96 * se) - np.exp(b)]], fmt="none", ecolor="k", elinewidth=0.7, capsize=1.5)
    ax.scatter(np.exp(b), yy, s=16, marker="s", c=c, zorder=3)
    ax.text(1.03, yy, f"{np.exp(b):.2f} ({np.exp(b - 1.96 * se):.2f}–{np.exp(b + 1.96 * se):.2f})", transform=ax.get_yaxis_transform(), va="center", fontsize=6.2)
    ax.text(1.50, yy, f"{p:.1e}".replace("e-0", "×10$^{-").replace("e-", "×10$^{-") + "}$" if p < 0.001 else f"{p:.3f}", transform=ax.get_yaxis_transform(), va="center", fontsize=6.2)
for lab, yh in heads:
    ax.text(-0.02, yh, lab, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.3, fontweight="bold")
y = ypos; ylab = [r[1] for r in rows]; rows_n = ymax
ax.text(1.03, rows_n - 0.3, "OR (95% CI)", transform=ax.get_yaxis_transform(), fontsize=6.3, fontweight="bold")
ax.text(1.50, rows_n - 0.3, "P", transform=ax.get_yaxis_transform(), fontsize=6.3, fontweight="bold")
ax.set_yticks(y); ax.set_yticklabels(ylab, fontsize=6.3); ax.set_ylim(-0.7, rows_n + 0.2)
ax.axvline(1, c="k", lw=0.6, ls="--"); ax.set_xscale("log")
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
ax.xaxis.set_major_locator(FixedLocator([0.25, 0.5, 1, 2, 4])); ax.xaxis.set_major_formatter(FixedFormatter(["0.25", "0.5", "1", "2", "4"])); ax.xaxis.set_minor_locator(NullLocator())
ax.set_xlim(0.18, 6); ax.set_xlabel("OR per SD of genetically predicted protein / expression")
ax.legend(handles=[Line2D([], [], ls="", marker="s", ms=4, mfc=c, mec=c, label=l) for c, l in [(C_DKD, f"P < {0.05/len(rows):.4f} (Bonferroni)"), ("#E69F00", "P < 0.05"), (GREY, "P ≥ 0.05")]],
          loc="upper center", bbox_to_anchor=(0.6, -0.15), ncol=3, frameon=False, fontsize=6.2, handletextpad=0.2, columnspacing=1.0)
fig.text(0.01, 0.48, "C", fontsize=10, fontweight="bold", va="bottom")
fig.savefig(OUT + "Figure3.pdf"); fig.savefig(OUT + "Figure3.png", dpi=300)
fig.savefig(OUT + "Figure3_tmp.png", dpi=600); Image.open(OUT + "Figure3_tmp.png").convert("RGB").save(OUT + "Figure3.tiff", compression="tiff_lzw", dpi=(600, 600)); os.remove(OUT + "Figure3_tmp.png"); plt.close(fig)

# ============================================================ Figure 4 (validation)
fig = plt.figure(figsize=(175 * MM, 150 * MM))
gs = gridspec.GridSpec(2, 3, figure=fig, height_ratios=[1.1, 1], hspace=0.5, wspace=0.6)
ax = fig.add_subplot(gs[0, :2]); img(ax, F + "umap_hires/UMAP_validation_GSE209781.png"); ax.set_anchor("W"); letter(ax, "A", x=-0.02, y=0.98)
pct = pd.read_csv(R + "S3_key_genes_pct_by_ct.csv", index_col=0); mean = pd.read_csv(R + "S3_key_genes_mean_by_ct.csv", index_col=0)
gl = ["LTB", "LTBR", "MAP3K14", "RELB", "NFKB2", "CCL19", "CCL21", "CCR7", "IL7R", "CD44"]
gl = [g for g in gl if g in pct.columns]
cts = ["T", "B", "NK", "Myeloid", "Mast", "FIB/MES", "EC", "PT", "TAL", "CD-PC"]
cts = [c for c in cts if c in pct.index]
ax = fig.add_subplot(gs[0, 2]); letter(ax, "B", x=-0.45)
Z = mean.loc[cts, gl]; Z = (Z - Z.min()) / (Z.max() - Z.min() + 1e-9)
for i, c in enumerate(cts):
    for j, gname in enumerate(gl):
        ax.scatter(j, i, s=pct.at[c, gname] / 100 * 60 + 0.5, c=[plt.cm.Reds(0.15 + 0.85 * Z.at[c, gname])], lw=0)
ax.set_xticks(range(len(gl))); ax.set_xticklabels(gl, rotation=90, style="italic"); ax.set_yticks(range(len(cts))); ax.set_yticklabels(cts)
ax.set_ylim(len(cts) - 0.5, -0.5)
for v in (10, 50, 100): ax.scatter([], [], s=v / 100 * 60, c="grey", label=f"{v}%")
ax.legend(title="% expressing", frameon=False, bbox_to_anchor=(1.0, 1.0), loc="upper left", fontsize=6, title_fontsize=6)
bb = ax.get_position(); cax = fig.add_axes([bb.x1 + 0.03, bb.y0 + 0.02, 0.012, bb.height * 0.32])
cb = fig.colorbar(plt.cm.ScalarMappable(cmap=matplotlib.colors.LinearSegmentedColormap.from_list("r", [plt.cm.Reds(0.15), plt.cm.Reds(1.0)]), norm=plt.Normalize(0, 1)), cax=cax, ticks=[0, 1])
cb.ax.set_yticklabels(["min", "max"], fontsize=6); cb.ax.tick_params(length=2); cb.outline.set_linewidth(0.5); cb.set_label("Scaled mean", fontsize=6, labelpad=1)
DL = pd.read_csv(R + "S6_donor_level_values.csv"); PRE = pd.read_csv(R + "S6_prespecified_axis_tests.csv")
panels = [("LTB", "T"), ("LTB", "B"), ("LTB", "NK"), ("MAP3K14", "FIB/MES"), ("RELB", "FIB/MES"), ("CCL21", "FIB/MES")]
ax = fig.add_subplot(gs[1, :2]); letter(ax, "C", x=-0.1)
for i, (g, c) in enumerate(panels):
    for off, cond_, col, mk in [(-0.17, "CTRL", C_CTRL, "o"), (0.17, "DKD", C_DKD, "s")]:
        v = DL[(DL.feature == g) & (DL.ct == c) & (DL.condition == cond_)]["mean"].values
        v = np.clip(v, 0.002, None)
        ax.scatter(np.full(len(v), i + off) + np.array([-0.05, 0, 0.05])[:len(v)], v, s=14, c=col, marker=mk, lw=0, zorder=3)
        ax.hlines(np.exp(np.mean(np.log(v))), i + off - 0.12, i + off + 0.12, color="k", lw=0.9, zorder=4)
    r = PRE[(PRE.gene == g) & (PRE.ct == c)].iloc[0]
    ax.text(i, 9, f"P = {r.P:.3f}\nFDR = {r.FDR_prespecified:.3f}", ha="center", fontsize=5.8, va="bottom")
ax.set_yscale("log"); ax.set_ylim(0.0015, 25)
ax.set_xticks(range(len(panels))); ax.set_xticklabels([f"$\\it{{{g}}}$\n{c}" for g, c in panels])
ax.set_ylabel("Donor-mean expression\n(log-normalized; log axis)")
from matplotlib.lines import Line2D as _L2
ax.legend(handles=[_L2([], [], ls="", marker="o", ms=4, mfc=C_CTRL, mec=C_CTRL, label="Control (n = 3)"), _L2([], [], ls="", marker="s", ms=4, mfc=C_DKD, mec=C_DKD, label="DKD (n = 3)")],
          frameon=False, loc="lower left", ncol=2, fontsize=6.3)
SG = pd.read_csv(R + "S6_signatures_noaxis.csv")
ax = fig.add_subplot(gs[1, 2]); letter(ax, "D", x=-0.45)
items = [("TLS_full", "T"), ("TLS_noaxis", "T"), ("TLS_full", "B"), ("TLS_noaxis", "B"), ("TLS_full", "NK"), ("TLS_noaxis", "NK"), ("Stromal_organiser", "FIB/MES")]
xs, cols_, labs = [], [], []
pos = [0, 0.8, 2, 2.8, 4, 4.8, 6.2]
for xpos, (sgn, c) in zip(pos, items):
    r = SG[(SG.signature == sgn) & (SG.ct == c)].iloc[0]
    fc = {"TLS_full": "#9E9E9E", "TLS_noaxis": C_DKD, "Stromal_organiser": "#0072B2"}[sgn]
    ax.bar(xpos, r["diff"], width=0.7, color=fc, edgecolor="k" if r.P_t < 0.05 else "none", lw=0.6)
    ax.text(xpos, max(r["diff"], 0) + 0.008, f"{r.P_t:.2f}", ha="center", fontsize=5.6, rotation=90, va="bottom")
ax.set_xticks([0.4, 2.4, 4.4, 6.2]); ax.set_xticklabels(["T", "B", "NK", "FIB/\nMES"]); ax.axhline(0, c="k", lw=0.6)
ax.set_ylabel("Module score, DKD − control\n(Welch P above bars)"); ax.set_ylim(-0.02, 0.46); ax.set_yticks([0, 0.1, 0.2, 0.3])
from matplotlib.patches import Patch as _P
ax.legend(handles=[_P(color="#9E9E9E", label="TLS score (all genes)"), _P(color=C_DKD, label="TLS score without\naxis genes"), _P(color="#0072B2", label="Stromal organiser")],
          frameon=False, fontsize=5.8, loc="upper left", handlelength=1.0, borderaxespad=0.1)
save(fig, "Figure4")
print("done")
