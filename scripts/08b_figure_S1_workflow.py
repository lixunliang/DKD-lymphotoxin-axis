# Additional file 2: Figure S1 - study workflow
import os, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper_figures") + os.sep
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42})
W, H = 174, 150
fig = plt.figure(figsize=(W / 25.4, H / 25.4)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
FILL = {"data": "#EAF2FA", "method": "#F2F2F2", "sc": "#FDEBDD", "mr": "#E6F4EA", "out": "#FFF4D6"}
def box(x, y, w, h, title, body, kind):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.5", fc=FILL[kind], ec="#404040", lw=0.6))
    n = body.count("\n") + 1
    ax.text(x + w / 2, y + h - 2.2, title, ha="center", va="top", fontsize=7, fontweight="bold")
    ax.text(x + w / 2, y + h - 7.0, body, ha="center", va="top", fontsize=6.3, linespacing=1.35)
def arr(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", lw=0.7, color="#404040", shrinkA=0, shrinkB=0, mutation_scale=7))
def line(xs, ys): ax.plot(xs, ys, c="#404040", lw=0.7, solid_capstyle="butt")
# row 1: data sources
y1, h1 = 128, 18
box(2, y1, 54, h1, "Glomerular microarrays", "GSE30122, GSE104948 (2 platforms),\nGSE96804\nDKD vs control, 4 comparisons", "data")
box(60, y1, 54, h1, "Tubulointerstitial microarrays", "GSE30122, GSE104954\nDKD vs control, 2 comparisons", "data")
box(118, y1, 54, h1, "Whole-biopsy RNA-seq", "GSE142025\nadvanced (n = 21) vs early (n = 6) DKD\nand vs control (n = 9)", "data")
# row 2: statistics
y2, h2 = 104, 12
box(2, y2, 112, h2, "Random-effects meta-analysis", "moderated t / Hedges' g per dataset; DerSimonian–Laird pooling; BH FDR", "method")
box(118, y2, 54, h2, "Differential expression", "moderated t, advanced vs early DKD", "method")
arr(29, y1, 29, y2 + h2); arr(87, y1, 87, y2 + h2); arr(145, y1, 145, y2 + h2)
# row 3: concordance
y3, h3 = 80, 12
box(32, y3, 110, h3, "Concordance filter across compartments and stage", "79 Tier-1 and 166 Tier-2 progression genes (Additional file 1)", "method")
arr(58, y2, 58, y3 + h3); line([145, 145], [y2, y3 + h3 / 2])
ax.annotate("", xy=(142, y3 + h3 / 2), xytext=(145, y3 + h3 / 2), arrowprops=dict(arrowstyle="-|>", lw=0.7, color="#404040", shrinkA=0, shrinkB=0, mutation_scale=7))
# bus to three branches
yb = 72
line([87, 87], [y3, yb]); line([29, 145], [yb, yb])
y4, h4 = 46, 20
for xc in (29, 87, 145): arr(xc, yb, xc, y4 + h4)
box(2, y4, 54, h4, "Functional annotation", "g:Profiler enrichment\n(Reactome, KEGG)\nHPA secretome annotation", "method")
box(60, y4, 54, h4, "Single-cell localisation", "Discovery snRNA-seq (CELLxGENE)\n5 DKD / 6 control; 39,176 nuclei\nValidation scRNA-seq (GSE209781)\n3 DKD / 3 control; 59,284 cells", "sc")
box(118, y4, 54, h4, "Protein Mendelian randomization", "cis-pQTL: UKB-PPP Olink (1,417)\ndeCODE SomaScan (1,389)\n→ FinnGen R12 diabetic nephropathy", "mr")
y5, h5 = 22, 12
box(60, y5, 54, h5, "In silico perturbation", "GRN virtual knockout / over-expression", "sc")
arr(87, y4, 87, y5 + h5)
# convergence
y6, h6 = 2, 12
box(2, y6, 170, h6, "Convergent axis", "lymphocyte LTB → stromal LTβR → NIK (MAP3K14) → RELB/NFKB2 → CCL19/CCL21 → tertiary lymphoid organisation", "out")
arr(29, y4, 29, y6 + h6); arr(87, y5, 87, y6 + h6); arr(145, y4, 145, y6 + h6)
fig.savefig(OUT + "FigureS1.pdf"); fig.savefig(OUT + "FigureS1.png", dpi=300)
fig.savefig(OUT + "FigureS1_tmp.png", dpi=600); Image.open(OUT + "FigureS1_tmp.png").convert("RGB").save(OUT + "FigureS1.tiff", compression="tiff_lzw", dpi=(600, 600)); os.remove(OUT + "FigureS1_tmp.png")
