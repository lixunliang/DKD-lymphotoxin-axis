# -*- coding: utf-8 -*-
"""High-resolution UMAP panels for Figure 2A (discovery snRNA-seq) and Figure 4A (validation scRNA-seq).
Reads the processed h5ad files written by scripts 04 and 05; writes PNG (600 dpi) + PDF to figures/umap_hires/."""
import os, sys
import numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
RAW = os.path.join(BASE, "raw"); OUT = os.path.join(BASE, "figures", "umap_hires"); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "pdf.fonttype": 42, "axes.linewidth": 0.8})
MM = 1 / 25.4
COND = {"CTRL": "#0072B2", "DKD": "#D55E00"}
# colour-blind-friendly categorical palette (Okabe-Ito + extensions)
PAL = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#882255", "#44AA99",
       "#117733", "#332288", "#AA4499", "#DDCC77", "#88CCEE"]

def place_labels(ax, xy, cats, cl):
    """Label each cluster at its median; push overlapping labels apart and add a leader line when a label moves."""
    rng = np.ptp(xy, axis=0); pos = {}; anc = {}
    for c in cl:
        m = np.asarray(cats) == c
        if m.sum() < 15: continue
        anc[c] = np.median(xy[m], axis=0); pos[c] = anc[c].copy()
        if m.mean() < 0.02:                     # rare cell type: move label sideways so the points stay visible
            side = -1.0 if anc[c][0] < np.median(xy[:, 0]) else 1.0
            lo, hi = xy[:, 0].min() + 0.07 * rng[0], xy[:, 0].max() - 0.07 * rng[0]
            if not lo <= anc[c][0] + side * 0.10 * rng[0] <= hi: side = -side   # keep the label inside the panel
            pos[c][0] += side * 0.10 * rng[0]
    keys = list(pos)
    for _ in range(300):                       # simple repulsion in axis-normalised units
        moved = False
        for i, a in enumerate(keys):
            for b in keys[i + 1:]:
                d = (pos[a] - pos[b]) / rng; need = np.array([0.14, 0.075])
                if abs(d[0]) < need[0] and abs(d[1]) < need[1]:
                    step = (np.sign(d[1]) if d[1] != 0 else 1.0) * 0.01 * rng[1]
                    pos[a][1] += step; pos[b][1] -= step; moved = True
        if not moved: break
    for c in keys:
        shift = np.linalg.norm((pos[c] - anc[c]) / rng)
        ax.annotate(c, anc[c], xytext=pos[c], textcoords="data", fontsize=7.5, fontweight="bold", ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8),
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="k", shrinkA=0, shrinkB=0) if shift > 0.02 else None)

def panel(ad, key, ax, title_txt, palette=None, label_on_data=True, size=None):
    xy = ad.obsm["X_umap"]; cats = ad.obs[key].astype("category"); cl = list(cats.cat.categories)
    pal = palette or {c: PAL[i % len(PAL)] for i, c in enumerate(cl)}
    lab = np.asarray(cats); frac = pd.Series(lab).value_counts(normalize=True)
    small = np.isin(lab, [c for c in cl if frac.get(c, 0) < 0.02]) if label_on_data else np.zeros(len(lab), bool)
    s = size or max(0.3, 60000 / len(xy) * 0.4)
    order = np.random.default_rng(0).permutation(np.where(~small)[0])      # random draw order so no group hides another
    ax.scatter(xy[order, 0], xy[order, 1], c=[pal[c] for c in lab[order]], s=s, lw=0, rasterized=True)
    if small.any():                                                         # rare cell types drawn on top and slightly larger
        ax.scatter(xy[small, 0], xy[small, 1], c=[pal[c] for c in lab[small]], s=s * 3, lw=0, rasterized=True)
    if label_on_data:
        place_labels(ax, xy, cats, cl)
    else:
        for c in cl: ax.scatter([], [], c=pal[c], s=20, label={"CTRL": "Control", "DKD": "DKD"}.get(c, c))
        ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=7.5, markerscale=1.2, handletextpad=0.2, borderaxespad=0.1)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.text(0.02, 0.98, title_txt, transform=ax.transAxes, fontsize=7, va="top", color="#404040")

def make(path, ct_key, cond_key, name, note):
    if not os.path.exists(path):
        print("missing", path); return
    ad = sc.read_h5ad(path)
    if "X_umap" not in ad.obsm:
        sc.pp.neighbors(ad); sc.tl.umap(ad)
    ad.obs[cond_key] = ad.obs[cond_key].astype(str)
    fig, axs = plt.subplots(1, 2, figsize=(120 * MM, 58 * MM))
    panel(ad, ct_key, axs[0], f"{ad.n_obs:,} {note}")
    panel(ad, cond_key, axs[1], "", palette=COND, label_on_data=False)
    fig.tight_layout(w_pad=1.5)
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=600); fig.savefig(os.path.join(OUT, name + ".pdf")); plt.close(fig)
    print("saved", name, ad.n_obs, "cells; groups:", ad.obs[ct_key].value_counts().to_dict())

make(os.path.join(RAW, "processed_discovery_snRNA.h5ad"), "ct", "condition", "UMAP_discovery_snRNA", "nuclei")
make(os.path.join(RAW, "processed_GSE209781.h5ad"), "ct", "condition", "UMAP_validation_GSE209781", "cells")
print("done")
