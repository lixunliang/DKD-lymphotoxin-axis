# -*- coding: utf-8 -*-
"""
In-depth single-nucleus analysis of human DKD kidney (discovery cohort, CELLxGENE h5ad)
  2a  re-processing (normalise, HVG, PCA, Harmony by donor, UMAP)
  2b  donor-level cell composition DKD vs control
  2c  bulk-derived progression signature scored at single-cell level
  2d  injury states: failed-repair PT (FR-PT), myofibroblast, leukocyte subtypes
  2e  ligand-receptor communication DKD vs control (LIANA consensus)
  2f  transcription-factor activity (decoupler + CollecTRI) on donor pseudobulk
Outputs: results/S2_*.csv, figures/Fig5_* (PNG 300 dpi + PDF), processed h5ad for later steps
"""
import os, sys, glob, time, warnings, traceback
warnings.filterwarnings("ignore")
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repository root; scripts/ is one level down
RAW, OUT, FIG = [os.path.join(BASE, d) for d in ("raw", "results", "figures")]
for d in (OUT, FIG): os.makedirs(d, exist_ok=True)
LOGF = open(os.path.join(BASE, "log_04_discovery.txt"), "a", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOGF.write(s + "\n"); LOGF.flush()

import numpy as np, pandas as pd
import scanpy as sc, anndata as ad
from scipy import stats, sparse
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
sc.settings.figdir = FIG; sc.settings.verbosity = 1

ABBR = {"epithelial cell of proximal tubule": "PT", "kidney loop of Henle thick ascending limb epithelial cell": "TAL",
        "kidney loop of Henle thin ascending limb epithelial cell": "tL", "kidney distal convoluted tubule epithelial cell": "DCT",
        "renal principal cell": "PC", "renal alpha-intercalated cell": "IC-A", "renal beta-intercalated cell": "IC-B",
        "endothelial cell": "EC", "fibroblast": "FIB", "mesangial cell": "MES", "parietal epithelial cell": "PEC",
        "podocyte": "POD", "leukocyte": "IMM"}
COND_COL = {"CTRL": "#4C72B0", "DKD": "#C44E52"}

def savefig(fig, name):
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight"); plt.close(fig)

def bh(p):
    p = np.asarray(p, float); out = np.full(len(p), np.nan); ok = np.isfinite(p)
    if ok.sum() == 0: return out
    q = p[ok]; o = np.argsort(q); r = np.empty(len(q)); r[o] = np.arange(1, len(q) + 1)
    v = np.minimum.accumulate((q * len(q) / r)[o][::-1])[::-1]; w = np.empty(len(q)); w[o] = np.minimum(v, 1)
    out[ok] = w; return out

# ------------------------------------------------------------------ load
def load():
    f = sorted(glob.glob(os.path.join(RAW, "sc_*.h5ad")), key=os.path.getsize)[-1]
    log("loading", f)
    a = sc.read_h5ad(f)
    if a.raw is not None:
        b = a.raw.to_adata(); b.obs = a.obs.copy()
    else:
        b = a.copy()
    sym = b.var["feature_name"].astype(str) if "feature_name" in b.var else b.var_names.to_series().astype(str)
    b.var["ensembl"] = b.var_names; b.var_names = sym.values; b.var_names_make_unique()
    b.obs["condition"] = np.where(b.obs["disease"].astype(str).str.contains("diabet", case=False), "DKD",
                         np.where(b.obs["disease"].astype(str) == "normal", "CTRL", "other"))
    b = b[b.obs.condition != "other"].copy()
    b.obs["ct"] = b.obs["cell_type"].astype(str).map(ABBR).fillna(b.obs["cell_type"].astype(str))
    b.obs["donor"] = b.obs["donor_id"].astype(str)
    b.obs["ct"] = b.obs["ct"].astype("category"); b.obs["condition"] = b.obs["condition"].astype("category")
    b.layers["counts"] = b.X.copy()
    log(f"  {b.n_obs} cells, {b.n_vars} genes; donors: {b.obs.groupby('condition').donor.nunique().to_dict()}")
    return b

def preprocess(b):
    sc.pp.normalize_total(b, target_sum=1e4); sc.pp.log1p(b)
    sc.pp.highly_variable_genes(b, n_top_genes=3000, batch_key="donor")
    sc.tl.pca(b, n_comps=40, use_highly_variable=True)
    try:
        sc.external.pp.harmony_integrate(b, "donor"); rep = "X_pca_harmony"
    except Exception as e:
        log("  harmony unavailable:", e); rep = "X_pca"
    sc.pp.neighbors(b, use_rep=rep, n_neighbors=20); sc.tl.umap(b, min_dist=0.3)
    return b

# ------------------------------------------------------------------ 2a/2b
def fig_umap(b):
    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    sc.pl.umap(b, color="ct", ax=ax[0], show=False, legend_loc="on data", legend_fontsize=6, frameon=False, title="Cell types")
    sc.pl.umap(b, color="condition", ax=ax[1], show=False, frameon=False, palette=COND_COL, title="Condition", size=2)
    savefig(fig, "Fig5A_umap")
    markers = {"PT": ["LRP2", "CUBN"], "tL": ["CRYAB"], "TAL": ["UMOD", "SLC12A1"], "DCT": ["SLC12A3"], "PC": ["AQP2"],
               "IC-A": ["SLC4A1"], "IC-B": ["SLC26A4"], "POD": ["NPHS1", "PODXL"], "PEC": ["CFH", "VCAM1"], "EC": ["PECAM1", "EMCN"],
               "MES": ["PDGFRB", "ITGA8"], "FIB": ["COL1A2", "DCN"], "IMM": ["PTPRC", "CD163"]}
    order = [c for c in markers if c in set(b.obs.ct)]
    gl = {k: [g for g in v if g in b.var_names] for k, v in markers.items() if k in order}
    dp = sc.pl.dotplot(b, gl, groupby="ct", categories_order=order, show=False, standard_scale="var", figsize=(8, 3.6))
    savefig(plt.gcf(), "Fig5B_marker_dotplot")

def composition(b):
    tab = pd.crosstab(b.obs.donor, b.obs.ct, normalize="index") * 100
    cond = b.obs.groupby("donor").condition.first()
    rows = []
    for c in tab.columns:
        x, y = tab.loc[cond == "DKD", c], tab.loc[cond == "CTRL", c]
        p = stats.mannwhitneyu(x, y, alternative="two-sided").pvalue if len(x) > 1 and len(y) > 1 else np.nan
        rows.append({"cell_type": c, "DKD_median_pct": x.median(), "CTRL_median_pct": y.median(),
                     "log2_ratio": np.log2((x.mean() + 0.01) / (y.mean() + 0.01)), "P_MWU": p})
    res = pd.DataFrame(rows); res["FDR"] = bh(res.P_MWU)
    res.to_csv(os.path.join(OUT, "S2_composition_DKDvsCTRL.csv"), index=False, encoding="utf-8-sig")
    tab.assign(condition=cond).to_csv(os.path.join(OUT, "S2_composition_per_donor.csv"), encoding="utf-8-sig")
    fig, ax = plt.subplots(figsize=(8, 3))
    cts = list(tab.columns); pos = np.arange(len(cts))
    for i, (g, off) in enumerate([("CTRL", -0.18), ("DKD", 0.18)]):
        d = tab.loc[cond == g]
        ax.boxplot([d[c] for c in cts], positions=pos + off, widths=0.3, showfliers=False,
                   patch_artist=True, boxprops=dict(facecolor=COND_COL[g], alpha=0.35), medianprops=dict(color="k"))
        for j, c in enumerate(cts):
            ax.scatter(np.full(len(d), pos[j] + off) + np.random.uniform(-0.05, 0.05, len(d)), d[c], s=6, c=COND_COL[g])
    for j, c in enumerate(cts):
        p = res.set_index("cell_type").at[c, "P_MWU"]
        if p < 0.05: ax.text(pos[j], tab[c].max() * 1.02, "*" if p >= 0.01 else "**", ha="center")
    ax.set_xticks(pos); ax.set_xticklabels(cts, rotation=45, ha="right"); ax.set_ylabel("% of nuclei per donor"); ax.set_yscale("symlog", linthresh=1)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=COND_COL[k], alpha=0.5) for k in ("CTRL", "DKD")], labels=["Control", "DKD"], frameon=False)
    savefig(fig, "Fig5C_composition")
    log("  composition:", res.sort_values("P_MWU").head(5).to_dict("records"))
    return res

# ------------------------------------------------------------------ helpers: donor-level comparison of a per-cell score
def donor_test(b, key, groupby="ct", min_cells=10):
    d = b.obs[[groupby, "donor", "condition", key]]
    m = d.groupby([groupby, "donor", "condition"])[key].agg(["mean", "size"]).reset_index()
    m = m[m["size"] >= min_cells]
    rows = []
    for c, g in m.groupby(groupby):
        x, y = g[g.condition == "DKD"]["mean"], g[g.condition == "CTRL"]["mean"]
        if len(x) >= 2 and len(y) >= 2:
            rows.append({groupby: c, "DKD": x.mean(), "CTRL": y.mean(), "diff": x.mean() - y.mean(),
                         "P_t": stats.ttest_ind(x, y, equal_var=False).pvalue, "nDKD": len(x), "nCTRL": len(y)})
    r = pd.DataFrame(rows)
    if len(r): r["FDR"] = bh(r.P_t)
    return r, m

# ------------------------------------------------------------------ 2c progression signature
def signature(b):
    cand = pd.read_csv(os.path.join(OUT, "progression_genes_tiered.csv"), index_col=0)
    t1 = cand[cand.tier.str.startswith("Tier1")]
    sets = {"Prog_up_T1": t1.index[t1.ADVvsEARLY_logFC > 0], "Prog_down_T1": t1.index[t1.ADVvsEARLY_logFC < 0],
            "Prog_up_T12": cand.index[cand.ADVvsEARLY_logFC > 0]}
    res = []
    for k, gl in sets.items():
        gl = [g for g in gl if g in b.var_names]
        sc.tl.score_genes(b, gl, score_name=k, use_raw=False)
        r, m = donor_test(b, k); r["signature"] = k; r["n_genes"] = len(gl); res.append(r)
        m.to_csv(os.path.join(OUT, f"S2_{k}_donor_means.csv"), index=False, encoding="utf-8-sig")
    R = pd.concat(res); R.to_csv(os.path.join(OUT, "S2_signature_by_celltype.csv"), index=False, encoding="utf-8-sig")
    # figure: mean score per cell type (bar) + DKD vs CTRL per cell type
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.2))
    mean_ct = b.obs.groupby("ct")["Prog_up_T1"].mean().sort_values(ascending=False)
    ax[0].bar(mean_ct.index, mean_ct.values, color="#8C8C8C"); ax[0].set_ylabel("Progression-up score (Tier 1)")
    ax[0].tick_params(axis="x", rotation=60)
    r = R[R.signature == "Prog_up_T1"].set_index("ct").reindex(mean_ct.index)
    ax[1].bar(r.index, r["diff"], color=[COND_COL["DKD"] if v > 0 else COND_COL["CTRL"] for v in r["diff"].fillna(0)])
    for i, (c, row) in enumerate(r.iterrows()):
        if row.P_t < 0.05: ax[1].text(i, row["diff"], "*", ha="center", va="bottom" if row["diff"] > 0 else "top")
    ax[1].axhline(0, c="k", lw=0.5); ax[1].set_ylabel("DKD - control (donor means)"); ax[1].tick_params(axis="x", rotation=60)
    savefig(fig, "Fig5D_progression_signature")
    sc.pl.umap(b, color=["Prog_up_T1", "Prog_down_T1"], cmap="Reds", frameon=False, show=False, size=3, vmax="p99")
    savefig(plt.gcf(), "Fig5E_signature_umap")
    log("  signature (up T1) highest in:", mean_ct.head(4).round(3).to_dict())
    return R

# ------------------------------------------------------------------ 2d injury states
STATES = {
    "FRPT": (["VCAM1", "DCDC2", "PROM1", "HAVCR1", "ITGB8", "CDH6", "TPM1", "SPP1", "CCL2"], "PT"),
    "MyoFIB": (["ACTA2", "POSTN", "COL1A1", "COL3A1", "FN1", "TAGLN", "COL1A2"], "FIB"),
}
IMM_SUB = {"Macrophage": ["CD163", "MRC1", "C1QA", "C1QB", "C1QC", "F13A1", "MS4A6A", "CD68"],
           "Monocyte": ["CD14", "FCN1", "S100A8", "S100A9", "LYZ", "VCAN"],
           "T cell": ["CD3E", "CD3D", "CD2", "IL7R", "CD247"],
           "NK": ["NKG7", "GNLY", "KLRD1", "PRF1"], "B cell": ["MS4A1", "CD79A", "BANK1"],
           "Plasma": ["JCHAIN", "MZB1", "IGKC"], "DC": ["CLEC10A", "FCER1A", "CD1C"], "Mast": ["CPA3", "TPSAB1", "KIT"]}

def states(b):
    out = {}
    for k, (gl, ct) in STATES.items():
        gl = [g for g in gl if g in b.var_names]
        sub = b[b.obs.ct == ct].copy()
        if sub.n_obs < 50: continue
        sc.tl.score_genes(sub, gl, score_name=k, use_raw=False)
        thr = np.percentile(sub.obs.loc[sub.obs.condition == "CTRL", k], 95)
        sub.obs[k + "_high"] = (sub.obs[k] > thr).astype(float) * 100
        r1, _ = donor_test(sub, k); r2, m2 = donor_test(sub, k + "_high", min_cells=20)
        r1["metric"] = k + " score"; r2["metric"] = k + " % high (>95th pct of control)"
        out[k] = pd.concat([r1, r2]); m2.to_csv(os.path.join(OUT, f"S2_{k}_pct_per_donor.csv"), index=False, encoding="utf-8-sig")
        b.obs.loc[sub.obs_names, k] = sub.obs[k]
    # immune subtypes by marker scores
    imm = b[b.obs.ct == "IMM"].copy()
    if imm.n_obs >= 30:
        for k, gl in IMM_SUB.items():
            gl = [g for g in gl if g in imm.var_names]
            if len(gl) >= 2: sc.tl.score_genes(imm, gl, score_name="s_" + k, use_raw=False)
        sc_cols = [c for c in imm.obs.columns if c.startswith("s_")]
        imm.obs["imm_sub"] = imm.obs[sc_cols].idxmax(axis=1).str[2:]
        imm.obs.loc[imm.obs[sc_cols].max(axis=1) <= 0, "imm_sub"] = "Unassigned"
        b.obs.loc[imm.obs_names, "imm_sub"] = imm.obs["imm_sub"]
        tab = pd.crosstab(imm.obs.imm_sub, imm.obs.condition); tab.to_csv(os.path.join(OUT, "S2_immune_subtypes.csv"), encoding="utf-8-sig")
        per = pd.crosstab([b.obs.donor, b.obs.condition], b.obs.imm_sub.fillna("non-immune"), normalize="index") * 100
        per.to_csv(os.path.join(OUT, "S2_immune_subtypes_pct_of_all_nuclei_per_donor.csv"), encoding="utf-8-sig")
        log("  immune subtypes:", tab.to_dict())
        gl = {k: [g for g in v if g in imm.var_names][:4] for k, v in IMM_SUB.items()}
        gl = {k: v for k, v in gl.items() if v}
        sc.pl.dotplot(imm, gl, groupby="imm_sub", show=False, standard_scale="var", figsize=(9, 2.8))
        savefig(plt.gcf(), "Fig5G_immune_subtypes")
    R = pd.concat(out.values()) if out else pd.DataFrame()
    R.to_csv(os.path.join(OUT, "S2_injury_states.csv"), index=False, encoding="utf-8-sig")
    # figure: per-donor % FR-PT and % MyoFIB
    fig, axs = plt.subplots(1, 2, figsize=(5, 2.6))
    for ax, k in zip(axs, ["FRPT", "MyoFIB"]):
        f = os.path.join(OUT, f"S2_{k}_pct_per_donor.csv")
        if not os.path.exists(f): continue
        m = pd.read_csv(f)
        for i, g in enumerate(["CTRL", "DKD"]):
            v = m[m.condition == g]["mean"]
            ax.boxplot(v, positions=[i], widths=0.5, showfliers=False)
            ax.scatter(np.full(len(v), i) + np.random.uniform(-0.08, 0.08, len(v)), v, c=COND_COL[g], s=12)
        rr = out[k]; rr = rr[rr.metric.str.contains("% high")]
        p = rr.P_t.iloc[0] if len(rr) else np.nan
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Control", "DKD"]); ax.set_title(f"{'FR-PT' if k=='FRPT' else 'Myofibroblast'}  P={p:.3g}")
        ax.set_ylabel("% of cell type")
    savefig(fig, "Fig5F_injury_states")
    return R

# ------------------------------------------------------------------ 2e ligand-receptor
def communication(b):
    try:
        import liana as li
    except Exception as e:
        log("  liana not available:", e); return None
    res = {}
    for g in ("CTRL", "DKD"):
        sub = b[b.obs.condition == g].copy()
        li.mt.rank_aggregate(sub, groupby="ct", resource_name="consensus", expr_prop=0.1, min_cells=10,
                             use_raw=False, n_perms=1000, verbose=False, key_added="liana")
        r = sub.uns["liana"].copy(); r["condition"] = g; res[g] = r
    allr = pd.concat(res.values()); allr.to_csv(os.path.join(OUT, "S2_liana_all.csv"), index=False, encoding="utf-8-sig")
    key = ["source", "target", "ligand_complex", "receptor_complex"]
    d = res["DKD"].set_index(key); c = res["CTRL"].set_index(key)
    cols = ["lr_means", "magnitude_rank", "specificity_rank"]
    m = d[cols].join(c[cols], lsuffix="_DKD", rsuffix="_CTRL", how="outer")
    m["lr_means_DKD"] = m["lr_means_DKD"].fillna(0); m["lr_means_CTRL"] = m["lr_means_CTRL"].fillna(0)
    m["delta_lr_means"] = m.lr_means_DKD - m.lr_means_CTRL
    m["sig_DKD"] = m.magnitude_rank_DKD < 0.05; m["sig_CTRL"] = m.magnitude_rank_CTRL < 0.05
    m = m.reset_index()
    m.to_csv(os.path.join(OUT, "S2_liana_DKDvsCTRL.csv"), index=False, encoding="utf-8-sig")
    gained = m[m.sig_DKD & ~m.sig_CTRL.fillna(False)].sort_values("delta_lr_means", ascending=False)
    gained.to_csv(os.path.join(OUT, "S2_liana_gained_in_DKD.csv"), index=False, encoding="utf-8-sig")
    # count of significant interactions per source-target (heatmaps)
    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    cts = sorted(set(m.source) | set(m.target))
    mats = []
    for g in ("CTRL", "DKD"):
        s = res[g][res[g].magnitude_rank < 0.05]
        M = pd.crosstab(s.source, s.target).reindex(index=cts, columns=cts).fillna(0); mats.append(M)
    for ax, M, t in zip(axs, mats + [mats[1] - mats[0]], ["Control", "DKD", "DKD - control"]):
        vmax = np.abs(M.values).max() or 1
        im = ax.imshow(M.values, cmap="RdBu_r" if "-" in t else "Reds", vmin=-vmax if "-" in t else 0, vmax=vmax)
        ax.set_xticks(range(len(cts))); ax.set_xticklabels(cts, rotation=90); ax.set_yticks(range(len(cts))); ax.set_yticklabels(cts)
        ax.set_title(t); ax.set_xlabel("receiver"); ax.set_ylabel("sender"); fig.colorbar(im, ax=ax, shrink=0.6)
    savefig(fig, "Fig6A_interaction_counts")
    # top gained interactions involving immune / fibroblast / PEC / PT
    focus = gained[gained.source.isin(["IMM", "FIB", "PEC", "PT", "EC"]) | gained.target.isin(["IMM", "FIB"])].head(30)
    if len(focus):
        fig, ax = plt.subplots(figsize=(5, 0.18 * len(focus) + 1))
        lab = focus.source + " → " + focus.target + ": " + focus.ligand_complex + "–" + focus.receptor_complex
        ax.barh(range(len(focus))[::-1], focus.delta_lr_means, color=COND_COL["DKD"])
        ax.set_yticks(range(len(focus))[::-1]); ax.set_yticklabels(lab, fontsize=6); ax.set_xlabel("Δ mean L–R expression (DKD − control)")
        savefig(fig, "Fig6B_gained_interactions")
    log(f"  L-R: significant DKD {int(m.sig_DKD.sum())}, control {int(m.sig_CTRL.fillna(False).sum())}, gained {len(gained)}")
    return m

# ------------------------------------------------------------------ 2f TF activity
def pseudobulk(b, min_cells=10):
    X = b.layers["counts"]; X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
    key = b.obs.ct.astype(str) + "|" + b.obs.donor.astype(str)
    rows, idx = [], []
    for k, ii in b.obs.groupby(key).indices.items():
        if len(ii) < min_cells: continue
        rows.append(np.asarray(X[ii].sum(0)).ravel()); idx.append(k)
    P = pd.DataFrame(np.vstack(rows), index=idx, columns=b.var_names)
    P = P.loc[:, (P > 0).mean() > 0.1]
    L = np.log1p(P.div(P.sum(1), axis=0) * 1e6)
    meta = pd.DataFrame({"ct": [i.split("|")[0] for i in idx], "donor": [i.split("|")[1] for i in idx]}, index=idx)
    meta["condition"] = meta.donor.map(b.obs.groupby("donor").condition.first())
    return L, meta

def tf_activity(b):
    try:
        import decoupler as dc
        net = dc.get_collectri(organism="human", split_complexes=False)
    except Exception as e:
        log("  decoupler/CollecTRI unavailable:", e); return None
    L, meta = pseudobulk(b)
    L.to_csv(os.path.join(OUT, "S2_pseudobulk_logCPM.csv.gz"), compression="gzip")
    meta.to_csv(os.path.join(OUT, "S2_pseudobulk_meta.csv"), encoding="utf-8-sig")
    est = dc.run_ulm(mat=L, net=net, source="source", target="target", weight="weight", verbose=False)
    acts = est[0] if isinstance(est, tuple) else est
    acts.to_csv(os.path.join(OUT, "S2_TF_activity_pseudobulk.csv.gz"), compression="gzip")
    rows = []
    for c in meta.ct.unique():
        mm = meta[meta.ct == c]
        dd, cc = mm.index[mm.condition == "DKD"], mm.index[mm.condition == "CTRL"]
        if len(dd) < 2 or len(cc) < 2: continue
        t, p = stats.ttest_ind(acts.loc[dd], acts.loc[cc], equal_var=False)
        rows.append(pd.DataFrame({"ct": c, "TF": acts.columns, "DKD": acts.loc[dd].mean().values, "CTRL": acts.loc[cc].mean().values,
                                  "t": t, "P": p}))
    R = pd.concat(rows); R["FDR"] = np.nan
    for c in R.ct.unique():
        i = R.ct == c; R.loc[i, "FDR"] = bh(R.loc[i, "P"].values)
    R.to_csv(os.path.join(OUT, "S2_TF_DKDvsCTRL.csv"), index=False, encoding="utf-8-sig")
    top = R[R.ct.isin(["PT", "FIB", "IMM", "PEC", "EC", "POD", "TAL"])].dropna(subset=["P"])
    sel = top.sort_values("P").groupby("ct").head(6).TF.unique()[:30]
    H = top[top.TF.isin(sel)].pivot_table(index="TF", columns="ct", values="t")
    fig, ax = plt.subplots(figsize=(0.45 * H.shape[1] + 2, 0.2 * H.shape[0] + 1))
    v = np.nanmax(np.abs(H.values)) if H.size else 1
    im = ax.imshow(H.values, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
    ax.set_xticks(range(H.shape[1])); ax.set_xticklabels(H.columns, rotation=90); ax.set_yticks(range(H.shape[0])); ax.set_yticklabels(H.index, fontsize=6)
    fig.colorbar(im, ax=ax, shrink=0.5, label="t (DKD vs control)"); savefig(fig, "Fig6C_TF_activity")
    log("  TF top hits:", R.sort_values("P").head(8)[["ct", "TF", "t", "P"]].round(4).to_dict("records"))
    return R

# ------------------------------------------------------------------ main
def main():
    log("=" * 70); log("discovery cohort analysis start", time.strftime("%Y-%m-%d %H:%M:%S"))
    log("versions: scanpy", sc.__version__, "anndata", ad.__version__)
    b = load(); b = preprocess(b)
    for name, fn in [("umap", fig_umap), ("composition", composition), ("signature", signature),
                     ("states", states), ("communication", communication), ("tf", tf_activity)]:
        try:
            log(f"-- {name}"); fn(b)
        except Exception as e:
            log(f"!! {name} failed: {e}"); log(traceback.format_exc())
    try:
        b.obs = b.obs.astype({c: str for c in b.obs.columns if b.obs[c].dtype == object})
        b.write_h5ad(os.path.join(RAW, "processed_discovery_snRNA.h5ad"), compression="gzip")
    except Exception as e:
        log("!! save failed:", e)
    log("discovery cohort analysis done", time.strftime("%Y-%m-%d %H:%M:%S"))

if __name__ == "__main__":
    main()
