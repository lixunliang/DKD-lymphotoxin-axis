# -*- coding: utf-8 -*-
"""
Single-cell robustness analyses
  R1  TLS evidence without circularity: module scores that EXCLUDE the lymphotoxin-axis genes
        TLS_noaxis   = TLS genes minus LTB, LTA, CCL19, CCL21, CXCL13, CCR7, IL7R      (lymphocytes)
        Stromal_organiser = VCAM1, ICAM1, MADCAM1, PDPN, CXCL12, TNFSF13B, IL7, CD34 (no CCL19/21)  (FIB/MES)
  R2  donor-level values (one row per donor x cell type) for key genes and signatures -> plotting every donor;
      BH correction across the pre-specified lymphotoxin-axis tests
  R3  stability of the virtual LTB knockout in validation-cohort lymphocytes:
        random seeds, network size (500 / 1,500 genes), principal components (2 / 5), cells per network (1,000),
        and an expression-matched random-knockout null
Inputs   raw/processed_GSE209781.h5ad (written by script 05), results/ALL_genes_integrated.csv
Outputs  results/S6_*.csv, log_10_robustness.txt     (about 60-90 min, mostly R3)
"""
import os, sys, time, importlib.util, traceback
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
spec = importlib.util.spec_from_file_location("st3", os.path.join(os.path.dirname(os.path.abspath(__file__)), "05_validation_cohort_and_virtual_knockout.py"))
st3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(st3)       # reuse QC'd functions (pcnet, propagate, donor_level ...)
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats, sparse
OUT = os.path.join(BASE, "results")
LOGF = open(os.path.join(BASE, "log_10_robustness.txt"), "a", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOGF.write(s + "\n"); LOGF.flush()

AXIS = ["LTB", "LTA", "CCL19", "CCL21", "CXCL13", "CCR7", "IL7R"]
TLS_NOAXIS = [g for g in st3.TLS if g not in AXIS]
STROMAL_ORG = ["VCAM1", "ICAM1", "MADCAM1", "PDPN", "CXCL12", "TNFSF13B", "IL7", "CD34"]
KEY = ["LTB", "LTA", "LTBR", "MAP3K14", "RELB", "NFKB2", "CCL19", "CCL21", "CXCL13", "CCR7", "IL7R", "CD44"]

def r1_r2(A):
    sets = {"TLS_full": st3.TLS, "TLS_noaxis": TLS_NOAXIS, "Stromal_organiser": STROMAL_ORG}
    sig_rows, donor_rows = [], []
    for k, gl in sets.items():
        gl = [g for g in gl if g in A.var_names]
        sc.tl.score_genes(A, gl, score_name=k, random_state=0)
        r, m = st3.donor_level(A, k); r["signature"] = k; r["genes"] = ",".join(gl); sig_rows.append(r)
        m["feature"] = k; donor_rows.append(m)
    gg = [g for g in KEY if g in A.var_names]
    X = A[:, gg].X; X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    kg = []
    for j, g in enumerate(gg):
        A.obs["g_" + g] = X[:, j]
        r, m = st3.donor_level(A, "g_" + g); r["gene"] = g; kg.append(r); m["feature"] = g; donor_rows.append(m)
    sig = pd.concat(sig_rows); sig.to_csv(os.path.join(OUT, "S6_signatures_noaxis.csv"), index=False, encoding="utf-8-sig")
    kg = pd.concat(kg)
    D = pd.concat(donor_rows); D.to_csv(os.path.join(OUT, "S6_donor_level_values.csv"), index=False, encoding="utf-8-sig")
    # pre-specified lymphotoxin-axis tests (ligand in lymphocytes, receptor-side genes in FIB/MES)
    pre = [("LTB", "T"), ("LTB", "B"), ("LTB", "NK"), ("LTBR", "FIB/MES"), ("MAP3K14", "FIB/MES"), ("RELB", "FIB/MES"),
           ("NFKB2", "FIB/MES"), ("CCL19", "FIB/MES"), ("CCL21", "FIB/MES")]
    rows = []
    for g, c in pre:
        x = kg[(kg.gene == g) & (kg.ct == c)]
        if len(x): rows.append(dict(gene=g, ct=c, DKD=x.DKD.iloc[0], CTRL=x.CTRL.iloc[0], P=x.P_t.iloc[0]))
    pre = pd.DataFrame(rows); pre["FDR_prespecified"] = st3.bh(pre.P.values)
    pre.to_csv(os.path.join(OUT, "S6_prespecified_axis_tests.csv"), index=False, encoding="utf-8-sig")
    log("R1 signatures (donor-level Welch):")
    for _, x in sig[sig.ct.isin(["T", "B", "NK", "Myeloid", "FIB/MES"])].iterrows():
        log(f"   {x.signature:<18} {x.ct:<8} diff={x['diff']:+.3f}  P={x.P_t:.3g}  (n {x.nDKD}/{x.nCTRL})")
    log("R2 pre-specified axis tests with BH FDR:"); log(pre.round(4).to_string(index=False))

def vko_variant(A, mask, gene, n_genes=1000, n_comp=3, n_cells=600, n_nets=3, seed=1, n_null=1000, matched_null=False):
    rng = np.random.default_rng(seed)
    sub = A[mask].copy(); sub.var_names_make_unique()
    if (sub.obs.condition == "DKD").sum() >= 300: sub = sub[sub.obs.condition == "DKD"].copy()
    vn = pd.Index(sub.var_names)
    good = ~(vn.str.startswith(("ENSG", "LINC", "MT-", "RPL", "RPS", "MIR", "SNOR")) | vn.str.contains(r"-AS\d|\.\d|^LOC|-DT$|^OR\d", regex=True))
    frac = np.asarray((sub.X > 0).mean(0)).ravel()
    always = [gene, "LTB", "LTBR", "RELB", "NFKB2", "MAP3K14", "CCL19", "CCL21", "CXCL13", "CD44", "IL7R"]
    sub = sub[:, (good & (frac >= 0.05)) | vn.isin(always)].copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=min(n_genes, sub.n_vars - 1))
    genes = list(sub.var_names[sub.var.highly_variable])
    genes = list(dict.fromkeys(genes + [g for g in always if g in sub.var_names and g not in genes]))
    X = sub[:, genes].X; X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    nets = []
    for k in range(n_nets):
        idx = rng.choice(X.shape[0], size=min(n_cells, X.shape[0]), replace=False)
        Xs = X[idx]; sd = Xs.std(0); sd[sd == 0] = 1; Xs = (Xs - Xs.mean(0)) / sd
        nets.append(st3.pcnet(Xs, n_comp=n_comp))
    W = np.mean(nets, axis=0); rad = np.max(np.abs(np.linalg.eigvals(W)))
    if rad >= 0.9: W = W * (0.9 / rad)
    g = genes.index(gene); s_ko = st3.propagate(W, g, -1.0)
    cand = [i for i in range(len(genes)) if i != g]
    if matched_null:      # random genes from the same mean-expression decile as the target
        mu = X.mean(0); dec = pd.qcut(pd.Series(mu).rank(method="first"), 10, labels=False).values
        cand = [i for i in cand if dec[i] == dec[g]]
    pick = rng.choice(cand, size=min(n_null, len(cand)), replace=False)
    null = np.vstack([st3.propagate(W, j, -1.0) for j in pick])
    fc = st3.load_progression_logfc().reindex(genes); ok = fc.notna().values & (np.arange(len(genes)) != g)
    rev = stats.spearmanr(s_ko[ok], fc.values[ok]).correlation
    nr = np.array([stats.spearmanr(n_[ok], fc.values[ok]).correlation for n_ in null])
    p = (np.sum(nr <= rev) + 1) / (len(nr) + 1)
    z = (s_ko - null.mean(0)) / (null.std(0) + 1e-12)
    order = np.argsort(z); top_down = [genes[i] for i in order[:12] if i != g][:10]
    rank = {gg: int(np.where(order == genes.index(gg))[0][0]) + 1 for gg in ("MAP3K14", "CCL19", "CCR7", "CD44") if gg in genes}
    return dict(rho=rev, P=p, n_cells=sub.n_obs, n_genes=len(genes), null_size=len(pick), top_down=",".join(top_down),
                **{f"rank_{k}": v for k, v in rank.items()})

def r3(A):
    lymph = A.obs.ct.isin(["T", "B", "Plasma", "NK"]).values
    grid = [("baseline (1,000 genes, 3 PCs, 600 cells, seed 1)", {}),
            ("seed 2", {"seed": 2}), ("seed 3", {"seed": 3}),
            ("500 genes", {"n_genes": 500}), ("1,500 genes", {"n_genes": 1500}),
            ("2 PCs", {"n_comp": 2}), ("5 PCs", {"n_comp": 5}),
            ("1,000 cells per network", {"n_cells": 1000}),
            ("expression-matched null", {"matched_null": True})]
    rows = []
    for lab, kw in grid:
        t0 = time.time()
        try:
            r = vko_variant(A, lymph, "LTB", **kw); r["setting"] = lab; rows.append(r)
            log(f"R3 {lab:<52} rho={r['rho']:+.3f}  P={r['P']:.3f}  ({(time.time()-t0)/60:.1f} min)  top down: {r['top_down']}")
        except Exception as e:
            log(f"R3 {lab}: failed {e}"); log(traceback.format_exc())
        pd.DataFrame(rows).to_csv(os.path.join(OUT, "S6_vKO_LTB_stability.csv"), index=False, encoding="utf-8-sig")

def main():
    log("=" * 78); log("single-cell robustness start", time.strftime("%Y-%m-%d %H:%M:%S"))
    A = sc.read_h5ad(os.path.join(BASE, "raw", "processed_GSE209781.h5ad"))
    log(f"  validation cohort: {A.n_obs} cells")
    try: r1_r2(A)
    except Exception as e: log("!! R1/R2 failed:", e); log(traceback.format_exc())
    try: r3(A)
    except Exception as e: log("!! R3 failed:", e); log(traceback.format_exc())
    log("single-cell robustness done", time.strftime("%Y-%m-%d %H:%M:%S"))

if __name__ == "__main__":
    main()
