# -*- coding: utf-8 -*-
"""
Validation cohort and virtual perturbation
  3A  Independent validation in GSE209781 (scRNA-seq, 3 DKD vs 3 control kidneys)
      download -> QC -> Harmony -> marker-based annotation -> composition, TLS/lymphotoxin axis,
      progression signature, CD44, RELB/NFKB2 TF activity, LTB-LTBR and CD44 communication
  3B  Virtual knockout / over-expression by gene-regulatory-network signal propagation
      (PC-regression GRN per cell type, CellOracle-style 3-step propagation, gene-level null from
      random-gene knockouts, and a test of whether the perturbation reverses the bulk progression signature)
      KO targets: LTB and RELB (lymphocytes/immune), LTBR (fibroblast/stromal), CD44 (immune)
Outputs: results/S3_*.csv, figures/Fig7_*, Fig8_*
"""
import os, sys, re, io, gzip, glob, json, ssl, time, tarfile, warnings, traceback, urllib.request
warnings.filterwarnings("ignore")
try:
    import certifi; ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())
except Exception:
    pass
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repository root; scripts/ is one level down
RAW, OUT, FIG = [os.path.join(BASE, d) for d in ("raw", "results", "figures")]
VDIR = os.path.join(RAW, "GSE209781")
for d in (OUT, FIG, VDIR): os.makedirs(d, exist_ok=True)
LOGF = open(os.path.join(BASE, "log_05_validation.txt"), "a", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOGF.write(s + "\n"); LOGF.flush()

import numpy as np, pandas as pd
import scanpy as sc, anndata as ad
from scipy import stats, sparse
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
COND_COL = {"CTRL": "#4C72B0", "DKD": "#C44E52"}
UA = {"User-Agent": "Mozilla/5.0 (research script)"}

def savefig(fig, name):
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=300, bbox_inches="tight"); fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight"); plt.close(fig)

def bh(p):
    p = np.asarray(p, float); out = np.full(len(p), np.nan); ok = np.isfinite(p)
    if ok.sum() == 0: return out
    q = p[ok]; o = np.argsort(q); r = np.empty(len(q)); r[o] = np.arange(1, len(q) + 1)
    v = np.minimum.accumulate((q * len(q) / r)[o][::-1])[::-1]; w = np.empty(len(q)); w[o] = np.minimum(v, 1); out[ok] = w; return out

def fetch(url, dest, tries=30):
    if os.path.exists(dest) and os.path.getsize(dest) > 0: return dest
    part = dest + ".part"
    for k in range(tries):
        try:
            have = os.path.getsize(part) if os.path.exists(part) else 0
            hdr = dict(UA)
            if have: hdr["Range"] = f"bytes={have}-"
            log(f"  downloading {url}" + (f" (resume {have/1e6:.0f} MB)" if have else ""))
            with urllib.request.urlopen(urllib.request.Request(url, headers=hdr), timeout=300) as r:
                if have and r.status != 206: have = 0
                t0 = time.time()
                with open(part, "ab" if have else "wb") as f:
                    while True:
                        b = r.read(4 << 20)
                        if not b: break
                        f.write(b)
                        if time.time() - t0 > 30: log(f"    {os.path.getsize(part)/1e6:.0f} MB"); t0 = time.time()
            os.replace(part, dest); return dest
        except urllib.error.HTTPError as e:
            if e.code == 416 and os.path.exists(part): os.replace(part, dest); return dest
            log(f"  retry {k+1}: {e}");
            if e.code in (403, 404): break
            time.sleep(10)
        except Exception as e:
            log(f"  retry {k+1}: {e}"); time.sleep(10)
    raise RuntimeError("download failed: " + url)

def listdir(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
        html = r.read().decode("utf-8", "replace")
    return sorted(set(re.findall(r'href="([^"?/][^"]*)"', html)))

# ===================================================================== 3A validation
ACC = "GSE209781"
GEO = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{ACC[:-3]}nnn/{ACC}/"

def sample_meta():
    fn = [f for f in listdir(GEO + "matrix/") if f.endswith("_series_matrix.txt.gz")][0]
    p = fetch(GEO + "matrix/" + fn, os.path.join(VDIR, fn))
    meta, chars = {}, []
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"): break
            if line.startswith("!Sample_"):
                parts = [x.strip('"') for x in line.rstrip("\n").split("\t")]
                k = parts[0][8:]
                if k == "characteristics_ch1": chars.append(parts[1:])
                else: meta.setdefault(k, parts[1:])
    m = pd.DataFrame({"gsm": meta["geo_accession"], "title": meta.get("title"), "source": meta.get("source_name_ch1")})
    for i, c in enumerate(chars): m[f"char{i+1}"] = c
    txt = m.drop(columns="gsm").astype(str).agg(" | ".join, axis=1).str.lower()
    title = m["title"].astype(str).str.upper()
    m["condition"] = np.where(title.str.match(r"^(NM|NC|CTRL|CON|NOR)"), "CTRL",
                     np.where(title.str.match(r"^(DKD|DN|DM)"), "DKD",
                     np.where(txt.str.contains(r"normal|control|healthy"), "CTRL", "other")))
    m.to_csv(os.path.join(OUT, "S3_GSE209781_samples.csv"), index=False, encoding="utf-8-sig")
    log("  samples:", m[["gsm", "title", "condition"]].values.tolist())
    return m

def read_member_table(raw, name):
    if name.endswith(".gz"): raw = gzip.decompress(raw)
    sep = "\t" if raw[:5000].count(b"\t") > raw[:5000].count(b",") else ","
    d = pd.read_csv(io.BytesIO(raw), sep=sep, index_col=0)
    return d

def build_counts(meta):
    supp = listdir(GEO + "suppl/")
    log("  supplementary:", supp)
    adatas = []
    tars = [f for f in supp if f.endswith(".tar")]
    others = [f for f in supp if f.startswith(ACC) and not f.endswith(".tar") and re.search(r"\.(h5ad|h5|rds|txt|csv|tsv|mtx)(\.gz)?$", f)]
    if tars:
        tp = fetch(GEO + "suppl/" + tars[0], os.path.join(VDIR, tars[0]))
        class _M:   # minimal member wrapper holding bytes
            def __init__(self, name, data): self.name, self.data = name, data
        class _TF:
            def extractfile(self, m): return io.BytesIO(m.data)
        with tarfile.open(tp) as outer:
            members = []
            for m in outer.getmembers():
                if not m.isfile(): continue
                data = outer.extractfile(m).read()
                g = re.match(r"(GSM\d+)", os.path.basename(m.name))
                if re.search(r"\.(tar\.gz|tgz|tar)$", m.name.lower()):
                    with tarfile.open(fileobj=io.BytesIO(data)) as inner:
                        for im in inner.getmembers():
                            if im.isfile() and "raw_feature" not in im.name.lower():
                                nm = os.path.basename(im.name)
                                if g and not nm.startswith("GSM"): nm = g.group(1) + "_" + nm
                                members.append(_M(nm, inner.extractfile(im).read()))
                else:
                    members.append(_M(os.path.basename(m.name), data))
        tf = _TF()
        if True:
            log("  files after unpacking (first 30):", [m.name for m in members[:30]])
            by = {}
            for m in members:
                g = re.match(r"(GSM\d+)", os.path.basename(m.name))
                if g: by.setdefault(g.group(1), []).append(m)
            for gsm, ms in by.items():
                names = {os.path.basename(m.name).lower(): m for m in ms}
                a = None
                h5 = [m for n, m in names.items() if n.endswith(".h5")]
                mtx = [m for n, m in names.items() if "matrix.mtx" in n]
                if h5:
                    tmp = os.path.join(VDIR, os.path.basename(h5[0].name))
                    if not os.path.exists(tmp):
                        with open(tmp, "wb") as f: f.write(tf.extractfile(h5[0]).read())
                    a = sc.read_10x_h5(tmp)
                elif mtx:
                    d = os.path.join(VDIR, gsm); os.makedirs(d, exist_ok=True)
                    for n, m in names.items():
                        tgt = None
                        if "matrix.mtx" in n: tgt = "matrix.mtx.gz"
                        elif "barcodes" in n: tgt = "barcodes.tsv.gz"
                        elif "features" in n or "genes" in n: tgt = "features.tsv.gz"
                        if tgt:
                            b = tf.extractfile(m).read()
                            if not n.endswith(".gz"): b = gzip.compress(b)
                            with open(os.path.join(d, tgt), "wb") as f: f.write(b)
                    try:
                        a = sc.read_10x_mtx(d, var_names="gene_symbols")
                    except Exception:
                        os.replace(os.path.join(d, "features.tsv.gz"), os.path.join(d, "genes.tsv.gz"))
                        a = sc.read_10x_mtx(d, var_names="gene_symbols")
                else:
                    tabs = [m for n, m in names.items() if re.search(r"\.(txt|csv|tsv)(\.gz)?$", n)]
                    if tabs:
                        d = read_member_table(tf.extractfile(tabs[0]).read(), tabs[0].name)
                        if d.shape[0] > d.shape[1]: d = d.T   # want cells x genes
                        a = ad.AnnData(sparse.csr_matrix(d.values.astype(np.float32)), obs=pd.DataFrame(index=d.index.astype(str)), var=pd.DataFrame(index=d.columns.astype(str)))
                if a is not None:
                    a.var_names_make_unique(); a.obs["gsm"] = gsm; a.obs_names = gsm + "_" + a.obs_names.astype(str)
                    adatas.append(a); log(f"  {gsm}: {a.n_obs} cells x {a.n_vars} genes")
    if not adatas and others:
        for f in others:
            p = fetch(GEO + "suppl/" + f, os.path.join(VDIR, f))
            if f.endswith(".h5ad"):
                a = sc.read_h5ad(p); adatas.append(a); log("  read", f, a.shape)
    if not adatas:
        raise RuntimeError("could not parse GSE209781 supplementary files; see log for file names")
    A = ad.concat(adatas, join="outer", fill_value=0) if len(adatas) > 1 else adatas[0]
    if "gsm" in A.obs:
        A.obs = A.obs.merge(meta.set_index("gsm")[["condition", "title"]], left_on="gsm", right_index=True, how="left")
    A.obs["donor"] = A.obs.get("gsm", "d").astype(str)
    return A

MAJOR = {"PT": ["LRP2", "CUBN", "SLC34A1", "ALDOB"], "TAL": ["UMOD", "SLC12A1"], "DCT": ["SLC12A3", "CALB1"],
         "CD-PC": ["AQP2", "GATA3"], "CD-IC": ["ATP6V1G3", "SLC4A1", "SLC26A4"], "POD": ["NPHS1", "NPHS2", "PODXL"],
         "EC": ["PECAM1", "EMCN", "PLVAP"], "FIB/MES": ["COL1A2", "DCN", "PDGFRB", "LUM"], "Myeloid": ["CD68", "C1QA", "CD14", "LYZ", "CD163"],
         "T": ["CD3E", "CD3D", "IL7R", "CD2"], "NK": ["NKG7", "GNLY", "KLRD1"], "B": ["MS4A1", "CD79A", "CD19"],
         "Plasma": ["JCHAIN", "MZB1", "IGHG1"], "Mast": ["TPSAB1", "CPA3"]}
IMMUNE = {"Myeloid", "T", "NK", "B", "Plasma", "Mast"}
TLS = ["CCL19", "CCL21", "CXCL13", "LTB", "CCR7", "SELL", "LAMP3", "MS4A1", "CD79A", "CD3E", "CXCR5", "IL7R"]

def process(A):
    A.var["mt"] = A.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    n0 = A.n_obs
    A = A[(A.obs.n_genes_by_counts >= 200) & (A.obs.n_genes_by_counts <= 6000) & (A.obs.pct_counts_mt < 30)].copy()
    sc.pp.filter_genes(A, min_cells=10)
    log(f"  QC: {n0} -> {A.n_obs} cells")
    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    sc.pp.highly_variable_genes(A, n_top_genes=3000, batch_key="donor")
    sc.tl.pca(A, n_comps=40, use_highly_variable=True)
    try: sc.external.pp.harmony_integrate(A, "donor"); rep = "X_pca_harmony"
    except Exception as e: log("  harmony:", e); rep = "X_pca"
    sc.pp.neighbors(A, use_rep=rep, n_neighbors=20); sc.tl.umap(A); sc.tl.leiden(A, resolution=1.0, key_added="leiden")
    for k, gl in MAJOR.items():
        gl = [g for g in gl if g in A.var_names]
        if len(gl) >= 1: sc.tl.score_genes(A, gl, score_name="m_" + k)
    cols = [c for c in A.obs if c.startswith("m_")]
    cl = A.obs.groupby("leiden")[cols].mean()
    lab = cl.idxmax(axis=1).str[2:]
    A.obs["ct"] = A.obs.leiden.map(lab).astype("category")
    cl.assign(label=lab).to_csv(os.path.join(OUT, "S3_cluster_annotation_scores.csv"), encoding="utf-8-sig")
    A.obs["is_immune"] = A.obs.ct.isin(IMMUNE)
    log("  cell types:", A.obs.groupby(["ct", "condition"]).size().unstack(fill_value=0).to_dict())
    return A

def donor_level(A, key, groupby="ct", min_cells=10):
    d = A.obs[[groupby, "donor", "condition", key]]
    m = d.groupby([groupby, "donor", "condition"], observed=True)[key].agg(["mean", "size"]).reset_index()
    m = m[m["size"] >= min_cells]; rows = []
    for c, g in m.groupby(groupby, observed=True):
        x, y = g[g.condition == "DKD"]["mean"], g[g.condition == "CTRL"]["mean"]
        if len(x) >= 2 and len(y) >= 2:
            rows.append({groupby: c, "DKD": x.mean(), "CTRL": y.mean(), "diff": x.mean() - y.mean(),
                         "P_t": stats.ttest_ind(x, y, equal_var=False).pvalue, "nDKD": len(x), "nCTRL": len(y)})
    return pd.DataFrame(rows), m

def validate(A):
    res = {}
    # composition
    tab = pd.crosstab(A.obs.donor, A.obs.ct, normalize="index") * 100
    cond = A.obs.groupby("donor").condition.first()
    comp = []
    for c in tab.columns:
        x, y = tab.loc[cond == "DKD", c], tab.loc[cond == "CTRL", c]
        comp.append({"ct": c, "DKD_mean_pct": x.mean(), "CTRL_mean_pct": y.mean(), "log2_ratio": np.log2((x.mean() + .01) / (y.mean() + .01)),
                     "P_t": stats.ttest_ind(x, y, equal_var=False).pvalue if len(x) > 1 and len(y) > 1 else np.nan,
                     "P_MWU": stats.mannwhitneyu(x, y).pvalue if len(x) > 1 and len(y) > 1 else np.nan})
    comp = pd.DataFrame(comp); comp.to_csv(os.path.join(OUT, "S3_composition.csv"), index=False, encoding="utf-8-sig")
    tab.assign(condition=cond).to_csv(os.path.join(OUT, "S3_composition_per_donor.csv"), encoding="utf-8-sig")
    log("  composition:", comp.round(3).to_dict("records"))
    imm_pct = A.obs.groupby("donor").is_immune.mean() * 100
    lym_pct = A.obs.groupby("donor").ct.apply(lambda s: s.isin(["T", "B", "Plasma", "NK"]).mean() * 100)
    # signatures
    cand = pd.read_csv(os.path.join(OUT, "progression_genes_tiered.csv"), index_col=0)
    t1 = cand[cand.tier.str.startswith("Tier1")]
    sets = {"Prog_up_T1": list(t1.index[t1.ADVvsEARLY_logFC > 0]), "Prog_down_T1": list(t1.index[t1.ADVvsEARLY_logFC < 0]),
            "TLS_score": TLS, "NIK_RELB_targets": ["NFKB2", "RELB", "CCL19", "CCL21", "CXCL13", "MADCAM1", "VCAM1", "ICAM1", "TNFAIP3", "BIRC3"]}
    sig = []
    for k, gl in sets.items():
        gl = [g for g in gl if g in A.var_names]
        sc.tl.score_genes(A, gl, score_name=k)
        r, m = donor_level(A, k); r["signature"] = k; sig.append(r)
    sig = pd.concat(sig); sig.to_csv(os.path.join(OUT, "S3_signatures_by_celltype.csv"), index=False, encoding="utf-8-sig")
    # key genes per cell type, donor-level
    genes = ["LTB", "LTA", "LTBR", "MAP3K14", "RELB", "NFKB2", "CCL19", "CCL21", "CXCL13", "CCR7", "IL7R", "CD44", "THBS1", "COL4A1", "EGF", "PPARA"]
    gg = [g for g in genes if g in A.var_names]
    X = A[:, gg].X; X = X.toarray() if sparse.issparse(X) else X
    for j, g in enumerate(gg): A.obs["g_" + g] = X[:, j]
    rows = []
    for g in gg:
        r, _ = donor_level(A, "g_" + g); r["gene"] = g; rows.append(r)
    kg = pd.concat(rows); kg.to_csv(os.path.join(OUT, "S3_key_genes_by_celltype.csv"), index=False, encoding="utf-8-sig")
    # where are LTB and LTBR expressed
    mean = pd.DataFrame({g: A.obs.groupby("ct", observed=True)["g_" + g].mean() for g in gg})
    pct = pd.DataFrame({g: A.obs.groupby("ct", observed=True)["g_" + g].apply(lambda v: (v > 0).mean() * 100) for g in gg})
    mean.to_csv(os.path.join(OUT, "S3_key_genes_mean_by_ct.csv"), encoding="utf-8-sig"); pct.to_csv(os.path.join(OUT, "S3_key_genes_pct_by_ct.csv"), encoding="utf-8-sig")
    # figures
    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    sc.pl.umap(A, color="ct", ax=ax[0], show=False, legend_loc="on data", legend_fontsize=6, frameon=False, title="GSE209781 cell types")
    sc.pl.umap(A, color="condition", ax=ax[1], show=False, frameon=False, palette=COND_COL, size=2, title="Condition")
    savefig(fig, "Fig7A_validation_umap")
    sc.pl.dotplot(A, [g for g in ["LTB", "LTBR", "MAP3K14", "RELB", "NFKB2", "CCL19", "CCL21", "CXCL13", "CCR7", "IL7R", "CD44"] if g in A.var_names],
                  groupby="ct", show=False, standard_scale="var", figsize=(6, 4))
    savefig(plt.gcf(), "Fig7B_lymphotoxin_axis_dotplot")
    fig, axs = plt.subplots(1, 3, figsize=(8, 2.6))
    for ax, (nm, s) in zip(axs, [("Immune cells (% of all)", imm_pct), ("Lymphocytes (% of all)", lym_pct),
                                 ("TLS score (all cells, donor mean)", A.obs.groupby("donor").TLS_score.mean())]):
        for i, g in enumerate(["CTRL", "DKD"]):
            v = s[cond[s.index] == g]
            ax.boxplot(v, positions=[i], widths=0.5, showfliers=False); ax.scatter(np.full(len(v), i), v, c=COND_COL[g], s=14)
        x, y = s[cond[s.index] == "DKD"], s[cond[s.index] == "CTRL"]
        p = stats.ttest_ind(x, y, equal_var=False).pvalue if len(x) > 1 and len(y) > 1 else np.nan
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Control", "DKD"]); ax.set_title(f"{nm}\nP={p:.3g}", fontsize=7)
    savefig(fig, "Fig7C_validation_immune_TLS")
    log("  immune % per donor:", imm_pct.round(2).to_dict()); log("  lymphocyte %:", lym_pct.round(2).to_dict())
    return A

def tf_and_lr(A):
    out = {}
    try:
        import decoupler as dc
        net = dc.get_collectri(organism="human", split_complexes=False)
        X = A.layers["counts"]; X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
        key = A.obs.ct.astype(str) + "|" + A.obs.donor.astype(str); rows, idx = [], []
        for k, ii in A.obs.groupby(key).indices.items():
            if len(ii) >= 10: rows.append(np.asarray(X[ii].sum(0)).ravel()); idx.append(k)
        P = pd.DataFrame(np.vstack(rows), index=idx, columns=A.var_names); P = P.loc[:, (P > 0).mean() > 0.1]
        L = np.log1p(P.div(P.sum(1), axis=0) * 1e6)
        est = dc.run_ulm(mat=L, net=net, verbose=False); acts = est[0] if isinstance(est, tuple) else est
        meta = pd.DataFrame({"ct": [i.split("|")[0] for i in idx], "donor": [i.split("|")[1] for i in idx]}, index=idx)
        meta["condition"] = meta.donor.map(A.obs.groupby("donor").condition.first())
        rr = []
        for c in meta.ct.unique():
            mm = meta[meta.ct == c]; dd, cc = mm.index[mm.condition == "DKD"], mm.index[mm.condition == "CTRL"]
            if len(dd) < 2 or len(cc) < 2: continue
            t, p = stats.ttest_ind(acts.loc[dd], acts.loc[cc], equal_var=False)
            rr.append(pd.DataFrame({"ct": c, "TF": acts.columns, "t": t, "P": p}))
        R = pd.concat(rr); R.to_csv(os.path.join(OUT, "S3_TF_DKDvsCTRL.csv"), index=False, encoding="utf-8-sig")
        log("  TF RELB/NFKB2/NFKB1/RELA:", R[R.TF.isin(["RELB", "NFKB2", "NFKB1", "RELA", "IRF1"])].round(4).to_dict("records"))
        out["tf"] = R
    except Exception as e:
        log("!! TF failed:", e)
    try:
        import liana as li
        res = []
        for g in ("CTRL", "DKD"):
            sub = A[A.obs.condition == g].copy()
            li.mt.rank_aggregate(sub, groupby="ct", resource_name="consensus", expr_prop=0.1, min_cells=10, use_raw=False, n_perms=1000, verbose=False, key_added="liana")
            r = sub.uns["liana"]; r["condition"] = g; res.append(r)
        R = pd.concat(res); R.to_csv(os.path.join(OUT, "S3_liana_all.csv"), index=False, encoding="utf-8-sig")
        foc = R[R.ligand_complex.isin(["LTB", "LTA", "CCL19", "CCL21", "CXCL13"]) | R.receptor_complex.str.contains("CD44|LTBR|CCR7", regex=True)]
        foc.to_csv(os.path.join(OUT, "S3_liana_lymphotoxin_CD44.csv"), index=False, encoding="utf-8-sig")
        log("  L-R focus rows:", len(foc))
        out["lr"] = R
    except Exception as e:
        log("!! liana failed:", e)
    return out

# ===================================================================== 3B virtual knockout
def pcnet(X, n_comp=3):
    """PC-regression gene network (Osorio et al.). X: cells x genes (scaled). returns genes x genes"""
    n, p = X.shape; W = np.zeros((p, p), dtype=np.float32)
    for i in range(p):
        y = X[:, i]; Z = np.delete(X, i, axis=1)
        U, S, Vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
        k = min(n_comp, len(S)); T = U[:, :k] * S[:k]
        coef = np.linalg.lstsq(T, y - y.mean(), rcond=None)[0]
        b = Vt[:k].T @ coef
        W[i, np.arange(p) != i] = b
    return W

def propagate(W, g, delta, steps=3):
    """CellOracle-style signal propagation on a linear GRN: W[i, j] = effect of gene j on gene i"""
    d = np.zeros(W.shape[0]); d[g] = delta; tot = d.copy()
    for _ in range(steps):
        d = W @ d; d[g] = 0; tot += d
    tot[g] = delta
    return tot

def load_progression_logfc():
    a = pd.read_csv(os.path.join(OUT, "ALL_genes_integrated.csv"), index_col=0)
    return a["ADVvsEARLY_logFC"].dropna()

def virtual_ko(A, ct_mask, gene, label, n_genes=1000, n_nets=3, n_cells=600, n_null=1000, seed=1):
    """virtual knockout / overexpression by GRN signal propagation (PC-regression network)"""
    rng = np.random.default_rng(seed)
    sub = A[ct_mask].copy(); sub.var_names_make_unique()
    if "condition" in sub.obs and (sub.obs.condition == "DKD").sum() >= 300:
        sub = sub[sub.obs.condition == "DKD"].copy()
    if gene not in sub.var_names: log(f"  vKO {label}: {gene} not detected"); return None
    if sub.n_obs < 80: log(f"  vKO {label}: only {sub.n_obs} cells - skipped"); return None
    vn = pd.Index(sub.var_names)
    good = ~(vn.str.startswith(("ENSG", "LINC", "MT-", "RPL", "RPS", "MIR", "SNOR")) | vn.str.contains(r"-AS\d|\.\d|^LOC|-DT$|^OR\d", regex=True))
    Xc = sub.X; frac = np.asarray((Xc > 0).mean(0)).ravel()
    keep = (good & (frac >= 0.05)) | vn.isin([gene, "LTB", "LTBR", "RELB", "NFKB2", "MAP3K14", "CCL19", "CCL21", "CXCL13", "CD44", "IL7R"])
    sub = sub[:, keep].copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=min(n_genes, sub.n_vars - 1))
    genes = list(sub.var_names[sub.var.highly_variable])
    extra = [g for g in [gene, "LTB", "LTBR", "RELB", "NFKB2", "MAP3K14", "CCL19", "CCL21", "CXCL13", "CD44", "IL7R"] if g in sub.var_names and g not in genes]
    genes = list(dict.fromkeys(genes + extra))   # no duplicated genes (a duplicated target would leak the knocked-out signal)
    X = sub[:, genes].X; X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    if (X[:, genes.index(gene)] > 0).mean() < 0.02:
        log(f"  vKO {label}: {gene} expressed in <2% of cells - skipped"); return None
    nets = []
    for k in range(n_nets):
        idx = rng.choice(X.shape[0], size=min(n_cells, X.shape[0]), replace=False)
        Xs = X[idx]; sd = Xs.std(0); sd[sd == 0] = 1; Xs = (Xs - Xs.mean(0)) / sd
        nets.append(pcnet(Xs)); log(f"    {label}: network {k+1}/{n_nets} ({sub.n_obs} cells, {len(genes)} genes)")
    W = np.mean(nets, axis=0)
    rad = np.max(np.abs(np.linalg.eigvals(W)))
    if rad >= 0.9: W = W * (0.9 / rad)                     # keep propagation stable
    g = genes.index(gene)
    s_ko = propagate(W, g, -1.0)                            # knock-out: -1 SD of the gene
    s_oe = propagate(W, g, +2.0)                            # over-expression: +2 SD
    r = pd.DataFrame({"gene": genes, "shift_KO": s_ko, "shift_OE": s_oe})
    # null for each gene's shift: knock out random genes
    null = np.vstack([propagate(W, j, -1.0) for j in rng.choice([i for i in range(len(genes)) if i != g], size=min(n_null, len(genes) - 1), replace=False)])
    z = (s_ko - null.mean(0)) / (null.std(0) + 1e-12)
    r["z_KO"] = z; r["P"] = 2 * stats.norm.sf(np.abs(z)); r.loc[r.gene == gene, "P"] = np.nan; r["FDR"] = bh(r.P.values)
    # does the perturbation reverse the bulk progression signature?
    fc = load_progression_logfc().reindex(genes)
    ok = fc.notna().values & (np.arange(len(genes)) != g)
    rev = stats.spearmanr(s_ko[ok], fc.values[ok]).correlation
    null_rev = np.array([stats.spearmanr(n_[ok], fc.values[ok]).correlation for n_ in null])
    p_rev = (np.sum(null_rev <= rev) + 1) / (len(null_rev) + 1)
    summ = {"label": label, "gene": gene, "cells": sub.n_obs, "genes_in_network": len(genes),
            "reversal_rho_KO": rev, "empirical_P_KO_reverses": p_rev, "reversal_rho_OE": -rev,
            "n_perturbed_FDR05": int((r.FDR < 0.05).sum())}
    r = r.sort_values("P"); r.to_csv(os.path.join(OUT, f"S3_vKO_{label}.csv"), index=False, encoding="utf-8-sig")
    log(f"  vKO {label}: rho(KO shift, progression logFC)={rev:.3f}, empirical P={p_rev:.3g}; "
        f"perturbed FDR<0.05: {summ['n_perturbed_FDR05']}; top down: {r[r.shift_KO < 0].head(10).gene.tolist()}")
    return r, summ, null_rev

def enrich(genes, background, label):
    try:
        body = json.dumps({"organism": "hsapiens", "query": list(genes), "sources": ["GO:BP", "KEGG", "REAC"], "user_threshold": 0.05,
                           "domain_scope": "custom", "background": list(background)}).encode()
        req = urllib.request.Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/", data=body, headers={"Content-Type": "application/json", **UA})
        js = json.loads(urllib.request.urlopen(req, timeout=300).read().decode())
        df = pd.DataFrame([{"source": x["source"], "term": x["name"], "P_adj": x["p_value"], "overlap": x["intersection_size"]} for x in js.get("result", [])])
        df.to_csv(os.path.join(OUT, f"S3_vKO_{label}_enrich.csv"), index=False, encoding="utf-8-sig")
        log(f"    enrichment {label}:", df.head(8).term.tolist() if len(df) else "none")
        return df
    except Exception as e:
        log("    enrichment failed:", e)

def plot_vko(results):
    rs = {k: v for k, v in results.items() if v is not None}
    if not rs: return
    fig, axs = plt.subplots(2, len(rs), figsize=(3.1 * len(rs), 5.6)); axs = np.array(axs).reshape(2, -1)
    for j, (k, (r, summ, null_rev)) in enumerate(rs.items()):
        ax = axs[0, j]; t = r.dropna(subset=["P"]).sort_values("shift_KO")
        top = pd.concat([t.head(8), t.tail(8)])
        ax.barh(range(len(top)), top.shift_KO, color=np.where(top.shift_KO < 0, "#4C72B0", "#C44E52"))
        ax.set_yticks(range(len(top))); ax.set_yticklabels(top.gene, fontsize=6, style="italic"); ax.axvline(0, c="k", lw=0.5)
        ax.set_title(f"virtual KO of {summ['gene']}\n({k.split('_',1)[1]})", fontsize=7); ax.set_xlabel("predicted shift")
        ax = axs[1, j]; ax.hist(null_rev, bins=30, color="#BBBBBB"); ax.axvline(summ["reversal_rho_KO"], c="#C44E52")
        ax.set_title(f"reversal of progression signature\nrho={summ['reversal_rho_KO']:.2f}, P={summ['empirical_P_KO_reverses']:.3g}", fontsize=7)
        ax.set_xlabel("Spearman rho (KO shift vs progression logFC)")
    fig.tight_layout(); savefig(fig, "Fig8_virtual_KO")

def run_vko(A, B):
    res = {}
    if A is not None:
        lymph = A.obs.ct.isin(["T", "B", "Plasma", "NK"]).values; imm = A.obs.is_immune.values
        dkdA = (A.obs.condition == "DKD").values
        res["LTB_lymphocytes"] = virtual_ko(A, lymph, "LTB", "LTB_lymphocytes")
        res["RELB_immune"] = virtual_ko(A, imm, "RELB", "RELB_immune")
        res["CD44_immune"] = virtual_ko(A, imm, "CD44", "CD44_immune")
        st = A.obs.ct.isin(["FIB/MES"]).values
        if st.sum() >= 80: res["LTBR_stroma"] = virtual_ko(A, st, "LTBR", "LTBR_stroma")
    if B is not None:
        res["LTBR_PT"] = virtual_ko(B, (B.obs.ct == "PT").values, "LTBR", "LTBR_PT_discovery")
        res["MAP3K14_PT"] = virtual_ko(B, (B.obs.ct == "PT").values, "MAP3K14", "MAP3K14_PT_discovery")
        if res.get("LTBR_stroma") is None:
            res["LTBR_stroma"] = virtual_ko(B, B.obs.ct.isin(["FIB", "MES"]).values, "LTBR", "LTBR_stroma_discovery")
    summ = pd.DataFrame([v[1] for v in res.values() if v is not None])
    summ.to_csv(os.path.join(OUT, "S3_vKO_summary.csv"), index=False, encoding="utf-8-sig"); log("  vKO summary:", summ.round(4).to_dict("records"))
    for k, v in res.items():
        if v is None: continue
        r = v[0]; sig = r[(r.FDR < 0.05)]
        if len(sig) >= 5:
            enrich(sig[sig.shift_KO < 0].gene, r.gene, k + "_down"); enrich(sig[sig.shift_KO > 0].gene, r.gene, k + "_up")
    plot_vko(res)

def main():
    log("=" * 70); log("validation analysis start", time.strftime("%Y-%m-%d %H:%M:%S"))
    A = None
    pa = os.path.join(RAW, "processed_GSE209781.h5ad")
    if "--vko-only" in sys.argv and os.path.exists(pa):
        A = sc.read_h5ad(pa); log("  vKO-only mode: loaded processed validation set")
        B = sc.read_h5ad(os.path.join(RAW, "processed_discovery_snRNA.h5ad")) if os.path.exists(os.path.join(RAW, "processed_discovery_snRNA.h5ad")) else None
        run_vko(A, B); log("virtual knockout done", time.strftime("%Y-%m-%d %H:%M:%S")); return
    try:
        if os.path.exists(pa):
            A = sc.read_h5ad(pa); log("  loaded processed validation set")
        else:
            meta = sample_meta(); A = build_counts(meta)
            A = A[A.obs.condition.isin(["DKD", "CTRL"])].copy()
            A = process(A)
        A = validate(A); tf_and_lr(A)
        A.obs = A.obs.astype({c: str for c in A.obs.columns if A.obs[c].dtype == object})
        A.write_h5ad(pa, compression="gzip")
    except Exception as e:
        log("!! validation failed:", e); log(traceback.format_exc())
    B = None
    try:
        B = sc.read_h5ad(os.path.join(RAW, "processed_discovery_snRNA.h5ad"))
    except Exception as e:
        log("  discovery h5ad not loaded:", e)
    try:
        run_vko(A, B)
    except Exception as e:
        log("!! virtual KO failed:", e); log(traceback.format_exc())
    log("validation analysis done", time.strftime("%Y-%m-%d %H:%M:%S"))

if __name__ == "__main__":
    main()
