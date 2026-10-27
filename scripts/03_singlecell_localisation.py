# -*- coding: utf-8 -*-
"""
Single-cell localisation of DKD progression candidate genes
  - finds human kidney single-cell/nucleus datasets on CELLxGENE that include diabetic samples
  - downloads the h5ad (resumable) and reads it with h5py (no scanpy/anndata needed)
  - per cell type: mean expression and % expressing; donor-level pseudobulk DKD vs control
  - dot plot + tables + HTML section
Input : results/progression_genes_tiered.csv (from the tissue analysis)
Output: results/SC_*.csv, figures/Fig4_*.png/pdf, 单细胞定位报告.html
"""
import os, sys, re, json, ssl, time, traceback, urllib.request
try:
    import certifi
    ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())
except Exception:
    pass
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repository root; scripts/ is one level down
RAW = os.path.join(BASE, "raw"); OUT = os.path.join(BASE, "results"); FIG = os.path.join(BASE, "figures")
for d in (RAW, OUT, FIG):
    os.makedirs(d, exist_ok=True)
LOGF = open(os.path.join(BASE, "log_03_localisation.txt"), "a", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOGF.write(s + "\n"); LOGF.flush()

UA = {"User-Agent": "Mozilla/5.0 (research script)"}
API = "https://api.cellxgene.cziscience.com/curation/v1"

def get_json(url, tries=5):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            log(f"  retry {k+1}: {e}"); time.sleep(5)
    raise RuntimeError("failed: " + url)

def fetch(url, dest, tries=30):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    part = dest + ".part"
    for k in range(tries):
        try:
            have = os.path.getsize(part) if os.path.exists(part) else 0
            hdr = dict(UA)
            if have: hdr["Range"] = f"bytes={have}-"
            log(f"  downloading {url}" + (f" (resume {have/1e9:.2f} GB)" if have else ""))
            with urllib.request.urlopen(urllib.request.Request(url, headers=hdr), timeout=300) as r:
                if have and r.status != 206: have = 0
                total = have + int(r.headers.get("Content-Length", 0) or 0)
                done, t0 = have, time.time()
                with open(part, "ab" if have else "wb") as f:
                    while True:
                        b = r.read(4 << 20)
                        if not b: break
                        f.write(b); done += len(b)
                        if time.time() - t0 > 30:
                            log(f"    {done/1e9:.2f} / {total/1e9:.2f} GB"); t0 = time.time()
            os.replace(part, dest); return dest
        except urllib.error.HTTPError as e:
            if e.code == 416 and os.path.exists(part):
                os.replace(part, dest); return dest
            log(f"  retry {k+1}: {e}"); time.sleep(10)
        except Exception as e:
            log(f"  retry {k+1}: {e}"); time.sleep(10)
    raise RuntimeError("download failed: " + url)

# ------------------------------------------------------------------ dataset discovery
def names(lst):
    return [x.get("label", "") if isinstance(x, dict) else str(x) for x in (lst or [])]

def discover():
    cache = os.path.join(OUT, "SC_candidate_datasets.csv")
    cols = get_json(f"{API}/collections")
    log(f"CELLxGENE collections: {len(cols)}")
    rows = []
    for c in cols:
        text = (c.get("name", "") + " " + (c.get("description") or "")).lower()
        if "kidney" not in text and "renal" not in text:
            continue
        try:
            full = get_json(f"{API}/collections/{c['collection_id']}")
        except Exception as e:
            log(f"  skip collection {c.get('name')}: {e}"); continue
        for ds in full.get("datasets", []):
            dis = names(ds.get("disease")); tis = names(ds.get("tissue")); org = names(ds.get("organism"))
            assay = names(ds.get("assay"))
            if "Homo sapiens" not in org: continue
            if not any("kidney" in t.lower() or "renal" in t.lower() or "nephron" in t.lower() for t in tis): continue
            if not any("diabet" in d.lower() for d in dis): continue
            if any("ATAC" in a or "Visium" in a or "Slide-seq" in a or "spatial" in a.lower() for a in assay): continue
            h5 = [a for a in ds.get("assets", []) if a.get("filetype") == "H5AD"]
            if not h5: continue
            rows.append({"collection": c.get("name"), "dataset_id": ds.get("dataset_id"), "title": ds.get("title"),
                         "cells": ds.get("cell_count"), "disease": "; ".join(dis), "assay": "; ".join(assay),
                         "tissue": "; ".join(tis)[:200], "url": h5[0]["url"], "GB": round(h5[0].get("filesize", 0) / 1e9, 2)})
    df = pd.DataFrame(rows).sort_values("cells", ascending=False)
    df.to_csv(cache, index=False, encoding="utf-8-sig")
    for _, r in df.iterrows():
        log(f"  * {r.title} | {r.cells} cells | {r.GB} GB | {r.assay} | {r.disease[:120]}")
    return df

def choose(df, max_gb=15):
    d = df[(df.GB <= max_gb) & (df.GB > 0)].copy()
    if not len(d): d = df.copy()
    # prefer integrated KPMP / atlas single-nucleus data with the most cells
    d["pref"] = d.title.str.contains("nucle|snRNA|atlas|KPMP|all", case=False).astype(int)
    d = d.sort_values(["pref", "cells"], ascending=False)
    return d.iloc[0]

# ------------------------------------------------------------------ h5ad reading via h5py
def _dec(a):
    a = np.asarray(a)
    if a.dtype.kind in ("S", "O"):
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a], dtype=object)
    return a

def read_col(f, grp, key):
    g = f[grp]
    if key not in g:
        return None
    node = g[key]
    if hasattr(node, "keys") and "codes" in node:          # anndata >= 0.8 categorical
        cats = _dec(node["categories"][()]); codes = node["codes"][()]
        return pd.Categorical.from_codes(np.where(codes < 0, -1, codes), categories=pd.Index(cats).astype(str)).astype(object)
    vals = node[()]
    cat_grp = f"{grp}/__categories"
    if cat_grp in f and key in f[cat_grp]:                 # anndata 0.7 categorical
        cats = _dec(f[cat_grp][key][()])
        return np.array([cats[c] if c >= 0 else None for c in vals], dtype=object)
    return _dec(vals)

def var_names(f, grp):
    for k in ("feature_name", "gene_symbols", "gene_symbol", "features"):
        v = read_col(f, grp, k)
        if v is not None:
            return np.asarray(v, dtype=object)
    idx = f[grp].attrs.get("_index", "_index")
    return _dec(f[grp][idx][()])

def matrix_info(f, path):
    g = f[path]
    enc = g.attrs.get("encoding-type", g.attrs.get("h5sparse_format", ""))
    if isinstance(enc, bytes): enc = enc.decode()
    shape = tuple(g.attrs.get("shape", g.attrs.get("h5sparse_shape", (0, 0))))
    return ("csc" if "csc" in str(enc) else "csr"), shape

def extract(f, xpath, var_grp, genes, chunk=20000):
    """returns dense (cells x genes) normalized log1p matrix for the target genes"""
    fmt, shape = matrix_info(f, xpath)
    vn = var_names(f, var_grp)
    pos = {g: i for i, g in enumerate(vn)}
    gi = np.array([pos[g] for g in genes if g in pos]); gsel = [g for g in genes if g in pos]
    n = shape[0]
    g = f[xpath]
    data, ind, ptr = g["data"], g["indices"], g["indptr"][()]
    out = np.zeros((n, len(gi)), dtype=np.float32)
    lib = np.zeros(n, dtype=np.float64)
    is_int = None
    lut = np.full(len(vn), -1, dtype=np.int64); lut[gi] = np.arange(len(gi))
    if fmt == "csr":
        for s in range(0, n, chunk):
            e = min(n, s + chunk)
            a, b = ptr[s], ptr[e]
            dd = data[a:b]; ii = ind[a:b]
            rows = np.repeat(np.arange(s, e), np.diff(ptr[s:e + 1]))
            if is_int is None and len(dd):
                is_int = np.allclose(dd[:5000], np.round(dd[:5000]))
            lib[s:e] = np.bincount(rows - s, weights=dd, minlength=e - s)
            cc = lut[ii]
            m = cc >= 0
            if m.any():
                out[rows[m], cc[m]] = dd[m]
            if (s // chunk) % 10 == 0:
                log(f"    read {e}/{n} cells")
    else:  # csc: columns are genes
        for j, c in enumerate(gi):
            a, b = ptr[c], ptr[c + 1]
            out[ind[a:b], j] = data[a:b]
        dd = data[:min(len(data), 5000)]
        is_int = np.allclose(dd, np.round(dd))
        if is_int:
            log("    computing library sizes (csc)")
            for s in range(0, len(ptr) - 1, 2000):
                a, b = ptr[s], ptr[min(len(ptr) - 1, s + 2000)]
                lib += np.bincount(ind[a:b], weights=data[a:b], minlength=n)
    if is_int:
        log("    raw counts -> log1p(CP10k)")
        out = np.log1p(out / np.maximum(lib, 1)[:, None] * 1e4).astype(np.float32)
    else:
        log("    matrix already normalized; used as is")
    return pd.DataFrame(out, columns=gsel)

# ------------------------------------------------------------------ analysis
def pick_label(f):
    for k in ("cell_type", "subclass.l1", "celltype", "cell_type_original"):
        v = read_col(f, "obs", k)
        if v is not None:
            return k, np.asarray(v, dtype=object)
    raise RuntimeError("no cell type column")

def main():
    log("=" * 70); log("single-cell start", time.strftime("%Y-%m-%d %H:%M:%S"))
    import h5py
    cand = pd.read_csv(os.path.join(OUT, "progression_genes_tiered.csv"), index_col=0)
    genes = cand.index.tolist()
    extra = ["EGF", "TNFRSF1A", "TNFRSF1B", "HAVCR1", "LCN2", "CCL2", "PLAUR", "TIMP1", "CD163", "C3", "C1QA", "TREM2", "TYROBP",
             # cell-type markers for sanity checking
             "NPHS1", "PODXL", "LRP2", "CUBN", "UMOD", "SLC12A1", "SLC12A3", "AQP2", "PECAM1", "EMCN", "PDGFRB", "ACTA2",
             "COL1A1", "PTPRC", "CD68", "CD3E", "MS4A1", "VCAM1", "HAVCR1"]
    genes = list(dict.fromkeys(genes + extra))
    df = discover()
    if not len(df):
        log("!! no diabetic kidney single-cell dataset found on CELLxGENE"); return
    ds = choose(df)
    log(f"CHOSEN: {ds.title} ({ds.cells} cells, {ds.GB} GB)")
    path = fetch(ds.url, os.path.join(RAW, f"sc_{ds.dataset_id}.h5ad"))
    with h5py.File(path, "r") as f:
        lab_key, ct = pick_label(f)
        dis = np.asarray(read_col(f, "obs", "disease"), dtype=object)
        donor = read_col(f, "obs", "donor_id")
        donor = np.asarray(donor if donor is not None else ["d"] * len(ct), dtype=object)
        log(f"  cell-type column: {lab_key}; diseases: {pd.Series(dis).value_counts().to_dict()}")
        if "raw" in f and "X" in f["raw"]:
            X = extract(f, "raw/X", "raw/var", genes)
        else:
            X = extract(f, "X", "var", genes)
    grp = np.where(pd.Series(dis).str.contains("diabet", case=False), "DKD",
                   np.where(pd.Series(dis).str.lower().eq("normal"), "CTRL", "other"))
    meta = pd.DataFrame({"cell_type": ct, "group": grp, "donor": donor})
    log(f"  groups: {meta.group.value_counts().to_dict()}; donors DKD={meta[meta.group=='DKD'].donor.nunique()}, "
        f"CTRL={meta[meta.group=='CTRL'].donor.nunique()}")
    keep = meta.group.isin(["DKD", "CTRL"]).values
    meta, X = meta[keep].reset_index(drop=True), X[keep].reset_index(drop=True)
    ctn = meta.cell_type.value_counts()
    big = ctn[ctn >= 100].index
    meta.loc[~meta.cell_type.isin(big), "cell_type"] = "other (<100 cells)"
    pd.DataFrame({"cells": meta.groupby(["cell_type", "group"]).size()}).to_csv(os.path.join(OUT, "SC_celltype_counts.csv"), encoding="utf-8-sig")
    # per cell type expression (all cells DKD+CTRL)
    mean = X.groupby(meta.cell_type.values).mean()
    pct = (X > 0).groupby(meta.cell_type.values).mean() * 100
    mean.T.to_csv(os.path.join(OUT, "SC_mean_expression_by_celltype.csv"), encoding="utf-8-sig")
    pct.T.to_csv(os.path.join(OUT, "SC_pct_expressing_by_celltype.csv"), encoding="utf-8-sig")
    # donor pseudobulk per cell type: DKD vs CTRL
    pb = X.groupby([meta.cell_type.values, meta.group.values, meta.donor.values]).mean()
    pbn = meta.groupby(["cell_type", "group", "donor"]).size()
    pb = pb[pbn.reindex(pb.index).values >= 10]
    rows = []
    for c in mean.index:
        if c not in pb.index.get_level_values(0): continue
        sub = pb.loc[c]
        if "DKD" not in sub.index.get_level_values(0) or "CTRL" not in sub.index.get_level_values(0): continue
        a, b = sub.loc["DKD"], sub.loc["CTRL"]
        if len(a) < 2 or len(b) < 2: continue
        for gname in X.columns:
            t, p = stats.ttest_ind(a[gname], b[gname], equal_var=False) if a[gname].std() + b[gname].std() > 0 else (np.nan, np.nan)
            rows.append({"gene": gname, "cell_type": c, "DKD_mean": a[gname].mean(), "CTRL_mean": b[gname].mean(),
                         "diff": a[gname].mean() - b[gname].mean(), "P": p, "n_DKD": len(a), "n_CTRL": len(b)})
    de = pd.DataFrame(rows)
    if len(de):
        de["FDR"] = np.nan
        ok = de.P.notna()
        from numpy import argsort
        p = de.loc[ok, "P"].values; o = np.argsort(p); r = np.empty(len(p)); r[o] = np.arange(1, len(p) + 1)
        q = np.minimum.accumulate((p * len(p) / r)[o][::-1])[::-1]; qq = np.empty(len(p)); qq[o] = np.minimum(q, 1)
        de.loc[ok, "FDR"] = qq
    de.to_csv(os.path.join(OUT, "SC_pseudobulk_DKDvsCTRL.csv"), index=False, encoding="utf-8-sig")
    # summary per gene
    summ = []
    for gname in X.columns:
        m = mean[gname]; top = m.idxmax(); srt = m.sort_values(ascending=False)
        spec = (srt.iloc[0] - srt.iloc[1]) / (srt.iloc[0] + 1e-9) if len(srt) > 1 else np.nan
        d = de[(de.gene == gname) & (de.cell_type == top)] if len(de) else pd.DataFrame()
        summ.append({"gene": gname, "top_cell_type": top, "top_mean": srt.iloc[0], "pct_in_top": pct.at[top, gname],
                     "second_cell_type": srt.index[1] if len(srt) > 1 else "", "specificity": spec,
                     "DKDvsCTRL_diff_in_top": d["diff"].iloc[0] if len(d) else np.nan,
                     "DKDvsCTRL_P_in_top": d["P"].iloc[0] if len(d) else np.nan,
                     "tier": cand["tier"].get(gname, "marker/extra") if "tier" in cand else "",
                     "blood_secreted": cand["blood_secreted"].get(gname, "") if "blood_secreted" in cand else ""})
    S = pd.DataFrame(summ).set_index("gene")
    S.to_csv(os.path.join(OUT, "SC_gene_localisation_summary.csv"), encoding="utf-8-sig")
    log("  top cell types of tier-1 genes: " + str(S[S.tier.astype(str).str.startswith("Tier1")].top_cell_type.value_counts().head(10).to_dict()))
    dotplot(mean, pct, S, cand)
    report(ds, meta, S)
    log("single-cell done", time.strftime("%Y-%m-%d %H:%M:%S"))

def dotplot(mean, pct, S, cand):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "pdf.fonttype": 42})
    t1 = [g for g in S.index if str(S.at[g, "tier"]).startswith("Tier1")]
    sec = [g for g in t1 if S.at[g, "blood_secreted"] in (True, "True")]
    focus = list(dict.fromkeys(sec + ["C3", "C1QA", "C1QB", "CD163", "MSR1", "MRC1", "TREM2", "TYROBP", "THBS2", "VCAN",
                                      "COL1A1", "COL3A1", "TGFBI", "CCL19", "CCL21", "TNFRSF1B", "TIMP1", "CCL2", "EGF", "HAVCR1", "VCAM1"]))
    focus = [g for g in focus if g in mean.columns]
    for name, gl in [("Fig4_dotplot_key_genes", focus), ("Fig4b_dotplot_tier1_all", [g for g in t1 if g in mean.columns])]:
        if not gl: continue
        M = mean[gl]; P = pct[gl]
        order_ct = M.max(axis=1).sort_values(ascending=False).index  # cell types
        gorder = sorted(gl, key=lambda g: (list(order_ct).index(M[g].idxmax()), -M[g].max()))
        M, P = M.loc[order_ct, gorder], P.loc[order_ct, gorder]
        Z = (M - M.min()) / (M.max() - M.min() + 1e-9)
        fig, ax = plt.subplots(figsize=(max(6, 0.22 * len(gorder) + 3), max(3.5, 0.25 * len(order_ct) + 1.5)))
        yy, xx = np.meshgrid(range(len(order_ct)), range(len(gorder)), indexing="ij")
        sc = ax.scatter(xx.ravel(), yy.ravel(), s=(P.values.ravel() / 100) * 120 + 0.5, c=Z.values.ravel(), cmap="Reds", vmin=0, vmax=1, lw=0)
        ax.set_xticks(range(len(gorder))); ax.set_xticklabels(gorder, rotation=90, fontsize=6.5, style="italic")
        ax.set_yticks(range(len(order_ct))); ax.set_yticklabels(order_ct, fontsize=7)
        ax.set_xlim(-0.7, len(gorder) - 0.3); ax.set_ylim(len(order_ct) - 0.5, -0.5)
        cb = fig.colorbar(sc, ax=ax, shrink=0.4, pad=0.01); cb.set_label("scaled mean expression")
        for v in (25, 50, 100):
            ax.scatter([], [], s=v / 100 * 120, c="grey", label=f"{v}%")
        ax.legend(title="% expressing", bbox_to_anchor=(1.12, 0.0), loc="lower left", frameon=False, fontsize=6)
        fig.tight_layout(); fig.savefig(os.path.join(FIG, name + ".png"), dpi=300); fig.savefig(os.path.join(FIG, name + ".pdf")); plt.close(fig)

def report(ds, meta, S):
    t = S.copy()
    t = t[t.tier.astype(str).str.startswith("Tier") | t.index.isin(["EGF", "TNFRSF1B", "C3", "CD163", "TIMP1", "CCL2"])]
    t = t.sort_values(["tier", "top_cell_type"])
    figs = [f for f in sorted(os.listdir(FIG)) if f.startswith("Fig4") and f.endswith(".png")]
    h = f"""<!doctype html><html><head><meta charset="utf-8"><title>单细胞定位</title>
<style>body{{font-family:"Microsoft YaHei",Arial;max-width:1100px;margin:24px auto;padding:0 16px}}table{{border-collapse:collapse;font-size:12px}}
td,th{{border-bottom:1px solid #eee;padding:3px 6px}}img{{max-width:100%}}</style></head><body>
<h1>DKD进展候选基因的单细胞定位</h1>
<p>数据集：{ds.title}（{ds.collection}；{ds.cells} 个细胞）。比较组：DKD vs 正常对照（disease=normal）。
DKD 供者 {meta[meta.group=='DKD'].donor.nunique()} 例，对照 {meta[meta.group=='CTRL'].donor.nunique()} 例。
差异检验为供者水平 pseudobulk Welch t 检验（每供者每细胞类型≥10个细胞）。</p>
{''.join(f'<p>{f}</p><img src="figures/{f}">' for f in figs)}
<h2>每个基因的主要表达细胞类型</h2>{t.to_html(float_format=lambda x: f"{x:.3g}")}
</body></html>"""
    open(os.path.join(BASE, "单细胞定位报告.html"), "w", encoding="utf-8").write(h)

if __name__ == "__main__":
    try:
        main()
    except Exception as ex:
        log("!! FATAL", ex); log(traceback.format_exc())
