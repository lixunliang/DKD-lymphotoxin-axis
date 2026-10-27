# -*- coding: utf-8 -*-
"""
DKD kidney-tissue public transcriptome analysis
  1) DKD vs control, glomeruli and tubulointerstitium (GEO microarrays) -> per-dataset moderated t + random-effects meta
  2) Progression: advanced vs early DKD (GSE142025 RNA-seq, whole biopsy)
  3) Integration -> progression-associated genes; HPA secretome annotation (blood-detectable candidates)
  4) Enrichment (g:Profiler), figures, HTML report
Run with:  python scripts/01_bulk_meta_analysis.py   (outputs are written under the repository root)
"""
import os, re, io, gzip, sys, json, time, tarfile, zipfile, urllib.request, urllib.error, traceback
import numpy as np
import pandas as pd
from scipy import stats, special

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repository root; scripts/ is one level down
RAW = os.path.join(BASE, "raw")
OUT = os.path.join(BASE, "results")
FIG = os.path.join(BASE, "figures")
for d in (RAW, OUT, FIG):
    os.makedirs(d, exist_ok=True)

LOGF = open(os.path.join(BASE, "run_log.txt"), "a", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOGF.write(s + "\n"); LOGF.flush()

# ---------------------------------------------------------------- download helpers
UA = {"User-Agent": "Mozilla/5.0 (research script)"}
def fetch(url, dest, tries=12):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    part = dest + ".part"
    for k in range(tries):
        try:
            have = os.path.getsize(part) if os.path.exists(part) else 0
            hdr = dict(UA)
            if have:
                hdr["Range"] = f"bytes={have}-"
            log(f"  downloading {url}" + (f" (resume from {have/1e6:.0f} MB)" if have else ""))
            req = urllib.request.Request(url, headers=hdr)
            with urllib.request.urlopen(req, timeout=300) as r:
                if have and r.status != 206:
                    have = 0  # server ignored Range -> restart
                with open(part, "ab" if have else "wb") as f:
                    while True:
                        b = r.read(1 << 20)
                        if not b:
                            break
                        f.write(b)
            os.replace(part, dest)
            return dest
        except urllib.error.HTTPError as e:
            if e.code == 416 and os.path.exists(part):  # already complete
                os.replace(part, dest)
                return dest
            log(f"  retry {k+1}/{tries}: {e}")
            if e.code in (404, 403, 416):
                break
            time.sleep(5)
        except Exception as e:
            log(f"  retry {k+1}/{tries}: {e}")
            time.sleep(5)
    raise RuntimeError("download failed: " + url)

def fetch_text(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            log(f"  retry {k+1}/{tries}: {e}")
            time.sleep(5 * (k + 1))
    raise RuntimeError("fetch failed: " + url)

def geo_dir(acc, kind):
    stub = acc[:-3] + "nnn"
    t = "series" if acc.startswith("GSE") else "platforms"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/{t}/{stub}/{acc}/{kind}/"

def list_dir(url):
    html = fetch_text(url)
    return sorted(set(re.findall(r'href="([^"?/][^"]*)"', html)))

# ---------------------------------------------------------------- parsing
def parse_series_matrix(path):
    meta, chars = {}, []
    rows, header, intable = [], None, False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                intable = True; continue
            if line.startswith("!series_matrix_table_end"):
                break
            if intable:
                parts = [p.strip('"') for p in line.split("\t")]
                if header is None:
                    header = parts
                else:
                    rows.append(parts)
                continue
            if line.startswith("!Sample_"):
                parts = [p.strip('"') for p in line.split("\t")]
                key = parts[0][8:]
                if key == "characteristics_ch1":
                    chars.append(parts[1:])
                else:
                    meta.setdefault(key, parts[1:])
    gsm = meta.get("geo_accession", [])
    info = pd.DataFrame({"gsm": gsm,
                         "title": meta.get("title", [""] * len(gsm)),
                         "source": meta.get("source_name_ch1", [""] * len(gsm))})
    for i, c in enumerate(chars):
        info[f"char{i+1}"] = c
    info["platform"] = meta.get("platform_id", [""] * len(gsm))
    expr = None
    if header and rows:
        df = pd.DataFrame(rows, columns=header).set_index(header[0])
        df = df.apply(pd.to_numeric, errors="coerce")
        df = df.dropna(how="all")
        if df.shape[0] > 100:
            expr = df
    return info, expr

def sample_text(info):
    cols = [c for c in info.columns if c not in ("gsm", "platform")]
    t = info[cols].astype(str).agg(" | ".join, axis=1).str.lower()
    return t.apply(lambda s: re.sub(r"[^a-z0-9]+", " ", s))

# ---------------------------------------------------------------- gene annotation
_GENEINFO = None
def gene_info():
    global _GENEINFO
    if _GENEINFO is None:
        p = fetch("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
                  os.path.join(RAW, "Homo_sapiens.gene_info.gz"))
        gi = pd.read_csv(p, sep="\t", usecols=["GeneID", "Symbol", "dbXrefs", "type_of_gene"], dtype=str)
        gi["ens"] = gi["dbXrefs"].str.extract(r"Ensembl:(ENSG\d+)")
        _GENEINFO = gi
    return _GENEINFO

def platform_map(gpl, ids):
    """probe id -> gene symbol (single-gene probes only)"""
    ids = pd.Index(ids).astype(str)
    # Brainarray ENTREZG custom CDF: "1234_at"
    if (ids.str.match(r"^\d+_at$")).mean() > 0.8:
        gi = gene_info()
        e2s = dict(zip(gi["GeneID"], gi["Symbol"]))
        return pd.Series([e2s.get(i.split("_")[0]) for i in ids], index=ids)
    # try GEO annot file
    try:
        p = fetch(geo_dir(gpl, "annot") + f"{gpl}.annot.gz", os.path.join(RAW, f"{gpl}.annot.gz"))
        with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
            lines = f.read().split("\n")
        b = [i for i, l in enumerate(lines) if l.startswith("!platform_table_begin")][0]
        e = [i for i, l in enumerate(lines) if l.startswith("!platform_table_end")][0]
        tab = pd.read_csv(io.StringIO("\n".join(lines[b + 1:e])), sep="\t", dtype=str)
        sym = tab.set_index("ID")["Gene symbol"]
    except Exception as ex:
        log(f"  annot not available for {gpl} ({ex}); using family SOFT")
        p = fetch(geo_dir(gpl, "soft") + f"{gpl}_family.soft.gz", os.path.join(RAW, f"{gpl}_family.soft.gz"))
        buf, intab = [], False
        with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
            for l in f:
                if l.startswith("!platform_table_begin"):
                    intab = True; continue
                if l.startswith("!platform_table_end"):
                    break
                if intab:
                    buf.append(l)
        tab = pd.read_csv(io.StringIO("".join(buf)), sep="\t", dtype=str, low_memory=False)
        cols = {c.lower(): c for c in tab.columns}
        if "gene symbol" in cols:
            sym = tab.set_index("ID")[cols["gene symbol"]]
        elif "gene_assignment" in cols:
            sym = tab.set_index("ID")[cols["gene_assignment"]].str.split(" // ").str[1]
        elif "symbol" in cols:
            sym = tab.set_index("ID")[cols["symbol"]]
        else:
            raise RuntimeError(f"no symbol column in {gpl}: {list(tab.columns)[:20]}")
    sym = sym.astype(str).str.strip()
    sym[sym.str.contains("///") | sym.isin(["", "nan", "---"])] = np.nan
    sym.index = sym.index.astype(str)
    return sym.reindex(ids)

def collapse(expr, sym):
    e = expr.copy()
    e["__g"] = sym.reindex(e.index.astype(str)).values
    e = e.dropna(subset=["__g"])
    m = e.drop(columns="__g").mean(axis=1)
    e = e.loc[m.groupby(e["__g"]).idxmax().values]
    e.index = e.pop("__g")
    return e

def ensure_log(e):
    v = e.values[np.isfinite(e.values)]
    if np.nanpercentile(v, 99) > 100:
        log("  data appear unlogged -> log2(x+1)")
        e = np.log2(e.clip(lower=0) + 1)
    return e

# ---------------------------------------------------------------- statistics
def trigamma_inv(x):
    x = np.asarray(x, float)
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.all(-dif / y < 1e-8):
            break
    return y

def moderated_t(e, g1, g2):
    """limma-style eBayes two-group comparison; g1 = case mask, g2 = reference mask"""
    a, b = e.loc[:, g1].values, e.loc[:, g2].values
    n1, n2 = a.shape[1], b.shape[1]
    m1, m2 = np.nanmean(a, 1), np.nanmean(b, 1)
    v1, v2 = np.nanvar(a, 1, ddof=1), np.nanvar(b, 1, ddof=1)
    d = n1 + n2 - 2
    s2 = ((n1 - 1) * v1 + (n2 - 1) * v2) / d
    ok = np.isfinite(s2) & (s2 > 1e-12)
    z = np.log(s2[ok])
    ee = z - special.digamma(d / 2) + np.log(d / 2)
    emean = ee.mean()
    evar = ee.var(ddof=1) - special.polygamma(1, d / 2)
    if evar > 0:
        d0 = 2 * trigamma_inv(evar)
        s0 = np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0, s0 = 1e6, np.exp(emean)
    s2p = (d0 * s0 + d * s2) / (d0 + d)
    lfc = m1 - m2
    t = lfc / np.sqrt(s2p * (1 / n1 + 1 / n2))
    p = 2 * stats.t.sf(np.abs(t), d + d0)
    sp = np.sqrt(s2)
    J = 1 - 3 / (4 * (n1 + n2) - 9)
    g = J * lfc / sp
    vg = (n1 + n2) / (n1 * n2) + g ** 2 / (2 * (n1 + n2))
    res = pd.DataFrame({"logFC": lfc, "AveExpr": (m1 * n1 + m2 * n2) / (n1 + n2), "t": t, "P": p,
                        "g": g, "vg": vg}, index=e.index)
    res = res[np.isfinite(res["t"])]
    res["FDR"] = bh(res["P"].values)
    return res.sort_values("P")

def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); r = np.empty(n); r[o] = np.arange(1, n + 1)
    q = p * n / r
    qs = np.minimum.accumulate(q[o][::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(qs, 1)
    return out

def meta_re(tabs, min_k=2):
    """DerSimonian-Laird random-effects meta of Hedges g across datasets"""
    g = pd.concat({k: v["g"] for k, v in tabs.items()}, axis=1)
    vg = pd.concat({k: v["vg"] for k, v in tabs.items()}, axis=1)
    k = g.notna().sum(1)
    g, vg = g[k >= min_k], vg[k >= min_k]
    w = 1 / vg
    fe = (w * g).sum(1) / w.sum(1)
    Q = (w * (g.sub(fe, axis=0)) ** 2).sum(1)
    kk = g.notna().sum(1)
    c = w.sum(1) - (w ** 2).sum(1) / w.sum(1)
    tau2 = ((Q - (kk - 1)) / c).clip(lower=0)
    wr = 1 / (vg.add(tau2, axis=0))
    re_ = (wr * g).sum(1) / wr.sum(1)
    se = np.sqrt(1 / wr.sum(1))
    z = re_ / se
    p = 2 * stats.norm.sf(np.abs(z))
    I2 = ((Q - (kk - 1)) / Q.replace(0, np.nan)).clip(lower=0).fillna(0) * 100
    sign_cons = (np.sign(g).eq(np.sign(re_), axis=0) & g.notna()).sum(1) / kk
    out = pd.DataFrame({"meta_g": re_, "meta_se": se, "meta_P": p, "k": kk, "I2": I2,
                        "dir_consistency": sign_cons})
    out["meta_FDR"] = bh(out["meta_P"].values)
    for col in g.columns:
        out[f"g_{col}"] = g[col]
    return out.sort_values("meta_P")

# ---------------------------------------------------------------- group rules
def assign_groups(acc, info):
    t = sample_text(info)
    grp = pd.Series("other", index=info.index)
    comp = pd.Series("", index=info.index)
    if acc == "GSE30122":
        comp[t.str.contains(r"glom")] = "glom"
        comp[t.str.contains(r"tubul")] = "tub"
        grp[t.str.contains(r"diabet|\bdkd\b|\bdn\b")] = "DKD"
        grp[t.str.contains(r"control|normal|healthy")] = "CTRL"
        # GSM758485-GSM758509 re-deposit the same control individuals ("Control Human Glomeruli NN") -> drop duplicates
        dup = info["title"].str.match(r"^Control Human (Glomeruli|Tubuli) \d+$")
        grp[dup.values] = "duplicate"
    elif acc in ("GSE104948", "GSE104954"):
        comp[:] = "glom" if acc == "GSE104948" else "tub"
        grp[t.str.contains(r"diabetic nephropathy|\bdn\b")] = "DKD"
        ld = t.str.contains(r"living donor|\bld\b")
        tn = t.str.contains(r"tumor nephrectomy|\btn\b")
        if ld.sum() >= 5:
            grp[ld & (grp == "other")] = "CTRL"
        else:
            grp[(ld | tn) & (grp == "other")] = "CTRL"
        # samples with mixed diagnoses (e.g. 'DN, FSGS') -> exclude
        mixed = t.str.contains(r"fsgs|mcd|igan|iga |lupus|sle|mgn|rpgn|hypertens") & (grp == "DKD")
        grp[mixed] = "other"
    elif acc == "GSE96804":
        comp[:] = "glom"
        grp[t.str.contains(r"diabet|\bdn\b|\bdkd\b")] = "DKD"
        grp[t.str.contains(r"control|normal|healthy|tumor|nephrectomy") & (grp == "other")] = "CTRL"
    elif acc == "GSE142025":
        comp[:] = "whole"
        grp[t.str.contains(r"early")] = "EARLY"
        grp[t.str.contains(r"advanced|late")] = "ADV"
        grp[t.str.contains(r"control|normal|healthy|nephrectomy") & (grp == "other")] = "CTRL"
    return grp, comp

# ---------------------------------------------------------------- microarray datasets
ARRAYS = ["GSE30122", "GSE104948", "GSE104954", "GSE96804"]

def load_array_series(acc):
    out = []
    files = [f for f in list_dir(geo_dir(acc, "matrix")) if f.endswith("_series_matrix.txt.gz")]
    log(f"{acc}: matrix files {files}")
    for fn in files:
        p = fetch(geo_dir(acc, "matrix") + fn, os.path.join(RAW, fn))
        info, expr = parse_series_matrix(p)
        if expr is None:
            log(f"  {fn}: no expression table, skipped"); continue
        gpl = info["platform"].iloc[0]
        grp, comp = assign_groups(acc, info)
        info["group"], info["compartment"] = grp.values, comp.values
        info.to_csv(os.path.join(OUT, f"groups_{fn.replace('_series_matrix.txt.gz','')}.csv"),
                    index=False, encoding="utf-8-sig")
        expr = ensure_log(expr)
        sym = platform_map(gpl, expr.index)
        e = collapse(expr, sym)
        log(f"  {fn}: platform {gpl}, {e.shape[0]} genes x {e.shape[1]} samples; groups "
            f"{info.groupby(['compartment','group']).size().to_dict()}")
        out.append((fn, gpl, info, e))
    return out

def run_arrays():
    tabs = {"glom": {}, "tub": {}}
    for acc in ARRAYS:
        try:
            for fn, gpl, info, e in load_array_series(acc):
                for c in ("glom", "tub"):
                    sel = info[info["compartment"] == c]
                    case = sel.loc[sel["group"] == "DKD", "gsm"].tolist()
                    ref = sel.loc[sel["group"] == "CTRL", "gsm"].tolist()
                    if len(case) >= 3 and len(ref) >= 3:
                        ee = e[case + ref]
                        ee = ee[ee.mean(1) > np.nanpercentile(ee.mean(1), 10)]
                        r = moderated_t(ee, case, ref)
                        name = f"{acc}_{gpl}_{c}"
                        r.to_csv(os.path.join(OUT, f"DE_{name}_DKDvsCTRL.csv"), encoding="utf-8-sig")
                        tabs[c][f"{acc}_{gpl}"] = r
                        log(f"  DE {name}: DKD n={len(case)} vs CTRL n={len(ref)}; "
                            f"FDR<0.05: {(r.FDR<0.05).sum()}")
                    elif len(case) or len(ref):
                        log(f"  {acc} {gpl} {c}: too few samples (DKD {len(case)}, CTRL {len(ref)}) - skipped")
        except Exception as ex:
            log(f"!! {acc} failed: {ex}"); log(traceback.format_exc())
    metas = {}
    for c in ("glom", "tub"):
        if len(tabs[c]) >= 2:
            m = meta_re(tabs[c])
            m.to_csv(os.path.join(OUT, f"META_{c}_DKDvsCTRL.csv"), encoding="utf-8-sig")
            metas[c] = m
            log(f"META {c}: {len(tabs[c])} datasets, {len(m)} genes, meta FDR<0.05 & |g|>0.8: "
                f"{((m.meta_FDR<0.05)&(m.meta_g.abs()>0.8)).sum()}")
        elif len(tabs[c]) == 1:
            k, r = list(tabs[c].items())[0]
            metas[c] = pd.DataFrame({"meta_g": r["g"], "meta_P": r["P"], "meta_FDR": r["FDR"], "k": 1})
            log(f"META {c}: only one dataset ({k}) - used directly")
    return tabs, metas

# ---------------------------------------------------------------- GSE142025 RNA-seq
def read_any_table(path):
    opener = gzip.open if path.endswith(".gz") else open
    name = path.lower()
    if name.endswith((".xlsx", ".xls")):
        return pd.read_excel(path)
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        txt = f.read()
    sep = "\t" if txt.count("\t") > txt.count(",") else ","
    return pd.read_csv(io.StringIO(txt), sep=sep)

def gene_ids_to_symbol(ids):
    ids = pd.Index(ids).astype(str).str.replace(r"\.\d+$", "", regex=True)
    gi = gene_info()
    if ids.str.startswith("ENSG").mean() > 0.5:
        m = dict(zip(gi["ens"].dropna(), gi.loc[gi["ens"].notna(), "Symbol"]))
        return [m.get(i) for i in ids]
    if ids.str.match(r"^\d+$").mean() > 0.5:
        m = dict(zip(gi["GeneID"], gi["Symbol"]))
        return [m.get(i) for i in ids]
    return list(ids)

def run_progression():
    acc = "GSE142025"
    files = [f for f in list_dir(geo_dir(acc, "matrix")) if f.endswith("_series_matrix.txt.gz")]
    p = fetch(geo_dir(acc, "matrix") + files[0], os.path.join(RAW, files[0]))
    info, _ = parse_series_matrix(p)
    grp, _ = assign_groups(acc, info)
    info["group"] = grp.values
    info.to_csv(os.path.join(OUT, f"groups_{acc}.csv"), index=False, encoding="utf-8-sig")
    log(f"{acc}: groups {info['group'].value_counts().to_dict()}")
    supp = list_dir(geo_dir(acc, "suppl"))
    log(f"{acc}: supplementary {supp}")
    cand = [f for f in supp if f.startswith(acc) and not f.endswith(".tar")
            and re.search(r"\.(txt|csv|tsv|xlsx|xls)(\.gz)?$", f)]
    mats = []
    for fn in cand:
        mats.append(read_any_table(fetch(geo_dir(acc, "suppl") + fn, os.path.join(RAW, fn))))
    if not mats:
        tars = [f for f in supp if f.endswith(".tar")]
        tp = fetch(geo_dir(acc, "suppl") + tars[0], os.path.join(RAW, tars[0]))
        cols = {}
        with tarfile.open(tp) as tf:
            for m in tf.getmembers():
                gsm = re.match(r"(GSM\d+)", os.path.basename(m.name))
                if not gsm: continue
                raw = tf.extractfile(m).read()
                if m.name.endswith(".gz"): raw = gzip.decompress(raw)
                d = pd.read_csv(io.BytesIO(raw), sep=None, engine="python", header=None)
                if pd.to_numeric(d.iloc[0, 1:], errors="coerce").isna().all():
                    d = d.iloc[1:]
                vals = d.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
                vcol = vals.notna().mean().idxmax()
                ser = pd.Series(vals[vcol].values, index=d.iloc[:, 0].astype(str)).dropna()
                cols[gsm.group(1)] = ser[~ser.index.duplicated()]
        mats = [pd.DataFrame(cols).reset_index()]
    d = mats[0]
    log(f"  count table shape {d.shape}; first columns {list(d.columns[:6])}")
    idcol = d.columns[0]
    num = d.drop(columns=[idcol]).apply(pd.to_numeric, errors="coerce")
    num = num.loc[:, num.notna().mean() > 0.9]
    num.index = d[idcol].astype(str)
    # map columns to samples
    title_norm = {re.sub(r"[^a-z0-9]", "", str(t).lower()): g for t, g in zip(info["title"], info["gsm"])}
    colmap = {}
    for c in num.columns:
        cn = re.sub(r"[^a-z0-9]", "", str(c).lower())
        m = re.search(r"GSM\d+", str(c))
        if m and m.group(0) in set(info["gsm"]): colmap[c] = m.group(0)
        elif cn in title_norm: colmap[c] = title_norm[cn]
    if len(colmap) >= 0.8 * num.shape[1]:
        num = num[list(colmap)].rename(columns=colmap)
        g = info.set_index("gsm")["group"].reindex(num.columns)
    else:
        log("  columns not matched to GSM/title; grouping by column names")
        cn = pd.Series([str(c).lower() for c in num.columns], index=num.columns)
        g = pd.Series("other", index=num.columns)
        g[cn.str.contains("early")] = "EARLY"
        g[cn.str.contains("adv|late")] = "ADV"
        g[cn.str.contains("control|normal|ctrl|^c\\d|^n\\d") & (g == "other")] = "CTRL"
    pd.DataFrame({"column": num.columns, "group": g.values}).to_csv(
        os.path.join(OUT, "groups_GSE142025_columns.csv"), index=False, encoding="utf-8-sig")
    log(f"  column groups {g.value_counts().to_dict()}")
    is_count = np.nanmax(num.values) > 1000 and np.allclose(num.fillna(0).values % 1, 0)
    sym = gene_ids_to_symbol(num.index)
    num.index = sym
    num = num[pd.notna(num.index)]
    num = num.groupby(level=0).sum() if is_count else num.groupby(level=0).mean()
    if is_count:
        cpm = num / num.sum(0) * 1e6
        keep = (cpm > 1).sum(1) >= max(3, int(0.25 * num.shape[1]))
        k = num[keep]
        # DESeq2 median-of-ratios size factors
        pos = k[(k > 0).all(1)]
        lg = np.log(pos)
        sf = np.exp((lg.sub(lg.mean(1), axis=0)).median(0))
        e = np.log2(k / sf + 1)
    else:
        e = ensure_log(num)
        e = e[e.mean(1) > np.nanpercentile(e.mean(1), 25)]
    log(f"  {'counts' if is_count else 'normalized values'} -> {e.shape[0]} genes")
    res = {}
    for name, a, b in [("ADVvsEARLY", "ADV", "EARLY"), ("ADVvsCTRL", "ADV", "CTRL"), ("EARLYvsCTRL", "EARLY", "CTRL")]:
        ca, cb = g.index[g == a].tolist(), g.index[g == b].tolist()
        if len(ca) >= 3 and len(cb) >= 3:
            r = moderated_t(e[ca + cb], ca, cb)
            r.to_csv(os.path.join(OUT, f"DE_GSE142025_{name}.csv"), encoding="utf-8-sig")
            res[name] = r
            log(f"  DE {name}: n={len(ca)} vs {len(cb)}; FDR<0.05: {(r.FDR<0.05).sum()}")
    return res, e, g

# ---------------------------------------------------------------- annotation & enrichment
def hpa():
    try:
        p = fetch("https://www.proteinatlas.org/download/proteinatlas.tsv.zip", os.path.join(RAW, "proteinatlas.tsv.zip"))
        with zipfile.ZipFile(p) as z:
            d = pd.read_csv(z.open(z.namelist()[0]), sep="\t", low_memory=False)
        keep = ["Gene"] + [c for c in d.columns if re.search(r"Secretome location|Protein class|Subcellular location|Conc\. blood", c)]
        return d[keep].drop_duplicates("Gene").set_index("Gene")
    except Exception as ex:
        log(f"!! HPA annotation unavailable: {ex}")
        return None

def gprofiler(genes, background, label):
    if len(genes) < 5:
        return None
    body = json.dumps({"organism": "hsapiens", "query": list(genes), "sources": ["GO:BP", "KEGG", "REAC", "WP"],
                       "user_threshold": 0.05, "domain_scope": "custom", "background": list(background),
                       "no_evidences": True}).encode()
    for k in range(3):
        try:
            req = urllib.request.Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/", data=body,
                                         headers={"Content-Type": "application/json", **UA})
            with urllib.request.urlopen(req, timeout=180) as r:
                js = json.loads(r.read().decode())
            rows = [{"source": x["source"], "term_id": x["native"], "term": x["name"], "P_adj": x["p_value"],
                     "term_size": x["term_size"], "overlap": x["intersection_size"],
                     "query_size": x["query_size"]} for x in js.get("result", [])]
            df = pd.DataFrame(rows).sort_values("P_adj") if rows else pd.DataFrame()
            df.to_csv(os.path.join(OUT, f"ENRICH_{label}.csv"), index=False, encoding="utf-8-sig")
            log(f"  enrichment {label}: {len(df)} terms")
            return df
        except Exception as ex:
            log(f"  g:Profiler retry {k+1}: {ex}"); time.sleep(10)
    return None

# ---------------------------------------------------------------- figures
def figures(prog, metas, integ, enr):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "Arial", "font.size": 9, "pdf.fonttype": 42})
    if "ADVvsEARLY" in prog:
        r = prog["ADVvsEARLY"]
        fig, ax = plt.subplots(figsize=(4.2, 4))
        sig = r.FDR < 0.05
        ax.scatter(r.logFC[~sig], -np.log10(r.P[~sig]), s=3, c="#BBBBBB", lw=0)
        ax.scatter(r.logFC[sig & (r.logFC > 0)], -np.log10(r.P[sig & (r.logFC > 0)]), s=4, c="#C0392B", lw=0)
        ax.scatter(r.logFC[sig & (r.logFC < 0)], -np.log10(r.P[sig & (r.logFC < 0)]), s=4, c="#2E86C1", lw=0)
        top = integ.head(12).index if integ is not None and len(integ) else r.head(12).index
        for gname in top:
            if gname in r.index:
                ax.annotate(gname, (r.at[gname, "logFC"], -np.log10(r.at[gname, "P"])), fontsize=6.5)
        ax.set_xlabel("log2 fold change (advanced vs early DKD)"); ax.set_ylabel("-log10 P")
        ax.set_title("GSE142025: advanced vs early DKD")
        fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig1_volcano_ADVvsEARLY.png"), dpi=300)
        fig.savefig(os.path.join(FIG, "Fig1_volcano_ADVvsEARLY.pdf")); plt.close(fig)
    if integ is not None and len(integ):
        top = integ.head(30)
        cols = [c for c in integ.columns if c.startswith("g_") or c in ("glom_meta_g", "tub_meta_g", "ADVvsEARLY_g")]
        cols = [c for c in ["glom_meta_g", "tub_meta_g", "ADVvsEARLY_g"] if c in integ.columns]
        M = top[cols].astype(float)
        fig, ax = plt.subplots(figsize=(3.2, 7))
        im = ax.imshow(M.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
        ax.set_yticks(range(len(M))); ax.set_yticklabels(M.index, fontsize=7)
        ax.set_xticks(range(len(cols))); ax.set_xticklabels(["Glom\nDKD vs Ctrl", "Tub\nDKD vs Ctrl", "Adv vs\nEarly DKD"][:len(cols)], fontsize=7)
        fig.colorbar(im, ax=ax, shrink=0.4, label="Hedges g")
        fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig2_heatmap_top30.png"), dpi=300)
        fig.savefig(os.path.join(FIG, "Fig2_heatmap_top30.pdf")); plt.close(fig)
    for lab, df in enr.items():
        if df is None or not len(df):
            continue
        d = df.head(15).iloc[::-1]
        fig, ax = plt.subplots(figsize=(6, 4.5))
        ax.barh(range(len(d)), -np.log10(d.P_adj), color="#C0392B" if "up" in lab else "#2E86C1")
        ax.set_yticks(range(len(d))); ax.set_yticklabels([f"{t[:55]} ({s})" for t, s in zip(d.term, d.source)], fontsize=7)
        ax.set_xlabel("-log10 adjusted P"); ax.set_title(lab)
        fig.tight_layout(); fig.savefig(os.path.join(FIG, f"Fig3_enrich_{lab}.png"), dpi=300); plt.close(fig)

# ---------------------------------------------------------------- report
PROPOSAL_MARKERS = ["TNFRSF1A", "TNFRSF1B", "PLAUR", "LCN2", "HAVCR1", "FABP1", "TGFB1", "CCL2", "SPP1", "COL1A1", "FN1"]

def html_report(prog, metas, integ, enr, prop):
    def tbl(df, n=40, cols=None):
        if df is None or not len(df): return "<p>（无结果）</p>"
        d = df.head(n) if cols is None else df[cols].head(n)
        return d.to_html(float_format=lambda x: f"{x:.3g}", border=0, classes="t")
    lines = open(os.path.join(BASE, "run_log.txt"), encoding="utf-8").read().splitlines()[-80:]
    figs = sorted(f for f in os.listdir(FIG) if f.endswith(".png"))
    h = f"""<!doctype html><html><head><meta charset="utf-8"><title>DKD公共组织转录组分析</title>
<style>body{{font-family:"Microsoft YaHei",Arial;max-width:1100px;margin:24px auto;padding:0 16px;color:#222}}
h1{{font-size:22px}}h2{{font-size:17px;border-bottom:1px solid #ccc;padding-bottom:4px;margin-top:28px}}
table.t{{border-collapse:collapse;font-size:12px}}table.t td,table.t th{{border-bottom:1px solid #eee;padding:3px 6px;text-align:right}}
img{{max-width:100%;border:1px solid #eee;margin:6px 0}}pre{{background:#f6f6f6;padding:8px;font-size:11px;overflow:auto}}</style></head><body>
<h1>DKD肾组织公共转录组：进展相关基因与血液可检测候选分子</h1>
<p>生成时间 {time.strftime('%Y-%m-%d %H:%M')}。所有数值表在 results 文件夹，图在 figures 文件夹（PNG 300 dpi + PDF）。</p>
<h2>1. 进展相关候选基因（组织 DKD vs 对照 与 晚期 vs 早期DKD 方向一致）</h2>
<p>筛选：至少一个部位（肾小球/小管间质）meta FDR&lt;0.05 且 |g|&gt;0.8；GSE142025 晚期 vs 早期 FDR&lt;0.05 且方向一致。按综合秩排序；secreted 列来自 HPA 分泌组注释。</p>
{tbl(integ, 60)}
<h2>2. 申请书中提及的标志物在组织中的表现</h2>{tbl(prop, 30)}
<h2>3. 通路富集</h2>"""
    for lab, df in enr.items():
        h += f"<h3>{lab}</h3>" + tbl(df, 20, ["source", "term_id", "term", "P_adj", "overlap"] if df is not None and len(df) else None)
    h += "<h2>4. 图</h2>" + "".join(f'<p>{f}</p><img src="figures/{f}">' for f in figs)
    h += "<h2>5. 运行日志（末80行）</h2><pre>" + "\n".join(lines) + "</pre></body></html>"
    open(os.path.join(BASE, "DKD公共数据库分析报告.html"), "w", encoding="utf-8").write(h)

# ---------------------------------------------------------------- main
def main():
    log("=" * 70); log("start", time.strftime("%Y-%m-%d %H:%M:%S"))
    tabs, metas = run_arrays()
    try:
        prog, e_prog, g_prog = run_progression()
    except Exception as ex:
        log(f"!! GSE142025 failed: {ex}"); log(traceback.format_exc()); prog = {}
    # integration
    integ = None
    parts = {}
    for c, m in metas.items():
        parts[f"{c}_meta_g"] = m["meta_g"]; parts[f"{c}_meta_FDR"] = m["meta_FDR"]
        if "k" in m: parts[f"{c}_k"] = m["k"]
        if "dir_consistency" in m: parts[f"{c}_dir"] = m["dir_consistency"]
    if "ADVvsEARLY" in prog:
        r = prog["ADVvsEARLY"]
        parts.update({"ADVvsEARLY_logFC": r["logFC"], "ADVvsEARLY_g": r["g"], "ADVvsEARLY_P": r["P"], "ADVvsEARLY_FDR": r["FDR"]})
    if "ADVvsCTRL" in prog:
        parts["ADVvsCTRL_logFC"] = prog["ADVvsCTRL"]["logFC"]; parts["ADVvsCTRL_FDR"] = prog["ADVvsCTRL"]["FDR"]
    allt = pd.DataFrame(parts)
    H = hpa()
    if H is not None:
        allt = allt.join(H, how="left")
        sec = [c for c in H.columns if "Secretome location" in c]
        if sec:
            allt["secreted"] = allt[sec[0]].fillna("")
    allt.to_csv(os.path.join(OUT, "ALL_genes_integrated.csv"), encoding="utf-8-sig")
    if "ADVvsEARLY_FDR" in allt:
        tissue_sig = pd.Series(False, index=allt.index); tissue_dir = pd.Series(0.0, index=allt.index)
        for c in ("glom", "tub"):
            if f"{c}_meta_g" in allt:
                s = (allt[f"{c}_meta_FDR"] < 0.05) & (allt[f"{c}_meta_g"].abs() > 0.8)
                tissue_dir = tissue_dir.where(~(s & (tissue_dir == 0)), np.sign(allt[f"{c}_meta_g"]))
                tissue_sig |= s.fillna(False)
        prog_sig = (allt["ADVvsEARLY_FDR"] < 0.05)
        same = np.sign(allt["ADVvsEARLY_logFC"]) == tissue_dir
        integ = allt[tissue_sig & prog_sig & same].copy()
        rk = allt["ADVvsEARLY_P"].rank()
        for c in ("glom", "tub"):
            if f"{c}_meta_g" in allt:
                rk = rk + allt[f"{c}_meta_FDR"].rank().fillna(len(allt))
        integ["rank_score"] = rk.reindex(integ.index)
        integ["direction"] = np.where(integ["ADVvsEARLY_logFC"] > 0, "up", "down")
        integ = integ.sort_values("rank_score")
        integ.to_csv(os.path.join(OUT, "PROGRESSION_candidates.csv"), encoding="utf-8-sig")
        log(f"INTEGRATED progression genes: {len(integ)} (up {(integ.direction=='up').sum()}, down {(integ.direction=='down').sum()})")
        if "secreted" in integ:
            b = integ[integ["secreted"].str.contains("blood", case=False, na=False)]
            b.to_csv(os.path.join(OUT, "PROGRESSION_candidates_secreted_to_blood.csv"), encoding="utf-8-sig")
            log(f"  of which secreted to blood (HPA): {len(b)}")
    enr = {}
    if integ is not None and len(integ):
        bg = allt.index[allt["ADVvsEARLY_P"].notna()].tolist()
        enr["progression_up"] = gprofiler(integ.index[integ.direction == "up"], bg, "progression_up")
        enr["progression_down"] = gprofiler(integ.index[integ.direction == "down"], bg, "progression_down")
    prop = allt.reindex([g for g in PROPOSAL_MARKERS if g in allt.index])
    prop.to_csv(os.path.join(OUT, "proposal_markers.csv"), encoding="utf-8-sig")
    try:
        figures(prog, metas, integ, enr)
    except Exception as ex:
        log(f"!! figures failed: {ex}"); log(traceback.format_exc())
    show = None
    if integ is not None:
        cols = [c for c in ["direction", "glom_meta_g", "glom_meta_FDR", "tub_meta_g", "tub_meta_FDR",
                            "ADVvsEARLY_logFC", "ADVvsEARLY_FDR", "secreted"] if c in integ.columns]
        show = integ[cols]
    pcols = [c for c in ["glom_meta_g", "glom_meta_FDR", "tub_meta_g", "tub_meta_FDR", "ADVvsEARLY_logFC", "ADVvsEARLY_FDR"] if c in prop.columns]
    html_report(prog, metas, show, enr, prop[pcols] if len(prop) else prop)
    log("done", time.strftime("%Y-%m-%d %H:%M:%S"))

if __name__ == "__main__":
    try:
        main()
    except Exception as ex:
        log("!! FATAL", ex); log(traceback.format_exc())
