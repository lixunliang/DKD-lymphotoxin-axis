# -*- coding: utf-8 -*-
"""
Genetic triangulation: is the LTB genetic signal for diabetic nephropathy specific, or an MHC / type-1-diabetes artefact?

Questions
  Q1  LTB lead cis-pQTL (rs2395470, ~300 kb from LTB, next to HLA-C) and the LTB-proximal conditional signal
      (rs3093553, inside the LTA-TNF-LTB cluster) -> diabetic nephropathy, T1D, T2D, T1D/T2D with renal complications,
      eGFR in people WITH diabetes (CKDGen, Winkler 2022), rapid eGFR decline / incident CKD (Gorski 2021), UACR.
  Q2  Do all MHC-region Olink proteins that "hit" diabetic nephropathy simply track their effect on type 1 diabetes?
      (Wald ratio on DN vs Wald ratio on T1D across all MHC proteins; LTB's position relative to that line.)
  Q3  Pathway genes OUTSIDE the MHC (RELB, NFKB2, MAP3K14, LTBR, CCL19, CCL21, CXCL13, CCR7): blood cis-eQTL (eQTLGen)
      Mendelian randomization on the same outcomes.

Data
  read-only from $MR_PROJ (folder with previously downloaded outcome GWAS and the Olink cis-instrument table)
      raw\\outcome\\finngen_r12_dm_nephropathy, finngen_r12_n14_chronkidneydis, ckdgen_winkler2022_dmstrat,
      ckdgen_gorski2021_rapid3, ckdgen_gorski2021_ckdi25, ckdgen_teumer2019_uacr_ea ; extracted\\olink_cis_instruments.tsv
  downloaded here into raw\\gwas  (resumable; about 3-4 GB)
      FinnGen R12: T1D, T2D, E4_DM1REN, E4_DM2REN ; eQTLGen cis-eQTL (FDR<0.05)
Outputs  results\\S5_MR_robustness_all.tsv, results\\S5_MR_MHC_vs_T1D.tsv, results\\S5_eQTL_instruments.tsv, log_09_genetics.txt
"""
import os, sys, gzip, math, time, urllib.request, ssl, traceback
try:
    import certifi; ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())
except Exception:
    pass
import numpy as np, pandas as pd

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
MRP = os.environ.get("MR_PROJ", os.path.join(BASE, "external_gwas"))   # folder holding raw/outcome/... and extracted/olink_cis_instruments.tsv
RAWO = os.path.join(MRP, "raw", "outcome")
GW = os.path.join(BASE, "raw", "gwas"); OUT = os.path.join(BASE, "results")
for d in (GW, OUT): os.makedirs(d, exist_ok=True)
LOGF = open(os.path.join(BASE, "log_09_genetics.txt"), "a", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOGF.write(s + "\n"); LOGF.flush()
UA = {"User-Agent": "Mozilla/5.0 (research script)"}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

# ---------------------------------------------------------------- downloads
FG = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_{}.gz"
NEW_FG = {"FG_T1D": ["T1D", "E4_DM1"], "FG_T2D": ["T2D", "E4_DM2"], "FG_T1D_REN": ["E4_DM1REN"], "FG_T2D_REN": ["E4_DM2REN"]}
EQTLGEN = ["https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"]

def head_ok(url):
    try:
        req = urllib.request.Request(url, headers={**UA, "Range": "bytes=0-1"})
        with urllib.request.urlopen(req, timeout=60) as r:
            cr = r.headers.get("Content-Range", "")
            return int(cr.split("/")[-1]) if "/" in cr else int(r.headers.get("Content-Length", 0) or 0)
    except Exception as e:
        log(f"    not available: {url} ({type(e).__name__}: {e})"); return None

def fetch(url, dest, total, tries=40):
    if os.path.exists(dest) and (not total or os.path.getsize(dest) == total): return True
    part = dest + ".part"
    for k in range(tries):
        try:
            have = os.path.getsize(part) if os.path.exists(part) else 0
            if total and have >= total: break
            hdr = dict(UA)
            if have: hdr["Range"] = f"bytes={have}-"
            with urllib.request.urlopen(urllib.request.Request(url, headers=hdr), timeout=300) as r:
                if have and r.status != 206: have = 0
                t0, last = time.time(), 0
                with open(part, "ab" if have else "wb") as f:
                    while True:
                        b = r.read(8 << 20)
                        if not b: break
                        f.write(b); have += len(b)
                        if time.time() - last > 20:
                            last = time.time()
                            log(f"    {os.path.basename(dest)}: {have/1e6:,.0f} / {total/1e6:,.0f} MB")
            if not total or have >= total: break
        except Exception as e:
            log(f"    interrupted ({type(e).__name__}); retry {k+1} in 10 s"); time.sleep(10)
    if os.path.exists(part) and (not total or os.path.getsize(part) >= total):
        os.replace(part, dest); return True
    return False

def downloads():
    paths = {}
    if os.environ.get("MR_SKIP_NET") == "1":
        for key, names in NEW_FG.items():
            for nm in names:
                p = os.path.join(GW, f"finngen_R12_{nm}.gz")
                if os.path.exists(p): paths[key] = p; break
        p = os.path.join(GW, os.path.basename(EQTLGEN[0]))
        if os.path.exists(p): paths["EQTLGEN"] = p
        return paths
    for key, names in NEW_FG.items():
        for nm in names:
            url = FG.format(nm); total = head_ok(url)
            if total:
                dest = os.path.join(GW, os.path.basename(url)); log(f"  {key}: {nm} ({total/1e6:,.0f} MB)")
                if fetch(url, dest, total): paths[key] = dest
                break
        if key not in paths: log(f"  !! {key}: no FinnGen endpoint downloaded")
    for url in EQTLGEN:
        total = head_ok(url)
        if total:
            dest = os.path.join(GW, os.path.basename(url)); log(f"  eQTLGen cis ({total/1e6:,.0f} MB)")
            if fetch(url, dest, total): paths["EQTLGEN"] = dest; break
    if "EQTLGEN" not in paths: log("  !! eQTLGen not downloaded - Q3 will be skipped")
    return paths

# ---------------------------------------------------------------- instruments
AXIS_NON_MHC = ["RELB", "NFKB2", "MAP3K14", "LTBR", "CCL19", "CCL21", "CXCL13", "CCR7", "NFKB1"]

def olink_instruments():
    p = os.path.join(MRP, "extracted", "olink_cis_instruments.tsv")
    d = pd.read_csv(p, sep="\t")
    d = d[d.EA.str.len().eq(1) & d.OA.str.len().eq(1)]
    pal = d.apply(lambda r: {r.EA, r.OA} in ({"A", "T"}, {"C", "G"}), axis=1)
    d = d[~pal]
    keep = (d.MHC == 1) & (d.source == "ST9_lead")                                    # every MHC lead (Q2)
    keep |= d.Gene.isin(["LTB", "LTA", "TNF"])                                         # all LTB/LTA/TNF signals (Q1)
    keep |= d.Gene.isin(["LTBR", "CCL19", "CCL21", "CXCL13", "IL7R", "NFKB1"]) & (d.source == "ST9_lead")
    d = d[keep].copy(); d["kind"] = "pQTL_Olink"
    d = d.rename(columns={"Beta": "bx", "SE": "sx"})[["kind", "ProteinID", "Gene", "source", "rsid", "CHR", "POS37", "EA", "OA", "EAF", "bx", "sx", "log10p", "MHC"]]
    log(f"  Olink instruments: {len(d)} ({d.Gene.nunique()} proteins; MHC leads {int(((d.MHC == 1) & (d.source == 'ST9_lead')).sum())})")
    return d

def eqtl_instruments(path):
    if not path: return pd.DataFrame()
    want = set(AXIS_NON_MHC + ["LTB"]); rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        h = f.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
        for line in f:
            v = line.rstrip("\n").split("\t")
            if v[ix["GeneSymbol"]] in want:
                rows.append({c: v[ix[c]] for c in ("SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele", "Zscore", "GeneSymbol", "NrSamples", "Pvalue", "FDR")})
    e = pd.DataFrame(rows)
    if e.empty: log("  !! no eQTLs found for axis genes"); return e
    for c in ("Zscore", "NrSamples", "Pvalue", "FDR"): e[c] = pd.to_numeric(e[c], errors="coerce")
    e = e[(e.FDR < 0.05) & e.AssessedAllele.str.len().eq(1) & e.OtherAllele.str.len().eq(1)]
    e = e[~e.apply(lambda r: {r.AssessedAllele, r.OtherAllele} in ({"A", "T"}, {"C", "G"}), axis=1)]
    e["absz"] = e.Zscore.abs(); lead = e.sort_values("absz", ascending=False).drop_duplicates("GeneSymbol")
    out = pd.DataFrame({"kind": "eQTL_eQTLGen", "ProteinID": lead.GeneSymbol, "Gene": lead.GeneSymbol, "source": "eQTLGen_lead",
                        "rsid": lead.SNP, "CHR": lead.SNPChr, "POS37": lead.SNPPos, "EA": lead.AssessedAllele, "OA": lead.OtherAllele,
                        "EAF": np.nan, "bx": lead.Zscore, "sx": np.nan, "log10p": -np.log10(lead.Pvalue.clip(lower=1e-300)), "MHC": 0, "N": lead.NrSamples})
    log(f"  eQTLGen lead cis-eQTLs: {out[['Gene', 'rsid']].values.tolist()}")
    return out

# ---------------------------------------------------------------- outcomes
def outcome_list(paths):
    L = [("FG_DN", "Diabetic nephropathy (FinnGen R12)", os.path.join(RAWO, "finngen_r12_dm_nephropathy", "finngen_R12_DM_NEPHROPATHY.gz"), "finngen", "logOR"),
         ("FG_CKD", "Chronic kidney disease (FinnGen R12)", os.path.join(RAWO, "finngen_r12_n14_chronkidneydis", "finngen_R12_N14_CHRONKIDNEYDIS.gz"), "finngen", "logOR")]
    lab = {"FG_T1D": "Type 1 diabetes (FinnGen R12)", "FG_T2D": "Type 2 diabetes (FinnGen R12)",
           "FG_T1D_REN": "T1D with renal complications (FinnGen R12)", "FG_T2D_REN": "T2D with renal complications (FinnGen R12)"}
    for k, l in lab.items():
        if k in paths: L.append((k, l, paths[k], "finngen", "logOR"))
    L += [("EGFR_DM", "log eGFR, people with diabetes (CKDGen 2022)", os.path.join(RAWO, "ckdgen_winkler2022_dmstrat", "meta_egfr_dmstrat_stage1plus2.txt.gz"), "winkler_dm", "beta"),
          ("EGFR_NODM", "log eGFR, people without diabetes (CKDGen 2022)", os.path.join(RAWO, "ckdgen_winkler2022_dmstrat", "meta_egfr_dmstrat_stage1plus2.txt.gz"), "winkler_nodm", "beta"),
          ("RAPID3", "Rapid eGFR decline (CKDGen 2021)", os.path.join(RAWO, "ckdgen_gorski2021_rapid3", "CKDGen_rapid3_overall.txt.gz"), "ckdgen", "logOR"),
          ("CKDI25", "Incident CKD / eGFR drop 25% (CKDGen 2021)", os.path.join(RAWO, "ckdgen_gorski2021_ckdi25", "CKDGen_ckdi25_overall.txt.gz"), "ckdgen", "logOR"),
          ("UACR", "log UACR (CKDGen 2019)", os.path.join(RAWO, "ckdgen_teumer2019_uacr_ea", "formatted_20180517-UACR_overall-EA-nstud_18-SumMac_400.tbl.rsid.gz"), "ckdgen", "beta")]
    return L

def opener(p):
    with open(p, "rb") as f: gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(p, "rt", encoding="utf-8", errors="replace") if gz else open(p, "r", encoding="utf-8", errors="replace")

def scan(path, fmt, rsids):
    """returns {rsid: (EA, OA, beta, se, eaf)}"""
    d = {}; n = 0; t0 = time.time()
    with opener(path) as f:
        head = f.readline().split(); low = [h.lower().lstrip("#") for h in head]; ix = {c: i for i, c in enumerate(low)}
        if fmt == "finngen":
            irs, iea, ioa, ib, ise, iaf = ix["rsids"], ix["alt"], ix["ref"], ix["beta"], ix["sebeta"], ix["af_alt"]
            for line in f:
                n += 1
                v = line.rstrip("\n").split("\t")
                if len(v) <= iaf: continue
                for rs in v[irs].split(","):
                    if rs in rsids:
                        try: d[rs] = (v[iea].upper(), v[ioa].upper(), float(v[ib]), float(v[ise]), float(v[iaf]))
                        except ValueError: pass
        elif fmt.startswith("winkler"):
            sfx = "dm" if fmt == "winkler_dm" else "nodm"
            irs, iea, ioa = ix["rsid"], ix["allele1"], ix["allele2"]
            ib, ise = ix[f"effect.{sfx}"], ix[f"stderr.{sfx}"]; iaf = ix.get(f"freq1.{sfx}")
            for line in f:
                n += 1
                v = line.split()
                if len(v) > ib and v[irs] in rsids:
                    try: d[v[irs]] = (v[iea].upper(), v[ioa].upper(), float(v[ib]), float(v[ise]), float(v[iaf]) if iaf is not None else np.nan)
                    except ValueError: pass
        else:
            irs, iea, ioa, ise = ix["rsid"], ix["allele1"], ix["allele2"], ix["stderr"]
            ib, ior, iaf = ix.get("effect"), ix.get("or"), ix.get("freq1")
            for line in f:
                n += 1
                v = line.split()
                if len(v) > ise and v[irs] in rsids:
                    try:
                        b = float(v[ib]) if ib is not None else math.log(float(v[ior]))
                        d[v[irs]] = (v[iea].upper(), v[ioa].upper(), b, float(v[ise]), float(v[iaf]) if iaf is not None else np.nan)
                    except ValueError: pass
    log(f"    scanned {n:,} rows in {(time.time()-t0)/60:.1f} min; matched {len(d)}/{len(rsids)} instruments")
    return d

def align(ea, oa, ea2, oa2):
    if (ea, oa) == (ea2, oa2): return 1
    if (ea, oa) == (oa2, ea2): return -1
    c = (COMP.get(ea2, "?"), COMP.get(oa2, "?"))
    if (ea, oa) == c: return 1
    if (ea, oa) == (c[1], c[0]): return -1
    return 0

def pz(z): return math.erfc(abs(z) / math.sqrt(2.0))

# ---------------------------------------------------------------- main
def main():
    log("=" * 78); log("genetic triangulation start", time.strftime("%Y-%m-%d %H:%M:%S"))
    log("A. downloads (resumable; rerun this bat if interrupted)")
    paths = downloads()
    log("B. instruments")
    ins = olink_instruments()
    eq = eqtl_instruments(paths.get("EQTLGEN"))
    allins = pd.concat([ins, eq], ignore_index=True)
    allins.to_csv(os.path.join(OUT, "S5_instruments.tsv"), sep="\t", index=False)
    rsids = set(allins.rsid.astype(str))
    log("C. outcomes")
    assoc = {}
    for key, lab, path, fmt, scale in outcome_list(paths):
        log(f"  {key}: {lab}")
        if not os.path.exists(path): log(f"    !! missing file {path}"); continue
        try: assoc[key] = (lab, scale, scan(path, fmt, rsids))
        except Exception as e: log(f"    !! scan failed: {type(e).__name__}: {e}"); log(traceback.format_exc())
    # eQTL effect sizes: z -> per-allele SD units using allele frequency from FinnGen DN (Finnish) or CKDGen
    if len(eq):
        for i, r in allins[allins.kind == "eQTL_eQTLGen"].iterrows():
            f = np.nan
            for key in ("FG_DN", "RAPID3", "EGFR_DM"):
                o = assoc.get(key, (None, None, {}))[2].get(str(r.rsid))
                if o is not None and np.isfinite(o[4]):
                    s = align(r.EA, r.OA, o[0], o[1]); f = o[4] if s == 1 else (1 - o[4] if s == -1 else np.nan)
                    if np.isfinite(f): break
            z, N = float(r.bx), float(r.N)
            if np.isfinite(f) and 0 < f < 1:
                b = z / math.sqrt(2 * f * (1 - f) * (N + z * z)); allins.loc[i, ["bx", "sx", "EAF"]] = [b, abs(b / z), f]
    log("D. Wald ratios")
    rows = []
    for _, r in allins.iterrows():
        if not (np.isfinite(r.bx) and np.isfinite(r.sx) and r.bx != 0): continue
        for key, (lab, scale, d) in assoc.items():
            o = d.get(str(r.rsid))
            if o is None: continue
            s = align(r.EA, r.OA, o[0], o[1])
            if s == 0: continue
            by, sy = s * o[2], o[3]
            w, wse = by / r.bx, abs(sy / r.bx)
            rows.append(dict(kind=r.kind, Gene=r.Gene, ProteinID=r.ProteinID, source=r.source, rsid=r.rsid, chr=r.CHR, pos37=r.POS37,
                             MHC=int(r.MHC), EA=r.EA, OA=r.OA, bx=r.bx, sx=r.sx, outcome=key, outcome_label=lab, scale=scale,
                             snp_beta_outcome=by, snp_se_outcome=sy, snp_p_outcome=pz(by / sy), wald=w, wald_se=wse, wald_p=pz(w / wse)))
    R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "S5_MR_robustness_all.tsv"), sep="\t", index=False)
    if R.empty: log("!! no results"); return
    # ---------------- Q1
    log(""); log("Q1  LTB instruments across outcomes  (Wald per SD protein; logOR or log eGFR)")
    ltb = R[R.Gene == "LTB"].copy()
    for (rs, src, pos), g in ltb.groupby(["rsid", "source", "pos37"], sort=False):
        log(f"  {rs} ({src}, chr6:{pos}, pQTL beta {g.bx.iloc[0]:+.3f})")
        for _, x in g.iterrows():
            eff = f"OR {math.exp(x.wald):.2f}" if x.scale == "logOR" else f"beta {x.wald:+.4f}"
            log(f"     {x.outcome:<11} {eff} (P={x.wald_p:.2g})  {x.outcome_label}")
    # IVW over LTB signals (lead + conditional signals; approximate - conditional pQTL effects with marginal outcome effects)
    log("  IVW across all LTB signals (approximate):")
    for key, g in ltb.groupby("outcome"):
        w = 1 / g.wald_se ** 2; b = np.sum(w * g.wald) / np.sum(w); se = 1 / math.sqrt(np.sum(w))
        log(f"     {key:<11} k={len(g)}  IVW {b:+.3f} (se {se:.3f}) P={pz(b/se):.2g}")
    # ---------------- Q2
    log(""); log("Q2  MHC lead pQTLs: effect on DN vs effect on T1D / T2D")
    mhc = R[(R.MHC == 1) & (R.source == "ST9_lead")]
    piv = mhc.pivot_table(index="Gene", columns="outcome", values="wald")
    pivp = mhc.pivot_table(index="Gene", columns="outcome", values="wald_p")
    tab = piv.join(pivp, rsuffix="_P"); tab.to_csv(os.path.join(OUT, "S5_MR_MHC_vs_T1D.tsv"), sep="\t")
    for oc in ("FG_T1D", "FG_T2D", "FG_T2D_REN", "EGFR_DM"):
        if "FG_DN" in piv and oc in piv:
            ok = piv[["FG_DN", oc]].dropna()
            if len(ok) >= 5:
                from scipy import stats
                rho = stats.spearmanr(ok.FG_DN, ok[oc]); sl = np.polyfit(ok[oc], ok.FG_DN, 1)
                log(f"  DN vs {oc}: {len(ok)} MHC proteins, Spearman rho={rho.correlation:.2f} (P={rho.pvalue:.2g}); slope={sl[0]:.2f}")
                if "LTB" in ok.index:
                    pred = np.polyval(sl, ok.loc["LTB", oc])
                    log(f"     LTB: DN {ok.loc['LTB','FG_DN']:+.3f} vs predicted from {oc} {pred:+.3f}")
    top = pivp.get("FG_DN", pd.Series(dtype=float)).sort_values().head(12).index
    cols = [c for c in ("FG_DN", "FG_T1D", "FG_T2D", "FG_T1D_REN", "FG_T2D_REN", "EGFR_DM") if c in piv]
    log("  top MHC proteins for DN (Wald / P):")
    for gname in top:
        log("     " + f"{gname:<10}" + "  ".join(f"{c}={piv.at[gname, c]:+.2f}(P={pivp.at[gname, c]:.1g})" if pd.notna(piv.at[gname, c]) else f"{c}=NA" for c in cols))
    # ---------------- Q3
    log(""); log("Q3  non-MHC pathway genes")
    for _, x in R[(R.MHC == 0) & R.Gene.isin(AXIS_NON_MHC + ["IL7R"])].sort_values(["Gene", "kind", "outcome"]).iterrows():
        eff = f"OR {math.exp(x.wald):.2f}" if x.scale == "logOR" else f"beta {x.wald:+.4f}"
        log(f"  {x.Gene:<8} {x.kind:<13} {x.rsid:<12} {x.outcome:<11} {eff} (P={x.wald_p:.2g})")
    log("genetic triangulation done", time.strftime("%Y-%m-%d %H:%M:%S"))

if __name__ == "__main__":
    try: main()
    except Exception as e:
        log("!! failed:", e); log(traceback.format_exc())
