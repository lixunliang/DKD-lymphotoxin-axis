# -*- coding: utf-8 -*-
# Pathway enrichment (g:Profiler) for tiered DKD progression genes
import os, sys, json, ssl, urllib.request
import pandas as pd
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(BASE, "results")
try:
    import certifi; CTX = ssl.create_default_context(cafile=certifi.where())
except Exception:
    CTX = ssl.create_default_context()
d = pd.read_csv(os.path.join(OUT, "progression_genes_tiered.csv"), index_col=0)
bg = pd.read_csv(os.path.join(OUT, "ALL_genes_integrated.csv"), index_col=0)
bg = bg.index[bg["ADVvsEARLY_P"].notna()].tolist()
def run(genes, label):
    body = json.dumps({"organism": "hsapiens", "query": list(genes), "sources": ["GO:BP", "KEGG", "REAC", "WP"],
                       "user_threshold": 0.05, "domain_scope": "custom", "background": bg, "no_evidences": False}).encode()
    req = urllib.request.Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/", data=body,
                                 headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"})
    js = json.loads(urllib.request.urlopen(req, timeout=300, context=CTX).read().decode())
    names = js.get("meta", {}).get("genes_metadata", {}).get("query", {}).get("query_1", {}).get("ensgs", [])
    rows = []
    for x in js.get("result", []):
        rows.append({"source": x["source"], "term_id": x["native"], "term": x["name"], "P_adj": x["p_value"],
                     "term_size": x["term_size"], "overlap": x["intersection_size"], "query_size": x["query_size"]})
    df = pd.DataFrame(rows).sort_values("P_adj") if rows else pd.DataFrame()
    df.to_csv(os.path.join(OUT, f"ENRICH_{label}.csv"), index=False, encoding="utf-8-sig")
    print(label, len(genes), "genes ->", len(df), "terms")
    print(df.head(15)[["source", "term", "P_adj", "overlap"]].to_string() if len(df) else "")
up = d.index[d["ADVvsEARLY_logFC"] > 0]; dn = d.index[d["ADVvsEARLY_logFC"] < 0]
run(up, "progression_up"); run(dn, "progression_down")
print("done")
