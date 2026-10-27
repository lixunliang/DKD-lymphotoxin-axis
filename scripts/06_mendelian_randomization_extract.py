# -*- coding: utf-8 -*-
"""Extract MR estimates for progression / lymphotoxin-pathway proteins from proteome-wide cis-pQTL MR results.
Inputs (from the authors' proteome-wide MR pipeline; Wald ratio / IVW, Benjamini-Hochberg FDR per platform and outcome):
  kidney_disease_protein_mr.tsv  (deCODE SomaScan; columns outcome, SeqId, Gene, beta, se, p, FDR, MHC)
  olink_kidney_disease_mr.tsv    (UKB-PPP Olink;   columns outcome, ProteinID, Gene, beta, se, p, FDR, MHC)
Outcomes include DKD_FG = FinnGen R12 DM_NEPHROPATHY."""
import os, sys, numpy as np, pandas as pd
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
MRDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(BASE, "mr_input")
R = os.path.join(BASE, "results")
cand = pd.read_csv(os.path.join(R, "progression_genes_tiered.csv"), index_col=0)
G = set(cand.index) | {"CD44", "RELB", "NFKB2", "CCL19", "CCL21", "LTB", "LTA", "LTBR", "CXCL13", "TNFRSF1B", "TIMP1", "CD163", "C3", "THBS1", "EGF", "IL7R", "MAP3K14"}
a = pd.read_csv(os.path.join(MRDIR, "kidney_disease_protein_mr.tsv"), sep="\t"); a["platform"] = "deCODE (SomaScan)"; a = a.rename(columns={"SeqId": "protein_id"})
b = pd.read_csv(os.path.join(MRDIR, "olink_kidney_disease_mr.tsv"), sep="\t"); b["platform"] = "UKB-PPP (Olink)"; b = b.rename(columns={"ProteinID": "protein_id"})
cols = ["platform", "outcome", "Gene", "protein_id", "beta", "se", "p", "FDR", "MHC"]
x = pd.concat([a[cols], b[cols]]); x = x[x.Gene.isin(G)].copy()
x["OR"] = np.exp(x.beta); x["OR_L"] = np.exp(x.beta - 1.96 * x.se); x["OR_U"] = np.exp(x.beta + 1.96 * x.se)
x["in_tier"] = x.Gene.map(cand.tier).fillna("pathway gene")
x.sort_values(["outcome", "p"]).to_csv(os.path.join(R, "S4_MR_candidates_kidney_disease.csv"), index=False, encoding="utf-8-sig")
print(x[x.outcome == "DKD_FG"].sort_values("p").head(10)[["platform", "Gene", "OR", "OR_L", "OR_U", "p", "FDR", "MHC"]])
