# -*- coding: utf-8 -*-
"""Define Tier-1 / Tier-2 progression genes from the outputs of 01_bulk_meta_analysis.py
Tier 1: meta FDR<0.05 & |g|>0.8 & full direction consistency in BOTH glomeruli and tubulointerstitium,
        advanced-vs-early FDR<0.01 & |log2FC|>=1 (same direction), advanced-vs-control FDR<0.05 (same direction),
        and no opposite change in early DKD vs control (FDR<0.05).
Tier 2: same criteria in ONE compartment."""
import os, numpy as np, pandas as pd
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
a = pd.read_csv(os.path.join(R, "ALL_genes_integrated.csv"), index_col=0)
e = pd.read_csv(os.path.join(R, "DE_GSE142025_EARLYvsCTRL.csv"), index_col=0)
a["E_lfc"], a["E_FDR"] = e.logFC, e.FDR
d = np.sign(a.ADVvsEARLY_logFC)
gs = (a.glom_meta_FDR < 0.05) & (a.glom_meta_g.abs() > 0.8) & (np.sign(a.glom_meta_g) == d) & (a.glom_dir == 1)
ts = (a.tub_meta_FDR < 0.05) & (a.tub_meta_g.abs() > 0.8) & (np.sign(a.tub_meta_g) == d) & (a.tub_dir == 1)
prog = (a.ADVvsEARLY_FDR < 0.01) & (a.ADVvsEARLY_logFC.abs() >= 1)
adv = (a.ADVvsCTRL_FDR < 0.05) & (np.sign(a.ADVvsCTRL_logFC) == d)
notopp = ~((a.E_FDR < 0.05) & (np.sign(a.E_lfc) == -d))
t1 = a[gs & ts & prog & adv & notopp].copy(); t1["tier"] = "Tier1_glom+tub"
t2 = a[(gs ^ ts) & prog & adv & notopp].copy(); t2["tier"] = "Tier2_single_compartment"
cols = ["tier", "glom_meta_g", "glom_meta_FDR", "tub_meta_g", "tub_meta_FDR", "E_lfc", "E_FDR", "ADVvsEARLY_logFC", "ADVvsEARLY_FDR",
        "ADVvsCTRL_logFC", "ADVvsCTRL_FDR", "secreted", "Protein class", "Blood concentration - Conc. blood MS [pg/L]"]
both = pd.concat([t1, t2]); out = both[[c for c in cols if c in both.columns]].rename(columns={"E_lfc": "EARLYvsCTRL_logFC", "E_FDR": "EARLYvsCTRL_FDR"})
out["blood_secreted"] = out.get("secreted", pd.Series("", index=out.index)).fillna("").str.contains("blood", case=False)
out.index.name = "gene"
out.to_csv(os.path.join(R, "progression_genes_tiered.csv"), encoding="utf-8-sig")
print(out.tier.value_counts())
