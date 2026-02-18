"""
14_gsea_risk_groups.py — Genome-wide GSEA between high and low risk groups

1. Rank ALL genes by Spearman correlation with risk score
2. Run preranked GSEA against Hallmark, KEGG, Reactome, GO BP
3. Visualize top enriched pathways
"""
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import json
import os
import warnings
warnings.filterwarnings('ignore')

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data", "tcga")
MODEL = os.path.join(BASE, "results", "model")
TABLES = os.path.join(BASE, "results", "tables")
FIGS = os.path.join(BASE, "results", "figures")

# ── Load data ───────────────────────────────────────────────────────────────
# Need FULL transcriptome for genome-wide GSEA
expr_full = pd.read_csv(os.path.join(DATA, "tcga_lihc_expression_full.csv"), index_col=0)
expr_full = expr_full.T  # samples as rows, genes as columns
expr_full.index.name = "patient_id"

clinical = pd.read_csv(os.path.join(DATA, "tcga_lihc_clinical.csv"))
pid_col = "submitter_id" if "submitter_id" in clinical.columns else "patientId"
clinical["match_id"] = clinical[pid_col].str[:12]
expr_full["match_id"] = [x[:12] for x in expr_full.index]

merged = expr_full.merge(clinical[["match_id", "OS_months", "OS_event"]].drop_duplicates(subset="match_id"),
                          on="match_id", how="inner")

with open(os.path.join(MODEL, "lasso_model.json")) as f:
    model = json.load(f)

selected_genes = model["genes"]
df = merged.dropna(subset=["OS_months", "OS_event"]).copy()
df = df[df["OS_months"] > 0]

# Compute risk score
risk = np.zeros(len(df))
for gene, coef in selected_genes.items():
    if gene in df.columns:
        risk += coef * df[gene].values
df["risk_score"] = risk

print(f"Patients: {len(df)}")

# ══════════════════════════════════════════════════════════════════════════════
# 1. RANK ALL GENES BY CORRELATION WITH RISK SCORE
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("1. RANKING GENES BY CORRELATION WITH RISK SCORE")
print("=" * 70)

# Get all gene columns (exclude clinical)
clinical_cols = {"patient_id", "sample_id", "OS_months", "OS_event", "age_years",
                 "gender", "pathologic_stage", "histologic_grade", "race", "ethnicity",
                 "risk_score", "match_id", "age_at_diagnosis", "tumor_stage"}
gene_cols = [c for c in df.columns if c not in clinical_cols]

# Compute Spearman correlation for each gene
gene_ranks = []
for gene in gene_cols:
    try:
        vals = df[gene].values.astype(float)
        if vals.std() < 1e-10:
            continue
        corr, pval = stats.spearmanr(df["risk_score"], vals)
        if np.isnan(corr):
            continue
        gene_ranks.append({"gene": gene, "correlation": corr, "p_value": pval})
    except (ValueError, TypeError):
        continue

rank_df = pd.DataFrame(gene_ranks).sort_values("correlation", ascending=False)
print(f"  Genes ranked: {len(rank_df)}")
print(f"  Top positively correlated: {rank_df.head(5)['gene'].tolist()}")
print(f"  Top negatively correlated: {rank_df.tail(5)['gene'].tolist()}")

# Create ranked gene list for GSEA (gene -> correlation)
rnk = rank_df.set_index("gene")["correlation"]

# ══════════════════════════════════════════════════════════════════════════════
# 2. PRERANKED GSEA
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("2. PRERANKED GSEA")
print("=" * 70)

gene_sets = [
    "MSigDB_Hallmark_2020",
    "KEGG_2021_Human",
    "Reactome_2022",
    "GO_Biological_Process_2023",
]

all_gsea = []
for gs_name in gene_sets:
    print(f"\n  Running GSEA: {gs_name}...")
    try:
        res = gp.prerank(
            rnk=rnk,
            gene_sets=gs_name,
            threads=4,
            min_size=15,
            max_size=500,
            permutation_num=1000,
            seed=42,
            no_plot=True,
        )
        result = res.res2d.copy()
        result["library"] = gs_name
        all_gsea.append(result)
        n_sig = (result["FDR q-val"].astype(float) < 0.25).sum()
        print(f"    Significant (FDR<0.25): {n_sig}/{len(result)}")
    except Exception as e:
        print(f"    Failed: {e}")

if all_gsea:
    gsea_combined = pd.concat(all_gsea, ignore_index=True)
    gsea_combined.to_csv(os.path.join(TABLES, "gsea_preranked_results.csv"), index=False)
    print(f"\n  Total pathways tested: {len(gsea_combined)}")
    print(f"  Significant (FDR<0.25): {(gsea_combined['FDR q-val'].astype(float) < 0.25).sum()}")

# ══════════════════════════════════════════════════════════════════════════════
# 3. VISUALIZATION
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("3. GSEA VISUALIZATION")
print("=" * 70)

if all_gsea:
    gsea_combined["NES"] = gsea_combined["NES"].astype(float)
    gsea_combined["FDR q-val"] = gsea_combined["FDR q-val"].astype(float)

    # ── Hallmark dotplot ──
    hallmark = gsea_combined[gsea_combined["library"] == "MSigDB_Hallmark_2020"].copy()
    if len(hallmark) > 0:
        hallmark = hallmark.sort_values("NES", ascending=True)
        # Take top 10 positive + top 10 negative NES
        top_pos = hallmark.tail(15)
        top_neg = hallmark.head(15)
        top = pd.concat([top_neg, top_pos]).drop_duplicates(subset=["Term"])
        top = top.sort_values("NES")

        fig, ax = plt.subplots(figsize=(10, 10))
        colors = ['#d62728' if nes > 0 else '#2166ac' for nes in top["NES"]]
        sizes = [max(20, min(200, -np.log10(fdr + 1e-10) * 30)) for fdr in top["FDR q-val"]]
        ax.scatter(top["NES"], range(len(top)), c=colors, s=sizes, edgecolors='black', linewidth=0.5, alpha=0.8)
        ax.set_yticks(range(len(top)))
        labels = [t.replace("_", " ").title() if len(t) < 50 else t[:47] + "..." for t in top["Term"]]
        ax.set_yticklabels(labels, fontsize=8)
        ax.axvline(0, color='grey', linestyle='--', linewidth=0.8)
        ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=11)
        ax.set_title("Hallmark Pathway Enrichment (GSEA Preranked)\nRed=enriched in high-risk, Blue=enriched in low-risk",
                      fontsize=12, fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(FIGS, "gsea_hallmark_dotplot.png"), dpi=200, bbox_inches='tight')
        print("Saved: gsea_hallmark_dotplot.png")

    # ── Combined top pathways across all libraries ──
    sig_paths = gsea_combined[gsea_combined["FDR q-val"] < 0.25].copy()
    if len(sig_paths) > 0:
        # Top 10 positive + top 10 negative across all libraries
        top_pos = sig_paths.nlargest(10, "NES")
        top_neg = sig_paths.nsmallest(10, "NES")
        top_all = pd.concat([top_neg, top_pos]).drop_duplicates().sort_values("NES")

        fig, ax = plt.subplots(figsize=(12, 10))
        colors = ['#d62728' if nes > 0 else '#2166ac' for nes in top_all["NES"]]
        bars = ax.barh(range(len(top_all)), top_all["NES"], color=colors, edgecolor='white', linewidth=0.5)
        ax.set_yticks(range(len(top_all)))
        labels = []
        for _, row in top_all.iterrows():
            term = row["Term"]
            lib = row["library"].split("_")[0]
            label = f"[{lib}] {term}"
            if len(label) > 60:
                label = label[:57] + "..."
            labels.append(label)
        ax.set_yticklabels(labels, fontsize=8)
        ax.axvline(0, color='black', linewidth=0.8)
        ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=11)
        ax.set_title("Top Enriched Pathways (All Libraries, FDR<0.25)\nRed=high-risk enriched, Blue=low-risk enriched",
                      fontsize=12, fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(FIGS, "gsea_top_pathways.png"), dpi=200, bbox_inches='tight')
        print("Saved: gsea_top_pathways.png")

    # ── KEGG pathways ──
    kegg = gsea_combined[gsea_combined["library"] == "KEGG_2021_Human"].copy()
    if len(kegg) > 0:
        kegg_sig = kegg[kegg["FDR q-val"] < 0.25].sort_values("NES")
        if len(kegg_sig) > 0:
            top_kegg = pd.concat([kegg_sig.head(10), kegg_sig.tail(10)]).drop_duplicates().sort_values("NES")
            fig, ax = plt.subplots(figsize=(10, 8))
            colors = ['#d62728' if nes > 0 else '#2166ac' for nes in top_kegg["NES"]]
            ax.barh(range(len(top_kegg)), top_kegg["NES"], color=colors, edgecolor='white')
            ax.set_yticks(range(len(top_kegg)))
            ax.set_yticklabels([t[:50] for t in top_kegg["Term"]], fontsize=9)
            ax.axvline(0, color='black', linewidth=0.8)
            ax.set_xlabel("Normalized Enrichment Score (NES)")
            ax.set_title("KEGG Pathway Enrichment (FDR<0.25)", fontsize=12, fontweight='bold')
            plt.tight_layout()
            plt.savefig(os.path.join(FIGS, "gsea_kegg_barplot.png"), dpi=200, bbox_inches='tight')
            print("Saved: gsea_kegg_barplot.png")

# ══════════════════════════════════════════════════════════════════════════════
# 4. SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("4. SUMMARY")
print("=" * 70)

if all_gsea:
    for gs_name in gene_sets:
        subset = gsea_combined[gsea_combined["library"] == gs_name]
        n_sig = (subset["FDR q-val"].astype(float) < 0.25).sum()
        top_pos_term = subset.nlargest(1, "NES")["Term"].values[0] if len(subset) > 0 else "N/A"
        top_neg_term = subset.nsmallest(1, "NES")["Term"].values[0] if len(subset) > 0 else "N/A"
        print(f"  {gs_name}: {n_sig} significant | Top+: {top_pos_term[:40]} | Top-: {top_neg_term[:40]}")

print("\n✓ Genome-wide GSEA complete.")
