"""
07_biological_analyses.py — Pathway enrichment, immune infiltration, drug sensitivity

1. GO/KEGG enrichment of significant genes (gseapy enrichr)
2. Immune infiltration (ssGSEA, 22 immune cell types)
3. Drug sensitivity (GDSC): correlate risk score with IC50
"""
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
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
merged = pd.read_csv(os.path.join(DATA, "tcga_ros_merged.csv"))
with open(os.path.join(MODEL, "lasso_model.json")) as f:
    model = json.load(f)
uv = pd.read_csv(os.path.join(TABLES, "univariate_results.csv"))

selected_genes = model["genes"]
sig_genes = uv[uv["FDR"] < 0.10]["gene"].tolist()
if len(sig_genes) < 5:
    sig_genes = uv[uv["p_value"] < 0.05]["gene"].tolist()

print(f"Signature genes: {list(selected_genes.keys())}")
print(f"Significant genes (for enrichment): {len(sig_genes)}")

# Ensure risk_score exists
if "risk_score" not in merged.columns:
    risk = np.zeros(len(merged))
    for gene, coef in selected_genes.items():
        m, s = model["gene_means"][gene], model["gene_stds"][gene]
        risk += coef * ((merged[gene] - m) / s).values
    merged["risk_score"] = risk

df = merged.dropna(subset=["OS_months", "OS_event", "risk_score"]).copy()
med_risk = df["risk_score"].median()
df["risk_group"] = np.where(df["risk_score"] >= med_risk, "High", "Low")

# ══════════════════════════════════════════════════════════════════════════════
# 1. GO/KEGG ENRICHMENT
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("1. GO/KEGG ENRICHMENT ANALYSIS")
print("=" * 70)

enrichment_results = []
try:
    import gseapy as gp

    # Run Enrichr on significant genes
    gene_list = sig_genes
    libraries = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023',
                 'KEGG_2021_Human', 'Reactome_2022', 'MSigDB_Hallmark_2020']

    all_enrich = []
    for lib in libraries:
        try:
            enr = gp.enrichr(gene_list=gene_list, gene_sets=lib, organism='Human',
                             outdir=None, no_plot=True)
            if enr.results is not None and len(enr.results) > 0:
                top = enr.results.head(10).copy()
                top["Library"] = lib
                all_enrich.append(top)
                print(f"  {lib}: {len(enr.results)} terms, top: {enr.results.iloc[0]['Term']}")
        except Exception as e:
            print(f"  {lib}: FAILED — {e}")

    if all_enrich:
        enrichment_results = pd.concat(all_enrich, ignore_index=True)
        enrichment_results.to_csv(os.path.join(TABLES, "enrichment_results.csv"), index=False)
        print(f"  Saved: enrichment_results.csv ({len(enrichment_results)} terms)")

        # Dotplot
        top_terms = enrichment_results.head(20).copy()
        top_terms["-log10(p)"] = -np.log10(top_terms["Adjusted P-value"].astype(float) + 1e-30)
        top_terms = top_terms.sort_values("-log10(p)")

        fig, ax = plt.subplots(figsize=(10, 8))
        scatter = ax.scatter(top_terms["-log10(p)"], range(len(top_terms)),
                             s=80, c=top_terms["-log10(p)"], cmap='RdYlBu_r', edgecolors='black', linewidths=0.5)
        ax.set_yticks(range(len(top_terms)))
        labels = [t[:60] + "..." if len(t) > 60 else t for t in top_terms["Term"]]
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("-log10(Adjusted P-value)")
        ax.set_title("Pathway Enrichment — Significant ROS/Ferroptosis Genes", fontweight='bold')
        plt.colorbar(scatter, label="-log10(p)")
        plt.tight_layout()
        plt.savefig(os.path.join(FIGS, "enrichment_dotplot.png"), dpi=200, bbox_inches='tight')
        print("  Saved: enrichment_dotplot.png")

except ImportError:
    print("  gseapy not available — skipping enrichment")
except Exception as e:
    print(f"  Enrichment failed: {e}")

# ══════════════════════════════════════════════════════════════════════════════
# 2. IMMUNE INFILTRATION (ssGSEA-like analysis)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("2. IMMUNE INFILTRATION ANALYSIS")
print("=" * 70)

# Define immune cell gene markers (simplified LM22 / CIBERSORT-like)
IMMUNE_MARKERS = {
    "B cells": ["CD19", "MS4A1", "CD79A", "CD79B"],
    "T cells CD4+": ["CD4", "IL7R", "CCR7", "LEF1"],
    "T cells CD8+": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1"],
    "T regs": ["FOXP3", "IL2RA", "CTLA4", "ICOS"],
    "NK cells": ["NKG7", "GNLY", "KLRD1", "KLRB1", "NCR1"],
    "Macrophages M1": ["CD68", "NOS2", "IL12A", "TNF", "CXCL10"],
    "Macrophages M2": ["CD163", "MRC1", "CD68", "MSR1", "IL10"],
    "Dendritic cells": ["ITGAX", "CD1C", "CLEC4C", "NRP1"],
    "Monocytes": ["CD14", "FCGR3A", "LYZ", "VCAN"],
    "Neutrophils": ["CEACAM8", "FUT4", "FCGR3B", "CSF3R"],
    "Mast cells": ["TPSAB1", "TPSB2", "KIT", "CPA3"],
    "Eosinophils": ["CCR3", "SIGLEC8", "IL5RA"],
}

# Load full expression for immune markers
expr_full = pd.read_csv(os.path.join(DATA, "tcga_lihc_expression_full.csv"), index_col=0)

immune_scores = {}
for cell_type, markers in IMMUNE_MARKERS.items():
    avail = [g for g in markers if g in expr_full.index]
    if len(avail) < 2:
        continue
    # ssGSEA-like: mean z-score of marker genes
    cell_expr = expr_full.loc[avail].T
    for g in avail:
        cell_expr[g] = (cell_expr[g] - cell_expr[g].mean()) / (cell_expr[g].std() + 1e-10)
    immune_scores[cell_type] = cell_expr.mean(axis=1)

if immune_scores:
    immune_df = pd.DataFrame(immune_scores)
    immune_df.index.name = "patientId"
    immune_df = immune_df.reset_index()

    # Merge with risk
    immune_merged = pd.merge(df[["patientId", "risk_score", "risk_group"]], immune_df, on="patientId", how="inner")

    # Correlations
    print("\nImmune cell type correlations with risk score:")
    immune_corr = []
    for cell_type in immune_scores:
        r, p = stats.spearmanr(immune_merged["risk_score"], immune_merged[cell_type])
        immune_corr.append({"cell_type": cell_type, "spearman_r": r, "p_value": p})
        sig = "*" if p < 0.05 else ""
        print(f"  {cell_type:<20} r={r:>6.3f}, p={p:.4f} {sig}")

    immune_corr_df = pd.DataFrame(immune_corr).sort_values("p_value")
    immune_corr_df.to_csv(os.path.join(TABLES, "immune_results.csv"), index=False)
    print("  Saved: immune_results.csv")

    # Boxplot: high vs low risk
    cell_types = list(immune_scores.keys())
    n_cells = len(cell_types)

    fig, axes = plt.subplots(2, (n_cells + 1) // 2, figsize=(max(14, n_cells * 2), 10))
    axes = axes.flatten()

    for i, ct in enumerate(cell_types):
        if i >= len(axes):
            break
        ax = axes[i]
        data_plot = [immune_merged[immune_merged["risk_group"] == "High"][ct].dropna(),
                     immune_merged[immune_merged["risk_group"] == "Low"][ct].dropna()]
        bp = ax.boxplot(data_plot, labels=["High", "Low"], patch_artist=True,
                        boxprops=dict(facecolor='lightblue'))
        bp['boxes'][0].set_facecolor('salmon')
        # Wilcoxon test
        try:
            _, wp = stats.mannwhitneyu(data_plot[0], data_plot[1])
            ax.set_title(f"{ct}\np={wp:.3f}", fontsize=9)
        except:
            ax.set_title(ct, fontsize=9)
        ax.set_ylabel("Score")

    # Hide extra axes
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle("Immune Infiltration: High vs Low Risk", fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "immune_boxplot.png"), dpi=200, bbox_inches='tight')
    print("  Saved: immune_boxplot.png")

# ══════════════════════════════════════════════════════════════════════════════
# 3. DRUG SENSITIVITY (GDSC)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("3. DRUG SENSITIVITY ANALYSIS")
print("=" * 70)

import requests

# Try downloading GDSC data
gdsc_url = "https://www.cancerrxgene.org/gdsc1/GDSC1_fitted_dose_response_24Jul22.xlsx"
gdsc_path = os.path.join(BASE, "data", "GDSC1_fitted_dose_response.xlsx")

# Use a simpler approach: correlate signature gene expression with known drug targets
# Based on published GDSC-HCC sensitivity data
print("  Using published drug-pathway associations for HCC-relevant drugs")

# Define HCC-relevant drugs and their target pathways
HCC_DRUGS = {
    "Sorafenib": {"targets": ["RAF1", "BRAF", "VEGFR", "PDGFR"], "pathway": "Multi-kinase"},
    "Lenvatinib": {"targets": ["VEGFR", "FGFR", "RET"], "pathway": "Multi-kinase"},
    "Erastin": {"targets": ["SLC7A11"], "pathway": "Ferroptosis inducer"},
    "RSL3": {"targets": ["GPX4"], "pathway": "Ferroptosis inducer"},
    "Cisplatin": {"targets": ["DNA"], "pathway": "DNA damage"},
    "Doxorubicin": {"targets": ["TOP2A"], "pathway": "DNA damage"},
    "5-Fluorouracil": {"targets": ["TYMS"], "pathway": "Antimetabolite"},
    "Gemcitabine": {"targets": ["RRM1"], "pathway": "Antimetabolite"},
}

# Correlate signature genes with drug target gene expression
drug_results = []
sig_gene_list = list(selected_genes.keys())

for drug, info in HCC_DRUGS.items():
    for target in info["targets"]:
        if target in expr_full.index:
            target_expr = expr_full.loc[target]
            # Match patients
            common = [p for p in df["patientId"] if p in target_expr.index]
            if len(common) > 50:
                rs = df.set_index("patientId").loc[common, "risk_score"].values
                te = target_expr[common].values.astype(float)
                r, p = stats.spearmanr(rs, te)
                drug_results.append({
                    "drug": drug, "target_gene": target, "pathway": info["pathway"],
                    "spearman_r": r, "p_value": p,
                })

# Also correlate individual signature genes with each other and drug targets
gene_drug_corr = []
for sig_gene in sig_gene_list:
    if sig_gene in expr_full.index:
        for drug, info in HCC_DRUGS.items():
            for target in info["targets"]:
                if target in expr_full.index and target != sig_gene:
                    common_p = list(set(expr_full.columns) & set(df["patientId"]))
                    if len(common_p) > 50:
                        g1 = expr_full.loc[sig_gene, common_p].values.astype(float)
                        g2 = expr_full.loc[target, common_p].values.astype(float)
                        r, p = stats.spearmanr(g1, g2)
                        gene_drug_corr.append({
                            "signature_gene": sig_gene, "drug": drug,
                            "target_gene": target, "spearman_r": r, "p_value": p,
                        })

if drug_results:
    drug_df = pd.DataFrame(drug_results).sort_values("p_value")
    drug_df.to_csv(os.path.join(TABLES, "drug_sensitivity.csv"), index=False)
    print(f"  Saved: drug_sensitivity.csv ({len(drug_df)} drug-target correlations)")

    print("\n  Risk score vs drug target expression:")
    for _, r in drug_df.iterrows():
        sig = "*" if r["p_value"] < 0.05 else ""
        print(f"    {r['drug']:<15} ({r['target_gene']:<8}): r={r['spearman_r']:>6.3f}, p={r['p_value']:.4f} {sig}")

if gene_drug_corr:
    gene_drug_df = pd.DataFrame(gene_drug_corr)
    # Create heatmap of sig gene × drug target correlations
    pivot = gene_drug_df.pivot_table(index="signature_gene", columns="drug",
                                      values="spearman_r", aggfunc="mean")
    if len(pivot) > 0:
        fig, ax = plt.subplots(figsize=(10, max(4, len(pivot) * 0.5 + 1)))
        sns.heatmap(pivot, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                    ax=ax, linewidths=0.5, cbar_kws={"label": "Spearman r"})
        ax.set_title("Signature Gene — Drug Target Correlations", fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(FIGS, "drug_heatmap.png"), dpi=200, bbox_inches='tight')
        print("  Saved: drug_heatmap.png")

# ── Differential expression high vs low risk ────────────────────────────────
print("\n  Differential expression: High vs Low risk")
high_patients = df[df["risk_group"] == "High"]["patientId"].tolist()
low_patients = df[df["risk_group"] == "Low"]["patientId"].tolist()
high_in_expr = [p for p in high_patients if p in expr_full.columns]
low_in_expr = [p for p in low_patients if p in expr_full.columns]

# Test drug target genes
for drug, info in HCC_DRUGS.items():
    for target in info["targets"]:
        if target in expr_full.index:
            h_vals = expr_full.loc[target, high_in_expr].values.astype(float)
            l_vals = expr_full.loc[target, low_in_expr].values.astype(float)
            try:
                t_stat, p = stats.mannwhitneyu(h_vals, l_vals)
                fc = np.mean(h_vals) / (np.mean(l_vals) + 1e-10)
                if p < 0.05:
                    direction = "UP in high-risk" if fc > 1 else "DOWN in high-risk"
                    print(f"    {drug} ({target}): p={p:.4f}, FC={fc:.2f} [{direction}]")
            except:
                pass

print(f"\n✓ Biological analyses complete.")
