"""
13_immune_checkpoints_tide.py — Immune checkpoint expression, TIDE dysfunction, ssGSEA

1. Compare immune checkpoint gene expression between high/low risk groups
2. Calculate TIDE-like immune dysfunction and exclusion scores
3. ssGSEA with 28 immune cell type gene sets
4. Correlate risk score with immunotherapy-relevant features
"""
import pandas as pd
import numpy as np
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
# Need FULL transcriptome (not just ROS genes) for checkpoint/immune analysis
expr_full = pd.read_csv(os.path.join(DATA, "tcga_lihc_expression_full.csv"), index_col=0)
expr_full = expr_full.T  # samples as rows, genes as columns
expr_full.index.name = "patient_id"

clinical = pd.read_csv(os.path.join(DATA, "tcga_lihc_clinical.csv"))
# Match patients: clinical patient_id to expression index (first 12 chars)
pid_col = "submitter_id" if "submitter_id" in clinical.columns else "patientId"
clinical["match_id"] = clinical[pid_col].str[:12]
expr_full["match_id"] = [x[:12] for x in expr_full.index]

# Select available clinical columns
clin_cols = ["match_id", "OS_months", "OS_event"]
for c in ["age_years", "age_at_diagnosis", "gender", "pathologic_stage", "tumor_stage"]:
    if c in clinical.columns:
        clin_cols.append(c)
merged_full = expr_full.merge(clinical[clin_cols].drop_duplicates(subset="match_id"),
                               on="match_id", how="inner")

with open(os.path.join(MODEL, "lasso_model.json")) as f:
    model = json.load(f)

selected_genes = model["genes"]
sig_genes = list(selected_genes.keys())

# Compute risk score
df = merged_full.dropna(subset=["OS_months", "OS_event"]).copy()
df = df[df["OS_months"] > 0]
risk = np.zeros(len(df))
for gene, coef in selected_genes.items():
    if gene in df.columns:
        risk += coef * df[gene].values
df["risk_score"] = risk
df["risk_group"] = (df["risk_score"] >= df["risk_score"].median()).map({True: "High", False: "Low"})

print(f"Patients: {len(df)} ({(df['risk_group']=='High').sum()} high, {(df['risk_group']=='Low').sum()} low)")

# ══════════════════════════════════════════════════════════════════════════════
# 1. IMMUNE CHECKPOINT EXPRESSION
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("1. IMMUNE CHECKPOINT EXPRESSION BY RISK GROUP")
print("=" * 70)

checkpoints = {
    "CD274": "PD-L1",
    "PDCD1": "PD-1",
    "CTLA4": "CTLA-4",
    "LAG3": "LAG-3",
    "HAVCR2": "TIM-3",
    "TIGIT": "TIGIT",
    "SIGLEC15": "SIGLEC15",
    "IDO1": "IDO1",
    "CD276": "B7-H3",
}

checkpoint_results = []
available_checkpoints = {}
for gene, alias in checkpoints.items():
    if gene in df.columns:
        available_checkpoints[gene] = alias
        high_vals = df.loc[df["risk_group"] == "High", gene].values
        low_vals = df.loc[df["risk_group"] == "Low", gene].values
        stat, pval = stats.mannwhitneyu(high_vals, low_vals, alternative='two-sided')
        fc = high_vals.mean() / (low_vals.mean() + 1e-10)
        direction = "Up in High-Risk" if fc > 1 else "Down in High-Risk"
        checkpoint_results.append({
            "gene": gene, "alias": alias,
            "high_mean": high_vals.mean(), "low_mean": low_vals.mean(),
            "fold_change": fc, "direction": direction,
            "p_value": pval
        })
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
        print(f"  {alias:10s} ({gene:10s}): High={high_vals.mean():.2f} vs Low={low_vals.mean():.2f}  FC={fc:.2f}  p={pval:.4f} {sig}")
    else:
        print(f"  {alias:10s} ({gene:10s}): not in expression data")

ckpt_df = pd.DataFrame(checkpoint_results)
ckpt_df.to_csv(os.path.join(TABLES, "checkpoint_expression.csv"), index=False)
print(f"\nSignificant checkpoints (p<0.05): {(ckpt_df['p_value']<0.05).sum()}/{len(ckpt_df)}")

# ── Checkpoint expression boxplot ──
if len(available_checkpoints) >= 3:
    fig, axes = plt.subplots(2, (len(available_checkpoints) + 1) // 2, figsize=(14, 8))
    axes = axes.flatten()
    for i, (gene, alias) in enumerate(available_checkpoints.items()):
        if i >= len(axes):
            break
        ax = axes[i]
        data_plot = df[["risk_group", gene]].copy()
        data_plot[gene] = np.log2(data_plot[gene] + 1)
        sns.boxplot(data=data_plot, x="risk_group", y=gene, ax=ax,
                    palette={"High": "#d62728", "Low": "#2ca02c"}, width=0.5,
                    order=["Low", "High"])
        pval = ckpt_df.loc[ckpt_df["gene"] == gene, "p_value"].values[0]
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
        ax.set_title(f"{alias}\n({sig}, p={pval:.3f})", fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("log2(TPM+1)")
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    plt.suptitle("Immune Checkpoint Expression by Risk Group", fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "checkpoint_expression.png"), dpi=200, bbox_inches='tight')
    print("Saved: checkpoint_expression.png")

# ══════════════════════════════════════════════════════════════════════════════
# 2. TIDE-LIKE IMMUNE DYSFUNCTION AND EXCLUSION SCORES
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("2. IMMUNE DYSFUNCTION AND EXCLUSION SCORING")
print("=" * 70)

# Dysfunction markers: high expression of checkpoints + Tregs + M2 macrophages
dysfunction_genes = ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "FOXP3", "IL10", "TGFB1"]
# Exclusion markers: stromal/EMT signals that exclude T cells
exclusion_genes = ["TWIST1", "TWIST2", "SNAI1", "SNAI2", "VIM", "CDH2", "ACTA2", "FAP", "TGFB1", "COL1A1"]
# Cytolytic activity (proxy for pre-existing anti-tumor immunity)
cytolytic_genes = ["GZMA", "GZMB", "PRF1", "GNLY", "NKG7", "IFNG"]

score_results = []
for name, gene_list in [("Dysfunction", dysfunction_genes), ("Exclusion", exclusion_genes), ("Cytolytic", cytolytic_genes)]:
    avail = [g for g in gene_list if g in df.columns]
    if len(avail) < 2:
        print(f"  {name}: insufficient genes ({len(avail)}/{len(gene_list)})")
        continue

    # Z-score each gene, average
    score = np.zeros(len(df))
    for g in avail:
        vals = df[g].values
        z = (vals - vals.mean()) / (vals.std() + 1e-10)
        score += z
    score /= len(avail)
    df[f"{name.lower()}_score"] = score

    high_score = score[df["risk_group"] == "High"]
    low_score = score[df["risk_group"] == "Low"]
    stat, pval = stats.mannwhitneyu(high_score, low_score, alternative='two-sided')
    corr, corr_p = stats.spearmanr(df["risk_score"], score)

    score_results.append({
        "score": name, "n_genes": len(avail),
        "high_mean": high_score.mean(), "low_mean": low_score.mean(),
        "mannwhitney_p": pval, "spearman_r": corr, "spearman_p": corr_p
    })
    sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
    print(f"  {name:12s}: {len(avail)} genes | High={high_score.mean():.3f} vs Low={low_score.mean():.3f} | p={pval:.4f} {sig} | r={corr:.3f}")

score_df = pd.DataFrame(score_results)
score_df.to_csv(os.path.join(TABLES, "immune_scores.csv"), index=False)

# ── Immune score comparison plot ──
score_cols = [c for c in ["dysfunction_score", "exclusion_score", "cytolytic_score"] if c in df.columns]
if len(score_cols) >= 2:
    fig, axes = plt.subplots(1, len(score_cols), figsize=(5 * len(score_cols), 5))
    if len(score_cols) == 1:
        axes = [axes]
    for i, col in enumerate(score_cols):
        ax = axes[i]
        sns.violinplot(data=df, x="risk_group", y=col, ax=ax,
                       palette={"High": "#d62728", "Low": "#2ca02c"},
                       order=["Low", "High"], inner="box", cut=0)
        pval = score_df.loc[score_df["score"] == col.replace("_score", "").capitalize(), "mannwhitney_p"].values
        if len(pval) > 0:
            sig = "***" if pval[0] < 0.001 else "**" if pval[0] < 0.01 else "*" if pval[0] < 0.05 else "ns"
            ax.set_title(f"{col.replace('_score', '').capitalize()} Score\n(p={pval[0]:.4f}, {sig})", fontsize=12, fontweight='bold')
        ax.set_xlabel("Risk Group")
        ax.set_ylabel("Score (z-score average)")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "immune_dysfunction_scores.png"), dpi=200, bbox_inches='tight')
    print("Saved: immune_dysfunction_scores.png")

# ══════════════════════════════════════════════════════════════════════════════
# 3. ssGSEA IMMUNE CELL PROFILING (28 cell types)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("3. ssGSEA IMMUNE CELL PROFILING")
print("=" * 70)

# Charoentong et al. immune cell markers (commonly used 28-type panel)
immune_signatures = {
    "aDC": ["CD1C", "ITGAX", "NRP1", "CD1E", "CCR7", "BATF3"],
    "B cells": ["CD19", "MS4A1", "CD79A", "CD79B", "BLK", "FCRL5"],
    "CD8 T cells": ["CD8A", "CD8B", "GZMB", "PRF1", "IFNG", "TBX21"],
    "Cytotoxic": ["GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY", "NKG7", "KLRK1"],
    "DC": ["HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA"],
    "Eosinophils": ["SIGLEC8", "CCR3", "IL5RA", "PRG2"],
    "iDC": ["CD1A", "CD1B", "CD1C", "CD1E", "CD209"],
    "Macrophages": ["CD68", "CD163", "MSR1", "MRC1", "CSF1R", "MARCO"],
    "Mast cells": ["TPSB2", "TPSAB1", "CPA3", "HDC", "KIT"],
    "Neutrophils": ["CEACAM8", "FPR1", "SIGLEC5", "CSF3R", "FCGR3B"],
    "NK cells": ["KLRC1", "KLRD1", "NCR1", "NCR3", "NKG7", "KLRF1"],
    "pDC": ["IL3RA", "CLEC4C", "NRP1", "LILRA4"],
    "T cells": ["CD3D", "CD3E", "CD3G", "CD2", "CD6"],
    "T helper": ["CD4", "IL2RA", "ICOS", "CD40LG"],
    "Tcm": ["CD27", "CCR7", "SELL", "IL7R", "LEF1"],
    "Tem": ["GZMK", "EOMES", "CXCR3"],
    "Tfh": ["CXCL13", "CXCR5", "BCL6", "PDCD1", "ICOS", "SH2D1A"],
    "Th1": ["TBX21", "IFNG", "TNF", "IL2", "STAT4"],
    "Th2": ["GATA3", "IL4", "IL5", "IL13", "STAT6"],
    "Th17": ["RORC", "IL17A", "IL17F", "IL22", "CCR6"],
    "TReg": ["FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2"],
    "Gamma_delta": ["TRDV1", "TRDV2", "TRGV9"],
    "M1 Macrophages": ["NOS2", "IL1B", "TNF", "IRF5", "CD80", "CD86"],
    "M2 Macrophages": ["CD163", "MRC1", "MSR1", "ARG1", "IL10", "TGFB1"],
}

# Compute ssGSEA-like score: rank-based enrichment
print(f"  Computing ssGSEA scores for {len(immune_signatures)} cell types...")

# Get all expression columns (genes)
all_genes = [c for c in df.columns if c not in ["patient_id", "OS_months", "OS_event",
             "risk_score", "risk_group", "sample_id", "age_years", "gender",
             "pathologic_stage", "histologic_grade", "race", "ethnicity",
             "dysfunction_score", "exclusion_score", "cytolytic_score",
             "match_id", "age_at_diagnosis", "tumor_stage"]]

ssgsea_results = []
ssgsea_scores = pd.DataFrame(index=df.index)

for cell_type, markers in immune_signatures.items():
    avail = [g for g in markers if g in df.columns]
    if len(avail) < 2:
        continue
    # Simple ssGSEA: z-score average of marker genes
    score = np.zeros(len(df))
    for g in avail:
        vals = df[g].values
        z = (vals - vals.mean()) / (vals.std() + 1e-10)
        score += z
    score /= len(avail)
    ssgsea_scores[cell_type] = score

    corr, pval = stats.spearmanr(df["risk_score"], score)
    ssgsea_results.append({
        "cell_type": cell_type, "n_markers": len(avail),
        "spearman_r": corr, "p_value": pval
    })

ssgsea_df = pd.DataFrame(ssgsea_results).sort_values("spearman_r")
ssgsea_df.to_csv(os.path.join(TABLES, "ssgsea_immune_profiles.csv"), index=False)
print(f"  Cell types with data: {len(ssgsea_df)}")
print(f"  Significantly correlated (p<0.05): {(ssgsea_df['p_value']<0.05).sum()}")

# ── ssGSEA correlation barplot ──
fig, ax = plt.subplots(figsize=(10, 8))
colors = ['#d62728' if r > 0 else '#2ca02c' for r in ssgsea_df["spearman_r"]]
sig_markers = ['*' if p < 0.05 else '' for p in ssgsea_df["p_value"]]
bars = ax.barh(range(len(ssgsea_df)), ssgsea_df["spearman_r"], color=colors, edgecolor='white', linewidth=0.5)
ax.set_yticks(range(len(ssgsea_df)))
ax.set_yticklabels([f"{ct} {sig}" for ct, sig in zip(ssgsea_df["cell_type"], sig_markers)], fontsize=9)
ax.set_xlabel("Spearman Correlation with Risk Score", fontsize=11)
ax.set_title("Immune Cell Infiltration vs Risk Score\n(ssGSEA, * p<0.05)", fontsize=13, fontweight='bold')
ax.axvline(0, color='black', linewidth=0.8)
ax.set_xlim(-0.4, 0.4)
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "ssgsea_immune_barplot.png"), dpi=200, bbox_inches='tight')
print("Saved: ssgsea_immune_barplot.png")

# ── Heatmap of immune scores by risk group ──
if len(ssgsea_scores.columns) > 3:
    fig, ax = plt.subplots(figsize=(14, 6))
    # Sort patients by risk score
    order = df["risk_score"].argsort().values
    plot_data = ssgsea_scores.iloc[order].T

    # Annotate risk group as color bar
    risk_colors = df.iloc[order]["risk_group"].map({"High": "#d62728", "Low": "#2ca02c"}).values

    sns.heatmap(plot_data, ax=ax, cmap="RdBu_r", center=0,
                xticklabels=False, yticklabels=True,
                cbar_kws={"label": "Enrichment score"})
    ax.set_title("Immune Cell Infiltration Landscape by Risk Score", fontsize=13, fontweight='bold')
    ax.set_xlabel(f"Patients (n={len(df)}, sorted by risk score →)")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "immune_landscape_heatmap.png"), dpi=200, bbox_inches='tight')
    print("Saved: immune_landscape_heatmap.png")

# ══════════════════════════════════════════════════════════════════════════════
# 4. IMMUNOTHERAPY RESPONSE MARKERS
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("4. IMMUNOTHERAPY RESPONSE MARKERS")
print("=" * 70)

# IFN-gamma signature (Ayers et al., JCI 2017) — predicts anti-PD-1 response
ifng_sig = ["IFNG", "STAT1", "CCR5", "CXCL9", "CXCL10", "CXCL11", "IDO1",
            "PRF1", "GZMA", "GZMB", "HLA-DRA", "HLA-DRB1", "HLA-E"]
# T-cell inflamed GEP (Ayers et al.)
tgep_sig = ["CD27", "CD274", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1",
            "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
            "PSMB10", "STAT1", "TIGIT"]
# Exclusion signature (Hugo et al., Cell 2016)
innate_anti_pd1 = ["AXL", "ROR2", "WNT5A", "TWIST2", "LOXL2", "FAP"]

marker_results = []
for name, genes in [("IFNg_signature", ifng_sig), ("T_cell_inflamed_GEP", tgep_sig),
                    ("Innate_anti_PD1", innate_anti_pd1)]:
    avail = [g for g in genes if g in df.columns]
    if len(avail) < 2:
        print(f"  {name}: insufficient genes ({len(avail)})")
        continue
    score = np.zeros(len(df))
    for g in avail:
        z = (df[g].values - df[g].values.mean()) / (df[g].values.std() + 1e-10)
        score += z
    score /= len(avail)

    corr, pval = stats.spearmanr(df["risk_score"], score)
    high_s = score[df["risk_group"] == "High"]
    low_s = score[df["risk_group"] == "Low"]
    _, mw_p = stats.mannwhitneyu(high_s, low_s, alternative='two-sided')

    marker_results.append({
        "signature": name, "n_genes": len(avail),
        "spearman_r": corr, "spearman_p": pval,
        "high_mean": high_s.mean(), "low_mean": low_s.mean(), "mannwhitney_p": mw_p
    })
    sig = "***" if mw_p < 0.001 else "**" if mw_p < 0.01 else "*" if mw_p < 0.05 else "ns"
    print(f"  {name:25s}: r={corr:.3f} (p={pval:.4f}) | High={high_s.mean():.3f} vs Low={low_s.mean():.3f} {sig}")

marker_df = pd.DataFrame(marker_results)
marker_df.to_csv(os.path.join(TABLES, "immunotherapy_markers.csv"), index=False)

# ══════════════════════════════════════════════════════════════════════════════
# 5. SUMMARY FIGURE: COMPREHENSIVE IMMUNE LANDSCAPE
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("5. SUMMARY")
print("=" * 70)

n_sig_checkpoints = (ckpt_df['p_value'] < 0.05).sum() if len(ckpt_df) > 0 else 0
n_sig_ssgsea = (ssgsea_df['p_value'] < 0.05).sum() if len(ssgsea_df) > 0 else 0
print(f"  Checkpoints differentially expressed: {n_sig_checkpoints}/{len(ckpt_df)}")
print(f"  Immune cell types correlated with risk: {n_sig_ssgsea}/{len(ssgsea_df)}")
if len(score_df) > 0:
    for _, row in score_df.iterrows():
        print(f"  {row['score']} score: r={row['spearman_r']:.3f}, p={row['spearman_p']:.4f}")

print("\n✓ Immune checkpoint and TIDE analysis complete.")
