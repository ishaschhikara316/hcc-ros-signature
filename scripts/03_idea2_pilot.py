"""
IDEA 2 PILOT: Oxidative Stress Gene Signature in HCC

Quick check:
- Run univariate Cox regression on oxidative stress / ROS pathway genes
- Test if oxidative-stress-high vs -low subgroups differ in survival
- Assess signal strength for a paper
"""
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
import os
import warnings
warnings.filterwarnings('ignore')

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Load merged data
merged = pd.read_csv(os.path.join(OUT_DIR, "tcga_lihc_merged.csv"))
print(f"Loaded merged data: {merged.shape}")

# Oxidative stress / ROS / Nrf2-KEAP1 pathway genes
oxidative_genes = [
    # Nrf2-KEAP1 pathway (master regulators)
    "NFE2L2", "KEAP1", "MAFG",
    # Superoxide dismutases
    "SOD1", "SOD2", "SOD3",
    # Catalase
    "CAT",
    # Glutathione peroxidases
    "GPX1", "GPX2", "GPX3", "GPX4",
    # Glutathione synthesis & recycling
    "GSR", "GCLC", "GCLM",
    # Thioredoxin system
    "TXNRD1", "TXN",
    # Peroxiredoxins
    "PRDX1", "PRDX2", "PRDX3",
    # Nrf2 target genes
    "NQO1", "HMOX1", "HMOX2", "FTH1", "FTL", "SLC7A11", "SRXN1", "SQSTM1",
    # Pentose phosphate pathway (NADPH generation)
    "G6PD", "PGD", "ME1",
    # Glutaredoxins
    "GLRX", "GLRX2", "MSRA",
    # NADPH oxidases (ROS producers)
    "NOX1", "NOX4", "CYBB", "NCF1", "NCF2", "RAC1",
    # Eicosanoid pathway (oxidative)
    "PTGS2", "ALOX5", "ALOX12", "ALOX15",
]

available = [g for g in oxidative_genes if g in merged.columns]
print(f"\nOxidative stress genes available: {len(available)}/{len(oxidative_genes)}")

# ============================================================
# UNIVARIATE COX REGRESSION ON EACH GENE
# ============================================================
print("\n" + "=" * 70)
print("UNIVARIATE COX REGRESSION — OXIDATIVE STRESS GENES")
print("=" * 70)

results = []
cph = CoxPHFitter()

for gene in available:
    df_gene = merged[["OS_months", "OS_event", gene]].dropna()
    df_gene[gene] = (df_gene[gene] - df_gene[gene].mean()) / df_gene[gene].std()

    try:
        cph.fit(df_gene, duration_col="OS_months", event_col="OS_event")
        summary = cph.summary
        hr = np.exp(summary["coef"].values[0])
        ci_low = np.exp(summary["coef lower 95%"].values[0])
        ci_high = np.exp(summary["coef upper 95%"].values[0])
        pval = summary["p"].values[0]
        results.append({
            "Gene": gene,
            "HR": hr,
            "CI_low": ci_low,
            "CI_high": ci_high,
            "p_value": pval,
            "significant_005": pval < 0.05,
            "significant_010": pval < 0.10,
        })
    except Exception as e:
        print(f"  {gene}: Cox regression failed - {e}")
        results.append({"Gene": gene, "HR": np.nan, "p_value": np.nan})

results_df = pd.DataFrame(results).sort_values("p_value")
results_df.to_csv(os.path.join(OUT_DIR, "idea2_cox_results.csv"), index=False)

print(f"\n{'Gene':<12} {'HR':>8} {'95% CI':>18} {'p-value':>10} {'Sig?':>6}")
print("-" * 60)
for _, row in results_df.iterrows():
    if pd.isna(row["HR"]):
        continue
    sig = "***" if row["p_value"] < 0.001 else ("**" if row["p_value"] < 0.01 else ("*" if row["p_value"] < 0.05 else ("~" if row["p_value"] < 0.10 else "")))
    print(f"{row['Gene']:<12} {row['HR']:>8.3f} ({row['CI_low']:.3f}-{row['CI_high']:.3f}) {row['p_value']:>10.4f} {sig:>6}")

n_sig_005 = results_df["significant_005"].sum()
n_sig_010 = results_df["significant_010"].sum()
print(f"\nSignificant at p<0.05: {n_sig_005}/{len(available)}")
print(f"Significant at p<0.10: {n_sig_010}/{len(available)}")

# ============================================================
# COMPOSITE OXIDATIVE STRESS SCORE
# ============================================================
print("\n" + "=" * 70)
print("COMPOSITE OXIDATIVE STRESS SCORE — SURVIVAL ANALYSIS")
print("=" * 70)

# Method 1: Simple mean z-score across all oxidative stress genes
df_score = merged[["patientId", "OS_months", "OS_event"] + available].dropna()
for gene in available:
    df_score[f"{gene}_z"] = (df_score[gene] - df_score[gene].mean()) / df_score[gene].std()

z_cols = [f"{g}_z" for g in available]
df_score["ox_stress_score"] = df_score[z_cols].mean(axis=1)

# Split by median
median_score = df_score["ox_stress_score"].median()
high_ox = df_score[df_score["ox_stress_score"] >= median_score]
low_ox = df_score[df_score["ox_stress_score"] < median_score]

lr = logrank_test(high_ox["OS_months"], low_ox["OS_months"],
                  event_observed_A=high_ox["OS_event"],
                  event_observed_B=low_ox["OS_event"])

print(f"Composite score (mean z-score of {len(available)} genes):")
print(f"  High oxidative stress: n={len(high_ox)}, deaths={high_ox['OS_event'].sum()}")
print(f"  Low oxidative stress:  n={len(low_ox)}, deaths={low_ox['OS_event'].sum()}")
print(f"  Log-rank p = {lr.p_value:.4f}")

# Method 2: Score using only significant genes (p < 0.10)
sig_genes = results_df[results_df["significant_010"] == True]["Gene"].tolist()
if len(sig_genes) >= 3:
    # Weight by direction: risk genes positive, protective genes negative
    sig_results = results_df[results_df["Gene"].isin(sig_genes)]
    z_sig_cols = []
    for _, row in sig_results.iterrows():
        gene = row["Gene"]
        col = f"{gene}_z"
        if col in df_score.columns:
            # If HR > 1 (risk), keep positive z-score; if HR < 1 (protective), flip
            if row["HR"] < 1:
                df_score[f"{gene}_weighted"] = -df_score[col]
            else:
                df_score[f"{gene}_weighted"] = df_score[col]
            z_sig_cols.append(f"{gene}_weighted")

    df_score["ox_risk_score"] = df_score[z_sig_cols].mean(axis=1)

    median_risk = df_score["ox_risk_score"].median()
    high_risk = df_score[df_score["ox_risk_score"] >= median_risk]
    low_risk = df_score[df_score["ox_risk_score"] < median_risk]

    lr2 = logrank_test(high_risk["OS_months"], low_risk["OS_months"],
                       event_observed_A=high_risk["OS_event"],
                       event_observed_B=low_risk["OS_event"])

    print(f"\nRefined score (weighted, {len(sig_genes)} significant genes):")
    print(f"  Genes used: {sig_genes}")
    print(f"  High risk: n={len(high_risk)}, deaths={high_risk['OS_event'].sum()}")
    print(f"  Low risk:  n={len(low_risk)}, deaths={low_risk['OS_event'].sum()}")
    print(f"  Log-rank p = {lr2.p_value:.6f}")

# ============================================================
# KM CURVES FOR TOP GENES
# ============================================================
top_genes = results_df.head(6)["Gene"].tolist()

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for idx, gene in enumerate(top_genes):
    ax = axes[idx]
    df_gene = merged[["OS_months", "OS_event", gene]].dropna()
    median_val = df_gene[gene].median()
    high = df_gene[df_gene[gene] >= median_val]
    low = df_gene[df_gene[gene] < median_val]
    lr_gene = logrank_test(high["OS_months"], low["OS_months"],
                           event_observed_A=high["OS_event"],
                           event_observed_B=low["OS_event"])

    kmf = KaplanMeierFitter()
    kmf.fit(high["OS_months"], high["OS_event"], label=f"High (n={len(high)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="red", alpha=0.7)
    kmf.fit(low["OS_months"], low["OS_event"], label=f"Low (n={len(low)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="blue", alpha=0.7)

    ax.set_title(f"{gene}\nLog-rank p = {lr_gene.p_value:.4f}", fontsize=11)
    ax.set_xlabel("Overall Survival (months)")
    ax.set_ylabel("Survival Probability")
    ax.legend(fontsize=9)
    ax.set_xlim(0, 100)

plt.suptitle("Idea 2: Oxidative Stress Genes — KM Survival in TCGA-LIHC",
             fontsize=14, fontweight="bold", y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "idea2_km_curves.png"), dpi=150, bbox_inches="tight")
print("\nSaved KM curves plot.")

# KM for composite score
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# All genes score
kmf = KaplanMeierFitter()
ax = axes[0]
kmf.fit(high_ox["OS_months"], high_ox["OS_event"], label=f"High OxStress (n={len(high_ox)})")
kmf.plot_survival_function(ax=ax, ci_show=True, color="red")
kmf.fit(low_ox["OS_months"], low_ox["OS_event"], label=f"Low OxStress (n={len(low_ox)})")
kmf.plot_survival_function(ax=ax, ci_show=True, color="blue")
ax.set_title(f"Composite Score (all {len(available)} genes)\nLog-rank p = {lr.p_value:.4f}")
ax.set_xlabel("Overall Survival (months)")
ax.set_ylabel("Survival Probability")
ax.legend()
ax.set_xlim(0, 100)

if len(sig_genes) >= 3:
    ax = axes[1]
    kmf.fit(high_risk["OS_months"], high_risk["OS_event"], label=f"High Risk (n={len(high_risk)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="red")
    kmf.fit(low_risk["OS_months"], low_risk["OS_event"], label=f"Low Risk (n={len(low_risk)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="blue")
    ax.set_title(f"Refined Score ({len(sig_genes)} sig genes)\nLog-rank p = {lr2.p_value:.6f}")
    ax.set_xlabel("Overall Survival (months)")
    ax.set_ylabel("Survival Probability")
    ax.legend()
    ax.set_xlim(0, 100)

plt.suptitle("Idea 2: Oxidative Stress Score — KM Survival in TCGA-LIHC",
             fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "idea2_composite_km.png"), dpi=150, bbox_inches="tight")
print("Saved composite KM plot.")

# Forest plot
fig, ax = plt.subplots(figsize=(10, 12))
valid = results_df.dropna(subset=["HR"]).sort_values("p_value", ascending=False)
y_pos = range(len(valid))
colors = ['red' if p < 0.05 else ('orange' if p < 0.10 else 'gray') for p in valid["p_value"]]
ax.scatter(valid["HR"], y_pos, c=colors, s=50, zorder=5)
for i, (_, row) in enumerate(valid.iterrows()):
    ax.plot([row["CI_low"], row["CI_high"]], [i, i], c=colors[i], linewidth=1.5)
ax.axvline(x=1.0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
ax.set_yticks(list(y_pos))
ax.set_yticklabels(valid["Gene"].tolist())
ax.set_xlabel("Hazard Ratio (95% CI)")
ax.set_title("Idea 2: Univariate Cox Regression — Oxidative Stress Genes in TCGA-LIHC\n(Red: p<0.05, Orange: p<0.10, Gray: NS)")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "idea2_forest_plot.png"), dpi=150, bbox_inches="tight")
print("Saved forest plot.")

# ============================================================
# CLUSTERING ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("HIERARCHICAL CLUSTERING — OXIDATIVE STRESS SUBGROUPS")
print("=" * 70)

# Cluster patients based on oxidative stress gene expression
expr_matrix = df_score[available].values
# Z-score normalize
expr_z = (expr_matrix - expr_matrix.mean(axis=0)) / expr_matrix.std(axis=0)

# Hierarchical clustering
dist = pdist(expr_z, metric='euclidean')
link = linkage(dist, method='ward')
clusters = fcluster(link, t=2, criterion='maxclust')

df_score["cluster"] = clusters
c1 = df_score[df_score["cluster"] == 1]
c2 = df_score[df_score["cluster"] == 2]

lr_clust = logrank_test(c1["OS_months"], c2["OS_months"],
                        event_observed_A=c1["OS_event"],
                        event_observed_B=c2["OS_event"])

# Determine which cluster is high/low oxidative stress
c1_mean = c1["ox_stress_score"].mean()
c2_mean = c2["ox_stress_score"].mean()
if c1_mean > c2_mean:
    high_label, low_label = "Cluster 1 (High OxStress)", "Cluster 2 (Low OxStress)"
else:
    high_label, low_label = "Cluster 1 (Low OxStress)", "Cluster 2 (High OxStress)"

print(f"Cluster 1: n={len(c1)}, mean OxStress score={c1_mean:.3f}, deaths={c1['OS_event'].sum()}")
print(f"Cluster 2: n={len(c2)}, mean OxStress score={c2_mean:.3f}, deaths={c2['OS_event'].sum()}")
print(f"Log-rank p = {lr_clust.p_value:.4f}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("IDEA 2 PILOT SUMMARY")
print("=" * 70)
print(f"Total oxidative stress genes tested: {len(available)}")
print(f"Significant at p < 0.05: {n_sig_005} ({n_sig_005/len(available)*100:.0f}%)")
print(f"Significant at p < 0.10: {n_sig_010} ({n_sig_010/len(available)*100:.0f}%)")
print(f"\nComposite score log-rank p: {lr.p_value:.4f}")
if len(sig_genes) >= 3:
    print(f"Refined score log-rank p: {lr2.p_value:.6f}")
print(f"Clustering log-rank p: {lr_clust.p_value:.4f}")
print(f"\nTop 5 genes by significance:")
for _, row in results_df.head(5).iterrows():
    direction = "risk" if row["HR"] > 1 else "protective"
    print(f"  {row['Gene']}: HR={row['HR']:.3f}, p={row['p_value']:.4f} ({direction})")

print(f"\nVERDICT: ", end="")
if n_sig_005 >= 8 and lr.p_value < 0.01:
    print("STRONG signal. Robust oxidative stress signature with clear survival separation.")
elif n_sig_005 >= 5:
    print("GOOD signal. Enough prognostic genes for a solid signature paper.")
elif n_sig_005 >= 3:
    print("MODERATE signal. Paper is feasible but needs careful framing.")
else:
    print("WEAK signal. May not be enough for a standalone paper.")
