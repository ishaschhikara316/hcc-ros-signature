"""
02_univariate_screening.py — Univariate Cox regression on all ROS/ferroptosis genes

For each gene: z-score → Cox PH → HR, 95% CI, p-value
Multiple testing correction: Benjamini-Hochberg FDR + Bonferroni
Generates: volcano plot, forest plot, KM curves for top 6 genes
"""
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data", "tcga")
TABLES = os.path.join(BASE, "results", "tables")
FIGS = os.path.join(BASE, "results", "figures")
os.makedirs(TABLES, exist_ok=True)
os.makedirs(FIGS, exist_ok=True)

# ── Load data ───────────────────────────────────────────────────────────────
merged = pd.read_csv(os.path.join(DATA, "tcga_ros_merged.csv"))
gene_info = pd.read_csv(os.path.join(DATA, "gene_set_info.csv"))
available_genes = gene_info[gene_info["available"]]["gene"].tolist()
available_genes = [g for g in available_genes if g in merged.columns]
print(f"Loaded {len(merged)} patients, {len(available_genes)} genes")

# ── Univariate Cox regression ───────────────────────────────────────────────
print("\n" + "=" * 70)
print("UNIVARIATE COX REGRESSION")
print("=" * 70)

results = []
cph = CoxPHFitter()

for gene in available_genes:
    df_gene = merged[["OS_months", "OS_event", gene]].dropna()
    if df_gene[gene].std() == 0:
        continue
    df_gene[gene] = (df_gene[gene] - df_gene[gene].mean()) / df_gene[gene].std()
    try:
        cph.fit(df_gene, duration_col="OS_months", event_col="OS_event")
        s = cph.summary
        results.append({
            "gene": gene,
            "pathway": gene_info[gene_info["gene"] == gene]["pathway"].values[0],
            "HR": np.exp(s["coef"].values[0]),
            "coef": s["coef"].values[0],
            "CI_low": np.exp(s["coef lower 95%"].values[0]),
            "CI_high": np.exp(s["coef upper 95%"].values[0]),
            "p_value": s["p"].values[0],
            "n": len(df_gene),
        })
    except Exception as e:
        print(f"  {gene}: FAILED — {e}")

uv = pd.DataFrame(results).sort_values("p_value")

# ── Multiple testing correction ─────────────────────────────────────────────
n_tests = len(uv)
# Benjamini-Hochberg FDR
ranked_idx = np.argsort(uv["p_value"].values)
fdr = np.zeros(n_tests)
raw_p = uv["p_value"].values
for i, idx in enumerate(ranked_idx):
    rank = i + 1
    fdr[idx] = raw_p[idx] * n_tests / rank
# Enforce monotonicity (step-up)
fdr_sorted_idx = np.argsort(-fdr)
running_min = 1.0
for idx in np.argsort(np.argsort(-raw_p)):
    fdr[idx] = min(fdr[idx], running_min)
    running_min = fdr[idx]
# Proper monotonicity: iterate from largest rank down
fdr2 = np.zeros(n_tests)
for i, idx in enumerate(ranked_idx):
    rank = i + 1
    fdr2[idx] = raw_p[idx] * n_tests / rank
running = 1.0
for i in range(n_tests - 1, -1, -1):
    idx = ranked_idx[i]
    fdr2[idx] = min(fdr2[idx], running)
    running = fdr2[idx]

uv["FDR"] = fdr2
uv["Bonferroni_p"] = np.minimum(uv["p_value"] * n_tests, 1.0)

# Summary
print(f"\n{'Gene':<12} {'HR':>7} {'95% CI':>18} {'p-value':>12} {'FDR':>10} {'Pathway'}")
print("-" * 85)
for _, r in uv.iterrows():
    sig = "***" if r["p_value"] < 0.001 else ("**" if r["p_value"] < 0.01 else ("*" if r["p_value"] < 0.05 else ""))
    print(f"{r['gene']:<12} {r['HR']:>7.3f} ({r['CI_low']:.3f}-{r['CI_high']:.3f}) {r['p_value']:>12.2e} {r['FDR']:>10.4f} {sig}  {r['pathway']}")

n_p05 = (uv["p_value"] < 0.05).sum()
n_fdr10 = (uv["FDR"] < 0.10).sum()
n_fdr05 = (uv["FDR"] < 0.05).sum()
n_bonf = (uv["Bonferroni_p"] < 0.05).sum()
print(f"\nSignificant genes:")
print(f"  p < 0.05:        {n_p05}/{n_tests}")
print(f"  FDR < 0.10:      {n_fdr10}/{n_tests}")
print(f"  FDR < 0.05:      {n_fdr05}/{n_tests}")
print(f"  Bonferroni<0.05: {n_bonf}/{n_tests}")

# Save
uv.to_csv(os.path.join(TABLES, "univariate_results.csv"), index=False)
print(f"\nSaved: univariate_results.csv")

# ── Volcano-style plot (HR vs -log10 p) ────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 7))
log_hr = np.log2(uv["HR"])
neg_logp = -np.log10(uv["p_value"])

# Color by significance
colors = []
for _, r in uv.iterrows():
    if r["FDR"] < 0.05:
        colors.append("red" if r["HR"] > 1 else "blue")
    elif r["FDR"] < 0.10:
        colors.append("salmon" if r["HR"] > 1 else "lightblue")
    else:
        colors.append("gray")

ax.scatter(log_hr, neg_logp, c=colors, s=60, edgecolors='black', linewidths=0.5, zorder=3)
ax.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, label='p=0.05')
ax.axvline(0, color='gray', linestyle='--', alpha=0.5)

# Label significant genes
for _, r in uv[uv["FDR"] < 0.10].iterrows():
    ax.annotate(r["gene"], (np.log2(r["HR"]), -np.log10(r["p_value"])),
                fontsize=8, ha='center', va='bottom', fontweight='bold')

ax.set_xlabel("log2(Hazard Ratio)", fontsize=12)
ax.set_ylabel("-log10(p-value)", fontsize=12)
ax.set_title("Univariate Cox Regression — ROS/Ferroptosis Genes in TCGA-LIHC", fontsize=13, fontweight='bold')
ax.legend(["FDR<0.05 risk", "FDR<0.05 protective", "p=0.05 threshold"], loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "univariate_volcano.png"), dpi=200, bbox_inches='tight')
print("Saved: univariate_volcano.png")

# ── Forest plot ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, max(8, len(uv) * 0.35)))
sorted_uv = uv.sort_values("p_value", ascending=False)
y_pos = range(len(sorted_uv))

fcolors = []
for _, r in sorted_uv.iterrows():
    if r["FDR"] < 0.05:
        fcolors.append("red")
    elif r["FDR"] < 0.10:
        fcolors.append("orange")
    elif r["p_value"] < 0.05:
        fcolors.append("goldenrod")
    else:
        fcolors.append("gray")

for i, (_, r) in enumerate(sorted_uv.iterrows()):
    ax.plot([r["CI_low"], r["CI_high"]], [i, i], color=fcolors[i], linewidth=1.5)
ax.scatter(sorted_uv["HR"], y_pos, c=fcolors, s=40, zorder=5, edgecolors='black', linewidths=0.3)
ax.axvline(1.0, color='black', linestyle='--', linewidth=0.8)
ax.set_yticks(list(y_pos))
ax.set_yticklabels(sorted_uv["gene"].tolist(), fontsize=9)
ax.set_xlabel("Hazard Ratio (95% CI)", fontsize=11)
ax.set_title("Forest Plot — Univariate Cox Regression\n(Red: FDR<0.05, Orange: FDR<0.10, Gold: p<0.05, Gray: NS)",
             fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "univariate_forest.png"), dpi=200, bbox_inches='tight')
print("Saved: univariate_forest.png")

# ── KM curves for top 6 genes ──────────────────────────────────────────────
top6 = uv.head(6)["gene"].tolist()
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
axes = axes.flatten()

for idx, gene in enumerate(top6):
    ax = axes[idx]
    df_g = merged[["OS_months", "OS_event", gene]].dropna()
    med = df_g[gene].median()
    high = df_g[df_g[gene] >= med]
    low = df_g[df_g[gene] < med]
    lr = logrank_test(high["OS_months"], low["OS_months"],
                      event_observed_A=high["OS_event"], event_observed_B=low["OS_event"])
    kmf = KaplanMeierFitter()
    kmf.fit(high["OS_months"], high["OS_event"], label=f"High (n={len(high)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="red", alpha=0.8)
    kmf.fit(low["OS_months"], low["OS_event"], label=f"Low (n={len(low)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="blue", alpha=0.8)

    hr_val = uv[uv["gene"] == gene]["HR"].values[0]
    fdr_val = uv[uv["gene"] == gene]["FDR"].values[0]
    ax.set_title(f"{gene}\nHR={hr_val:.2f}, log-rank p={lr.p_value:.2e}", fontsize=10, fontweight='bold')
    ax.set_xlabel("Months")
    ax.set_ylabel("Survival Probability")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 100)

plt.suptitle("Top 6 ROS/Ferroptosis Genes — Kaplan-Meier Survival", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "top_genes_km.png"), dpi=200, bbox_inches='tight')
print("Saved: top_genes_km.png")

print("\n✓ Univariate screening complete.")
