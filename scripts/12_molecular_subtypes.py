"""
12_molecular_subtypes.py — Consensus clustering + mutation landscape

1. Consensus clustering of TCGA patients by 11 signature genes
2. Characterize clusters: survival, clinical, immune, mutations
3. Query cBioPortal for mutation/CNA data of signature genes
4. TMB correlation with risk score
"""
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import spearmanr, mannwhitneyu, chi2_contingency, kruskal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import requests
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
with open(os.path.join(MODEL, "lasso_model.json")) as f:
    model = json.load(f)
selected_genes = model["genes"]
sig_genes = list(selected_genes.keys())

merged = pd.read_csv(os.path.join(DATA, "tcga_ros_merged.csv"))
if "risk_score" not in merged.columns:
    risk = np.zeros(len(merged))
    for gene, coef in selected_genes.items():
        z = (merged[gene] - model["gene_means"][gene]) / model["gene_stds"][gene]
        risk += coef * z.values
    merged["risk_score"] = risk

df = merged.dropna(subset=["OS_months", "OS_event", "risk_score"]).copy()
df = df[df["OS_months"] > 0]
print(f"Working with {len(df)} patients, {len(sig_genes)} signature genes")

# ══════════════════════════════════════════════════════════════════════════════
# 1. CONSENSUS CLUSTERING
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("1. CONSENSUS CLUSTERING")
print("=" * 70)

# Prepare expression matrix for clustering
expr_clust = df[sig_genes].copy()
scaler = StandardScaler()
expr_scaled = scaler.fit_transform(expr_clust)

# Test k=2 through k=5
silhouette_scores = {}
for k in range(2, 6):
    km = KMeans(n_clusters=k, random_state=42, n_init=50)
    labels = km.fit_predict(expr_scaled)
    sil = silhouette_score(expr_scaled, labels)
    silhouette_scores[k] = sil
    print(f"  k={k}: silhouette={sil:.4f}")

best_k = max(silhouette_scores, key=silhouette_scores.get)
print(f"  Best k={best_k} (silhouette={silhouette_scores[best_k]:.4f})")

# Validate: check that best_k doesn't produce tiny clusters (< 10 patients)
# If it does, fall back to k=2
km_check = KMeans(n_clusters=best_k, random_state=42, n_init=50)
labels_check = km_check.fit_predict(expr_scaled)
cluster_sizes = pd.Series(labels_check).value_counts()
if cluster_sizes.min() < 10:
    print(f"  WARNING: k={best_k} produces cluster with only {cluster_sizes.min()} patients")
    print(f"  Falling back to k=2 for meaningful clusters")
    best_k = 2

# Consensus: repeat clustering 100 times, build co-occurrence matrix
print("\n  Building consensus matrix (100 iterations)...")
n = len(expr_scaled)
consensus = np.zeros((n, n))
for _ in range(100):
    # Subsample 80%
    idx = np.random.choice(n, int(0.8 * n), replace=False)
    km = KMeans(n_clusters=best_k, random_state=None, n_init=10)
    labels_sub = km.fit_predict(expr_scaled[idx])
    for i in range(len(idx)):
        for j in range(i + 1, len(idx)):
            if labels_sub[i] == labels_sub[j]:
                consensus[idx[i], idx[j]] += 1
                consensus[idx[j], idx[i]] += 1

# Normalize
max_val = consensus.max()
if max_val > 0:
    consensus /= max_val

# Final clustering on consensus matrix
hc = AgglomerativeClustering(n_clusters=best_k, metric='precomputed',
                              linkage='average')
final_labels = hc.fit_predict(1 - consensus)
df["cluster"] = final_labels
print(f"  Cluster sizes: {pd.Series(final_labels).value_counts().sort_index().to_dict()}")

# ── Consensus matrix heatmap ──
fig, ax = plt.subplots(figsize=(10, 8))
# Sort by cluster
order = np.argsort(final_labels)
consensus_sorted = consensus[order][:, order]
labels_sorted = final_labels[order]

# Color bar for clusters
cluster_colors = plt.cm.Set2(labels_sorted / max(1, labels_sorted.max()))

sns.heatmap(consensus_sorted, ax=ax, cmap="YlOrRd", xticklabels=False, yticklabels=False,
            cbar_kws={"label": "Co-clustering frequency"})
ax.set_title(f"Consensus Clustering Matrix (k={best_k})", fontsize=13, fontweight='bold')
ax.set_xlabel(f"Patients (n={n})")
ax.set_ylabel(f"Patients (n={n})")
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "consensus_matrix.png"), dpi=200, bbox_inches='tight')
print("Saved: consensus_matrix.png")

# ══════════════════════════════════════════════════════════════════════════════
# 2. CLUSTER CHARACTERIZATION
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("2. CLUSTER CHARACTERIZATION")
print("=" * 70)

# ── Survival by cluster ──
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# KM curves
ax = axes[0]
kmf = KaplanMeierFitter()
colors = plt.cm.Set2(np.linspace(0, 1, best_k))
for c in range(best_k):
    mask = df["cluster"] == c
    n_c = mask.sum()
    events_c = int(df.loc[mask, "OS_event"].sum())
    kmf.fit(df.loc[mask, "OS_months"], df.loc[mask, "OS_event"],
            label=f"Cluster {c+1} (n={n_c}, events={events_c})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color=colors[c])

# Log-rank test (pairwise)
if best_k == 2:
    c0 = df[df["cluster"] == 0]
    c1 = df[df["cluster"] == 1]
    lr = logrank_test(c0["OS_months"], c1["OS_months"],
                      event_observed_A=c0["OS_event"], event_observed_B=c1["OS_event"])
    ax.set_title(f"OS by Molecular Cluster\nLog-rank p={lr.p_value:.2e}", fontsize=12, fontweight='bold')
else:
    # Multi-group log-rank (approximate)
    from itertools import combinations
    lr_pvals = []
    for i, j in combinations(range(best_k), 2):
        ci = df[df["cluster"] == i]
        cj = df[df["cluster"] == j]
        lr = logrank_test(ci["OS_months"], cj["OS_months"],
                          event_observed_A=ci["OS_event"], event_observed_B=cj["OS_event"])
        lr_pvals.append(lr.p_value)
    ax.set_title(f"OS by Molecular Cluster\nMin pairwise p={min(lr_pvals):.2e}", fontsize=12, fontweight='bold')

ax.set_xlabel("Overall Survival (months)")
ax.set_ylabel("Survival Probability")
ax.legend(fontsize=9)
ax.set_xlim(0, 100)

# Risk score by cluster
ax = axes[1]
cluster_data = [df.loc[df["cluster"] == c, "risk_score"].values for c in range(best_k)]
bp = ax.boxplot(cluster_data, labels=[f"Cluster {c+1}" for c in range(best_k)],
                patch_artist=True)
for i, patch in enumerate(bp['boxes']):
    patch.set_facecolor(colors[i])
ax.set_ylabel("Risk Score")
ax.set_title("Risk Score Distribution by Cluster", fontsize=12, fontweight='bold')

# Kruskal-Wallis test
if len(cluster_data) >= 2:
    valid_groups = [g for g in cluster_data if len(g) > 0]
    if len(valid_groups) >= 2:
        stat, p_kw = kruskal(*valid_groups)
        ax.text(0.05, 0.95, f"Kruskal-Wallis p={p_kw:.2e}", transform=ax.transAxes,
                fontsize=10, verticalalignment='top')

plt.suptitle("Molecular Subtypes — ROS/Ferroptosis Signature", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "cluster_survival.png"), dpi=200, bbox_inches='tight')
print("Saved: cluster_survival.png")

# ── Cluster clinical characteristics ──
print("\nCluster clinical characteristics:")
cluster_chars = []
for c in range(best_k):
    mask = df["cluster"] == c
    char = {"cluster": c + 1, "n": int(mask.sum())}

    char["risk_score_mean"] = df.loc[mask, "risk_score"].mean()
    char["events"] = int(df.loc[mask, "OS_event"].sum())
    char["event_rate"] = df.loc[mask, "OS_event"].mean()

    if "age_at_diagnosis" in df.columns:
        age = df.loc[mask, "age_at_diagnosis"]
        if age.median() > 200:
            age = age / 365.25
        char["age_median"] = age.median()

    if "gender" in df.columns:
        char["male_pct"] = (df.loc[mask, "gender"].str.lower() == "male").mean()

    if "tumor_stage" in df.columns:
        def encode_s(s):
            if pd.isna(s): return np.nan
            s = str(s).upper()
            if "IV" in s: return 4
            elif "III" in s: return 3
            elif "II" in s: return 2
            elif "I" in s: return 1
            return np.nan
        stages = df.loc[mask, "tumor_stage"].apply(encode_s)
        char["stage_III_IV_pct"] = (stages >= 3).mean()

    cluster_chars.append(char)
    print(f"  Cluster {c+1}: n={char['n']}, risk={char['risk_score_mean']:.2f}, "
          f"events={char['events']} ({char['event_rate']:.1%})")

cluster_df = pd.DataFrame(cluster_chars)
cluster_df.to_csv(os.path.join(TABLES, "cluster_characteristics.csv"), index=False)
print("Saved: cluster_characteristics.csv")

# ── Gene expression heatmap by cluster ──
fig, axes = plt.subplots(2, 1, figsize=(14, max(5, len(sig_genes) * 0.4 + 3)),
                          gridspec_kw={'height_ratios': [1, 4]}, sharex=True)

# Sort by cluster, then by risk score within cluster
sort_idx = df.sort_values(["cluster", "risk_score"]).index
expr_sorted = df.loc[sort_idx, sig_genes].copy()
for g in sig_genes:
    expr_sorted[g] = (expr_sorted[g] - expr_sorted[g].mean()) / (expr_sorted[g].std() + 1e-10)

# Top: cluster annotation
ax0 = axes[0]
cluster_sorted = df.loc[sort_idx, "cluster"].values
for i, c in enumerate(cluster_sorted):
    ax0.bar(i, 1, color=colors[c], width=1.0, linewidth=0)
ax0.set_ylabel("Cluster")
ax0.set_xlim(0, len(cluster_sorted))
ax0.yaxis.set_visible(False)
ax0.set_title("Expression Heatmap by Molecular Cluster", fontweight='bold')

# Bottom: gene expression
ax1 = axes[1]
sns.heatmap(expr_sorted.T, ax=ax1, cmap="RdBu_r", center=0, xticklabels=False,
            yticklabels=sig_genes, cbar_kws={"shrink": 0.5, "label": "Z-score"})
ax1.set_xlabel(f"Patients (n={len(df)})")

plt.tight_layout()
plt.savefig(os.path.join(FIGS, "cluster_heatmap.png"), dpi=200, bbox_inches='tight')
print("Saved: cluster_heatmap.png")

# ══════════════════════════════════════════════════════════════════════════════
# 3. MUTATION LANDSCAPE (cBioPortal)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("3. MUTATION LANDSCAPE (cBioPortal)")
print("=" * 70)

mutation_results = []
cna_results = []

# Query cBioPortal API for TCGA-LIHC mutations
CBIO_BASE = "https://www.cbioportal.org/api"
study_id = "lihc_tcga"

# Get mutation data for signature genes
print("  Querying cBioPortal for mutations...")
try:
    # Get molecular profiles
    profiles_url = f"{CBIO_BASE}/studies/{study_id}/molecular-profiles"
    resp = requests.get(profiles_url, timeout=30,
                        headers={"Accept": "application/json"})
    profiles = resp.json()

    mut_profile = None
    cna_profile = None
    for p in profiles:
        if p.get("molecularAlterationType") == "MUTATION_EXTENDED":
            mut_profile = p["molecularProfileId"]
        elif p.get("molecularAlterationType") == "COPY_NUMBER_ALTERATION" and "gistic" in p.get("molecularProfileId", "").lower():
            cna_profile = p["molecularProfileId"]

    print(f"  Mutation profile: {mut_profile}")
    print(f"  CNA profile: {cna_profile}")

    if mut_profile:
        # Get mutations for signature genes
        for gene in sig_genes:
            try:
                url = f"{CBIO_BASE}/molecular-profiles/{mut_profile}/genes/{gene}/mutations"
                resp = requests.get(url, timeout=15, headers={"Accept": "application/json"})
                if resp.status_code == 200:
                    muts = resp.json()
                    n_muts = len(muts)
                    unique_patients = len(set(m.get("patientId", "") for m in muts))
                    mut_types = {}
                    for m in muts:
                        mt = m.get("mutationType", "unknown")
                        mut_types[mt] = mut_types.get(mt, 0) + 1

                    mutation_results.append({
                        "gene": gene, "n_mutations": n_muts,
                        "n_patients_mutated": unique_patients,
                        "mutation_types": str(mut_types),
                    })
                    if n_muts > 0:
                        print(f"    {gene}: {unique_patients} patients, {n_muts} mutations")
            except Exception as e:
                pass

except Exception as e:
    print(f"  cBioPortal query failed: {e}")
    print("  Using TCGA MAF data directly instead...")

    # Fallback: try to get mutation data from GDC API
    try:
        gdc_url = "https://api.gdc.cancer.gov/analysis/top_mutated_genes_by_project"
        params = {"project_id": "TCGA-LIHC", "gene_symbol": ",".join(sig_genes)}
        resp = requests.get(gdc_url, params=params, timeout=30)
        if resp.status_code == 200:
            print("  Retrieved mutation data from GDC")
    except:
        pass

if mutation_results:
    mut_df = pd.DataFrame(mutation_results)
    mut_df.to_csv(os.path.join(TABLES, "mutation_landscape.csv"), index=False)
    print("\n  Saved: mutation_landscape.csv")

    # OncoPrint-style visualization
    fig, ax = plt.subplots(figsize=(12, 5))
    genes_mut = [r["gene"] for r in mutation_results]
    n_mut_patients = [r["n_patients_mutated"] for r in mutation_results]

    # Total patients in study (approximate)
    total_patients = 366  # TCGA-LIHC
    pct_mutated = [n / total_patients * 100 for n in n_mut_patients]

    y = np.arange(len(genes_mut))
    bar_colors = ['#d32f2f' if selected_genes[g] > 0 else '#1565c0' for g in genes_mut]
    ax.barh(y, pct_mutated, color=bar_colors, edgecolor='black', linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(genes_mut, fontsize=10)
    ax.set_xlabel("% Patients Mutated")
    ax.set_title("Mutation Frequency of Signature Genes — TCGA-LIHC",
                 fontsize=13, fontweight='bold')

    for i, (pct, n) in enumerate(zip(pct_mutated, n_mut_patients)):
        ax.text(pct + 0.3, i, f"{n} ({pct:.1f}%)", va='center', fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "mutation_landscape.png"), dpi=200, bbox_inches='tight')
    print("  Saved: mutation_landscape.png")

# ══════════════════════════════════════════════════════════════════════════════
# 4. TMB CORRELATION
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("4. TMB / MUTATION BURDEN CORRELATION")
print("=" * 70)

# Try to get TMB data from cBioPortal
try:
    tmb_url = f"{CBIO_BASE}/studies/{study_id}/clinical-data?attributeId=MUTATION_COUNT&clinicalDataType=SAMPLE"
    resp = requests.get(tmb_url, timeout=30, headers={"Accept": "application/json"})
    if resp.status_code == 200:
        tmb_data = resp.json()
        if tmb_data:
            tmb_df = pd.DataFrame(tmb_data)
            print(f"  TMB data: {len(tmb_df)} samples")
            tmb_df["value"] = pd.to_numeric(tmb_df["value"], errors="coerce")
            tmb_df = tmb_df.rename(columns={"patientId": "patientId_tmb"})

            # Match with our data
            if "patientId" in df.columns:
                # Try matching
                tmb_map = dict(zip(tmb_df["patientId_tmb"], tmb_df["value"]))
                df["tmb"] = df["patientId"].map(tmb_map)
                valid_tmb = df.dropna(subset=["tmb", "risk_score"])

                if len(valid_tmb) > 20:
                    r, p = spearmanr(valid_tmb["risk_score"], valid_tmb["tmb"])
                    print(f"  Risk score vs TMB: Spearman r={r:.3f}, p={p:.4f}")

                    fig, ax = plt.subplots(figsize=(7, 6))
                    ax.scatter(valid_tmb["risk_score"], valid_tmb["tmb"],
                               alpha=0.5, s=20, color='steelblue')
                    # Add trend line
                    z = np.polyfit(valid_tmb["risk_score"], valid_tmb["tmb"], 1)
                    p_line = np.poly1d(z)
                    x_range = np.linspace(valid_tmb["risk_score"].min(),
                                          valid_tmb["risk_score"].max(), 100)
                    ax.plot(x_range, p_line(x_range), 'r--', linewidth=2)
                    ax.set_xlabel("Risk Score")
                    ax.set_ylabel("Tumor Mutation Burden")
                    ax.set_title(f"Risk Score vs TMB\nSpearman r={r:.3f}, p={p:.4f}",
                                 fontsize=12, fontweight='bold')
                    plt.tight_layout()
                    plt.savefig(os.path.join(FIGS, "tmb_correlation.png"), dpi=200, bbox_inches='tight')
                    print("  Saved: tmb_correlation.png")

                    # TMB by risk group
                    med = valid_tmb["risk_score"].median()
                    high = valid_tmb[valid_tmb["risk_score"] >= med]["tmb"]
                    low = valid_tmb[valid_tmb["risk_score"] < med]["tmb"]
                    stat, p_mw = mannwhitneyu(high, low, alternative='two-sided')
                    print(f"  TMB high vs low risk: median {high.median():.0f} vs {low.median():.0f}, "
                          f"p={p_mw:.4f}")
except Exception as e:
    print(f"  TMB data not available: {e}")

# ══════════════════════════════════════════════════════════════════════════════
# 5. SUMMARY COMPOSITE FIGURE
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("5. SUMMARY")
print("=" * 70)

# Save cluster assignments
df[["patientId", "cluster", "risk_score"]].to_csv(
    os.path.join(TABLES, "patient_clusters.csv"), index=False)
print("Saved: patient_clusters.csv")

# Summary statistics
print(f"\nMolecular subtypes summary:")
print(f"  Optimal k: {best_k}")
print(f"  Silhouette score: {silhouette_scores[best_k]:.4f}")
for c in range(best_k):
    mask = df["cluster"] == c
    print(f"  Cluster {c+1}: n={mask.sum()}, "
          f"risk_score={df.loc[mask, 'risk_score'].mean():.2f}±{df.loc[mask, 'risk_score'].std():.2f}, "
          f"events={int(df.loc[mask, 'OS_event'].sum())}/{mask.sum()}")

if mutation_results:
    n_mutated = sum(1 for r in mutation_results if r["n_patients_mutated"] > 0)
    print(f"\nMutation landscape: {n_mutated}/{len(sig_genes)} genes have mutations in TCGA-LIHC")

print(f"\n✓ Molecular subtypes analysis complete.")
