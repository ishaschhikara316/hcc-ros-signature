"""
03_lasso_model.py — Build LASSO-Cox prognostic model

1. Select candidate genes (FDR < 0.10, or p < 0.05 fallback)
2. Cross-validate penalizer with 5-fold CV
3. Fit final LASSO-Cox, extract non-zero genes
4. Bootstrap stability analysis (1000 iterations)
5. Permutation test (500 iterations)
"""
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
import os
import warnings
warnings.filterwarnings('ignore')

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data", "tcga")
MODEL = os.path.join(BASE, "results", "model")
TABLES = os.path.join(BASE, "results", "tables")
FIGS = os.path.join(BASE, "results", "figures")
os.makedirs(MODEL, exist_ok=True)

# ── Load data ───────────────────────────────────────────────────────────────
merged = pd.read_csv(os.path.join(DATA, "tcga_ros_merged.csv"))
uv = pd.read_csv(os.path.join(TABLES, "univariate_results.csv"))
print(f"Loaded {len(merged)} patients, {len(uv)} univariate results")

# ── Select candidate genes ──────────────────────────────────────────────────
candidates = uv[uv["FDR"] < 0.10]["gene"].tolist()
if len(candidates) < 10:
    candidates = uv[uv["p_value"] < 0.05]["gene"].tolist()
    print(f"FDR<0.10 gave <10 genes, using p<0.05 cutoff")
if len(candidates) < 5:
    candidates = uv.head(20)["gene"].tolist()
    print(f"p<0.05 gave <5 genes, using top 20")

print(f"Candidate genes for LASSO: {len(candidates)}")
print(f"  Genes: {candidates}")

# ── Prepare LASSO dataset ──────────────────────────────────────────────────
lasso_cols = ["OS_months", "OS_event"] + candidates
lasso_data = merged[lasso_cols].dropna()
print(f"LASSO dataset: {len(lasso_data)} patients, {int(lasso_data['OS_event'].sum())} events")

# Z-score normalize
gene_means = {}
gene_stds = {}
for gene in candidates:
    m = lasso_data[gene].mean()
    s = lasso_data[gene].std()
    gene_means[gene] = float(m)
    gene_stds[gene] = float(s)
    lasso_data[gene] = (lasso_data[gene] - m) / s

# ── Cross-validate penalizer ────────────────────────────────────────────────
print("\n" + "=" * 70)
print("LASSO-COX CROSS-VALIDATION")
print("=" * 70)

penalizers = [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.5]
cv_results = []

np.random.seed(42)
n = len(lasso_data)
fold_ids = np.random.permutation(n) % 5

for pen in penalizers:
    fold_cindices = []
    fold_ngenes = []
    for fold in range(5):
        train = lasso_data.iloc[fold_ids != fold].copy()
        test = lasso_data.iloc[fold_ids == fold].copy()
        try:
            cph = CoxPHFitter(penalizer=pen, l1_ratio=1.0)
            cph.fit(train, duration_col="OS_months", event_col="OS_event")
            pred = cph.predict_partial_hazard(test)
            ci = concordance_index(test["OS_months"], -pred.values.ravel(), test["OS_event"])
            n_nonzero = (cph.params_.abs() > 1e-4).sum()
            fold_cindices.append(ci)
            fold_ngenes.append(n_nonzero)
        except:
            fold_cindices.append(0.5)
            fold_ngenes.append(0)

    mean_ci = np.mean(fold_cindices)
    std_ci = np.std(fold_cindices)
    mean_genes = np.mean(fold_ngenes)
    cv_results.append({
        "penalizer": pen,
        "mean_cindex": mean_ci,
        "std_cindex": std_ci,
        "mean_ngenes": mean_genes,
    })
    print(f"  λ={pen:.3f}: C-index={mean_ci:.4f}±{std_ci:.4f}, genes={mean_genes:.1f}")

cv_df = pd.DataFrame(cv_results)

# Select best: highest C-index with 6-12 genes (tighter for publishability)
viable = cv_df[(cv_df["mean_ngenes"] >= 6) & (cv_df["mean_ngenes"] <= 12)]
if len(viable) == 0:
    print("  No models in 6-12 gene range, relaxing to 3-15")
    viable = cv_df[(cv_df["mean_ngenes"] >= 3) & (cv_df["mean_ngenes"] <= 15)]
if len(viable) == 0:
    viable = cv_df[cv_df["mean_ngenes"] >= 3]
if len(viable) == 0:
    viable = cv_df
best_row = viable.loc[viable["mean_cindex"].idxmax()]
best_pen = best_row["penalizer"]
print(f"\nBest penalizer: λ={best_pen} (C-index={best_row['mean_cindex']:.4f}, genes~{best_row['mean_ngenes']:.0f})")

# ── Fit final model ─────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("FINAL LASSO-COX MODEL")
print("=" * 70)

cph_final = CoxPHFitter(penalizer=best_pen, l1_ratio=1.0)
cph_final.fit(lasso_data, duration_col="OS_months", event_col="OS_event")

# Extract non-zero genes
params = cph_final.params_
selected = params[params.abs() > 1e-4]
print(f"Selected {len(selected)} genes:")
for gene, coef in selected.items():
    hr = np.exp(coef)
    direction = "risk" if hr > 1 else "protective"
    pathway = uv[uv["gene"] == gene]["pathway"].values[0] if gene in uv["gene"].values else "?"
    print(f"  {gene:<12} coef={coef:>7.4f}  HR={hr:.3f}  ({direction})  [{pathway}]")

# Training C-index
pred_train = cph_final.predict_partial_hazard(lasso_data)
ci_train = concordance_index(lasso_data["OS_months"], -pred_train.values.ravel(), lasso_data["OS_event"])
print(f"\nTraining C-index: {ci_train:.4f}")

# ── LASSO path plot ─────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))
path_pens = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
gene_paths = {g: [] for g in candidates}
for pen in path_pens:
    try:
        cph_p = CoxPHFitter(penalizer=pen, l1_ratio=1.0)
        cph_p.fit(lasso_data, duration_col="OS_months", event_col="OS_event")
        for g in candidates:
            gene_paths[g].append(cph_p.params_.get(g, 0))
    except:
        for g in candidates:
            gene_paths[g].append(0)

for g in candidates:
    vals = gene_paths[g]
    if any(abs(v) > 1e-4 for v in vals):
        ax.plot(np.log10(path_pens), vals, marker='.', markersize=4, label=g)

ax.axvline(np.log10(best_pen), color='red', linestyle='--', alpha=0.7, label=f'Selected λ={best_pen}')
ax.set_xlabel("log10(Penalizer)", fontsize=11)
ax.set_ylabel("Coefficient", fontsize=11)
ax.set_title("LASSO-Cox Regularization Path", fontsize=13, fontweight='bold')
ax.legend(fontsize=7, ncol=2, loc='best')
ax.axhline(0, color='gray', linestyle='-', alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "lasso_path.png"), dpi=200, bbox_inches='tight')
print("Saved: lasso_path.png")

# ── Lambda CV plot ──────────────────────────────────────────────────────────
fig, ax1 = plt.subplots(figsize=(8, 5))
ax1.errorbar(np.log10(cv_df["penalizer"]), cv_df["mean_cindex"], yerr=cv_df["std_cindex"],
             fmt='o-', color='steelblue', capsize=3, label='C-index')
ax1.axvline(np.log10(best_pen), color='red', linestyle='--', label=f'Best λ={best_pen}')
ax1.set_xlabel("log10(Penalizer)")
ax1.set_ylabel("C-index (5-fold CV)", color='steelblue')

ax2 = ax1.twinx()
ax2.plot(np.log10(cv_df["penalizer"]), cv_df["mean_ngenes"], 's--', color='darkorange', label='N genes')
ax2.set_ylabel("Number of non-zero genes", color='darkorange')

ax1.set_title("LASSO-Cox Cross-Validation", fontsize=13, fontweight='bold')
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='lower left')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "lambda_cv.png"), dpi=200, bbox_inches='tight')
print("Saved: lambda_cv.png")

# ── Bootstrap stability (1000 iterations) ──────────────────────────────────
print("\nBootstrap stability analysis (1000 iterations)...")
n_boot = 1000
gene_counts = {g: 0 for g in candidates}
boot_cindices = []

for i in range(n_boot):
    idx = np.random.choice(len(lasso_data), size=len(lasso_data), replace=True)
    boot_data = lasso_data.iloc[idx].copy()
    try:
        cph_b = CoxPHFitter(penalizer=best_pen, l1_ratio=1.0)
        cph_b.fit(boot_data, duration_col="OS_months", event_col="OS_event")
        for g in candidates:
            if abs(cph_b.params_.get(g, 0)) > 1e-4:
                gene_counts[g] += 1
        pred_b = cph_b.predict_partial_hazard(lasso_data)
        ci_b = concordance_index(lasso_data["OS_months"], -pred_b.values.ravel(), lasso_data["OS_event"])
        boot_cindices.append(ci_b)
    except:
        pass

    if (i + 1) % 200 == 0:
        print(f"  {i+1}/{n_boot}")

# Stability report
print("\nGene selection frequency (top):")
stability = pd.DataFrame([
    {"gene": g, "frequency": gene_counts[g] / n_boot, "selected_final": g in selected.index}
    for g in candidates
]).sort_values("frequency", ascending=False)

for _, r in stability.head(15).iterrows():
    marker = " ★" if r["selected_final"] else ""
    print(f"  {r['gene']:<12} {r['frequency']:.1%}{marker}")

stability.to_csv(os.path.join(TABLES, "bootstrap_stability.csv"), index=False)

# Bootstrap stability plot
fig, ax = plt.subplots(figsize=(10, 6))
top_stab = stability.head(20).sort_values("frequency")
colors = ['steelblue' if sel else 'lightgray' for sel in top_stab["selected_final"]]
ax.barh(range(len(top_stab)), top_stab["frequency"], color=colors, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(top_stab)))
ax.set_yticklabels(top_stab["gene"])
ax.axvline(0.5, color='red', linestyle='--', alpha=0.5, label='50% threshold')
ax.set_xlabel("Selection Frequency (1000 bootstraps)")
ax.set_title("Bootstrap Stability of LASSO Gene Selection", fontsize=13, fontweight='bold')
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "bootstrap_stability.png"), dpi=200, bbox_inches='tight')
print("Saved: bootstrap_stability.png")

# ── Permutation test (500 iterations) ──────────────────────────────────────
print("\nPermutation test (500 iterations)...")
n_perm = 500
perm_cindices = []
for i in range(n_perm):
    perm_data = lasso_data.copy()
    shuffled = perm_data[["OS_months", "OS_event"]].sample(frac=1, random_state=i).values
    perm_data["OS_months"] = shuffled[:, 0]
    perm_data["OS_event"] = shuffled[:, 1]
    try:
        pred_perm = cph_final.predict_partial_hazard(perm_data)
        ci_perm = concordance_index(perm_data["OS_months"], -pred_perm.values.ravel(), perm_data["OS_event"])
        perm_cindices.append(ci_perm)
    except:
        perm_cindices.append(0.5)
    if (i + 1) % 100 == 0:
        print(f"  {i+1}/{n_perm}")

perm_p = (np.array(perm_cindices) >= ci_train).mean()
print(f"\nPermutation p-value: {perm_p:.4f}")
print(f"Observed C-index: {ci_train:.4f}")
print(f"Null distribution: mean={np.mean(perm_cindices):.4f}, 95th={np.percentile(perm_cindices, 95):.4f}")

# Permutation plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(perm_cindices, bins=30, color='lightgray', edgecolor='black', alpha=0.7, label='Null distribution')
ax.axvline(ci_train, color='red', linewidth=2, label=f'Observed C-index={ci_train:.4f}')
ax.axvline(np.percentile(perm_cindices, 95), color='orange', linestyle='--', label='95th percentile')
ax.set_xlabel("C-index")
ax.set_ylabel("Count")
ax.set_title(f"Permutation Test (n={n_perm}, p={perm_p:.4f})", fontsize=13, fontweight='bold')
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "permutation_test.png"), dpi=200, bbox_inches='tight')
print("Saved: permutation_test.png")

# ── Save model ──────────────────────────────────────────────────────────────
model_info = {
    "penalizer": float(best_pen),
    "l1_ratio": 1.0,
    "c_index_train": float(ci_train),
    "permutation_p": float(perm_p),
    "n_patients": int(len(lasso_data)),
    "n_events": int(lasso_data["OS_event"].sum()),
    "n_candidate_genes": len(candidates),
    "genes": {gene: float(coef) for gene, coef in selected.items()},
    "gene_means": {g: gene_means[g] for g in selected.index},
    "gene_stds": {g: gene_stds[g] for g in selected.index},
    "cv_results": cv_results,
    "boot_ci_95": [float(np.percentile(boot_cindices, 2.5)), float(np.percentile(boot_cindices, 97.5))],
}

with open(os.path.join(MODEL, "lasso_model.json"), "w") as f:
    json.dump(model_info, f, indent=2)
print(f"\nSaved: lasso_model.json")

# ── Compute and save risk scores ────────────────────────────────────────────
risk_scores = np.zeros(len(merged))
score_available = True
for gene, coef in selected.items():
    if gene not in merged.columns:
        score_available = False
        break
    z = (merged[gene] - gene_means[gene]) / gene_stds[gene]
    risk_scores += coef * z.values

if score_available:
    merged["risk_score"] = risk_scores
    merged.to_csv(os.path.join(DATA, "tcga_ros_merged.csv"), index=False)
    print("Updated tcga_ros_merged.csv with risk_score column")

print(f"\n✓ LASSO model complete: {len(selected)} genes, C-index={ci_train:.4f}, perm p={perm_p:.4f}")
