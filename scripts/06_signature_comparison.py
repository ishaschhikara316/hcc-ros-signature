"""
06_signature_comparison.py — Head-to-head comparison with published signatures

Compares our ROS/ferroptosis signature against:
1. Hong 8-gene oxidative stress signature
2. Buffa 15-gene hypoxia signature
3. MKI67 4-gene proliferation signature
"""
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
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
GEO = os.path.join(BASE, "data", "geo_cohorts")
MODEL = os.path.join(BASE, "results", "model")
TABLES = os.path.join(BASE, "results", "tables")
FIGS = os.path.join(BASE, "results", "figures")

# ── Load model ──────────────────────────────────────────────────────────────
with open(os.path.join(MODEL, "lasso_model.json")) as f:
    model = json.load(f)
our_genes = model["genes"]
our_means = model["gene_means"]
our_stds = model["gene_stds"]

# ── Define competing signatures ─────────────────────────────────────────────
COMPETING = {
    "Hong 8-gene OS": {
        "genes": ["G6PD", "MT3", "CBX2", "CDKN2B", "CCNA2", "MAPT", "EZH2", "SLC7A11"],
        "reference": "Hong et al. (2021)"
    },
    "Buffa hypoxia (15-gene)": {
        "genes": ["VEGFA", "SLC2A1", "PGAM1", "ENO1", "LDHA", "TPI1", "ALDOA",
                   "PGK1", "GPI", "P4HA1", "MRPS17", "CDKN3", "ADM", "NDRG1", "TUBB6"],
        "reference": "Buffa et al. (2010)"
    },
    "MKI67 proliferation (4-gene)": {
        "genes": ["MKI67", "TOP2A", "PCNA", "MCM2"],
        "reference": "Proliferation markers"
    },
}


def evaluate_signature(name, genes, expr_df, os_months, os_event, is_our=False, our_coefs=None):
    """Evaluate a gene signature on a cohort. Returns dict of metrics."""
    available = [g for g in genes if g in expr_df.columns]
    if len(available) < 2:
        print(f"  {name}: Only {len(available)}/{len(genes)} genes available — skipping")
        return None

    df = expr_df[available + ["OS_months", "OS_event"]].dropna()
    df = df[df["OS_months"] > 0]

    if len(df) < 30:
        return None

    # Compute score
    if is_our and our_coefs:
        score = np.zeros(len(df))
        for g in available:
            if g in our_coefs:
                # Use cohort-level z-scores (works across platforms)
                z = (df[g] - df[g].mean()) / (df[g].std() + 1e-10)
                score += our_coefs[g] * z
    else:
        # Unweighted z-score mean for competing signatures
        score = np.zeros(len(df))
        for g in available:
            z = (df[g] - df[g].mean()) / (df[g].std() + 1e-10)
            score += z
        score /= len(available)

    # C-index
    ci = concordance_index(df["OS_months"], -score, df["OS_event"])

    # Bootstrap CI
    boot = []
    for _ in range(500):
        idx = np.random.choice(len(df), len(df), replace=True)
        bd = df.iloc[idx]
        try:
            boot.append(concordance_index(bd["OS_months"], -score[idx], bd["OS_event"]))
        except:
            pass
    ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5]) if boot else (np.nan, np.nan)

    # KM
    med = np.median(score)
    high_idx = score >= med
    h = df[high_idx]
    l = df[~high_idx]
    lr = logrank_test(h["OS_months"], l["OS_months"],
                      event_observed_A=h["OS_event"], event_observed_B=l["OS_event"])

    # Cox
    cox_df = df[["OS_months", "OS_event"]].copy()
    cox_df["score"] = score
    try:
        cph = CoxPHFitter()
        cph.fit(cox_df, duration_col="OS_months", event_col="OS_event")
        hr = np.exp(cph.params_["score"])
        hr_ci = np.exp(cph.confidence_intervals_.values[0])
    except:
        hr, hr_ci = np.nan, [np.nan, np.nan]

    return {
        "signature": name, "n_genes_total": len(genes), "n_genes_available": len(available),
        "n": len(df), "events": int(df["OS_event"].sum()),
        "c_index": ci, "ci_lo": ci_lo, "ci_hi": ci_hi,
        "HR": hr, "HR_lower": hr_ci[0], "HR_upper": hr_ci[1],
        "logrank_p": lr.p_value,
    }


# ── Load TCGA expression ───────────────────────────────────────────────────
merged = pd.read_csv(os.path.join(DATA, "tcga_ros_merged.csv"))
# Also load full expression to get competing signature genes
expr_full = pd.read_csv(os.path.join(DATA, "tcga_lihc_expression_full.csv"), index_col=0)

# Add any missing competing genes to merged
all_competing_genes = set()
for sig in COMPETING.values():
    all_competing_genes.update(sig["genes"])

for gene in all_competing_genes:
    if gene not in merged.columns and gene in expr_full.index:
        merged[gene] = expr_full.loc[gene, merged["patientId"]].values

# ══════════════════════════════════════════════════════════════════════════════
# EVALUATE ON TCGA (TRAINING)
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("SIGNATURE COMPARISON — TCGA-LIHC (Training)")
print("=" * 70)

results_train = []

# Our signature
our_result = evaluate_signature(
    f"Our ROS/Ferroptosis ({len(our_genes)}-gene)", list(our_genes.keys()),
    merged, merged["OS_months"], merged["OS_event"],
    is_our=True, our_coefs=our_genes
)
if our_result:
    results_train.append(our_result)
    print(f"  Our signature: C={our_result['c_index']:.3f}, HR={our_result['HR']:.2f}, p={our_result['logrank_p']:.2e}")

# Competing signatures
for name, info in COMPETING.items():
    result = evaluate_signature(name, info["genes"], merged, merged["OS_months"], merged["OS_event"])
    if result:
        results_train.append(result)
        print(f"  {name}: C={result['c_index']:.3f}, HR={result['HR']:.2f}, p={result['logrank_p']:.2e}")

if results_train:
    comp_df = pd.DataFrame(results_train)
    comp_df.to_csv(os.path.join(TABLES, "signature_comparison.csv"), index=False)
    print("\nSaved: signature_comparison.csv")
    print(comp_df[["signature", "c_index", "HR", "logrank_p"]].to_string(index=False))


# ══════════════════════════════════════════════════════════════════════════════
# EVALUATE ON GSE14520 (VALIDATION)
# ══════════════════════════════════════════════════════════════════════════════
gse14520_cache_path = os.path.join(GEO, "gse14520_gene_expr_cache.csv")
suppl_path = os.path.join(GEO, "GSE14520_Extra_Supplement.txt.gz")
results_val = []

if os.path.exists(gse14520_cache_path) and os.path.exists(suppl_path):
    print("\n" + "=" * 70)
    print("SIGNATURE COMPARISON — GSE14520 (Validation)")
    print("=" * 70)

    # Load cached GSE14520 gene expression (gene x sample)
    expr14_genes = pd.read_csv(gse14520_cache_path, index_col=0)

    # Load supplement for survival
    suppl14 = pd.read_csv(suppl_path, sep='\t', compression='gzip')
    surv_time_col = surv_event_col = gsm_col = None
    for col in suppl14.columns:
        cl = col.lower()
        if 'survival' in cl and 'month' in cl: surv_time_col = col
        elif 'survival' in cl and 'status' in cl: surv_event_col = col
        elif cl == 'affy_gsm': gsm_col = col
    if gsm_col is None:
        for col in suppl14.columns:
            if 'gsm' in col.lower():
                gsm_col = col
                break

    if surv_time_col and surv_event_col and gsm_col:
        suppl14 = suppl14.dropna(subset=[surv_time_col, surv_event_col, gsm_col])
        suppl14[gsm_col] = suppl14[gsm_col].astype(str)
        suppl14 = suppl14.set_index(gsm_col)
        matched14 = [s for s in suppl14.index if s in expr14_genes.columns]

        if len(matched14) >= 30:
            # Build patient-level DataFrame for GSE14520
            val_expr = pd.DataFrame(index=matched14)
            for gene in expr14_genes.index:
                if gene in expr14_genes.index:
                    val_expr[gene] = expr14_genes.loc[gene, matched14].values
            val_expr["OS_months"] = suppl14.loc[matched14, surv_time_col].values.astype(float)
            val_expr["OS_event"] = suppl14.loc[matched14, surv_event_col].values.astype(float)

            # Our signature on validation
            our_val = evaluate_signature(
                f"Our ROS/Ferroptosis ({len(our_genes)}-gene)", list(our_genes.keys()),
                val_expr, val_expr["OS_months"], val_expr["OS_event"],
                is_our=True, our_coefs=our_genes
            )
            if our_val:
                results_val.append(our_val)
                print(f"  Our signature: C={our_val['c_index']:.3f}, HR={our_val['HR']:.2f}, p={our_val['logrank_p']:.2e}")

            # Competing signatures on validation
            for name, info in COMPETING.items():
                result = evaluate_signature(name, info["genes"], val_expr,
                                            val_expr["OS_months"], val_expr["OS_event"])
                if result:
                    results_val.append(result)
                    print(f"  {name}: C={result['c_index']:.3f}, HR={result['HR']:.2f}, p={result['logrank_p']:.2e}")

            if results_val:
                comp_val_df = pd.DataFrame(results_val)
                comp_val_df.to_csv(os.path.join(TABLES, "signature_comparison_validation.csv"), index=False)
                print("\nSaved: signature_comparison_validation.csv")
                print(comp_val_df[["signature", "c_index", "HR", "logrank_p"]].to_string(index=False))
else:
    print("\n  GSE14520 expression cache not found — skipping validation comparison")
    print("  (Run 05_external_validation.py first to generate cache)")

# ══════════════════════════════════════════════════════════════════════════════
# COMPARISON VISUALIZATION
# ══════════════════════════════════════════════════════════════════════════════
if results_train:
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # 1. C-index comparison
    ax = axes[0]
    names = [r["signature"].split("(")[0].strip() for r in results_train]
    cidx = [r["c_index"] for r in results_train]
    ci_err = [[r["c_index"] - r["ci_lo"] for r in results_train],
              [r["ci_hi"] - r["c_index"] for r in results_train]]
    colors = ['steelblue'] + ['lightcoral'] * (len(results_train) - 1)
    ax.barh(range(len(names)), cidx, xerr=ci_err, color=colors, edgecolor='black', capsize=3)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("C-index")
    ax.set_title("C-index Comparison", fontweight='bold')
    ax.axvline(0.5, color='gray', linestyle='--', alpha=0.5)

    # 2. HR comparison
    ax = axes[1]
    hrs = [r["HR"] for r in results_train]
    hr_err = [[r["HR"] - r["HR_lower"] for r in results_train],
              [r["HR_upper"] - r["HR"] for r in results_train]]
    ax.barh(range(len(names)), hrs, xerr=hr_err, color=colors, edgecolor='black', capsize=3)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("Hazard Ratio")
    ax.set_title("Hazard Ratio Comparison", fontweight='bold')
    ax.axvline(1.0, color='black', linestyle='--', alpha=0.5)

    # 3. -log10(p) comparison
    ax = axes[2]
    logp = [-np.log10(r["logrank_p"] + 1e-30) for r in results_train]
    ax.barh(range(len(names)), logp, color=colors, edgecolor='black')
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("-log10(p-value)")
    ax.set_title("Log-rank Significance", fontweight='bold')
    ax.axvline(-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
    ax.legend()

    plt.suptitle("Signature Comparison — TCGA-LIHC (Training)", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "signature_comparison.png"), dpi=200, bbox_inches='tight')
    print("\nSaved: signature_comparison.png")

# Validation comparison plot
if results_val:
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    names_v = [r["signature"].split("(")[0].strip() for r in results_val]
    cidx_v = [r["c_index"] for r in results_val]
    ci_err_v = [[r["c_index"] - r["ci_lo"] for r in results_val],
                [r["ci_hi"] - r["c_index"] for r in results_val]]
    colors_v = ['steelblue'] + ['lightcoral'] * (len(results_val) - 1)

    ax = axes[0]
    ax.barh(range(len(names_v)), cidx_v, xerr=ci_err_v, color=colors_v, edgecolor='black', capsize=3)
    ax.set_yticks(range(len(names_v)))
    ax.set_yticklabels(names_v, fontsize=9)
    ax.set_xlabel("C-index")
    ax.set_title("C-index Comparison", fontweight='bold')
    ax.axvline(0.5, color='gray', linestyle='--', alpha=0.5)

    ax = axes[1]
    hrs_v = [r["HR"] for r in results_val]
    hr_err_v = [[r["HR"] - r["HR_lower"] for r in results_val],
                [r["HR_upper"] - r["HR"] for r in results_val]]
    ax.barh(range(len(names_v)), hrs_v, xerr=hr_err_v, color=colors_v, edgecolor='black', capsize=3)
    ax.set_yticks(range(len(names_v)))
    ax.set_yticklabels(names_v, fontsize=9)
    ax.set_xlabel("Hazard Ratio")
    ax.set_title("Hazard Ratio Comparison", fontweight='bold')
    ax.axvline(1.0, color='black', linestyle='--', alpha=0.5)

    ax = axes[2]
    logp_v = [-np.log10(r["logrank_p"] + 1e-30) for r in results_val]
    ax.barh(range(len(names_v)), logp_v, color=colors_v, edgecolor='black')
    ax.set_yticks(range(len(names_v)))
    ax.set_yticklabels(names_v, fontsize=9)
    ax.set_xlabel("-log10(p-value)")
    ax.set_title("Log-rank Significance", fontweight='bold')
    ax.axvline(-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
    ax.legend()

    plt.suptitle("Signature Comparison — GSE14520 (Validation)", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "signature_comparison_validation.png"), dpi=200, bbox_inches='tight')
    print("Saved: signature_comparison_validation.png")

print(f"\n✓ Signature comparison complete.")
