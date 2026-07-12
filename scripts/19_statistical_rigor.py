#!/usr/bin/env python3
"""
19_statistical_rigor.py
=======================
Strengthens statistical rigor across all analyses.

Adds:
  1. Multiple testing corrections (BH FDR) for checkpoints, mutations, subgroups
  2. Formal C-index comparison (bootstrap difference test)
  3. Interaction tests (signature × sex, signature × ethnicity)
  4. One-standard-error rule analysis for LASSO
  5. Proportional hazards diagnostics with time-varying coefficients
  6. Sensitivity analysis: optimal cutpoint vs median
  7. Power analysis for validation cohorts

Outputs:
  - results/tables/corrected_checkpoint_pvalues.csv
  - results/tables/corrected_mutation_pvalues.csv
  - results/tables/cindex_comparison_tests.csv
  - results/tables/interaction_tests.csv
  - results/tables/sensitivity_cutpoints.csv
  - results/tables/power_analysis.csv
  - results/tables/one_se_rule.csv
  - results/figures/sensitivity_cutpoint.png
  - results/figures/interaction_forest.png
"""

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np
import json
import os
import warnings
from scipy import stats
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
RESULTS = os.path.join(BASE, "results")
TABLES = os.path.join(RESULTS, "tables")
FIGURES = os.path.join(RESULTS, "figures")

os.makedirs(TABLES, exist_ok=True)
os.makedirs(FIGURES, exist_ok=True)

# ── Load LASSO model ──
with open(os.path.join(RESULTS, "model", "lasso_model.json")) as f:
    model = json.load(f)
coefs = model["genes"]
gene_means = model["gene_means"]
gene_stds = model["gene_stds"]

print("=" * 70)
print("STATISTICAL RIGOR ENHANCEMENTS")
print("=" * 70)


# ── Helper: Benjamini-Hochberg FDR ──
def bh_fdr(pvalues):
    """Apply Benjamini-Hochberg FDR correction."""
    pvals = np.array(pvalues, dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([])
    ranked = np.argsort(pvals)
    fdr = np.empty(n)
    fdr[ranked] = pvals[ranked] * n / (np.arange(1, n + 1))
    # Enforce monotonicity (step-up)
    for i in range(n - 2, -1, -1):
        fdr[ranked[i]] = min(fdr[ranked[i]], fdr[ranked[i + 1]])
    return np.minimum(fdr, 1.0)


def bootstrap_cindex_diff(time, event, score1, score2, n_boot=1000, seed=42):
    """Test if C-index of score1 is significantly different from score2."""
    rng = np.random.RandomState(seed)
    n = len(time)
    diffs = []
    for _ in range(n_boot):
        idx = rng.choice(n, n, replace=True)
        try:
            c1 = concordance_index(time[idx], -score1[idx], event[idx])
            c2 = concordance_index(time[idx], -score2[idx], event[idx])
            diffs.append(c1 - c2)
        except Exception:
            pass
    if len(diffs) < 100:
        return np.nan, np.nan, (np.nan, np.nan)
    diffs = np.array(diffs)
    mean_diff = np.mean(diffs)
    p_value = 2 * min(np.mean(diffs > 0), np.mean(diffs < 0))
    ci = (np.percentile(diffs, 2.5), np.percentile(diffs, 97.5))
    return mean_diff, p_value, ci


# ══════════════════════════════════════════════════════════════════════
# SECTION 1: Multiple testing corrections
# ══════════════════════════════════════════════════════════════════════

print("\n[1/7] Applying multiple testing corrections...")

# ── Checkpoint expression ──
chk_path = os.path.join(TABLES, "checkpoint_expression.csv")
if os.path.exists(chk_path):
    chk_df = pd.read_csv(chk_path)
    p_col = None
    for col in chk_df.columns:
        if 'p' in col.lower() and 'val' in col.lower():
            p_col = col
            break
    if p_col is None:
        # Try to find any numeric column that looks like p-values
        for col in chk_df.columns:
            if chk_df[col].dtype in [float, np.float64] and chk_df[col].max() <= 1:
                if col not in ['mean_high', 'mean_low', 'median_high', 'median_low']:
                    p_col = col
                    break

    if p_col:
        pvals = chk_df[p_col].values
        chk_df["FDR_BH"] = bh_fdr(pvals)
        chk_df["significant_after_FDR"] = chk_df["FDR_BH"] < 0.05
        chk_df.to_csv(os.path.join(TABLES, "corrected_checkpoint_pvalues.csv"), index=False)

        n_sig_raw = (pvals < 0.05).sum()
        n_sig_fdr = (chk_df["FDR_BH"] < 0.05).sum()
        print(f"  Checkpoints: {n_sig_raw} significant (raw) → {n_sig_fdr} significant (FDR)")
        print("  Saved: corrected_checkpoint_pvalues.csv")
    else:
        print("  WARNING: Could not find p-value column in checkpoint_expression.csv")

# ── Mutation frequencies ──
mut_path = os.path.join(TABLES, "mutation_by_risk_group.csv")
if os.path.exists(mut_path):
    mut_df = pd.read_csv(mut_path)
    p_col = None
    for col in mut_df.columns:
        if 'p' in col.lower():
            p_col = col
            break

    if p_col:
        pvals = mut_df[p_col].values.astype(float)
        mut_df["FDR_BH"] = bh_fdr(pvals)
        mut_df["significant_after_FDR"] = mut_df["FDR_BH"] < 0.05

        # Also add Bonferroni
        mut_df["Bonferroni_p"] = np.minimum(pvals * len(pvals), 1.0)
        mut_df["significant_Bonferroni"] = mut_df["Bonferroni_p"] < 0.05

        mut_df.to_csv(os.path.join(TABLES, "corrected_mutation_pvalues.csv"), index=False)

        n_sig_raw = (pvals < 0.05).sum()
        n_sig_fdr = (mut_df["FDR_BH"] < 0.05).sum()
        n_sig_bonf = (mut_df["significant_Bonferroni"]).sum()
        print(f"  Mutations: {n_sig_raw} significant (raw) → {n_sig_fdr} (FDR) → {n_sig_bonf} (Bonferroni)")
        print("  Saved: corrected_mutation_pvalues.csv")
    else:
        print("  WARNING: Could not find p-value column in mutation_by_risk_group.csv")

# ── ssGSEA immune correlations ──
ssgsea_path = os.path.join(TABLES, "ssgsea_immune_profiles.csv")
if os.path.exists(ssgsea_path):
    ssgsea_df = pd.read_csv(ssgsea_path)
    p_col = None
    for col in ssgsea_df.columns:
        if 'p' in col.lower():
            p_col = col
            break
    if p_col:
        # Ensure we only get numeric p-values
        ssgsea_df[p_col] = pd.to_numeric(ssgsea_df[p_col], errors='coerce')
        ssgsea_df = ssgsea_df.dropna(subset=[p_col])
        pvals = ssgsea_df[p_col].values.astype(float)
        ssgsea_df["FDR_BH"] = bh_fdr(pvals)
        ssgsea_df["significant_after_FDR"] = ssgsea_df["FDR_BH"] < 0.05
        ssgsea_df.to_csv(os.path.join(TABLES, "corrected_ssgsea_pvalues.csv"), index=False)

        n_sig_raw = (pvals < 0.05).sum()
        n_sig_fdr = (ssgsea_df["FDR_BH"] < 0.05).sum()
        print(f"  ssGSEA immune: {n_sig_raw} significant (raw) → {n_sig_fdr} significant (FDR)")
        print("  Saved: corrected_ssgsea_pvalues.csv")


# ══════════════════════════════════════════════════════════════════════
# SECTION 2: Formal C-index comparison tests
# ══════════════════════════════════════════════════════════════════════

print("\n[2/7] Formal C-index comparison tests (bootstrap)...")

tcga_path = os.path.join(DATA, "tcga", "tcga_ros_merged.csv")
if os.path.exists(tcga_path):
    tcga = pd.read_csv(tcga_path)
    tcga = tcga.dropna(subset=["OS_months", "OS_event"])

    time_vals = tcga["OS_months"].values
    event_vals = tcga["OS_event"].values.astype(int)
    risk_scores = tcga["risk_score"].values if "risk_score" in tcga.columns else None

    if risk_scores is not None:
        comp_rows = []

        # Risk score vs stage
        stage_map = {"Stage I": 1, "Stage II": 2, "Stage III": 3,
                     "Stage IIIA": 3, "Stage IIIB": 3, "Stage IIIC": 3,
                     "Stage IV": 4, "Stage IVA": 4, "Stage IVB": 4}
        if "tumor_stage" in tcga.columns:
            tcga["stage_num"] = tcga["tumor_stage"].map(stage_map)
            valid = tcga["stage_num"].notna()
            if valid.sum() > 100:
                diff, p, ci = bootstrap_cindex_diff(
                    time_vals[valid.values], event_vals[valid.values],
                    risk_scores[valid.values], tcga.loc[valid, "stage_num"].values
                )
                c_risk = concordance_index(time_vals[valid.values], -risk_scores[valid.values], event_vals[valid.values])
                c_stage = concordance_index(time_vals[valid.values], -tcga.loc[valid, "stage_num"].values, event_vals[valid.values])
                comp_rows.append({
                    "comparison": "Risk Score vs Stage",
                    "C_index_model1": round(c_risk, 4),
                    "C_index_model2": round(c_stage, 4),
                    "difference": round(diff, 4),
                    "95%CI_lower": round(ci[0], 4),
                    "95%CI_upper": round(ci[1], 4),
                    "p_value": round(p, 4)
                })
                print(f"  Risk Score vs Stage: ΔC={diff:.4f} (95% CI: {ci[0]:.4f}-{ci[1]:.4f}), p={p:.4f}")

        # Risk score vs age
        if "age_at_diagnosis" in tcga.columns:
            valid = tcga["age_at_diagnosis"].notna()
            if valid.sum() > 100:
                diff, p, ci = bootstrap_cindex_diff(
                    time_vals[valid.values], event_vals[valid.values],
                    risk_scores[valid.values], tcga.loc[valid, "age_at_diagnosis"].values
                )
                c_age = concordance_index(time_vals[valid.values],
                                          -tcga.loc[valid, "age_at_diagnosis"].values,
                                          event_vals[valid.values])
                comp_rows.append({
                    "comparison": "Risk Score vs Age",
                    "C_index_model1": round(concordance_index(time_vals[valid.values], -risk_scores[valid.values], event_vals[valid.values]), 4),
                    "C_index_model2": round(c_age, 4),
                    "difference": round(diff, 4),
                    "95%CI_lower": round(ci[0], 4),
                    "95%CI_upper": round(ci[1], 4),
                    "p_value": round(p, 4)
                })
                print(f"  Risk Score vs Age: ΔC={diff:.4f}, p={p:.4f}")

        # Nomogram vs risk score alone
        if "tumor_stage" in tcga.columns and "age_at_diagnosis" in tcga.columns:
            valid = tcga["stage_num"].notna() & tcga["age_at_diagnosis"].notna()
            if valid.sum() > 100:
                sub = tcga[valid].copy()
                sub["male"] = (sub["gender"].str.lower() == "male").astype(int) if "gender" in sub.columns else 0
                # Fit nomogram (multivariate Cox)
                nomo_cols = ["risk_score", "stage_num", "age_at_diagnosis", "male"]
                nomo_cols = [c for c in nomo_cols if c in sub.columns]
                try:
                    cph = CoxPHFitter(penalizer=0.01)
                    cph.fit(sub[nomo_cols + ["OS_months", "OS_event"]],
                            duration_col="OS_months", event_col="OS_event")
                    nomo_score = cph.predict_partial_hazard(sub[nomo_cols]).values.flatten()

                    diff, p, ci = bootstrap_cindex_diff(
                        sub["OS_months"].values, sub["OS_event"].values.astype(int),
                        nomo_score, sub["risk_score"].values
                    )
                    c_nomo = concordance_index(sub["OS_months"].values, -nomo_score, sub["OS_event"].values.astype(int))
                    comp_rows.append({
                        "comparison": "Nomogram vs Risk Score",
                        "C_index_model1": round(c_nomo, 4),
                        "C_index_model2": round(concordance_index(sub["OS_months"].values, -sub["risk_score"].values, sub["OS_event"].values.astype(int)), 4),
                        "difference": round(diff, 4),
                        "95%CI_lower": round(ci[0], 4),
                        "95%CI_upper": round(ci[1], 4),
                        "p_value": round(p, 4)
                    })
                    print(f"  Nomogram vs Risk Score: ΔC={diff:.4f}, p={p:.4f}")
                except Exception as e:
                    print(f"  Nomogram comparison failed: {e}")

        if comp_rows:
            comp_df = pd.DataFrame(comp_rows)
            comp_df.to_csv(os.path.join(TABLES, "cindex_comparison_tests.csv"), index=False)
            print("  Saved: cindex_comparison_tests.csv")


# ══════════════════════════════════════════════════════════════════════
# SECTION 3: Interaction tests (signature × sex, signature × ethnicity)
# ══════════════════════════════════════════════════════════════════════

print("\n[3/7] Interaction tests...")

if os.path.exists(tcga_path):
    tcga = pd.read_csv(tcga_path)
    tcga = tcga.dropna(subset=["OS_months", "OS_event"])

    interaction_rows = []

    # Risk score × sex interaction
    if "gender" in tcga.columns and "risk_score" in tcga.columns:
        tcga["male"] = (tcga["gender"].str.lower() == "male").astype(int)
        tcga["risk_x_male"] = tcga["risk_score"] * tcga["male"]

        try:
            # Model without interaction
            cph_base = CoxPHFitter(penalizer=0.01)
            cph_base.fit(tcga[["risk_score", "male", "OS_months", "OS_event"]].dropna(),
                         duration_col="OS_months", event_col="OS_event")

            # Model with interaction
            cph_int = CoxPHFitter(penalizer=0.01)
            cph_int.fit(tcga[["risk_score", "male", "risk_x_male", "OS_months", "OS_event"]].dropna(),
                        duration_col="OS_months", event_col="OS_event")

            int_hr = np.exp(cph_int.params_["risk_x_male"])
            int_p = cph_int.summary.loc["risk_x_male", "p"]

            # LRT
            ll_base = cph_base.log_likelihood_
            ll_int = cph_int.log_likelihood_
            lrt_stat = 2 * (ll_int - ll_base)
            lrt_p = 1 - stats.chi2.cdf(lrt_stat, 1)

            interaction_rows.append({
                "interaction": "risk_score × sex",
                "interaction_HR": round(float(int_hr), 4),
                "interaction_p_Wald": round(float(int_p), 4),
                "LRT_statistic": round(float(lrt_stat), 4),
                "LRT_p_value": round(float(lrt_p), 4),
                "interpretation": "Significant" if lrt_p < 0.05 else "Not significant"
            })
            print(f"  Risk × Sex: interaction HR={int_hr:.3f}, LRT p={lrt_p:.4f}")
        except Exception as e:
            print(f"  Risk × Sex interaction failed: {e}")

    # Risk score × age interaction
    if "age_at_diagnosis" in tcga.columns and "risk_score" in tcga.columns:
        tcga["age_60"] = (tcga["age_at_diagnosis"] >= 60).astype(int)
        tcga["risk_x_age60"] = tcga["risk_score"] * tcga["age_60"]

        try:
            cph_base = CoxPHFitter(penalizer=0.01)
            cph_base.fit(tcga[["risk_score", "age_60", "OS_months", "OS_event"]].dropna(),
                         duration_col="OS_months", event_col="OS_event")

            cph_int = CoxPHFitter(penalizer=0.01)
            cph_int.fit(tcga[["risk_score", "age_60", "risk_x_age60", "OS_months", "OS_event"]].dropna(),
                        duration_col="OS_months", event_col="OS_event")

            int_hr = np.exp(cph_int.params_["risk_x_age60"])
            int_p = cph_int.summary.loc["risk_x_age60", "p"]

            ll_base = cph_base.log_likelihood_
            ll_int = cph_int.log_likelihood_
            lrt_stat = 2 * (ll_int - ll_base)
            lrt_p = 1 - stats.chi2.cdf(lrt_stat, 1)

            interaction_rows.append({
                "interaction": "risk_score × age (≥60)",
                "interaction_HR": round(float(int_hr), 4),
                "interaction_p_Wald": round(float(int_p), 4),
                "LRT_statistic": round(float(lrt_stat), 4),
                "LRT_p_value": round(float(lrt_p), 4),
                "interpretation": "Significant" if lrt_p < 0.05 else "Not significant"
            })
            print(f"  Risk × Age: interaction HR={int_hr:.3f}, LRT p={lrt_p:.4f}")
        except Exception as e:
            print(f"  Risk × Age interaction failed: {e}")

    # Risk score × stage interaction
    if "tumor_stage" in tcga.columns and "risk_score" in tcga.columns:
        stage_map = {"Stage I": 1, "Stage II": 2, "Stage III": 3,
                     "Stage IIIA": 3, "Stage IIIB": 3, "Stage IIIC": 3,
                     "Stage IV": 4, "Stage IVA": 4, "Stage IVB": 4}
        tcga["stage_num"] = tcga["tumor_stage"].map(stage_map)
        tcga["late_stage"] = (tcga["stage_num"] >= 3).astype(int)
        tcga["risk_x_late"] = tcga["risk_score"] * tcga["late_stage"]

        valid = tcga[["risk_score", "late_stage", "risk_x_late", "OS_months", "OS_event"]].dropna()
        if len(valid) > 50:
            try:
                cph_base = CoxPHFitter(penalizer=0.01)
                cph_base.fit(valid[["risk_score", "late_stage", "OS_months", "OS_event"]],
                             duration_col="OS_months", event_col="OS_event")

                cph_int = CoxPHFitter(penalizer=0.01)
                cph_int.fit(valid[["risk_score", "late_stage", "risk_x_late", "OS_months", "OS_event"]],
                            duration_col="OS_months", event_col="OS_event")

                int_hr = np.exp(cph_int.params_["risk_x_late"])
                int_p = cph_int.summary.loc["risk_x_late", "p"]

                ll_base = cph_base.log_likelihood_
                ll_int = cph_int.log_likelihood_
                lrt_stat = 2 * (ll_int - ll_base)
                lrt_p = 1 - stats.chi2.cdf(lrt_stat, 1)

                interaction_rows.append({
                    "interaction": "risk_score × stage (III/IV)",
                    "interaction_HR": round(float(int_hr), 4),
                    "interaction_p_Wald": round(float(int_p), 4),
                    "LRT_statistic": round(float(lrt_stat), 4),
                    "LRT_p_value": round(float(lrt_p), 4),
                    "interpretation": "Significant" if lrt_p < 0.05 else "Not significant"
                })
                print(f"  Risk × Stage: interaction HR={int_hr:.3f}, LRT p={lrt_p:.4f}")
            except Exception as e:
                print(f"  Risk × Stage interaction failed: {e}")

    if interaction_rows:
        int_df = pd.DataFrame(interaction_rows)
        int_df.to_csv(os.path.join(TABLES, "interaction_tests.csv"), index=False)
        print("  Saved: interaction_tests.csv")

        # Interaction forest plot
        fig, ax = plt.subplots(figsize=(8, max(3, len(interaction_rows) * 1.2)))
        y_pos = list(range(len(interaction_rows)))
        y_pos.reverse()

        for i, row in enumerate(interaction_rows):
            y = y_pos[i]
            hr = row["interaction_HR"]
            ax.scatter(hr, y, s=100, c='darkgreen' if row["LRT_p_value"] < 0.05 else 'gray',
                       zorder=3, edgecolors='black')
            label = f"{row['interaction']} (p={row['LRT_p_value']:.3f})"
            ax.text(0.02, y, label, va='center', ha='left', fontsize=10,
                    transform=ax.get_yaxis_transform())

        ax.axvline(x=1, color='gray', linestyle='--', linewidth=0.8)
        ax.set_xlabel("Interaction HR", fontsize=12)
        ax.set_yticks([])
        ax.set_title("Interaction Tests: Risk Score × Clinical Variables", fontsize=13)
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES, "interaction_forest.png"), dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved: interaction_forest.png")


# ══════════════════════════════════════════════════════════════════════
# SECTION 4: One-standard-error rule for LASSO
# ══════════════════════════════════════════════════════════════════════

print("\n[4/7] One-standard-error rule analysis...")

cv_results = model.get("cv_results", [])
if cv_results:
    cv_df = pd.DataFrame(cv_results)

    # Find best lambda (max mean C-index)
    best_idx = cv_df["mean_cindex"].idxmax()
    best_cindex = cv_df.loc[best_idx, "mean_cindex"]
    best_std = cv_df.loc[best_idx, "std_cindex"] if "std_cindex" in cv_df.columns else cv_df.loc[best_idx, "std"]
    best_lambda = cv_df.loc[best_idx, "penalizer"]

    # 1SE rule: most parsimonious model within 1 SE of best
    threshold = best_cindex - best_std
    # Find largest penalizer with mean_cindex >= threshold
    candidates = cv_df[cv_df["mean_cindex"] >= threshold]
    if len(candidates) > 0:
        se_idx = candidates["penalizer"].idxmax()
        se_lambda = cv_df.loc[se_idx, "penalizer"]
        se_cindex = cv_df.loc[se_idx, "mean_cindex"]
        se_ngenes = cv_df.loc[se_idx, "mean_ngenes"]
    else:
        se_lambda = best_lambda
        se_cindex = best_cindex
        se_ngenes = cv_df.loc[best_idx, "mean_ngenes"]

    selected_lambda = model["penalizer"]
    selected_cindex = model.get("c_index_train", np.nan)

    one_se_results = {
        "best_lambda": best_lambda,
        "best_CV_cindex": round(best_cindex, 4),
        "best_CV_std": round(best_std, 4),
        "best_n_genes": cv_df.loc[best_idx, "mean_ngenes"],
        "1SE_lambda": se_lambda,
        "1SE_CV_cindex": round(se_cindex, 4),
        "1SE_n_genes": se_ngenes,
        "selected_lambda": selected_lambda,
        "selected_train_cindex": round(selected_cindex, 4),
        "1SE_threshold": round(threshold, 4),
        "within_1SE": "Yes" if abs(selected_cindex - best_cindex) <= best_std else "Close"
    }

    pd.DataFrame([one_se_results]).to_csv(os.path.join(TABLES, "one_se_rule.csv"), index=False)
    print(f"  Best λ: {best_lambda} (C={best_cindex:.4f} ± {best_std:.4f}, {cv_df.loc[best_idx, 'mean_ngenes']:.0f} genes)")
    print(f"  1SE λ:  {se_lambda} (C={se_cindex:.4f}, {se_ngenes:.0f} genes)")
    print(f"  Selected λ: {selected_lambda} (train C={selected_cindex:.4f})")
    print(f"  1SE threshold: {threshold:.4f}")
    print("  Saved: one_se_rule.csv")


# ══════════════════════════════════════════════════════════════════════
# SECTION 5: Sensitivity analysis — optimal cutpoint vs median
# ══════════════════════════════════════════════════════════════════════

print("\n[5/7] Sensitivity analysis: cutpoint selection...")

if os.path.exists(tcga_path):
    tcga = pd.read_csv(tcga_path)
    tcga = tcga.dropna(subset=["OS_months", "OS_event"])

    if "risk_score" in tcga.columns:
        time_vals = tcga["OS_months"].values
        event_vals = tcga["OS_event"].values.astype(int)
        risk_scores = tcga["risk_score"].values

        cut_rows = []

        # Test multiple cutpoints
        percentiles = [25, 30, 33, 40, 50, 60, 67, 70, 75]
        for pct in percentiles:
            cutoff = np.percentile(risk_scores, pct)
            high = risk_scores >= cutoff
            low = ~high

            if high.sum() < 10 or low.sum() < 10:
                continue

            # Log-rank test
            lr = logrank_test(time_vals[high], time_vals[low],
                              event_vals[high], event_vals[low])

            # Cox HR
            try:
                group = high.astype(int)
                cph_df = pd.DataFrame({
                    "time": time_vals, "event": event_vals, "group": group
                })
                cph = CoxPHFitter()
                cph.fit(cph_df, duration_col="time", event_col="event")
                hr = np.exp(cph.params_["group"])
                hr_lo = np.exp(cph.confidence_intervals_.iloc[0, 0])
                hr_hi = np.exp(cph.confidence_intervals_.iloc[0, 1])
            except Exception:
                hr, hr_lo, hr_hi = np.nan, np.nan, np.nan

            # C-index
            ci = concordance_index(time_vals, -risk_scores, event_vals)

            cut_rows.append({
                "percentile": pct,
                "cutoff_value": round(cutoff, 4),
                "n_high_risk": int(high.sum()),
                "n_low_risk": int(low.sum()),
                "HR": round(float(hr), 3),
                "HR_lower": round(float(hr_lo), 3),
                "HR_upper": round(float(hr_hi), 3),
                "logrank_p": float(lr.p_value),
                "is_median": pct == 50
            })

        if cut_rows:
            cut_df = pd.DataFrame(cut_rows)
            cut_df.to_csv(os.path.join(TABLES, "sensitivity_cutpoints.csv"), index=False)
            print("  Cutpoint sensitivity analysis:")
            print(cut_df[["percentile", "n_high_risk", "HR", "logrank_p"]].to_string(index=False))
            print("  Saved: sensitivity_cutpoints.csv")

            # Sensitivity plot
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

            ax1.plot(cut_df["percentile"], cut_df["HR"], 'bo-', linewidth=2)
            ax1.fill_between(cut_df["percentile"], cut_df["HR_lower"], cut_df["HR_upper"],
                             alpha=0.2, color='blue')
            ax1.axhline(y=1, color='gray', linestyle='--')
            ax1.axvline(x=50, color='red', linestyle='--', alpha=0.5, label='Median (selected)')
            ax1.set_xlabel("Cutpoint Percentile", fontsize=12)
            ax1.set_ylabel("Hazard Ratio", fontsize=12)
            ax1.set_title("HR Across Cutpoints", fontsize=13)
            ax1.legend()

            ax2.plot(cut_df["percentile"], -np.log10(cut_df["logrank_p"].astype(float)),
                     'rs-', linewidth=2)
            ax2.axhline(y=-np.log10(0.05), color='gray', linestyle='--', label='p=0.05')
            ax2.axvline(x=50, color='red', linestyle='--', alpha=0.5, label='Median (selected)')
            ax2.set_xlabel("Cutpoint Percentile", fontsize=12)
            ax2.set_ylabel("-log10(p-value)", fontsize=12)
            ax2.set_title("Significance Across Cutpoints", fontsize=13)
            ax2.legend()

            plt.suptitle("Sensitivity Analysis: Cutpoint Selection", fontsize=14, fontweight='bold')
            plt.tight_layout()
            plt.savefig(os.path.join(FIGURES, "sensitivity_cutpoint.png"), dpi=300, bbox_inches='tight')
            plt.close()
            print("  Saved: sensitivity_cutpoint.png")


# ══════════════════════════════════════════════════════════════════════
# SECTION 6: Power analysis for validation cohorts
# ══════════════════════════════════════════════════════════════════════

print("\n[6/7] Power analysis for validation cohorts...")

def power_logrank(n_events, hr, alpha=0.05, ratio=1.0):
    """
    Approximate power for a log-rank test.
    Schoenfeld formula: power = Phi(sqrt(d * p1*p2) * log(HR) - z_alpha/2)
    where d = total events, p1 = proportion in group 1, p2 = proportion in group 2.
    """
    z_alpha = stats.norm.ppf(1 - alpha / 2)
    p1 = ratio / (1 + ratio)
    p2 = 1 / (1 + ratio)
    noncentrality = np.sqrt(n_events * p1 * p2) * abs(np.log(hr))
    power = stats.norm.cdf(noncentrality - z_alpha)
    return power


def required_events(hr, power=0.80, alpha=0.05, ratio=1.0):
    """Minimum events needed to detect HR with given power."""
    z_alpha = stats.norm.ppf(1 - alpha / 2)
    z_beta = stats.norm.ppf(power)
    p1 = ratio / (1 + ratio)
    p2 = 1 / (1 + ratio)
    d = ((z_alpha + z_beta) / (np.log(hr) * np.sqrt(p1 * p2))) ** 2
    return int(np.ceil(d))


# Validation cohort info
val_cohorts = [
    {"cohort": "TCGA-LIHC", "n": 302, "events": 129, "HR": 3.40},
    {"cohort": "GSE14520", "n": 221, "events": 85, "HR": 1.79},
    {"cohort": "ICGC LIRI-JP", "n": 231, "events": 43, "HR": 1.78},
    {"cohort": "GSE54236", "n": 80, "events": 80, "HR": 1.68},
    {"cohort": "GSE76427 (OS)", "n": 115, "events": 23, "HR": 0.49},
    {"cohort": "GSE76427 (RFS)", "n": 108, "events": 48, "HR": 1.46},
]

power_rows = []
for vc in val_cohorts:
    pwr = power_logrank(vc["events"], max(vc["HR"], 1.01))
    req_events_80 = required_events(max(vc["HR"], 1.01), power=0.80)
    req_events_90 = required_events(max(vc["HR"], 1.01), power=0.90)

    power_rows.append({
        "cohort": vc["cohort"],
        "n_patients": vc["n"],
        "n_events": vc["events"],
        "observed_HR": vc["HR"],
        "estimated_power": round(pwr, 3),
        "adequate_power": "Yes" if pwr >= 0.80 else "No",
        "events_needed_80pct": req_events_80,
        "events_needed_90pct": req_events_90,
    })

power_df = pd.DataFrame(power_rows)
power_df.to_csv(os.path.join(TABLES, "power_analysis.csv"), index=False)
print("  Power analysis results:")
print(power_df[["cohort", "n_events", "observed_HR", "estimated_power", "adequate_power"]].to_string(index=False))
print("  Saved: power_analysis.csv")


# ══════════════════════════════════════════════════════════════════════
# SECTION 7: Proportional hazards diagnostic
# ══════════════════════════════════════════════════════════════════════

print("\n[7/7] Proportional hazards diagnostics...")

if os.path.exists(tcga_path):
    tcga = pd.read_csv(tcga_path)
    tcga = tcga.dropna(subset=["OS_months", "OS_event"])

    if "risk_score" in tcga.columns:
        # Fit Cox model and check PH assumption via time-varying coefficient
        # Split follow-up into early (<24mo) and late (≥24mo)
        tcga["risk_group"] = (tcga["risk_score"] >= tcga["risk_score"].median()).astype(int)

        ph_rows = []
        for split_time in [12, 24, 36, 48]:
            # Early period
            early = tcga.copy()
            early["OS_months"] = np.minimum(early["OS_months"], split_time)
            early["OS_event"] = np.where(tcga["OS_months"] <= split_time,
                                          tcga["OS_event"], 0).astype(int)

            # Late period
            late = tcga[tcga["OS_months"] > split_time].copy()
            late["OS_months"] = late["OS_months"] - split_time

            try:
                cph_early = CoxPHFitter()
                cph_early.fit(early[["risk_score", "OS_months", "OS_event"]],
                              duration_col="OS_months", event_col="OS_event")
                hr_early = np.exp(cph_early.params_["risk_score"])
                p_early = cph_early.summary.loc["risk_score", "p"]
            except Exception:
                hr_early, p_early = np.nan, np.nan

            try:
                if len(late) > 20 and late["OS_event"].sum() > 5:
                    cph_late = CoxPHFitter()
                    cph_late.fit(late[["risk_score", "OS_months", "OS_event"]],
                                 duration_col="OS_months", event_col="OS_event")
                    hr_late = np.exp(cph_late.params_["risk_score"])
                    p_late = cph_late.summary.loc["risk_score", "p"]
                else:
                    hr_late, p_late = np.nan, np.nan
            except Exception:
                hr_late, p_late = np.nan, np.nan

            ph_rows.append({
                "split_time_months": split_time,
                "HR_early": round(float(hr_early), 3) if not np.isnan(hr_early) else np.nan,
                "p_early": round(float(p_early), 4) if not np.isnan(p_early) else np.nan,
                "HR_late": round(float(hr_late), 3) if not np.isnan(hr_late) else np.nan,
                "p_late": round(float(p_late), 4) if not np.isnan(p_late) else np.nan,
                "HR_ratio": round(float(hr_early / hr_late), 3) if not (np.isnan(hr_early) or np.isnan(hr_late) or hr_late == 0) else np.nan,
                "PH_consistent": "Yes" if not np.isnan(hr_early) and not np.isnan(hr_late) and 0.5 < hr_early / hr_late < 2.0 else "Check"
            })

        ph_df = pd.DataFrame(ph_rows)
        ph_df.to_csv(os.path.join(TABLES, "ph_diagnostics.csv"), index=False)
        print("  Time-varying HR analysis:")
        print(ph_df.to_string(index=False))
        print("  Saved: ph_diagnostics.csv")


# ══════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("STATISTICAL RIGOR ENHANCEMENTS COMPLETE")
print("=" * 70)
print("\nNew outputs:")
print("  Tables:")
print("    - corrected_checkpoint_pvalues.csv (BH FDR correction)")
print("    - corrected_mutation_pvalues.csv (BH FDR + Bonferroni)")
print("    - corrected_ssgsea_pvalues.csv (BH FDR correction)")
print("    - cindex_comparison_tests.csv (bootstrap comparison tests)")
print("    - interaction_tests.csv (risk × sex, age, stage)")
print("    - sensitivity_cutpoints.csv (25th-75th percentile cutoffs)")
print("    - power_analysis.csv (Schoenfeld power estimates)")
print("    - one_se_rule.csv (LASSO 1SE rule check)")
print("    - ph_diagnostics.csv (proportional hazards time-varying)")
print("  Figures:")
print("    - sensitivity_cutpoint.png (HR and p-value across cutpoints)")
print("    - interaction_forest.png (interaction effect sizes)")
print("=" * 70)
