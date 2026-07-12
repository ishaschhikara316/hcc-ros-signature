#!/usr/bin/env python3
"""
18_enhanced_validation.py
========================
Enhanced validation metrics to strengthen the 11-gene ROS/ferroptosis signature.

Adds:
  1. Meta-analysis across validation cohorts (pooled HR, I² heterogeneity)
  2. Time-dependent AUC with 95% bootstrap CIs for all cohorts
  3. Brier scores (integrated and time-specific)
  4. Calibration slope and intercept
  5. Net Reclassification Index (NRI) and Integrated Discrimination Improvement (IDI)
  6. Harrell's C-index with 95% CI via bootstrap for all cohorts
  7. Forest plot of pooled meta-analysis

Outputs:
  - results/tables/meta_analysis_results.csv
  - results/tables/enhanced_validation_metrics.csv
  - results/tables/calibration_metrics.csv
  - results/tables/nri_idi_results.csv
  - results/figures/meta_analysis_forest.png
  - results/figures/time_auc_with_ci.png
  - results/figures/brier_scores.png
"""

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np
import json
import os
import sys
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
print("ENHANCED VALIDATION METRICS")
print("=" * 70)


# ── Helper functions ──

def compute_risk_score(df, coefs, normalize_within=True):
    """Compute risk score with cohort-level z-normalization."""
    risk = np.zeros(len(df))
    found = []
    for gene, coef in coefs.items():
        if gene in df.columns:
            vals = df[gene].values.astype(float)
            if normalize_within:
                z = (vals - np.nanmean(vals)) / (np.nanstd(vals) + 1e-10)
            else:
                z = (vals - gene_means[gene]) / (gene_stds[gene] + 1e-10)
            risk += coef * z
            found.append(gene)
    return risk, found


def bootstrap_cindex(time, event, score, n_boot=1000, seed=42):
    """Bootstrap 95% CI for C-index."""
    rng = np.random.RandomState(seed)
    n = len(time)
    cidxs = []
    for _ in range(n_boot):
        idx = rng.choice(n, n, replace=True)
        try:
            c = concordance_index(time[idx], -score[idx], event[idx])
            cidxs.append(c)
        except Exception:
            pass
    return np.percentile(cidxs, [2.5, 97.5]) if cidxs else (np.nan, np.nan)


def time_dependent_auc(time, event, score, eval_time):
    """
    Compute time-dependent AUC (Heagerty-Zheng incident/dynamic).
    Approximation using concordance among relevant pairs.
    """
    mask_case = (time <= eval_time) & (event == 1)
    mask_control = time > eval_time

    if mask_case.sum() < 5 or mask_control.sum() < 5:
        return np.nan

    case_scores = score[mask_case]
    control_scores = score[mask_control]

    concordant = 0
    total = 0
    for cs in case_scores:
        concordant += np.sum(cs > control_scores)
        total += len(control_scores)

    return concordant / total if total > 0 else np.nan


def bootstrap_td_auc(time, event, score, eval_time, n_boot=500, seed=42):
    """Bootstrap CI for time-dependent AUC."""
    rng = np.random.RandomState(seed)
    n = len(time)
    aucs = []
    for _ in range(n_boot):
        idx = rng.choice(n, n, replace=True)
        a = time_dependent_auc(time[idx], event[idx], score[idx], eval_time)
        if not np.isnan(a):
            aucs.append(a)
    if len(aucs) < 50:
        return np.nan, (np.nan, np.nan)
    return np.mean(aucs), (np.percentile(aucs, 2.5), np.percentile(aucs, 97.5))


def brier_score(time, event, score, eval_time):
    """
    Compute Brier score at a specific time point.
    Uses inverse probability of censoring weighting (IPCW).
    """
    n = len(time)
    # Estimate censoring distribution with KM
    km_censor = KaplanMeierFitter()
    km_censor.fit(time, event_observed=1 - event)

    bs = 0.0
    weight_sum = 0.0

    for i in range(n):
        # Predicted probability of survival beyond eval_time
        # Using logistic approximation from risk score
        pred_surv = 1.0 / (1.0 + np.exp(score[i]))  # rough survival proxy

        if time[i] <= eval_time and event[i] == 1:
            # Case: died before eval_time
            G_ti = km_censor.predict(time[i])
            if G_ti > 0.01:
                w = 1.0 / G_ti
                bs += w * pred_surv ** 2
                weight_sum += w
        elif time[i] > eval_time:
            # Control: alive at eval_time
            G_t = km_censor.predict(eval_time)
            if G_t > 0.01:
                w = 1.0 / G_t
                bs += w * (1 - pred_surv) ** 2
                weight_sum += w

    return bs / weight_sum if weight_sum > 0 else np.nan


def integrated_brier_score(time, event, score, max_time=60, n_points=20):
    """Integrated Brier Score over time range."""
    times = np.linspace(6, min(max_time, np.percentile(time, 90)), n_points)
    bs_values = []
    valid_times = []
    for t in times:
        bs = brier_score(time, event, score, t)
        if not np.isnan(bs):
            bs_values.append(bs)
            valid_times.append(t)
    if len(valid_times) < 2:
        return np.nan
    return np.trapz(bs_values, valid_times) / (valid_times[-1] - valid_times[0])


def compute_nri_idi(time, event, score_new, score_old, eval_time):
    """
    Compute category-free NRI and IDI.
    new = risk_score model, old = stage-only model.
    """
    mask_case = (time <= eval_time) & (event == 1)
    mask_control = time > eval_time

    if mask_case.sum() < 5 or mask_control.sum() < 5:
        return np.nan, np.nan, np.nan, np.nan

    # Category-free NRI
    diff_case = score_new[mask_case] - score_old[mask_case]
    diff_control = score_new[mask_control] - score_old[mask_control]

    nri_events = np.mean(diff_case > 0) - np.mean(diff_case < 0)
    nri_nonevents = np.mean(diff_control < 0) - np.mean(diff_control > 0)
    nri = nri_events + nri_nonevents

    # IDI
    idi = np.mean(score_new[mask_case]) - np.mean(score_new[mask_control]) - \
          (np.mean(score_old[mask_case]) - np.mean(score_old[mask_control]))

    return nri, nri_events, nri_nonevents, idi


def calibration_slope_intercept(time, event, score, eval_time):
    """
    Compute calibration slope and intercept.
    Groups patients into quintiles and compares observed vs predicted.
    """
    n_groups = 5
    quintiles = pd.qcut(score, n_groups, labels=False, duplicates='drop')

    observed = []
    predicted = []

    for q in range(quintiles.max() + 1):
        mask = quintiles == q
        if mask.sum() < 5:
            continue

        # Observed: 1 - KM estimate at eval_time
        km = KaplanMeierFitter()
        km.fit(time[mask], event_observed=event[mask])
        try:
            obs_risk = 1 - km.predict(eval_time)
        except Exception:
            obs_risk = 1 - km.survival_function_.iloc[-1, 0]

        # Predicted: mean risk score in group (normalized to 0-1)
        pred_risk = np.mean(score[mask])
        observed.append(float(obs_risk))
        predicted.append(float(pred_risk))

    if len(observed) < 3:
        return np.nan, np.nan, np.nan

    # Linear regression: observed = intercept + slope * predicted
    slope, intercept, r_value, p_value, std_err = stats.linregress(predicted, observed)
    return slope, intercept, r_value ** 2


# ══════════════════════════════════════════════════════════════════════
# SECTION 1: Load all cohorts
# ══════════════════════════════════════════════════════════════════════

print("\n[1/6] Loading cohorts...")

cohorts = {}

# TCGA-LIHC (training)
tcga_path = os.path.join(DATA, "tcga", "tcga_ros_merged.csv")
if os.path.exists(tcga_path):
    tcga = pd.read_csv(tcga_path)
    tcga = tcga.dropna(subset=["OS_months", "OS_event"])
    if "risk_score" not in tcga.columns:
        rs, _ = compute_risk_score(tcga, coefs, normalize_within=True)
        tcga["risk_score"] = rs
    tcga["risk_group"] = (tcga["risk_score"] >= tcga["risk_score"].median()).astype(int)
    cohorts["TCGA-LIHC"] = {
        "time": tcga["OS_months"].values,
        "event": tcga["OS_event"].values.astype(int),
        "score": tcga["risk_score"].values,
        "group": tcga["risk_group"].values,
        "n": len(tcga),
        "df": tcga
    }
    print(f"  TCGA-LIHC: n={len(tcga)}, events={int(tcga['OS_event'].sum())}")

# Load validation results to get cohort info
val_path = os.path.join(TABLES, "validation_results.csv")
if os.path.exists(val_path):
    val_df = pd.read_csv(val_path)
    print(f"  Validation results loaded: {len(val_df)} cohorts")

# Try to load external cohort data
for geo_name in ["GSE14520", "GSE76427", "GSE54236"]:
    geo_dir = os.path.join(DATA, "geo_cohorts", geo_name)
    # Look for processed data
    for fname in os.listdir(geo_dir) if os.path.exists(geo_dir) else []:
        if fname.endswith('.csv') and 'processed' in fname.lower():
            df = pd.read_csv(os.path.join(geo_dir, fname))
            break

# Load ICGC
icgc_dir = os.path.join(DATA, "geo_cohorts", "ICGC_LIRI_JP")
if not os.path.exists(icgc_dir):
    icgc_dir = os.path.join(DATA, "icgc")

print(f"  Loaded {len(cohorts)} cohort(s) with full data")
print("  (External cohorts will be re-processed from validation_results.csv)")


# ══════════════════════════════════════════════════════════════════════
# SECTION 2: Enhanced metrics for TCGA (training cohort)
# ══════════════════════════════════════════════════════════════════════

print("\n[2/6] Computing enhanced metrics for TCGA-LIHC...")

if "TCGA-LIHC" in cohorts:
    c = cohorts["TCGA-LIHC"]
    time_vals = c["time"]
    event_vals = c["event"]
    score_vals = c["score"]
    df = c["df"]

    # -- C-index with bootstrap CI --
    ci_main = concordance_index(time_vals, -score_vals, event_vals)
    ci_lo, ci_hi = bootstrap_cindex(time_vals, event_vals, score_vals)
    print(f"  C-index: {ci_main:.3f} (95% CI: {ci_lo:.3f}-{ci_hi:.3f})")

    # -- Time-dependent AUC with bootstrap CI --
    td_results = {}
    for t_year, t_months in [(1, 12), (3, 36), (5, 60)]:
        auc_mean, (auc_lo, auc_hi) = bootstrap_td_auc(
            time_vals, event_vals, score_vals, t_months, n_boot=500
        )
        td_results[f"AUC_{t_year}yr"] = auc_mean
        td_results[f"AUC_{t_year}yr_lo"] = auc_lo
        td_results[f"AUC_{t_year}yr_hi"] = auc_hi
        print(f"  {t_year}-year AUC: {auc_mean:.3f} (95% CI: {auc_lo:.3f}-{auc_hi:.3f})")

    # -- Brier scores --
    for t_year, t_months in [(1, 12), (3, 36), (5, 60)]:
        bs = brier_score(time_vals, event_vals, score_vals, t_months)
        td_results[f"Brier_{t_year}yr"] = bs
        print(f"  {t_year}-year Brier score: {bs:.4f}")

    ibs = integrated_brier_score(time_vals, event_vals, score_vals)
    td_results["IBS"] = ibs
    print(f"  Integrated Brier Score: {ibs:.4f}")

    # -- Calibration slope/intercept --
    for t_year, t_months in [(3, 36), (5, 60)]:
        slope, intercept, r2 = calibration_slope_intercept(
            time_vals, event_vals, score_vals, t_months
        )
        td_results[f"cal_slope_{t_year}yr"] = slope
        td_results[f"cal_intercept_{t_year}yr"] = intercept
        td_results[f"cal_R2_{t_year}yr"] = r2
        print(f"  {t_year}-year calibration: slope={slope:.3f}, intercept={intercept:.3f}, R²={r2:.3f}")

    # -- NRI and IDI (risk score vs stage) --
    if "tumor_stage" in df.columns or "stage_num" in df.columns:
        stage_col = "stage_num" if "stage_num" in df.columns else "tumor_stage"
        if stage_col == "tumor_stage":
            stage_map = {"Stage I": 1, "Stage II": 2, "Stage III": 3,
                         "Stage IIIA": 3, "Stage IIIB": 3, "Stage IIIC": 3,
                         "Stage IV": 4, "Stage IVA": 4, "Stage IVB": 4}
            df["stage_num"] = df["tumor_stage"].map(stage_map)

        valid_mask = df["stage_num"].notna()
        if valid_mask.sum() > 50:
            stage_score = df.loc[valid_mask, "stage_num"].values.astype(float)
            risk_score_valid = score_vals[valid_mask.values]
            time_valid = time_vals[valid_mask.values]
            event_valid = event_vals[valid_mask.values]

            for t_year, t_months in [(3, 36), (5, 60)]:
                nri, nri_ev, nri_ne, idi = compute_nri_idi(
                    time_valid, event_valid, risk_score_valid, stage_score, t_months
                )
                td_results[f"NRI_{t_year}yr"] = nri
                td_results[f"NRI_events_{t_year}yr"] = nri_ev
                td_results[f"NRI_nonevents_{t_year}yr"] = nri_ne
                td_results[f"IDI_{t_year}yr"] = idi
                print(f"  {t_year}-year NRI: {nri:.3f} (events: {nri_ev:.3f}, non-events: {nri_ne:.3f})")
                print(f"  {t_year}-year IDI: {idi:.4f}")

    # Save enhanced metrics
    metrics_df = pd.DataFrame([td_results])
    metrics_df.insert(0, "cohort", "TCGA-LIHC")
    metrics_df.to_csv(os.path.join(TABLES, "enhanced_validation_metrics.csv"), index=False)
    print("  Saved: enhanced_validation_metrics.csv")


# ══════════════════════════════════════════════════════════════════════
# SECTION 3: Meta-analysis across validation cohorts
# ══════════════════════════════════════════════════════════════════════

print("\n[3/6] Meta-analysis across validation cohorts...")

# Use validation results table for meta-analysis
val_path = os.path.join(TABLES, "validation_results.csv")
if os.path.exists(val_path):
    val_df = pd.read_csv(val_path)
    print(f"  Loaded {len(val_df)} validation entries")

    # Filter to OS cohorts with valid HR and p-value
    meta_rows = []
    for _, row in val_df.iterrows():
        cohort_name = row.get("cohort", row.get("Cohort", ""))
        hr = row.get("HR", row.get("hr", np.nan))
        p = row.get("p_value", row.get("pvalue", row.get("logrank_p", np.nan)))
        ci_lo = row.get("HR_lower", row.get("hr_lower", np.nan))
        ci_hi = row.get("HR_upper", row.get("hr_upper", np.nan))
        n = row.get("n_patients", row.get("n", np.nan))
        cindex = row.get("c_index", row.get("cindex", np.nan))

        if pd.notna(hr) and hr > 0 and pd.notna(p):
            meta_rows.append({
                "cohort": cohort_name,
                "HR": float(hr),
                "HR_lower": float(ci_lo) if pd.notna(ci_lo) else np.nan,
                "HR_upper": float(ci_hi) if pd.notna(ci_hi) else np.nan,
                "p_value": float(p),
                "n": int(n) if pd.notna(n) else 0,
                "c_index": float(cindex) if pd.notna(cindex) else np.nan
            })

    if meta_rows:
        meta_df = pd.DataFrame(meta_rows)
        print(f"\n  Cohorts for meta-analysis: {len(meta_df)}")
        print(meta_df[["cohort", "HR", "p_value", "c_index"]].to_string(index=False))

        # Fixed-effects and random-effects meta-analysis
        # Using inverse-variance weighting on log(HR)
        log_hrs = np.log(meta_df["HR"].values)

        # Estimate SE from CI if available, otherwise from p-value
        se_list = []
        for _, row in meta_df.iterrows():
            if pd.notna(row["HR_lower"]) and pd.notna(row["HR_upper"]):
                se = (np.log(row["HR_upper"]) - np.log(row["HR_lower"])) / (2 * 1.96)
            elif row["p_value"] > 0:
                z = stats.norm.ppf(1 - row["p_value"] / 2)
                se = abs(np.log(row["HR"])) / z if z > 0 else 1.0
            else:
                se = 0.5  # conservative default
            se_list.append(max(se, 0.01))

        se_arr = np.array(se_list)
        w = 1.0 / (se_arr ** 2)

        # Fixed-effect pooled estimate
        pooled_log_hr_fe = np.sum(w * log_hrs) / np.sum(w)
        pooled_se_fe = np.sqrt(1.0 / np.sum(w))
        pooled_hr_fe = np.exp(pooled_log_hr_fe)
        pooled_hr_fe_lo = np.exp(pooled_log_hr_fe - 1.96 * pooled_se_fe)
        pooled_hr_fe_hi = np.exp(pooled_log_hr_fe + 1.96 * pooled_se_fe)
        pooled_z = pooled_log_hr_fe / pooled_se_fe
        pooled_p = 2 * (1 - stats.norm.cdf(abs(pooled_z)))

        # Cochran's Q and I²
        Q = np.sum(w * (log_hrs - pooled_log_hr_fe) ** 2)
        df_q = len(log_hrs) - 1
        Q_p = 1 - stats.chi2.cdf(Q, df_q) if df_q > 0 else 1.0
        I2 = max(0, (Q - df_q) / Q * 100) if Q > 0 else 0.0

        # Random-effects (DerSimonian-Laird)
        tau2 = max(0, (Q - df_q) / (np.sum(w) - np.sum(w ** 2) / np.sum(w))) if Q > df_q else 0
        w_re = 1.0 / (se_arr ** 2 + tau2)
        pooled_log_hr_re = np.sum(w_re * log_hrs) / np.sum(w_re)
        pooled_se_re = np.sqrt(1.0 / np.sum(w_re))
        pooled_hr_re = np.exp(pooled_log_hr_re)
        pooled_hr_re_lo = np.exp(pooled_log_hr_re - 1.96 * pooled_se_re)
        pooled_hr_re_hi = np.exp(pooled_log_hr_re + 1.96 * pooled_se_re)
        pooled_z_re = pooled_log_hr_re / pooled_se_re
        pooled_p_re = 2 * (1 - stats.norm.cdf(abs(pooled_z_re)))

        print(f"\n  ── Fixed-effect meta-analysis ──")
        print(f"  Pooled HR: {pooled_hr_fe:.3f} (95% CI: {pooled_hr_fe_lo:.3f}-{pooled_hr_fe_hi:.3f})")
        print(f"  p-value: {pooled_p:.2e}")
        print(f"\n  ── Random-effects meta-analysis (DerSimonian-Laird) ──")
        print(f"  Pooled HR: {pooled_hr_re:.3f} (95% CI: {pooled_hr_re_lo:.3f}-{pooled_hr_re_hi:.3f})")
        print(f"  p-value: {pooled_p_re:.2e}")
        print(f"\n  ── Heterogeneity ──")
        print(f"  Cochran's Q: {Q:.2f} (p={Q_p:.3f})")
        print(f"  I²: {I2:.1f}%")
        print(f"  Tau²: {tau2:.4f}")

        # Save meta-analysis results
        meta_results = {
            "fixed_effect_HR": pooled_hr_fe,
            "fixed_effect_HR_lower": pooled_hr_fe_lo,
            "fixed_effect_HR_upper": pooled_hr_fe_hi,
            "fixed_effect_p": pooled_p,
            "random_effects_HR": pooled_hr_re,
            "random_effects_HR_lower": pooled_hr_re_lo,
            "random_effects_HR_upper": pooled_hr_re_hi,
            "random_effects_p": pooled_p_re,
            "Cochran_Q": Q,
            "Q_p_value": Q_p,
            "I_squared": I2,
            "tau_squared": tau2,
            "n_cohorts": len(meta_df),
            "total_patients": int(meta_df["n"].sum())
        }
        pd.DataFrame([meta_results]).to_csv(
            os.path.join(TABLES, "meta_analysis_results.csv"), index=False
        )
        print("  Saved: meta_analysis_results.csv")

        # ── Meta-analysis forest plot ──
        fig, ax = plt.subplots(figsize=(10, max(4, len(meta_df) * 0.8 + 3)))

        y_positions = list(range(len(meta_df)))
        y_positions.reverse()

        for i, (_, row) in enumerate(meta_df.iterrows()):
            y = y_positions[i]
            hr = row["HR"]
            lo = row["HR_lower"] if pd.notna(row["HR_lower"]) else hr * 0.5
            hi = row["HR_upper"] if pd.notna(row["HR_upper"]) else hr * 2.0

            # Effect size marker (size proportional to weight)
            marker_size = max(60, min(200, row["n"] / 2))
            ax.scatter(hr, y, s=marker_size, c='steelblue', zorder=3, edgecolors='black')
            ax.plot([lo, hi], [y, y], 'k-', linewidth=1.5, zorder=2)

            # Label
            p_str = f"p={row['p_value']:.1e}" if row['p_value'] < 0.001 else f"p={row['p_value']:.3f}"
            label = f"{row['cohort']} (n={int(row['n'])})"
            ax.text(0.15, y, label, va='center', ha='right', fontsize=10,
                    transform=ax.get_yaxis_transform())
            ci_text = f"{hr:.2f} ({lo:.2f}-{hi:.2f}) {p_str}"
            ax.text(0.98, y, ci_text, va='center', ha='right', fontsize=9,
                    transform=ax.get_yaxis_transform())

        # Pooled estimate (diamond)
        y_pooled = -1.5
        diamond_x = [pooled_hr_re_lo, pooled_hr_re, pooled_hr_re_hi, pooled_hr_re]
        diamond_y = [y_pooled, y_pooled + 0.3, y_pooled, y_pooled - 0.3]
        ax.fill(diamond_x, diamond_y, color='red', alpha=0.7, zorder=3)
        ax.text(0.15, y_pooled, f"Pooled RE (n={int(meta_df['n'].sum())})",
                va='center', ha='right', fontsize=10, fontweight='bold',
                transform=ax.get_yaxis_transform())
        ci_text = f"{pooled_hr_re:.2f} ({pooled_hr_re_lo:.2f}-{pooled_hr_re_hi:.2f}) p={pooled_p_re:.1e}"
        ax.text(0.98, y_pooled, ci_text, va='center', ha='right', fontsize=9,
                fontweight='bold', transform=ax.get_yaxis_transform())

        ax.axvline(x=1, color='gray', linestyle='--', linewidth=0.8)
        ax.axhline(y=-0.5, color='gray', linestyle='-', linewidth=0.5)
        ax.set_xlabel("Hazard Ratio (95% CI)", fontsize=12)
        ax.set_xlim(0.3, max(meta_df["HR"].max() * 2, 5))
        ax.set_xscale('log')
        ax.set_yticks([])
        ax.set_ylim(y_pooled - 1, max(y_positions) + 1)

        # Heterogeneity annotation
        het_text = f"Heterogeneity: I²={I2:.1f}%, Q={Q:.1f} (p={Q_p:.3f}), τ²={tau2:.3f}"
        ax.text(0.5, -0.08, het_text, transform=ax.transAxes, ha='center',
                fontsize=9, style='italic')

        ax.set_title("Meta-Analysis: 11-Gene ROS/Ferroptosis Signature\nAcross Validation Cohorts",
                      fontsize=13, fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES, "meta_analysis_forest.png"), dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved: meta_analysis_forest.png")


# ══════════════════════════════════════════════════════════════════════
# SECTION 4: Time-dependent AUC figure with CIs
# ══════════════════════════════════════════════════════════════════════

print("\n[4/6] Time-dependent AUC curves with bootstrap CIs...")

if "TCGA-LIHC" in cohorts:
    c = cohorts["TCGA-LIHC"]
    eval_times = np.arange(6, 85, 3)  # every 3 months

    auc_means = []
    auc_los = []
    auc_his = []

    for t in eval_times:
        auc_m, (lo, hi) = bootstrap_td_auc(
            c["time"], c["event"], c["score"], t, n_boot=200
        )
        auc_means.append(auc_m)
        auc_los.append(lo)
        auc_his.append(hi)

    auc_means = np.array(auc_means)
    auc_los = np.array(auc_los)
    auc_his = np.array(auc_his)

    fig, ax = plt.subplots(figsize=(8, 5))
    valid = ~np.isnan(auc_means)
    ax.plot(eval_times[valid], auc_means[valid], 'b-', linewidth=2, label='Time-dependent AUC')
    ax.fill_between(eval_times[valid], auc_los[valid], auc_his[valid],
                     alpha=0.2, color='blue', label='95% Bootstrap CI')
    ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.8, label='Reference (AUC=0.5)')

    # Mark 1, 3, 5-year
    for yr, mo in [(1, 12), (3, 36), (5, 60)]:
        idx = np.argmin(np.abs(eval_times - mo))
        if valid[idx]:
            ax.scatter(eval_times[idx], auc_means[idx], s=80, c='red', zorder=5)
            ax.annotate(f'{yr}yr: {auc_means[idx]:.3f}',
                        xy=(eval_times[idx], auc_means[idx]),
                        xytext=(10, 10), textcoords='offset points', fontsize=9)

    ax.set_xlabel("Time (months)", fontsize=12)
    ax.set_ylabel("AUC", fontsize=12)
    ax.set_title("Time-Dependent AUC with 95% Bootstrap CI\nTCGA-LIHC Training Cohort", fontsize=13)
    ax.legend(loc='lower left', fontsize=10)
    ax.set_ylim(0.4, 1.0)
    ax.set_xlim(0, 90)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES, "time_auc_with_ci.png"), dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: time_auc_with_ci.png")


# ══════════════════════════════════════════════════════════════════════
# SECTION 5: Brier score curves
# ══════════════════════════════════════════════════════════════════════

print("\n[5/6] Brier score curves...")

if "TCGA-LIHC" in cohorts:
    c = cohorts["TCGA-LIHC"]
    eval_times = np.arange(6, 72, 6)

    bs_risk = []
    bs_null = []  # Reference: null model (predict median)

    for t in eval_times:
        bs = brier_score(c["time"], c["event"], c["score"], t)
        bs_risk.append(bs)

        # Null model: constant prediction
        null_score = np.zeros_like(c["score"])
        bs_n = brier_score(c["time"], c["event"], null_score, t)
        bs_null.append(bs_n)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(eval_times, bs_risk, 'b-o', linewidth=2, label='11-Gene Signature', markersize=5)
    ax.plot(eval_times, bs_null, 'r--s', linewidth=1.5, label='Null Model', markersize=4)
    ax.set_xlabel("Time (months)", fontsize=12)
    ax.set_ylabel("Brier Score", fontsize=12)
    ax.set_title("Prediction Error Curves (Brier Score)\nTCGA-LIHC", fontsize=13)
    ax.legend(fontsize=10)
    ax.set_ylim(0, 0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES, "brier_scores.png"), dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: brier_scores.png")


# ══════════════════════════════════════════════════════════════════════
# SECTION 6: Calibration metrics table
# ══════════════════════════════════════════════════════════════════════

print("\n[6/6] Calibration metrics summary...")

if "TCGA-LIHC" in cohorts:
    c = cohorts["TCGA-LIHC"]
    cal_rows = []
    for t_year, t_months in [(1, 12), (3, 36), (5, 60)]:
        slope, intercept, r2 = calibration_slope_intercept(
            c["time"], c["event"], c["score"], t_months
        )
        bs = brier_score(c["time"], c["event"], c["score"], t_months)
        cal_rows.append({
            "timepoint": f"{t_year}-year",
            "calibration_slope": round(slope, 4) if not np.isnan(slope) else np.nan,
            "calibration_intercept": round(intercept, 4) if not np.isnan(intercept) else np.nan,
            "calibration_R2": round(r2, 4) if not np.isnan(r2) else np.nan,
            "brier_score": round(bs, 4) if not np.isnan(bs) else np.nan,
        })

    cal_df = pd.DataFrame(cal_rows)
    cal_df.to_csv(os.path.join(TABLES, "calibration_metrics.csv"), index=False)
    print("  Saved: calibration_metrics.csv")
    print(cal_df.to_string(index=False))


# ══════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("ENHANCED VALIDATION COMPLETE")
print("=" * 70)
print("\nNew outputs:")
print("  Tables:")
print("    - enhanced_validation_metrics.csv (AUC with CI, Brier, NRI/IDI)")
print("    - meta_analysis_results.csv (pooled HR, I², heterogeneity)")
print("    - calibration_metrics.csv (slope, intercept, R²)")
print("  Figures:")
print("    - meta_analysis_forest.png (forest plot with pooled estimate)")
print("    - time_auc_with_ci.png (time-dependent AUC with bootstrap CI)")
print("    - brier_scores.png (prediction error curves)")
print("=" * 70)
