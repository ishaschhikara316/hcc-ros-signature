"""
11_nomogram_dca.py — True nomogram + Decision Curve Analysis

1. Build points-based nomogram (risk_score + age + stage)
2. Predict 1/3/5-year survival probabilities
3. Decision curve analysis showing net benefit
4. Calibration with nomogram predictions
"""
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

merged = pd.read_csv(os.path.join(DATA, "tcga_ros_merged.csv"))

# Compute risk score if needed
if "risk_score" not in merged.columns:
    risk = np.zeros(len(merged))
    for gene, coef in selected_genes.items():
        z = (merged[gene] - model["gene_means"][gene]) / model["gene_stds"][gene]
        risk += coef * z.values
    merged["risk_score"] = risk

# Prepare clinical variables
def encode_stage(s):
    if pd.isna(s): return np.nan
    s = str(s).upper().strip()
    if "IV" in s: return 4
    if "III" in s: return 3
    if "II" in s: return 2
    if "I" in s: return 1
    return np.nan

merged["stage_num"] = merged["tumor_stage"].apply(encode_stage) if "tumor_stage" in merged.columns else np.nan
merged["male"] = (merged["gender"].str.lower() == "male").astype(int) if "gender" in merged.columns else np.nan

if "age_at_diagnosis" in merged.columns:
    age = merged["age_at_diagnosis"]
    merged["age_years"] = np.where(age > 200, age / 365.25, age)
elif "age" in merged.columns:
    merged["age_years"] = merged["age"]

df = merged.dropna(subset=["OS_months", "OS_event", "risk_score"]).copy()
df = df[df["OS_months"] > 0]
print(f"Working with {len(df)} patients ({int(df['OS_event'].sum())} events)")
print(f"Signature: {len(selected_genes)} genes")

# ══════════════════════════════════════════════════════════════════════════════
# 1. FIT NOMOGRAM COX MODEL
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("1. NOMOGRAM COX MODEL")
print("=" * 70)

# Select variables for nomogram
nomo_vars = ["risk_score"]
if "age_years" in df.columns and df["age_years"].notna().sum() > 100:
    nomo_vars.append("age_years")
if "stage_num" in df.columns and df["stage_num"].notna().sum() > 100:
    nomo_vars.append("stage_num")
if "male" in df.columns and df["male"].notna().sum() > 100:
    nomo_vars.append("male")

nomo_cols = ["OS_months", "OS_event"] + nomo_vars
nomo_df = df[nomo_cols].dropna()
nomo_df = nomo_df[nomo_df["OS_months"] > 0]
print(f"Nomogram variables: {nomo_vars}")
print(f"Patients for nomogram: {len(nomo_df)}")

cph = CoxPHFitter()
cph.fit(nomo_df, duration_col="OS_months", event_col="OS_event")
print("\nCox model summary:")
print(cph.summary[["coef", "exp(coef)", "p"]].to_string())

# ══════════════════════════════════════════════════════════════════════════════
# 2. POINTS-BASED NOMOGRAM
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("2. POINTS-BASED NOMOGRAM")
print("=" * 70)

# Build nomogram: convert each variable to a points scale (0-100)
var_ranges = {}
var_points_func = {}
max_total_beta = 0

for var in nomo_vars:
    coef = cph.params_[var]
    v_min = nomo_df[var].quantile(0.01)
    v_max = nomo_df[var].quantile(0.99)
    var_ranges[var] = (v_min, v_max)
    beta_range = abs(coef * (v_max - v_min))
    max_total_beta += beta_range

# Assign points proportional to beta contribution
points_per_unit_beta = 100 / (max_total_beta + 1e-10)

fig = plt.figure(figsize=(14, max(6, 2 + len(nomo_vars) * 1.8 + 4)))
gs = gridspec.GridSpec(len(nomo_vars) + 3, 1,
                        height_ratios=[1] * (len(nomo_vars) + 1) + [0.5, 1.5])

# Points scale (top)
ax_points = fig.add_subplot(gs[0])
ax_points.set_xlim(0, 100)
ax_points.set_xticks(np.arange(0, 101, 10))
ax_points.set_xlabel("")
ax_points.set_title("Points", fontsize=11, fontweight='bold')
ax_points.yaxis.set_visible(False)
ax_points.spines['left'].set_visible(False)
ax_points.spines['right'].set_visible(False)

# Variable scales
total_points_data = np.zeros(len(nomo_df))

for i, var in enumerate(nomo_vars):
    ax = fig.add_subplot(gs[i + 1])
    coef = cph.params_[var]
    v_min, v_max = var_ranges[var]

    # Map variable values to points
    beta_contrib = coef * (nomo_df[var].values - v_min)
    points = beta_contrib * points_per_unit_beta
    # Ensure all positive
    if coef < 0:
        points = -coef * (v_max - nomo_df[var].values) * points_per_unit_beta
    total_points_data += points

    ax.set_xlim(0, 100)

    # Create tick marks mapping points → variable values
    n_ticks = 8
    point_ticks = np.linspace(0, abs(coef * (v_max - v_min)) * points_per_unit_beta, n_ticks)
    if coef > 0:
        val_ticks = v_min + point_ticks / (abs(coef) * points_per_unit_beta + 1e-10)
    else:
        val_ticks = v_max - point_ticks / (abs(coef) * points_per_unit_beta + 1e-10)

    ax.set_xticks(point_ticks)
    # Format labels
    if var == "age_years":
        ax.set_xticklabels([f"{v:.0f}" for v in val_ticks], fontsize=8)
        label = "Age (years)"
    elif var == "stage_num":
        ax.set_xticklabels([f"{v:.0f}" for v in val_ticks], fontsize=8)
        label = "TNM Stage"
    elif var == "risk_score":
        ax.set_xticklabels([f"{v:.1f}" for v in val_ticks], fontsize=8)
        label = "Risk Score"
    elif var == "male":
        ax.set_xticks([0, abs(coef) * points_per_unit_beta])
        ax.set_xticklabels(["Female", "Male"], fontsize=9)
        label = "Sex"
    else:
        ax.set_xticklabels([f"{v:.1f}" for v in val_ticks], fontsize=8)
        label = var

    ax.set_title(label, fontsize=11, fontweight='bold', loc='left')
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Total points
ax_total = fig.add_subplot(gs[len(nomo_vars) + 1])
total_range = (total_points_data.min(), total_points_data.max())
ax_total.set_xlim(0, 100)  # Normalized
ax_total.set_title("Total Points", fontsize=11, fontweight='bold', loc='left')
total_ticks = np.linspace(total_range[0], total_range[1], 10)
tick_positions = (total_ticks - total_range[0]) / (total_range[1] - total_range[0] + 1e-10) * 100
ax_total.set_xticks(tick_positions)
ax_total.set_xticklabels([f"{v:.0f}" for v in total_ticks], fontsize=8)
ax_total.yaxis.set_visible(False)
ax_total.spines['left'].set_visible(False)
ax_total.spines['right'].set_visible(False)

# Survival probability scales
ax_surv = fig.add_subplot(gs[len(nomo_vars) + 2])

# Predict survival at time points
surv_func = cph.predict_survival_function(nomo_df)
for tp, label, color in [(12, "1-year", "green"), (36, "3-year", "blue"), (60, "5-year", "red")]:
    if tp <= surv_func.index.max():
        surv_at_tp = surv_func.loc[tp] if tp in surv_func.index else surv_func.iloc[
            np.abs(surv_func.index - tp).argmin()]
        # Map survival to total points position
        ax_surv.scatter(
            (total_points_data - total_range[0]) / (total_range[1] - total_range[0] + 1e-10) * 100,
            surv_at_tp, s=2, alpha=0.3, color=color, label=f"{label} Survival"
        )

ax_surv.set_xlim(0, 100)
ax_surv.set_ylim(0, 1.05)
ax_surv.set_xlabel("Total Points (normalized)")
ax_surv.set_ylabel("Survival Prob.")
ax_surv.legend(fontsize=8, loc='upper right')
ax_surv.spines['top'].set_visible(False)
ax_surv.spines['right'].set_visible(False)

plt.suptitle("Prognostic Nomogram — ROS/Ferroptosis Signature",
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "nomogram_full.png"), dpi=200, bbox_inches='tight')
print("Saved: nomogram_full.png")

# ══════════════════════════════════════════════════════════════════════════════
# 3. PREDICTED SURVIVAL PROBABILITIES
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("3. PREDICTED SURVIVAL PROBABILITIES")
print("=" * 70)

surv_preds = cph.predict_survival_function(nomo_df)
for tp_name, tp_months in [("1-year", 12), ("3-year", 36), ("5-year", 60)]:
    if tp_months <= surv_preds.index.max():
        closest_idx = surv_preds.index[np.abs(surv_preds.index - tp_months).argmin()]
        probs = surv_preds.loc[closest_idx]
        print(f"  {tp_name} survival: median={probs.median():.3f}, "
              f"range={probs.min():.3f}-{probs.max():.3f}")

# Risk group predicted survival
med_rs = nomo_df["risk_score"].median()
high_risk = nomo_df[nomo_df["risk_score"] >= med_rs]
low_risk = nomo_df[nomo_df["risk_score"] < med_rs]

print("\n  High-risk group:")
surv_high = cph.predict_survival_function(high_risk)
for tp_name, tp_months in [("1-year", 12), ("3-year", 36), ("5-year", 60)]:
    if tp_months <= surv_high.index.max():
        closest_idx = surv_high.index[np.abs(surv_high.index - tp_months).argmin()]
        print(f"    {tp_name}: {surv_high.loc[closest_idx].median():.3f}")

print("  Low-risk group:")
surv_low = cph.predict_survival_function(low_risk)
for tp_name, tp_months in [("1-year", 12), ("3-year", 36), ("5-year", 60)]:
    if tp_months <= surv_low.index.max():
        closest_idx = surv_low.index[np.abs(surv_low.index - tp_months).argmin()]
        print(f"    {tp_name}: {surv_low.loc[closest_idx].median():.3f}")

# ══════════════════════════════════════════════════════════════════════════════
# 4. DECISION CURVE ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("4. DECISION CURVE ANALYSIS")
print("=" * 70)

def decision_curve_analysis(time, event, predicted_prob, timepoint, thresholds=None):
    """
    Compute net benefit for a model at a given timepoint.

    Net benefit = TP/N - FP/N * (pt / (1-pt))
    where pt is the threshold probability.
    """
    if thresholds is None:
        thresholds = np.arange(0.01, 0.99, 0.01)

    # Binary outcome at timepoint: event happened before timepoint
    outcome = ((time <= timepoint) & (event == 1)).astype(int)
    # Predicted probability of event = 1 - survival probability
    pred_event = 1 - predicted_prob

    n = len(outcome)
    prevalence = outcome.mean()

    net_benefits_model = []
    net_benefits_all = []
    net_benefits_none = []

    for pt in thresholds:
        # Treat all
        nb_all = prevalence - (1 - prevalence) * pt / (1 - pt + 1e-10)
        net_benefits_all.append(nb_all)

        # Treat none
        net_benefits_none.append(0)

        # Model
        predicted_positive = pred_event >= pt
        tp = ((predicted_positive) & (outcome == 1)).sum()
        fp = ((predicted_positive) & (outcome == 0)).sum()
        nb_model = tp / n - fp / n * pt / (1 - pt + 1e-10)
        net_benefits_model.append(nb_model)

    return thresholds, net_benefits_model, net_benefits_all, net_benefits_none


# Compute DCA for each timepoint
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
timepoints = {"1-year": 12, "3-year": 36, "5-year": 60}

dca_results = []

for ax, (tp_name, tp_months) in zip(axes, timepoints.items()):
    # Get predicted survival probability at timepoint
    if tp_months <= surv_preds.index.max():
        closest_idx = surv_preds.index[np.abs(surv_preds.index - tp_months).argmin()]
        pred_surv = surv_preds.loc[closest_idx].values

        thresholds, nb_model, nb_all, nb_none = decision_curve_analysis(
            nomo_df["OS_months"].values, nomo_df["OS_event"].values,
            pred_surv, tp_months
        )

        ax.plot(thresholds, nb_model, color='steelblue', linewidth=2,
                label='ROS/Ferroptosis Nomogram')
        ax.plot(thresholds, nb_all, color='gray', linewidth=1.5,
                linestyle='--', label='Treat All')
        ax.plot(thresholds, nb_none, color='black', linewidth=1.5,
                linestyle=':', label='Treat None')

        # Also compute for risk_score alone (simpler model)
        # Use Cox with risk_score only
        simple_cph = CoxPHFitter()
        simple_df = nomo_df[["OS_months", "OS_event", "risk_score"]].copy()
        simple_cph.fit(simple_df, duration_col="OS_months", event_col="OS_event")
        simple_surv = simple_cph.predict_survival_function(simple_df)
        simple_pred = simple_surv.loc[closest_idx].values if closest_idx in simple_surv.index else \
            simple_surv.iloc[(simple_surv.index - tp_months).abs().argmin()].values

        _, nb_simple, _, _ = decision_curve_analysis(
            nomo_df["OS_months"].values, nomo_df["OS_event"].values,
            simple_pred, tp_months
        )
        ax.plot(thresholds, nb_simple, color='darkorange', linewidth=1.5,
                linestyle='-', alpha=0.8, label='Risk Score Only')

        # Stage alone (if available)
        if "stage_num" in nomo_vars:
            stage_cph = CoxPHFitter()
            stage_df = nomo_df[["OS_months", "OS_event", "stage_num"]].dropna()
            if len(stage_df) > 50:
                stage_cph.fit(stage_df, duration_col="OS_months", event_col="OS_event")
                stage_surv = stage_cph.predict_survival_function(stage_df)
                stage_pred_idx = stage_surv.index[np.abs(stage_surv.index - tp_months).argmin()]
                stage_pred = stage_surv.loc[stage_pred_idx].values

                _, nb_stage, _, _ = decision_curve_analysis(
                    stage_df["OS_months"].values, stage_df["OS_event"].values,
                    stage_pred, tp_months
                )
                # Interpolate to match thresholds
                ax.plot(thresholds[:len(nb_stage)], nb_stage, color='green',
                        linewidth=1.5, linestyle='-.', alpha=0.8, label='Stage Only')

        ax.set_xlabel("Threshold Probability")
        ax.set_ylabel("Net Benefit")
        ax.set_title(f"{tp_name} Survival", fontsize=12, fontweight='bold')
        ax.legend(fontsize=8, loc='upper right')
        ax.set_xlim(0, 0.8)
        ax.set_ylim(-0.05, max(0.3, max(nb_model) * 1.2))
        ax.axhline(0, color='black', linewidth=0.5, alpha=0.3)

        # Find range where model has positive net benefit
        pos_range = [t for t, nb in zip(thresholds, nb_model) if nb > 0 and nb > 0]
        if pos_range:
            dca_results.append({
                "timepoint": tp_name,
                "positive_nb_range": f"{min(pos_range):.2f}-{max(pos_range):.2f}",
                "max_nb": max(nb_model),
            })
            print(f"  {tp_name}: Positive net benefit at thresholds {min(pos_range):.2f}-{max(pos_range):.2f}")

plt.suptitle("Decision Curve Analysis — ROS/Ferroptosis Nomogram",
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "dca_curves.png"), dpi=200, bbox_inches='tight')
print("Saved: dca_curves.png")

if dca_results:
    pd.DataFrame(dca_results).to_csv(os.path.join(TABLES, "dca_results.csv"), index=False)
    print("Saved: dca_results.csv")

# ══════════════════════════════════════════════════════════════════════════════
# 5. CALIBRATION WITH PREDICTED PROBABILITIES
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("5. NOMOGRAM CALIBRATION")
print("=" * 70)

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for ax, (tp_name, tp_months) in zip(axes, timepoints.items()):
    if tp_months <= surv_preds.index.max():
        closest_idx = surv_preds.index[np.abs(surv_preds.index - tp_months).argmin()]
        pred_surv = surv_preds.loc[closest_idx].values
        pred_event = 1 - pred_surv

        # Bin patients by predicted risk (quintiles)
        try:
            groups = pd.qcut(pred_event, q=5, labels=False, duplicates='drop')
        except:
            groups = pd.qcut(pred_event, q=4, labels=False, duplicates='drop')

        predicted_means = []
        observed_means = []
        ci_lower = []
        ci_upper = []

        for g in sorted(pd.Series(groups).dropna().unique()):
            mask = groups == g
            grp_time = nomo_df["OS_months"].values[mask]
            grp_event = nomo_df["OS_event"].values[mask]
            grp_pred = pred_event[mask]

            predicted_means.append(grp_pred.mean())

            # Observed event rate using KM
            kmf = KaplanMeierFitter()
            kmf.fit(grp_time, grp_event)
            surv_km = kmf.predict(tp_months)
            obs_event = 1 - surv_km
            observed_means.append(obs_event)

            # Bootstrap CI
            boot_obs = []
            for _ in range(200):
                idx = np.random.choice(len(grp_time), len(grp_time), replace=True)
                kmf_b = KaplanMeierFitter()
                try:
                    kmf_b.fit(grp_time[idx], grp_event[idx])
                    boot_obs.append(1 - kmf_b.predict(tp_months))
                except:
                    pass
            if boot_obs:
                ci_lower.append(np.percentile(boot_obs, 2.5))
                ci_upper.append(np.percentile(boot_obs, 97.5))
            else:
                ci_lower.append(obs_event)
                ci_upper.append(obs_event)

        # Plot
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Ideal')
        ax.errorbar(predicted_means, observed_means,
                     yerr=[np.array(observed_means) - np.array(ci_lower),
                           np.array(ci_upper) - np.array(observed_means)],
                     fmt='o-', color='steelblue', capsize=4, markersize=8,
                     label='Nomogram')
        ax.set_xlabel("Predicted Event Probability")
        ax.set_ylabel("Observed Event Rate")
        ax.set_title(f"{tp_name} Calibration", fontsize=12, fontweight='bold')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.legend(fontsize=9)
        ax.set_aspect('equal')

plt.suptitle("Nomogram Calibration Curves", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "nomogram_calibration.png"), dpi=200, bbox_inches='tight')
print("Saved: nomogram_calibration.png")

# ══════════════════════════════════════════════════════════════════════════════
# 6. C-INDEX COMPARISON (FORMAL BOOTSTRAP TEST)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("6. FORMAL C-INDEX COMPARISON (paired bootstrap)")
print("=" * 70)

# Compare: risk_score_only vs nomogram vs stage_only
models_to_compare = {
    "Risk Score Only": nomo_df[["OS_months", "OS_event", "risk_score"]].dropna(),
    "Nomogram (full)": nomo_df[["OS_months", "OS_event"] + nomo_vars].dropna(),
}
if "stage_num" in nomo_vars:
    models_to_compare["Stage Only"] = nomo_df[["OS_months", "OS_event", "stage_num"]].dropna()

boot_cindices = {name: [] for name in models_to_compare}
n_boot = 1000

for b in range(n_boot):
    idx = np.random.choice(len(nomo_df), len(nomo_df), replace=True)
    boot_data = nomo_df.iloc[idx]

    for name, cols_df in models_to_compare.items():
        try:
            cph_b = CoxPHFitter()
            pred_vars = [c for c in cols_df.columns if c not in ["OS_months", "OS_event"]]
            boot_sub = boot_data[["OS_months", "OS_event"] + pred_vars].dropna()
            cph_b.fit(boot_sub, duration_col="OS_months", event_col="OS_event")
            pred = cph_b.predict_partial_hazard(boot_sub)
            ci = concordance_index(boot_sub["OS_months"], -pred.values.ravel(), boot_sub["OS_event"])
            boot_cindices[name].append(ci)
        except:
            pass

    if (b + 1) % 200 == 0:
        print(f"  Bootstrap {b+1}/{n_boot}")

print("\nModel comparison (bootstrap C-index):")
comparison_results = []
for name, cis in boot_cindices.items():
    if cis:
        mean_ci = np.mean(cis)
        lo, hi = np.percentile(cis, [2.5, 97.5])
        comparison_results.append({
            "model": name, "c_index_mean": mean_ci,
            "c_index_lo": lo, "c_index_hi": hi,
        })
        print(f"  {name:<25} C-index: {mean_ci:.4f} ({lo:.4f}-{hi:.4f})")

# Paired difference test
if "Risk Score Only" in boot_cindices and "Nomogram (full)" in boot_cindices:
    n_paired = min(len(boot_cindices["Risk Score Only"]), len(boot_cindices["Nomogram (full)"]))
    diff = np.array(boot_cindices["Nomogram (full)"][:n_paired]) - \
           np.array(boot_cindices["Risk Score Only"][:n_paired])
    p_diff = (diff <= 0).mean()  # One-sided: nomogram better?
    print(f"\n  Nomogram vs Risk Score Only: delta C-index = {np.mean(diff):.4f}, "
          f"p = {p_diff:.4f} (one-sided)")

if comparison_results:
    pd.DataFrame(comparison_results).to_csv(os.path.join(TABLES, "cindex_comparison.csv"), index=False)
    print("Saved: cindex_comparison.csv")

print(f"\n✓ Nomogram and DCA complete.")
