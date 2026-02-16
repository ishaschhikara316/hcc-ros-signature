"""
04_training_evaluation.py — Full evaluation on TCGA training cohort

C-index (bootstrapped CI), KM curves, time-dependent ROC, multivariate Cox,
subgroup analysis, internal train-test split, heatmap, nomogram, calibration.
"""
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
import json
import os
import warnings
warnings.filterwarnings('ignore')

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data", "tcga")
MODEL = os.path.join(BASE, "results", "model")
TABLES = os.path.join(BASE, "results", "tables")
FIGS = os.path.join(BASE, "results", "figures")

# ── Load ────────────────────────────────────────────────────────────────────
merged = pd.read_csv(os.path.join(DATA, "tcga_ros_merged.csv"))
with open(os.path.join(MODEL, "lasso_model.json")) as f:
    model = json.load(f)

selected_genes = model["genes"]
gene_means = model["gene_means"]
gene_stds = model["gene_stds"]
print(f"Model: {len(selected_genes)} genes, training C-index={model['c_index_train']:.4f}")

# Compute risk score if not present
if "risk_score" not in merged.columns:
    risk = np.zeros(len(merged))
    for gene, coef in selected_genes.items():
        z = (merged[gene] - gene_means[gene]) / gene_stds[gene]
        risk += coef * z.values
    merged["risk_score"] = risk

df = merged.dropna(subset=["OS_months", "OS_event", "risk_score"]).copy()
df = df[df["OS_months"] > 0]
print(f"Evaluating on {len(df)} patients ({int(df['OS_event'].sum())} events)")

# ── 1. C-index with bootstrap CI ───────────────────────────────────────────
print("\n" + "=" * 60)
print("1. C-INDEX")
ci = concordance_index(df["OS_months"], -df["risk_score"], df["OS_event"])
boot_cis = []
for _ in range(1000):
    idx = np.random.choice(len(df), len(df), replace=True)
    bd = df.iloc[idx]
    try:
        boot_cis.append(concordance_index(bd["OS_months"], -bd["risk_score"], bd["OS_event"]))
    except:
        pass
ci_lo, ci_hi = np.percentile(boot_cis, [2.5, 97.5])
print(f"C-index: {ci:.4f} (95% CI: {ci_lo:.4f}-{ci_hi:.4f})")

# ── 2. KM curves (median split) ────────────────────────────────────────────
print("\n" + "=" * 60)
print("2. KAPLAN-MEIER ANALYSIS")
med = df["risk_score"].median()
high = df[df["risk_score"] >= med]
low = df[df["risk_score"] < med]
lr = logrank_test(high["OS_months"], low["OS_months"],
                  event_observed_A=high["OS_event"], event_observed_B=low["OS_event"])

cph_rs = CoxPHFitter()
cph_rs.fit(df[["OS_months", "OS_event", "risk_score"]], duration_col="OS_months", event_col="OS_event")
hr = np.exp(cph_rs.params_["risk_score"])
hr_ci = np.exp(cph_rs.confidence_intervals_.values[0])

print(f"High risk: n={len(high)}, events={int(high['OS_event'].sum())}")
print(f"Low risk:  n={len(low)}, events={int(low['OS_event'].sum())}")
print(f"Log-rank p: {lr.p_value:.2e}")
print(f"HR: {hr:.3f} (95% CI: {hr_ci[0]:.3f}-{hr_ci[1]:.3f})")

fig, ax = plt.subplots(figsize=(8, 6))
kmf = KaplanMeierFitter()
kmf.fit(high["OS_months"], high["OS_event"], label=f"High Risk (n={len(high)})")
kmf.plot_survival_function(ax=ax, ci_show=True, color="red")
kmf.fit(low["OS_months"], low["OS_event"], label=f"Low Risk (n={len(low)})")
kmf.plot_survival_function(ax=ax, ci_show=True, color="blue")
ax.set_title(f"TCGA-LIHC Training Cohort\nHR={hr:.2f}, Log-rank p={lr.p_value:.2e}", fontsize=13, fontweight='bold')
ax.set_xlabel("Overall Survival (months)")
ax.set_ylabel("Survival Probability")
ax.legend(fontsize=11)
ax.set_xlim(0, 100)
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "training_km.png"), dpi=200, bbox_inches='tight')
print("Saved: training_km.png")

# ── 3. Time-dependent ROC ──────────────────────────────────────────────────
print("\n" + "=" * 60)
print("3. TIME-DEPENDENT ROC")

def time_dependent_roc(time, event, score, timepoint, n_thresh=200):
    cases = (time <= timepoint) & (event == 1)
    controls = time > timepoint
    if cases.sum() == 0 or controls.sum() == 0:
        return 0.5, [0, 1], [0, 1]
    thresholds = np.percentile(score, np.linspace(0, 100, n_thresh))
    tpr_list, fpr_list = [], []
    for t in thresholds:
        pred_pos = score >= t
        tpr_list.append((pred_pos & cases).sum() / cases.sum())
        fpr_list.append((pred_pos & controls).sum() / controls.sum())
    # Sort by FPR
    order = np.argsort(fpr_list)
    fpr_s = np.array(fpr_list)[order]
    tpr_s = np.array(tpr_list)[order]
    auc = np.trapezoid(tpr_s, fpr_s)
    return auc, fpr_s, tpr_s

timepoints = {"1-year": 12, "3-year": 36, "5-year": 60}
fig, ax = plt.subplots(figsize=(7, 6))
colors = {"1-year": "green", "3-year": "blue", "5-year": "red"}

roc_results = {}
for name, tp in timepoints.items():
    auc, fpr, tpr = time_dependent_roc(df["OS_months"].values, df["OS_event"].values,
                                        df["risk_score"].values, tp)
    ax.plot(fpr, tpr, color=colors[name], linewidth=2, label=f"{name} AUC={auc:.3f}")
    roc_results[name] = auc
    print(f"  {name}: AUC = {auc:.3f}")

ax.plot([0, 1], [0, 1], 'k--', alpha=0.5)
ax.set_xlabel("False Positive Rate", fontsize=11)
ax.set_ylabel("True Positive Rate", fontsize=11)
ax.set_title("Time-Dependent ROC — TCGA Training", fontsize=13, fontweight='bold')
ax.legend(fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "roc_curves.png"), dpi=200, bbox_inches='tight')
print("Saved: roc_curves.png")

# ── 4. Multivariate Cox ────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("4. MULTIVARIATE COX")

def encode_stage(s):
    if pd.isna(s):
        return np.nan
    s = str(s).upper().strip()
    if "IV" in s: return 4
    if "III" in s: return 3
    if "II" in s: return 2
    if "I" in s: return 1
    return np.nan

df["stage_num"] = df["tumor_stage"].apply(encode_stage) if "tumor_stage" in df.columns else np.nan
df["male"] = (df["gender"].str.lower() == "male").astype(int) if "gender" in df.columns else np.nan

# Age in years
if "age_at_diagnosis" in df.columns:
    age = df["age_at_diagnosis"]
    df["age_years"] = np.where(age > 200, age / 365.25, age)
elif "age" in df.columns:
    df["age_years"] = df["age"]

# Grade
grade_map = {"G1": 1, "G2": 2, "G3": 3, "G4": 4}
if "tumor_grade" in df.columns:
    df["grade_num"] = df["tumor_grade"].map(grade_map)

mv_results = []

# Model 1: risk_score alone
cph1 = CoxPHFitter()
d1 = df[["OS_months", "OS_event", "risk_score"]].dropna()
cph1.fit(d1, duration_col="OS_months", event_col="OS_event")
print("\nModel 1: Risk score only")
print(cph1.summary[["coef", "exp(coef)", "p"]])

# Model 2: risk_score + age + sex + stage
mv_cols = ["OS_months", "OS_event", "risk_score"]
labels = ["risk_score"]
for col, name in [("age_years", "age"), ("male", "sex"), ("stage_num", "stage")]:
    if col in df.columns and df[col].notna().sum() > 100:
        mv_cols.append(col)
        labels.append(name)

if "grade_num" in df.columns and df["grade_num"].notna().sum() > 100:
    mv_cols.append("grade_num")
    labels.append("grade")

d2 = df[mv_cols].dropna()
if len(d2) >= 50:
    cph2 = CoxPHFitter()
    cph2.fit(d2, duration_col="OS_months", event_col="OS_event")
    print(f"\nModel 2: Risk score + {', '.join(labels[1:])} (n={len(d2)})")
    print(cph2.summary[["coef", "exp(coef)", "p"]])

    # Save multivariate table
    mv_table = cph2.summary[["coef", "exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].copy()
    mv_table.columns = ["coef", "HR", "HR_lower", "HR_upper", "p_value"]
    mv_table.to_csv(os.path.join(TABLES, "multivariate_cox.csv"))
    print("Saved: multivariate_cox.csv")

    # Forest plot for multivariate
    fig, ax = plt.subplots(figsize=(8, max(3, len(mv_table) * 0.6 + 1)))
    y_pos = range(len(mv_table))
    for i, (var, row) in enumerate(mv_table.iterrows()):
        color = 'red' if row['p_value'] < 0.05 else 'gray'
        ax.plot([row['HR_lower'], row['HR_upper']], [i, i], color=color, linewidth=2)
        ax.scatter(row['HR'], i, color=color, s=60, zorder=5)
    ax.axvline(1.0, color='black', linestyle='--')
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(mv_table.index)
    ax.set_xlabel("Hazard Ratio (95% CI)")
    ax.set_title("Multivariate Cox Regression — TCGA-LIHC", fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "multivariate_forest.png"), dpi=200, bbox_inches='tight')
    print("Saved: multivariate_forest.png")

# ── 5. Subgroup analysis ───────────────────────────────────────────────────
print("\n" + "=" * 60)
print("5. SUBGROUP ANALYSIS")

subgroups = {}
if "stage_num" in df.columns:
    subgroups["Stage I-II"] = df[df["stage_num"].isin([1, 2])]
    subgroups["Stage III-IV"] = df[df["stage_num"].isin([3, 4])]
if "age_years" in df.columns:
    subgroups["Age ≤60"] = df[df["age_years"] <= 60]
    subgroups["Age >60"] = df[df["age_years"] > 60]
if "male" in df.columns:
    subgroups["Male"] = df[df["male"] == 1]
    subgroups["Female"] = df[df["male"] == 0]

sub_results = []
for name, sub_df in subgroups.items():
    sub_df = sub_df.dropna(subset=["OS_months", "OS_event", "risk_score"])
    if len(sub_df) < 20 or sub_df["OS_event"].sum() < 5:
        continue
    ci_sub = concordance_index(sub_df["OS_months"], -sub_df["risk_score"], sub_df["OS_event"])
    med_sub = sub_df["risk_score"].median()
    h = sub_df[sub_df["risk_score"] >= med_sub]
    l = sub_df[sub_df["risk_score"] < med_sub]
    if len(h) > 5 and len(l) > 5 and h["OS_event"].sum() > 0 and l["OS_event"].sum() > 0:
        lr_sub = logrank_test(h["OS_months"], l["OS_months"],
                              event_observed_A=h["OS_event"], event_observed_B=l["OS_event"])
        try:
            cph_sub = CoxPHFitter()
            cph_sub.fit(sub_df[["OS_months", "OS_event", "risk_score"]],
                        duration_col="OS_months", event_col="OS_event")
            hr_sub = np.exp(cph_sub.params_["risk_score"])
            hr_ci_sub = np.exp(cph_sub.confidence_intervals_.values[0])
        except:
            hr_sub, hr_ci_sub = np.nan, [np.nan, np.nan]

        sub_results.append({
            "subgroup": name, "n": len(sub_df), "events": int(sub_df["OS_event"].sum()),
            "c_index": ci_sub, "HR": hr_sub, "HR_lower": hr_ci_sub[0], "HR_upper": hr_ci_sub[1],
            "logrank_p": lr_sub.p_value,
        })
        print(f"  {name}: n={len(sub_df)}, C-index={ci_sub:.3f}, HR={hr_sub:.2f}, p={lr_sub.p_value:.4f}")

if sub_results:
    sub_df_out = pd.DataFrame(sub_results)
    sub_df_out.to_csv(os.path.join(TABLES, "subgroup_analysis.csv"), index=False)
    print("Saved: subgroup_analysis.csv")

# ── 6. Internal train-test split (70/30, stratified by event) ──────────────
print("\n" + "=" * 60)
print("6. INTERNAL TRAIN-TEST SPLIT (70/30, stratified)")

train_df, test_df = train_test_split(df, test_size=0.3, random_state=42,
                                      stratify=df["OS_event"])

# Re-compute risk scores on test set using training set normalization
sig_genes_list = list(selected_genes.keys())
train_means = {g: train_df[g].mean() for g in sig_genes_list}
train_stds = {g: train_df[g].std() for g in sig_genes_list}

test_risk = np.zeros(len(test_df))
for gene, coef in selected_genes.items():
    z = (test_df[gene].values - train_means[gene]) / (train_stds[gene] + 1e-10)
    test_risk += coef * z
test_df = test_df.copy()
test_df["risk_score_norm"] = test_risk

ci_trn = concordance_index(train_df["OS_months"], -train_df["risk_score"], train_df["OS_event"])
ci_tst = concordance_index(test_df["OS_months"], -test_df["risk_score_norm"], test_df["OS_event"])
print(f"Training (70%): n={len(train_df)}, events={int(train_df['OS_event'].sum())}, C-index={ci_trn:.4f}")
print(f"Testing  (30%): n={len(test_df)}, events={int(test_df['OS_event'].sum())}, C-index={ci_tst:.4f}")
print(f"  (Test set normalized using training set mean/std)")

# KM for train-test
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax, (name, sub, rs_col) in zip(axes, [("Training", train_df, "risk_score"),
                                            ("Testing", test_df, "risk_score_norm")]):
    m = sub[rs_col].median()
    h = sub[sub[rs_col] >= m]
    l = sub[sub[rs_col] < m]
    lr_s = logrank_test(h["OS_months"], l["OS_months"],
                        event_observed_A=h["OS_event"], event_observed_B=l["OS_event"])
    kmf = KaplanMeierFitter()
    kmf.fit(h["OS_months"], h["OS_event"], label=f"High (n={len(h)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="red")
    kmf.fit(l["OS_months"], l["OS_event"], label=f"Low (n={len(l)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="blue")
    ci_s = concordance_index(sub["OS_months"], -sub[rs_col], sub["OS_event"])
    ax.set_title(f"{name} Set\nC-index={ci_s:.3f}, p={lr_s.p_value:.2e}", fontsize=12, fontweight='bold')
    ax.set_xlabel("Months")
    ax.set_ylabel("Survival Probability")
    ax.legend()
    ax.set_xlim(0, 100)

plt.suptitle("Internal Validation (70/30 Split)", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "train_test_km.png"), dpi=200, bbox_inches='tight')
print("Saved: train_test_km.png")

# ── 7. Heatmap of signature genes ──────────────────────────────────────────
print("\n" + "=" * 60)
print("7. SIGNATURE GENE HEATMAP")

sig_genes = list(selected_genes.keys())
heatmap_data = df[sig_genes].copy()
for g in sig_genes:
    heatmap_data[g] = (heatmap_data[g] - heatmap_data[g].mean()) / heatmap_data[g].std()

# Sort by risk score
order = df["risk_score"].argsort().values
heatmap_data = heatmap_data.iloc[order]

fig, axes = plt.subplots(2, 1, figsize=(14, max(4, len(sig_genes) * 0.5 + 2)),
                          gridspec_kw={'height_ratios': [1, 4]}, sharex=True)

# Risk score bar
ax0 = axes[0]
rs_sorted = df["risk_score"].values[order]
ax0.bar(range(len(rs_sorted)), rs_sorted, color=np.where(rs_sorted >= np.median(rs_sorted), 'red', 'blue'),
        width=1.0, linewidth=0)
ax0.set_ylabel("Risk Score")
ax0.set_title("Signature Gene Expression Heatmap (ordered by risk score)", fontweight='bold')

# Expression heatmap
ax1 = axes[1]
sns.heatmap(heatmap_data.T, ax=ax1, cmap="RdBu_r", center=0, xticklabels=False,
            yticklabels=sig_genes, cbar_kws={"shrink": 0.5, "label": "Z-score"})
ax1.set_xlabel(f"Patients (n={len(df)})")

plt.tight_layout()
plt.savefig(os.path.join(FIGS, "heatmap.png"), dpi=200, bbox_inches='tight')
print("Saved: heatmap.png")

# ── 8. Nomogram (simplified as coefficient plot) ───────────────────────────
print("\n" + "=" * 60)
print("8. NOMOGRAM (coefficient contribution plot)")

# Build multivariate model with risk_score + clinical for nomogram
nomo_cols = ["OS_months", "OS_event", "risk_score"]
if "age_years" in df.columns: nomo_cols.append("age_years")
if "male" in df.columns: nomo_cols.append("male")
if "stage_num" in df.columns: nomo_cols.append("stage_num")

d_nomo = df[nomo_cols].dropna()
if len(d_nomo) > 50:
    cph_nomo = CoxPHFitter()
    cph_nomo.fit(d_nomo, duration_col="OS_months", event_col="OS_event")

    fig, ax = plt.subplots(figsize=(8, 4))
    coefs = cph_nomo.params_.sort_values()
    colors_n = ['red' if v > 0 else 'blue' for v in coefs]
    ax.barh(range(len(coefs)), coefs, color=colors_n, edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(coefs)))
    ax.set_yticklabels(coefs.index)
    ax.set_xlabel("Cox Coefficient")
    ax.set_title("Nomogram — Variable Contribution", fontsize=12, fontweight='bold')
    ax.axvline(0, color='black', linestyle='-', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "nomogram.png"), dpi=200, bbox_inches='tight')
    print("Saved: nomogram.png")

# ── 9. Calibration plot ────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("9. CALIBRATION CURVES")

fig, ax = plt.subplots(figsize=(7, 6))
for tp_name, tp_months in [("1-year", 12), ("3-year", 36), ("5-year", 60)]:
    # Bin patients by predicted risk (quintiles)
    df_cal = df[["OS_months", "OS_event", "risk_score"]].dropna()
    df_cal["risk_group"] = pd.qcut(df_cal["risk_score"], q=5, labels=False, duplicates='drop')

    predicted = []
    observed = []
    for g in sorted(df_cal["risk_group"].unique()):
        grp = df_cal[df_cal["risk_group"] == g]
        # Predicted: mean risk score (higher = worse)
        pred_risk = grp["risk_score"].mean()
        # Observed: 1 - KM survival at timepoint
        kmf_cal = KaplanMeierFitter()
        kmf_cal.fit(grp["OS_months"], grp["OS_event"])
        surv = kmf_cal.predict(tp_months)
        obs_event_rate = 1 - surv
        predicted.append(pred_risk)
        observed.append(obs_event_rate)

    # Normalize predicted to 0-1 range
    pred_norm = (np.array(predicted) - min(predicted)) / (max(predicted) - min(predicted) + 1e-10)
    ax.plot(pred_norm, observed, 'o-', label=f"{tp_name}", markersize=8)

ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Ideal')
ax.set_xlabel("Predicted Risk (normalized)")
ax.set_ylabel("Observed Event Rate")
ax.set_title("Calibration Curves", fontsize=13, fontweight='bold')
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "calibration.png"), dpi=200, bbox_inches='tight')
print("Saved: calibration.png")

# ── Summary table ───────────────────────────────────────────────────────────
perf = pd.DataFrame([{
    "cohort": "TCGA-LIHC (training)",
    "n": len(df),
    "events": int(df["OS_event"].sum()),
    "c_index": ci,
    "c_index_95ci_lo": ci_lo,
    "c_index_95ci_hi": ci_hi,
    "HR": hr,
    "HR_ci_lo": hr_ci[0],
    "HR_ci_hi": hr_ci[1],
    "logrank_p": lr.p_value,
    "AUC_1yr": roc_results.get("1-year", np.nan),
    "AUC_3yr": roc_results.get("3-year", np.nan),
    "AUC_5yr": roc_results.get("5-year", np.nan),
}])
perf.to_csv(os.path.join(TABLES, "training_performance.csv"), index=False)
print("\nSaved: training_performance.csv")
print(f"\n✓ Training evaluation complete.")
