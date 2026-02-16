"""
08_address_limitations.py — Pre-emptively address reviewer concerns

1. European/White subset validation (TCGA race via GDC API)
2. Recurrence-free survival (RFS) on GSE76427 if available
3. PH assumption test (Schoenfeld residuals)
4. Demographics table (TCGA vs GSE14520 vs LIRI-JP)
"""
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import requests
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

# ── Load data ───────────────────────────────────────────────────────────────
merged = pd.read_csv(os.path.join(DATA, "tcga_ros_merged.csv"))
with open(os.path.join(MODEL, "lasso_model.json")) as f:
    model = json.load(f)

selected_genes = model["genes"]
gene_means = model["gene_means"]
gene_stds = model["gene_stds"]

# Ensure risk_score
if "risk_score" not in merged.columns:
    risk = np.zeros(len(merged))
    for gene, coef in selected_genes.items():
        risk += coef * ((merged[gene] - gene_means[gene]) / gene_stds[gene]).values
    merged["risk_score"] = risk

df = merged.dropna(subset=["OS_months", "OS_event", "risk_score"]).copy()
df = df[df["OS_months"] > 0]


# ══════════════════════════════════════════════════════════════════════════════
# 1. EUROPEAN/WHITE SUBSET VALIDATION
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("1. EUROPEAN/WHITE SUBSET VALIDATION")
print("=" * 70)

# Download race data from GDC API
race_data = None
try:
    endpoint = "https://api.gdc.cancer.gov/cases"
    fields = [
        "submitter_id",
        "demographic.race",
        "demographic.ethnicity",
    ]
    params = {
        "filters": json.dumps({
            "op": "in",
            "content": {
                "field": "project.project_id",
                "value": ["TCGA-LIHC"]
            }
        }),
        "fields": ",".join(fields),
        "size": 500,
        "format": "JSON"
    }

    resp = requests.get(endpoint, params=params, timeout=30)
    if resp.status_code == 200:
        hits = resp.json()["data"]["hits"]
        race_records = []
        for h in hits:
            pid = h.get("submitter_id", "")
            demo = h.get("demographic", {})
            race_records.append({
                "patientId": pid,
                "race": demo.get("race", "unknown"),
                "ethnicity": demo.get("ethnicity", "unknown"),
            })
        race_df = pd.DataFrame(race_records)
        print(f"  Downloaded race data for {len(race_df)} patients")
        print(f"  Race distribution:\n{race_df['race'].value_counts().to_string()}")

        # Merge
        df_race = pd.merge(df, race_df, on="patientId", how="left")

        # European/White subset
        white = df_race[df_race["race"].str.lower().str.contains("white", na=False)]
        asian = df_race[df_race["race"].str.lower().str.contains("asian", na=False)]

        european_results = []
        for name, subset in [("White/European", white), ("Asian", asian)]:
            if len(subset) >= 20 and subset["OS_event"].sum() >= 5:
                ci = concordance_index(subset["OS_months"], -subset["risk_score"], subset["OS_event"])

                # Bootstrap CI
                boot = []
                for _ in range(500):
                    idx = np.random.choice(len(subset), len(subset), replace=True)
                    bd = subset.iloc[idx]
                    try:
                        boot.append(concordance_index(bd["OS_months"], -bd["risk_score"], bd["OS_event"]))
                    except:
                        pass
                ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5]) if boot else (np.nan, np.nan)

                # KM
                med = subset["risk_score"].median()
                h = subset[subset["risk_score"] >= med]
                l = subset[subset["risk_score"] < med]
                lr = logrank_test(h["OS_months"], l["OS_months"],
                                  event_observed_A=h["OS_event"], event_observed_B=l["OS_event"])

                try:
                    cph = CoxPHFitter()
                    cph.fit(subset[["OS_months", "OS_event", "risk_score"]],
                            duration_col="OS_months", event_col="OS_event")
                    hr = np.exp(cph.params_["risk_score"])
                except:
                    hr = np.nan

                european_results.append({
                    "subset": name, "n": len(subset), "events": int(subset["OS_event"].sum()),
                    "c_index": ci, "ci_lo": ci_lo, "ci_hi": ci_hi,
                    "HR": hr, "logrank_p": lr.p_value,
                })
                print(f"  {name}: n={len(subset)}, events={int(subset['OS_event'].sum())}, "
                      f"C-index={ci:.3f} ({ci_lo:.3f}-{ci_hi:.3f}), HR={hr:.2f}, p={lr.p_value:.4f}")

        if european_results:
            eur_df = pd.DataFrame(european_results)
            eur_df.to_csv(os.path.join(TABLES, "european_validation.csv"), index=False)
            print("  Saved: european_validation.csv")

            # KM plot for European subset
            if len(white) >= 20:
                fig, ax = plt.subplots(figsize=(8, 6))
                med = white["risk_score"].median()
                h = white[white["risk_score"] >= med]
                l = white[white["risk_score"] < med]
                lr_w = logrank_test(h["OS_months"], l["OS_months"],
                                    event_observed_A=h["OS_event"], event_observed_B=l["OS_event"])
                kmf = KaplanMeierFitter()
                kmf.fit(h["OS_months"], h["OS_event"], label=f"High Risk (n={len(h)})")
                kmf.plot_survival_function(ax=ax, ci_show=True, color="red")
                kmf.fit(l["OS_months"], l["OS_event"], label=f"Low Risk (n={len(l)})")
                kmf.plot_survival_function(ax=ax, ci_show=True, color="blue")
                ci_w = concordance_index(white["OS_months"], -white["risk_score"], white["OS_event"])
                ax.set_title(f"White/European Subset (n={len(white)})\nC-index={ci_w:.3f}, p={lr_w.p_value:.4f}",
                             fontsize=12, fontweight='bold')
                ax.set_xlabel("Overall Survival (months)")
                ax.set_ylabel("Survival Probability")
                ax.legend()
                ax.set_xlim(0, 100)
                plt.tight_layout()
                plt.savefig(os.path.join(FIGS, "european_km.png"), dpi=200, bbox_inches='tight')
                print("  Saved: european_km.png")

        race_data = race_df
except Exception as e:
    print(f"  GDC API request failed: {e}")
    print("  Proceeding without race-stratified analysis")


# ══════════════════════════════════════════════════════════════════════════════
# 2. RECURRENCE-FREE SURVIVAL (GSE76427)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("2. RECURRENCE-FREE SURVIVAL")
print("=" * 70)

# GSE76427 may have RFS data in additional columns
gse76_path = os.path.join(GEO, "GSE76427_series_matrix.txt.gz")
if os.path.exists(gse76_path):
    import gzip
    content = gzip.open(gse76_path, 'rt', errors='replace').read()
    lines = content.split('\n')

    # Check for RFS-related characteristics
    rfs_found = False
    for line in lines:
        if line.startswith('!Sample_characteristics_ch1'):
            lower = line.lower()
            if 'recurrence' in lower or 'rfs' in lower or 'dfs' in lower or 'relapse' in lower:
                rfs_found = True
                print(f"  Found RFS-related data: {line[:200]}")
                break

    if not rfs_found:
        print("  No RFS data found in GSE76427 series matrix")
else:
    print("  GSE76427 data not available locally")


# ══════════════════════════════════════════════════════════════════════════════
# 3. PROPORTIONAL HAZARDS ASSUMPTION TEST
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("3. PROPORTIONAL HAZARDS ASSUMPTION")
print("=" * 70)

# Test PH assumption using Schoenfeld residuals
try:
    cph_ph = CoxPHFitter()
    cph_ph.fit(df[["OS_months", "OS_event", "risk_score"]], duration_col="OS_months", event_col="OS_event")

    # Test PH assumption
    ph_test = cph_ph.check_assumptions(df[["OS_months", "OS_event", "risk_score"]],
                                        p_value_threshold=0.05, show_plots=False)
    print(f"  PH assumption test completed")

    # Manual Schoenfeld residual test
    # Compute scaled Schoenfeld residuals correlation with time
    schoenfeld = cph_ph.compute_residuals(df[["OS_months", "OS_event", "risk_score"]], kind="schoenfeld")
    if len(schoenfeld) > 0:
        rho, p_sch = stats.spearmanr(schoenfeld.index, schoenfeld["risk_score"])
        print(f"  Schoenfeld residual test: rho={rho:.4f}, p={p_sch:.4f}")
        if p_sch > 0.05:
            print("  → PH assumption NOT violated (p > 0.05)")
        else:
            print("  → PH assumption may be violated (p < 0.05)")

        # Plot
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.scatter(schoenfeld.index, schoenfeld["risk_score"], alpha=0.4, s=20, color='steelblue')
        # LOWESS smoothing
        z = np.polyfit(schoenfeld.index.astype(float), schoenfeld["risk_score"].values, 2)
        p_fit = np.poly1d(z)
        x_smooth = np.linspace(schoenfeld.index.min(), schoenfeld.index.max(), 100)
        ax.plot(x_smooth, p_fit(x_smooth), 'r-', linewidth=2, label='Polynomial fit')
        ax.axhline(0, color='black', linestyle='--', alpha=0.5)
        ax.set_xlabel("Time (months)")
        ax.set_ylabel("Scaled Schoenfeld Residual")
        ax.set_title(f"Proportional Hazards Test\nSchoenfeld rho={rho:.3f}, p={p_sch:.3f}", fontweight='bold')
        ax.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(FIGS, "schoenfeld_residuals.png"), dpi=200, bbox_inches='tight')
        print("  Saved: schoenfeld_residuals.png")

except Exception as e:
    print(f"  PH test failed: {e}")
    # Fallback: simple log-log plot
    fig, ax = plt.subplots(figsize=(8, 6))
    med = df["risk_score"].median()
    for label, subset, color in [("High Risk", df[df["risk_score"] >= med], "red"),
                                  ("Low Risk", df[df["risk_score"] < med], "blue")]:
        kmf = KaplanMeierFitter()
        kmf.fit(subset["OS_months"], subset["OS_event"])
        sf = kmf.survival_function_
        sf = sf[sf.values > 0]
        ax.plot(np.log(sf.index + 1), np.log(-np.log(sf.values.flatten())), color=color, label=label)
    ax.set_xlabel("log(time)")
    ax.set_ylabel("log(-log(S(t)))")
    ax.set_title("Log-Log Plot — PH Assumption Check", fontweight='bold')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "loglog_plot.png"), dpi=200, bbox_inches='tight')
    print("  Saved: loglog_plot.png")


# ══════════════════════════════════════════════════════════════════════════════
# 4. DEMOGRAPHICS TABLE
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("4. DEMOGRAPHICS TABLE")
print("=" * 70)

def encode_stage(s):
    if pd.isna(s): return np.nan
    s = str(s).upper().strip()
    if "IV" in s: return 4
    if "III" in s: return 3
    if "II" in s: return 2
    if "I" in s: return 1
    return np.nan

def compute_demographics(name, data, os_col="OS_months", event_col="OS_event",
                         age_col=None, sex_col=None, stage_col=None):
    """Compute cohort demographics."""
    demo = {"Cohort": name, "N": len(data)}

    # Events
    if event_col in data.columns:
        demo["Events"] = int(data[event_col].sum())
        demo["Event_rate"] = f"{data[event_col].mean():.1%}"

    # Median follow-up
    if os_col in data.columns:
        demo["Median_OS_months"] = f"{data[os_col].median():.1f}"

    # Age
    if age_col and age_col in data.columns:
        age = data[age_col]
        if age.median() > 200:
            age = age / 365.25
        demo["Age_median"] = f"{age.median():.0f}"
        demo["Age_range"] = f"{age.min():.0f}-{age.max():.0f}"

    # Sex
    if sex_col and sex_col in data.columns:
        sex = data[sex_col].astype(str).str.lower()
        n_male = sex.isin(["male", "m", "1"]).sum()
        demo["Male_n"] = n_male
        demo["Male_pct"] = f"{n_male / len(data) * 100:.0f}%"

    # Stage
    if stage_col and stage_col in data.columns:
        stages = data[stage_col].apply(encode_stage)
        for s in [1, 2, 3, 4]:
            n_s = (stages == s).sum()
            if n_s > 0:
                demo[f"Stage_{s}"] = n_s

    return demo

# TCGA demographics
demos = []
demos.append(compute_demographics(
    "TCGA-LIHC", df, age_col="age_at_diagnosis",
    sex_col="gender", stage_col="tumor_stage"
))

# If validation data available, add their demographics with clinical info
val_results = os.path.join(TABLES, "validation_results.csv")
if os.path.exists(val_results):
    val_df = pd.read_csv(val_results)

    # Try to load GSE14520 supplement for demographics
    gse14520_suppl_path = os.path.join(GEO, "GSE14520_Extra_Supplement.txt.gz")
    if os.path.exists(gse14520_suppl_path):
        try:
            suppl14 = pd.read_csv(gse14520_suppl_path, sep='\t', compression='gzip')
            # Filter to tumor samples with Affy GSM
            suppl14 = suppl14.dropna(subset=["Affy_GSM", "Survival months", "Survival status"])
            n14 = len(suppl14)
            demo14 = {"Cohort": "GSE14520", "N": n14,
                       "Events": int(suppl14["Survival status"].sum()),
                       "Event_rate": f"{suppl14['Survival status'].mean():.1%}",
                       "Median_OS_months": f"{suppl14['Survival months'].median():.1f}"}
            if "Age" in suppl14.columns:
                age14 = pd.to_numeric(suppl14["Age"], errors="coerce")
                demo14["Age_median"] = f"{age14.median():.0f}"
                demo14["Age_range"] = f"{age14.min():.0f}-{age14.max():.0f}"
            if "Gender" in suppl14.columns:
                n_male14 = (suppl14["Gender"].str.upper().str.strip() == "M").sum()
                demo14["Male_n"] = n_male14
                demo14["Male_pct"] = f"{n_male14 / n14 * 100:.0f}%"
            if "TNM staging" in suppl14.columns:
                tnm = suppl14["TNM staging"].apply(encode_stage)
                for s in [1, 2, 3, 4]:
                    n_s = (tnm == s).sum()
                    if n_s > 0:
                        demo14[f"Stage_{s}"] = n_s
            demos.append(demo14)
        except Exception as e:
            print(f"  Could not parse GSE14520 demographics: {e}")

    # Try to load ICGC LIRI-JP cached donors for demographics
    icgc_donor_cache = os.path.join(GEO, "liri_jp_donors.csv")
    if os.path.exists(icgc_donor_cache):
        try:
            icgc_d = pd.read_csv(icgc_donor_cache, index_col=0)
            n_icgc = len(icgc_d)
            demo_icgc = {"Cohort": "ICGC LIRI-JP", "N": n_icgc}
            if "os_event" in icgc_d.columns:
                demo_icgc["Events"] = int(icgc_d["os_event"].sum())
                demo_icgc["Event_rate"] = f"{icgc_d['os_event'].mean():.1%}"
            if "os_months" in icgc_d.columns:
                demo_icgc["Median_OS_months"] = f"{icgc_d['os_months'].median():.1f}"
            if "age_years" in icgc_d.columns:
                age_icgc = icgc_d["age_years"].dropna()
                demo_icgc["Age_median"] = f"{age_icgc.median():.0f}"
                demo_icgc["Age_range"] = f"{age_icgc.min():.0f}-{age_icgc.max():.0f}"
            if "sex_binary" in icgc_d.columns:
                n_male_icgc = int(icgc_d["sex_binary"].sum())
                demo_icgc["Male_n"] = n_male_icgc
                demo_icgc["Male_pct"] = f"{n_male_icgc / n_icgc * 100:.0f}%"
            if "stage_numeric" in icgc_d.columns:
                for s in [1, 2, 3, 4]:
                    n_s = int((icgc_d["stage_numeric"] == s).sum())
                    if n_s > 0:
                        demo_icgc[f"Stage_{s}"] = n_s
            demos.append(demo_icgc)
        except Exception as e:
            print(f"  Could not parse ICGC demographics: {e}")

    # Add remaining cohorts without detailed demographics
    added_cohorts = {d["Cohort"] for d in demos}
    for _, row in val_df.iterrows():
        if row["cohort"] not in added_cohorts:
            demos.append({
                "Cohort": row["cohort"],
                "N": int(row["n"]),
                "Events": int(row["events"]),
                "Event_rate": f"{row['events'] / row['n']:.1%}",
            })

demo_df = pd.DataFrame(demos)
demo_df.to_csv(os.path.join(TABLES, "demographics.csv"), index=False)
print(demo_df.to_string(index=False))
print("\nSaved: demographics.csv")


# ══════════════════════════════════════════════════════════════════════════════
# 5. ADDITIONAL ROBUSTNESS: TERTILE STRATIFICATION
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("5. TERTILE STRATIFICATION")
print("=" * 70)

tertiles = np.percentile(df["risk_score"], [33.33, 66.67])
df["risk_tertile"] = pd.cut(df["risk_score"],
                             bins=[-np.inf, tertiles[0], tertiles[1], np.inf],
                             labels=["Low", "Medium", "High"])

fig, ax = plt.subplots(figsize=(8, 6))
colors_t = {"Low": "blue", "Medium": "green", "High": "red"}
kmf = KaplanMeierFitter()

for tert in ["High", "Medium", "Low"]:
    sub = df[df["risk_tertile"] == tert]
    kmf.fit(sub["OS_months"], sub["OS_event"], label=f"{tert} (n={len(sub)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color=colors_t[tert])

# Overall log-rank across 3 groups
from lifelines.statistics import multivariate_logrank_test
try:
    mlr = multivariate_logrank_test(df["OS_months"], df["risk_tertile"], df["OS_event"])
    p_tert = mlr.p_value
except:
    p_tert = np.nan

ax.set_title(f"Tertile Stratification — TCGA-LIHC\nMultivariate log-rank p={p_tert:.2e}", fontweight='bold')
ax.set_xlabel("Overall Survival (months)")
ax.set_ylabel("Survival Probability")
ax.legend()
ax.set_xlim(0, 100)
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "tertile_km.png"), dpi=200, bbox_inches='tight')
print(f"Tertile log-rank p = {p_tert:.2e}")
print("Saved: tertile_km.png")

print(f"\n✓ Limitations analyses complete.")
