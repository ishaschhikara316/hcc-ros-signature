"""
17_pan_cancer.py — Pan-cancer analysis of ROS/ferroptosis signature

Test the signature across multiple TCGA cancer types where
oxidative stress / ferroptosis biology is implicated:
1. TCGA-KIRC (kidney renal clear cell carcinoma)
2. TCGA-LUSC (lung squamous cell carcinoma)
3. TCGA-STAD (stomach adenocarcinoma)
4. TCGA-PAAD (pancreatic adenocarcinoma)
"""
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
from lifelines.statistics import logrank_test
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import json
import os
import requests
import gzip
import io
import warnings
warnings.filterwarnings('ignore')

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_BASE = os.path.join(BASE, "data")
MODEL = os.path.join(BASE, "results", "model")
TABLES = os.path.join(BASE, "results", "tables")
FIGS = os.path.join(BASE, "results", "figures")

PAN_DATA = os.path.join(DATA_BASE, "pan_cancer")
os.makedirs(PAN_DATA, exist_ok=True)

# ── Load signature ──────────────────────────────────────────────────────────
with open(os.path.join(MODEL, "lasso_model.json")) as f:
    model = json.load(f)
selected_genes = model["genes"]
sig_genes = list(selected_genes.keys())
print(f"Signature genes ({len(sig_genes)}): {sig_genes}")

# TCGA cancer types to test
CANCERS = {
    "KIRC": {"project": "TCGA-KIRC", "name": "Kidney Renal Clear Cell Carcinoma"},
    "LUSC": {"project": "TCGA-LUSC", "name": "Lung Squamous Cell Carcinoma"},
    "STAD": {"project": "TCGA-STAD", "name": "Stomach Adenocarcinoma"},
    "PAAD": {"project": "TCGA-PAAD", "name": "Pancreatic Adenocarcinoma"},
}


def download_tcga_cancer(project_id, cancer_abbr):
    """Download TCGA expression + clinical data via GDC API."""
    expr_cache = os.path.join(PAN_DATA, f"{cancer_abbr}_expression.csv")
    clin_cache = os.path.join(PAN_DATA, f"{cancer_abbr}_clinical.csv")

    if os.path.exists(expr_cache) and os.path.exists(clin_cache):
        print(f"  Loading cached {cancer_abbr}...")
        expr = pd.read_csv(expr_cache, index_col=0)
        clin = pd.read_csv(clin_cache, index_col=0)
        return expr, clin

    gdc = "https://api.gdc.cancer.gov"

    # ── Clinical data ──
    print(f"  Downloading {project_id} clinical data...")
    clin_records = []
    size = 500
    offset = 0
    while True:
        payload = {
            "filters": {
                "op": "and",
                "content": [
                    {"op": "=", "content": {"field": "project.project_id", "value": project_id}},
                    {"op": "=", "content": {"field": "demographic.vital_status", "value": ["alive", "dead"]}},
                ]
            },
            "fields": "submitter_id,demographic.vital_status,demographic.days_to_death,"
                       "diagnoses.days_to_last_follow_up,demographic.gender,demographic.age_at_index,"
                       "diagnoses.tumor_stage",
            "size": size,
            "from": offset,
        }
        resp = requests.post(f"{gdc}/cases", json=payload, timeout=60)
        data = resp.json()
        hits = data.get("data", {}).get("hits", [])
        if not hits:
            break
        clin_records.extend(hits)
        if len(hits) < size:
            break
        offset += size

    # Parse clinical
    clin_rows = []
    for case in clin_records:
        pid = case.get("submitter_id", "")
        demo = case.get("demographic", {})
        diag = case.get("diagnoses", [{}])[0] if case.get("diagnoses") else {}

        vital = demo.get("vital_status", "")
        days_death = demo.get("days_to_death")
        days_fu = diag.get("days_to_last_follow_up")

        if vital.lower() == "dead" and days_death is not None:
            os_months = float(days_death) / 30.44
            os_event = 1
        elif days_fu is not None:
            os_months = float(days_fu) / 30.44
            os_event = 0
        else:
            continue

        if os_months <= 0:
            continue

        clin_rows.append({
            "patient_id": pid,
            "OS_months": os_months,
            "OS_event": os_event,
            "age_years": demo.get("age_at_index"),
            "gender": demo.get("gender"),
            "stage": diag.get("tumor_stage", ""),
        })

    clin = pd.DataFrame(clin_rows).set_index("patient_id")
    print(f"  Clinical: {len(clin)} patients, {int(clin['OS_event'].sum())} events")

    # ── Expression data via cBioPortal (pre-aggregated, fast) ──
    # cBioPortal study ID: e.g. TCGA-KIRC -> kirc_tcga
    cancer_code = project_id.split("-")[1].lower()
    cbio_study = f"{cancer_code}_tcga"
    cbio_url = "https://www.cbioportal.org/api"
    print(f"  Downloading {project_id} expression via cBioPortal ({cbio_study})...")
    try:
        # Get RNA-seq v2 profile (raw RSEM values)
        rna_profile = f"{cbio_study}_rna_seq_v2_mrna"
        # Verify it exists
        resp_check = requests.get(f"{cbio_url}/molecular-profiles/{rna_profile}", timeout=10)
        if resp_check.status_code != 200:
            # Fallback: try other profile names
            profiles = requests.get(f"{cbio_url}/studies/{cbio_study}/molecular-profiles", timeout=30).json()
            rna_profile = None
            for p in profiles:
                pid = p["molecularProfileId"]
                if "rna_seq" in pid and "Zscores" not in pid and p["molecularAlterationType"] == "MRNA_EXPRESSION":
                    rna_profile = pid
                    break
            if rna_profile is None:
                for p in profiles:
                    if p["molecularAlterationType"] == "MRNA_EXPRESSION" and "Zscores" not in p["molecularProfileId"]:
                        rna_profile = p["molecularProfileId"]
                        break
            if rna_profile is None:
                raise Exception(f"No RNA expression profile found")
        print(f"  RNA profile: {rna_profile}")

        # Get entrez gene IDs for signature genes
        gene_ids = []
        gene_id_map = {}
        for gene in sig_genes:
            try:
                resp_g = requests.get(f"{cbio_url}/genes/{gene}", timeout=10).json()
                eid = resp_g.get("entrezGeneId")
                if eid:
                    gene_ids.append(eid)
                    gene_id_map[eid] = gene
            except:
                pass
        print(f"  Gene IDs resolved: {len(gene_ids)}/{len(sig_genes)}")

        # Fetch expression for signature genes
        payload = {
            "entrezGeneIds": gene_ids,
            "sampleListId": f"{cbio_study}_all"
        }
        resp_expr = requests.post(
            f"{cbio_url}/molecular-profiles/{rna_profile}/molecular-data/fetch",
            json=payload, timeout=120,
            headers={"Content-Type": "application/json", "Accept": "application/json"}
        )
        if resp_expr.status_code != 200:
            raise Exception(f"Expression fetch returned {resp_expr.status_code}: {resp_expr.text[:200]}")

        expr_data = resp_expr.json()
        print(f"  Fetched {len(expr_data)} expression values")

        # Build expression matrix
        rows = {}
        for entry in expr_data:
            # entry has: entrezGeneId, value, sampleId (flat keys, not nested)
            eid = entry.get("entrezGeneId")
            gene_name = gene_id_map.get(eid, "")
            sample_id = entry.get("sampleId", "")[:12]
            value = entry.get("value")
            if gene_name and sample_id and value is not None and not np.isnan(value):
                if gene_name not in rows:
                    rows[gene_name] = {}
                rows[gene_name][sample_id] = value

        if rows:
            expr = pd.DataFrame(rows).T  # genes as rows, samples as columns
            print(f"  Expression matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")
        else:
            expr = pd.DataFrame()

    except Exception as e:
        print(f"  cBioPortal download failed: {e}")
        expr = pd.DataFrame()

    if len(expr) > 0 and len(clin) > 0:
        expr.to_csv(expr_cache)
        clin.to_csv(clin_cache)

    return expr, clin


# ══════════════════════════════════════════════════════════════════════════════
# RUN PAN-CANCER ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════
all_results = []
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
axes = axes.flatten()

for idx, (abbr, info) in enumerate(CANCERS.items()):
    print("\n" + "=" * 70)
    print(f"{idx+1}. {info['project']} — {info['name']}")
    print("=" * 70)

    try:
        expr, clin = download_tcga_cancer(info["project"], abbr)

        if len(expr) == 0 or len(clin) == 0:
            print(f"  No data for {abbr}")
            continue

        # Check signature gene availability
        avail = [g for g in sig_genes if g in expr.index]
        print(f"  Signature genes found: {len(avail)}/{len(sig_genes)}: {avail}")

        if len(avail) < 5:
            print(f"  Insufficient genes — skipping")
            continue

        # Match patients
        common = list(set(expr.columns) & set(clin.index))
        print(f"  Matched patients: {len(common)}")

        if len(common) < 30:
            print(f"  Too few patients — skipping")
            continue

        # Extract expression
        expr_sub = expr.loc[avail, common].T

        # Compute risk score
        risk = np.zeros(len(expr_sub))
        for gene in avail:
            if gene in selected_genes:
                z = (expr_sub[gene] - expr_sub[gene].mean()) / (expr_sub[gene].std() + 1e-10)
                risk += selected_genes[gene] * z.values

        clin_sub = clin.loc[common].copy()
        clin_sub["risk_score"] = risk
        clin_sub = clin_sub.dropna(subset=["OS_months", "OS_event"])
        clin_sub = clin_sub[clin_sub["OS_months"] > 0]

        median_risk = clin_sub["risk_score"].median()
        clin_sub["risk_group"] = (clin_sub["risk_score"] >= median_risk).map({True: "High", False: "Low"})

        n = len(clin_sub)
        events = int(clin_sub["OS_event"].sum())
        print(f"  Final: n={n}, events={events}")

        if events < 10:
            print(f"  Too few events — skipping")
            continue

        # C-index
        c_idx = concordance_index(clin_sub["OS_months"], -clin_sub["risk_score"], clin_sub["OS_event"])

        # Cox
        cox_df = clin_sub[["risk_score", "OS_months", "OS_event"]].copy()
        cph = CoxPHFitter()
        cph.fit(cox_df, duration_col="OS_months", event_col="OS_event")
        hr = np.exp(cph.params_["risk_score"])
        ci_lo = np.exp(cph.confidence_intervals_.iloc[0, 0])
        ci_hi = np.exp(cph.confidence_intervals_.iloc[0, 1])
        cox_p = cph.summary["p"]["risk_score"]

        # Log-rank
        high = clin_sub[clin_sub["risk_group"] == "High"]
        low = clin_sub[clin_sub["risk_group"] == "Low"]
        lr = logrank_test(high["OS_months"], low["OS_months"], high["OS_event"], low["OS_event"])
        lr_p = lr.p_value

        print(f"  C-index: {c_idx:.3f}")
        print(f"  HR: {hr:.2f} ({ci_lo:.2f}-{ci_hi:.2f}), p={cox_p:.4f}")
        print(f"  Log-rank p: {lr_p:.4f}")

        all_results.append({
            "cancer": abbr, "full_name": info["name"],
            "n": n, "events": events, "genes_available": len(avail),
            "c_index": c_idx, "HR": hr, "HR_lower": ci_lo, "HR_upper": ci_hi,
            "logrank_p": lr_p, "cox_p": cox_p
        })

        # KM plot
        ax = axes[idx]
        kmf = KaplanMeierFitter()
        for group, color in [("High", "#d62728"), ("Low", "#2ca02c")]:
            subset = clin_sub[clin_sub["risk_group"] == group]
            kmf.fit(subset["OS_months"], subset["OS_event"], label=f"{group} (n={len(subset)})")
            kmf.plot_survival_function(ax=ax, color=color, ci_show=True, ci_alpha=0.15)
        sig = "***" if lr_p < 0.001 else "**" if lr_p < 0.01 else "*" if lr_p < 0.05 else "ns"
        ax.set_title(f"{info['project']}\nC={c_idx:.3f}, HR={hr:.2f}, p={lr_p:.2e} {sig}",
                     fontsize=11, fontweight='bold')
        ax.set_xlabel("Time (months)")
        ax.set_ylabel("Survival probability")
        ax.legend(fontsize=9)

    except Exception as e:
        print(f"  {abbr} failed: {e}")
        import traceback
        traceback.print_exc()
        axes[idx].set_title(f"{info['project']}\n(failed)", fontsize=11)
        axes[idx].set_visible(False)

plt.suptitle("Pan-Cancer Validation of ROS/Ferroptosis Signature", fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIGS, "pan_cancer_km.png"), dpi=200, bbox_inches='tight')
print("\nSaved: pan_cancer_km.png")

# ══════════════════════════════════════════════════════════════════════════════
# FOREST PLOT
# ══════════════════════════════════════════════════════════════════════════════
if all_results:
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(os.path.join(TABLES, "pan_cancer_results.csv"), index=False)

    print("\n" + "=" * 70)
    print("PAN-CANCER SUMMARY")
    print("=" * 70)
    print(results_df[["cancer", "n", "events", "c_index", "HR", "logrank_p"]].to_string(index=False))

    # Forest plot
    fig, ax = plt.subplots(figsize=(10, 5))
    for i, (_, row) in enumerate(results_df.iterrows()):
        color = '#d62728' if row["HR"] > 1 else '#2166ac'
        ax.plot([row["HR_lower"], row["HR_upper"]], [i, i], color=color, linewidth=2.5)
        ax.plot(row["HR"], i, 'D', color=color, markersize=8)
    ax.axvline(1.0, color='grey', linestyle='--', linewidth=0.8)
    ax.set_yticks(range(len(results_df)))
    labels = []
    for _, row in results_df.iterrows():
        sig = "***" if row["logrank_p"] < 0.001 else "**" if row["logrank_p"] < 0.01 else "*" if row["logrank_p"] < 0.05 else ""
        labels.append(f"{row['cancer']} (n={row['n']}) {sig}")
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlabel("Hazard Ratio (95% CI)", fontsize=11)
    ax.set_title("Pan-Cancer Forest Plot: ROS/Ferroptosis Signature", fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "pan_cancer_forest.png"), dpi=200, bbox_inches='tight')
    print("Saved: pan_cancer_forest.png")

    n_sig = (results_df["logrank_p"] < 0.05).sum()
    print(f"\n  Significant in {n_sig}/{len(results_df)} cancer types")

print("\n✓ Pan-cancer analysis complete.")
