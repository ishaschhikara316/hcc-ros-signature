"""
15_mutation_landscape.py — Somatic mutation landscape by risk group

1. Query cBioPortal for top mutated genes in TCGA-LIHC
2. Compare mutation rates between high/low risk groups (Fisher exact)
3. Map risk groups to established molecular classes (TP53 vs CTNNB1)
4. Oncoplot-style visualization
"""
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import json
import os
import requests
import warnings
warnings.filterwarnings('ignore')

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data", "tcga")
MODEL = os.path.join(BASE, "results", "model")
TABLES = os.path.join(BASE, "results", "tables")
FIGS = os.path.join(BASE, "results", "figures")

# ── Load data ───────────────────────────────────────────────────────────────
merged = pd.read_csv(os.path.join(DATA, "tcga_ros_merged.csv"))
with open(os.path.join(MODEL, "lasso_model.json")) as f:
    model = json.load(f)

selected_genes = model["genes"]
df = merged.dropna(subset=["OS_months", "OS_event"]).copy()
df = df[df["OS_months"] > 0]

risk = np.zeros(len(df))
for gene, coef in selected_genes.items():
    if gene in df.columns:
        risk += coef * df[gene].values
df["risk_score"] = risk
df["risk_group"] = (df["risk_score"] >= df["risk_score"].median()).map({True: "High", False: "Low"})

# Standardize sample IDs to TCGA barcode format (first 12 chars)
if "patientId" in df.columns:
    df["tcga_id"] = df["patientId"].str[:12]
elif "patient_id" in df.columns:
    df["tcga_id"] = df["patient_id"].str[:12]
elif "sample_id" in df.columns:
    df["tcga_id"] = df["sample_id"].str[:12]
else:
    # Try index
    df["tcga_id"] = df.index.astype(str).str[:12]

print(f"Patients: {len(df)} ({(df['risk_group']=='High').sum()} high, {(df['risk_group']=='Low').sum()} low)")

# ══════════════════════════════════════════════════════════════════════════════
# 1. FETCH MUTATION DATA FROM cBioPortal
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("1. FETCHING MUTATION DATA FROM cBioPortal")
print("=" * 70)

# Top recurrently mutated genes in HCC
target_genes = ["TP53", "CTNNB1", "ARID1A", "AXIN1", "ALB", "ARID2",
                "BAP1", "RB1", "PIK3CA", "NFE2L2", "KEAP1", "TERT",
                "CDKN2A", "RPS6KA3", "APOB", "TSC1", "TSC2", "LZTR1",
                "BRD7", "ATM"]

# Also include our signature genes
target_genes += list(selected_genes.keys())
target_genes = list(set(target_genes))

cbio_url = "https://www.cbioportal.org/api"
study_id = "lihc_tcga"

# Get mutation profile
print("  Fetching mutation profile...")
try:
    profiles = requests.get(f"{cbio_url}/studies/{study_id}/molecular-profiles", timeout=30).json()
    mut_profile = None
    for p in profiles:
        if p["molecularAlterationType"] == "MUTATION_EXTENDED":
            mut_profile = p["molecularProfileId"]
            break
    print(f"  Mutation profile: {mut_profile}")
except Exception as e:
    print(f"  Failed to get profiles: {e}")
    mut_profile = "lihc_tcga_mutations"

# Fetch mutations — use per-gene query with proper entrezGeneId
print(f"  Querying mutations for {len(target_genes)} genes...")

# First get entrezGeneIds for our target genes
print("  Resolving gene IDs...")
gene_to_entrez = {}
try:
    for gene in target_genes:
        url = f"{cbio_url}/genes/{gene}"
        resp = requests.get(url, timeout=10, headers={"Accept": "application/json"})
        if resp.status_code == 200:
            data = resp.json()
            gene_to_entrez[gene] = data.get("entrezGeneId")
except Exception as e:
    print(f"  Gene ID resolution failed: {e}")

print(f"  Resolved {len(gene_to_entrez)} gene IDs")

# Fetch mutations per gene
mutation_matrix = {}
for gene, entrez_id in gene_to_entrez.items():
    if entrez_id is None:
        continue
    try:
        url = f"{cbio_url}/molecular-profiles/{mut_profile}/mutations"
        resp = requests.get(url, params={
            "sampleListId": f"{study_id}_all",
            "entrezGeneId": entrez_id,
            "projection": "SUMMARY"
        }, timeout=30, headers={"Accept": "application/json"})
        if resp.status_code == 200:
            muts = resp.json()
            mutated_samples = set(m.get("sampleId", "")[:12] for m in muts)
            if mutated_samples:
                mutation_matrix[gene] = mutated_samples
    except:
        pass

# If per-gene approach doesn't work, try POST with all entrez IDs at once
if len(mutation_matrix) < 5:
    print("  Per-gene fetch incomplete, trying batch POST...")
    try:
        entrez_ids = [eid for eid in gene_to_entrez.values() if eid is not None]
        payload = {
            "entrezGeneIds": entrez_ids,
            "sampleListId": f"{study_id}_all"
        }
        resp = requests.post(
            f"{cbio_url}/molecular-profiles/{mut_profile}/mutations/fetch",
            json=payload, timeout=120,
            params={"projection": "SUMMARY"},
            headers={"Content-Type": "application/json", "Accept": "application/json"}
        )
        if resp.status_code == 200:
            all_muts = resp.json()
            print(f"  Fetched {len(all_muts)} mutations via POST")
            for m in all_muts:
                gene = m.get("gene", {}).get("hugoGeneSymbol", "")
                if gene in target_genes:
                    sid = m.get("sampleId", "")[:12]
                    if gene not in mutation_matrix:
                        mutation_matrix[gene] = set()
                    mutation_matrix[gene].add(sid)
        else:
            print(f"  POST returned {resp.status_code}: {resp.text[:200]}")
    except Exception as e:
        print(f"  POST fetch failed: {e}")

print(f"  Genes with mutation data: {len(mutation_matrix)}")
for gene in ["TP53", "CTNNB1", "ARID1A", "AXIN1", "NFE2L2", "KEAP1"]:
    if gene in mutation_matrix:
        print(f"    {gene}: {len(mutation_matrix[gene])} mutated samples")

# ══════════════════════════════════════════════════════════════════════════════
# 2. MUTATION RATES BY RISK GROUP
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("2. MUTATION RATES BY RISK GROUP")
print("=" * 70)

if len(mutation_matrix) > 0:
    mut_results = []
    for gene, mutated in sorted(mutation_matrix.items(), key=lambda x: -len(x[1])):
        if len(mutated) < 3:
            continue
        # Match to our patients
        high_ids = set(df.loc[df["risk_group"] == "High", "tcga_id"].values)
        low_ids = set(df.loc[df["risk_group"] == "Low", "tcga_id"].values)

        high_mut = len(mutated & high_ids)
        high_wt = len(high_ids) - high_mut
        low_mut = len(mutated & low_ids)
        low_wt = len(low_ids) - low_mut

        # Fisher exact test
        odds, pval = stats.fisher_exact([[high_mut, high_wt], [low_mut, low_wt]])

        high_rate = high_mut / len(high_ids) * 100 if len(high_ids) > 0 else 0
        low_rate = low_mut / len(low_ids) * 100 if len(low_ids) > 0 else 0

        mut_results.append({
            "gene": gene,
            "total_mutated": len(mutated),
            "high_risk_mutated": high_mut, "high_risk_rate": high_rate,
            "low_risk_mutated": low_mut, "low_risk_rate": low_rate,
            "odds_ratio": odds, "fisher_p": pval
        })

    mut_df = pd.DataFrame(mut_results).sort_values("total_mutated", ascending=False)
    mut_df.to_csv(os.path.join(TABLES, "mutation_by_risk_group.csv"), index=False)

    # Print top results
    for _, row in mut_df.head(15).iterrows():
        sig = "*" if row["fisher_p"] < 0.05 else ""
        enriched = "HIGH" if row["high_risk_rate"] > row["low_risk_rate"] else "LOW"
        print(f"  {row['gene']:10s}: High={row['high_risk_rate']:5.1f}%  Low={row['low_risk_rate']:5.1f}%  "
              f"OR={row['odds_ratio']:.2f}  p={row['fisher_p']:.4f} {sig} → enriched in {enriched}")

    # ══════════════════════════════════════════════════════════════════════════
    # 3. MOLECULAR CLASS MAPPING
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("3. MOLECULAR CLASS MAPPING (TP53 vs CTNNB1)")
    print("=" * 70)

    if "TP53" in mutation_matrix and "CTNNB1" in mutation_matrix:
        tp53_mut = mutation_matrix["TP53"]
        ctnnb1_mut = mutation_matrix["CTNNB1"]

        # Classify patients
        df["TP53_mut"] = df["tcga_id"].isin(tp53_mut).astype(int)
        df["CTNNB1_mut"] = df["tcga_id"].isin(ctnnb1_mut).astype(int)

        def mol_class(row):
            if row["TP53_mut"] == 1 and row["CTNNB1_mut"] == 0:
                return "TP53-mutant (proliferation)"
            elif row["CTNNB1_mut"] == 1 and row["TP53_mut"] == 0:
                return "CTNNB1-mutant (Wnt)"
            elif row["TP53_mut"] == 1 and row["CTNNB1_mut"] == 1:
                return "Both"
            else:
                return "Neither"

        df["mol_class"] = df.apply(mol_class, axis=1)

        # Cross-tab with risk group
        ct = pd.crosstab(df["mol_class"], df["risk_group"])
        print(ct)
        print()

        # Risk score by molecular class
        for mc in df["mol_class"].unique():
            subset = df[df["mol_class"] == mc]
            print(f"  {mc:35s}: n={len(subset)}, mean risk={subset['risk_score'].mean():.3f}, "
                  f"high-risk={100*(subset['risk_group']=='High').mean():.1f}%")

        # Chi-squared test for association
        chi2, chi_p, dof, expected = stats.chi2_contingency(ct)
        print(f"\n  Chi-squared test (mol class vs risk group): chi2={chi2:.2f}, p={chi_p:.4f}")

        # ── Molecular class barplot ──
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Left: stacked bar of risk group by molecular class
        ct_norm = ct.div(ct.sum(axis=1), axis=0) * 100
        ct_norm.plot(kind='barh', stacked=True, ax=ax1, color={"High": "#d62728", "Low": "#2ca02c"})
        ax1.set_xlabel("Percentage (%)")
        ax1.set_title("Risk Group Distribution by Molecular Class", fontweight='bold')
        ax1.legend(title="Risk Group")

        # Right: risk score distribution by molecular class
        order = df.groupby("mol_class")["risk_score"].mean().sort_values().index.tolist()
        sns.boxplot(data=df, y="mol_class", x="risk_score", ax=ax2, order=order,
                    palette="Set2", orient='h')
        ax2.axvline(df["risk_score"].median(), color='grey', linestyle='--', linewidth=0.8, label='Median')
        ax2.set_xlabel("Risk Score")
        ax2.set_title("Risk Score by Molecular Class", fontweight='bold')
        ax2.set_ylabel("")

        plt.tight_layout()
        plt.savefig(os.path.join(FIGS, "molecular_class_risk.png"), dpi=200, bbox_inches='tight')
        print("Saved: molecular_class_risk.png")

    # ══════════════════════════════════════════════════════════════════════════
    # 4. ONCOPLOT-STYLE VISUALIZATION
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("4. ONCOPLOT")
    print("=" * 70)

    # Top 15 mutated genes
    top_genes = mut_df.head(15)["gene"].tolist()

    # Build binary matrix
    all_patient_ids = df.sort_values("risk_score", ascending=False)["tcga_id"].values
    onco_matrix = pd.DataFrame(0, index=top_genes, columns=all_patient_ids)
    for gene in top_genes:
        if gene in mutation_matrix:
            for pid in mutation_matrix[gene]:
                if pid in onco_matrix.columns:
                    onco_matrix.loc[gene, pid] = 1

    fig, (ax_risk, ax_onco) = plt.subplots(2, 1, figsize=(16, 8),
                                            gridspec_kw={"height_ratios": [1, 6]},
                                            sharex=True)

    # Top bar: risk group
    risk_vals = df.set_index("tcga_id").loc[all_patient_ids, "risk_score"].values
    colors_risk = ['#d62728' if r >= df["risk_score"].median() else '#2ca02c' for r in risk_vals]
    ax_risk.bar(range(len(all_patient_ids)), risk_vals, color=colors_risk, width=1.0, linewidth=0)
    ax_risk.set_ylabel("Risk\nScore", fontsize=9)
    ax_risk.set_xlim(-0.5, len(all_patient_ids) - 0.5)
    ax_risk.axhline(df["risk_score"].median(), color='black', linewidth=0.5, linestyle='--')

    # Oncoplot
    sns.heatmap(onco_matrix.values.astype(float), ax=ax_onco,
                cmap=["#f0f0f0", "#333333"], cbar=False,
                xticklabels=False, yticklabels=top_genes,
                linewidths=0.1, linecolor='white')
    ax_onco.set_xlabel(f"Patients (n={len(all_patient_ids)}, sorted by risk score)")

    # Add mutation rates on right
    for i, gene in enumerate(top_genes):
        rate = len(mutation_matrix.get(gene, set())) / len(all_patient_ids) * 100
        ax_onco.text(len(all_patient_ids) + 1, i + 0.5, f"{rate:.0f}%", va='center', fontsize=8)

    plt.suptitle("Somatic Mutation Landscape by Risk Score", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "oncoplot_risk.png"), dpi=200, bbox_inches='tight')
    print("Saved: oncoplot_risk.png")

    # ── NRF2/KEAP1 mutation analysis ──
    print("\n  NRF2 pathway mutations:")
    for gene in ["NFE2L2", "KEAP1"]:
        if gene in mutation_matrix:
            mutated = mutation_matrix[gene]
            high_ids = set(df.loc[df["risk_group"] == "High", "tcga_id"].values)
            low_ids = set(df.loc[df["risk_group"] == "Low", "tcga_id"].values)
            print(f"    {gene}: {len(mutated & high_ids)} in high-risk, {len(mutated & low_ids)} in low-risk")

else:
    print("  No mutation data available — skipping visualization")

# ══════════════════════════════════════════════════════════════════════════════
# 5. SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("5. SUMMARY")
print("=" * 70)

if len(mutation_matrix) > 0:
    sig_genes_mut = mut_df[mut_df["fisher_p"] < 0.05]
    print(f"  Genes with differential mutation rates (p<0.05): {len(sig_genes_mut)}")
    for _, row in sig_genes_mut.iterrows():
        enriched = "high-risk" if row["high_risk_rate"] > row["low_risk_rate"] else "low-risk"
        print(f"    {row['gene']}: enriched in {enriched} (OR={row['odds_ratio']:.2f}, p={row['fisher_p']:.4f})")

print("\n✓ Mutation landscape analysis complete.")
