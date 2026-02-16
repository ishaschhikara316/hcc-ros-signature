"""
09_protein_validation.py — Protein-level validation using Human Protein Atlas

1. Download HPA normal tissue + pathology (cancer) protein expression
2. Compare signature gene protein levels: normal liver vs HCC
3. Show IHC staining evidence summary
4. Correlate mRNA expression with protein levels
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import requests
import json
import os
import warnings
warnings.filterwarnings('ignore')

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
HPA = os.path.join(DATA, "hpa")
MODEL = os.path.join(BASE, "results", "model")
TABLES = os.path.join(BASE, "results", "tables")
FIGS = os.path.join(BASE, "results", "figures")
os.makedirs(HPA, exist_ok=True)

# ── Load model ──────────────────────────────────────────────────────────────
with open(os.path.join(MODEL, "lasso_model.json")) as f:
    model = json.load(f)
sig_genes = list(model["genes"].keys())
sig_coefs = model["genes"]
print(f"Signature genes ({len(sig_genes)}): {sig_genes}")

# ══════════════════════════════════════════════════════════════════════════════
# 1. DOWNLOAD HPA DATA
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("1. DOWNLOADING HUMAN PROTEIN ATLAS DATA")
print("=" * 70)

# HPA v25 download — single comprehensive TSV (contains tissue + pathology data)
hpa_url = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
hpa_path = os.path.join(HPA, "proteinatlas.tsv.zip")
if not os.path.exists(hpa_path):
    print("  Downloading HPA comprehensive data (proteinatlas.tsv.zip)...")
    try:
        resp = requests.get(hpa_url, timeout=300)
        resp.raise_for_status()
        with open(hpa_path, 'wb') as f:
            f.write(resp.content)
        print(f"  Saved ({len(resp.content) / 1e6:.1f} MB)")
    except Exception as e:
        print(f"  Failed to download HPA: {e}")

# Also try direct normal_tissue and pathology endpoints (v24 format)
normal_path = os.path.join(HPA, "normal_tissue.tsv.zip")
pathology_path = os.path.join(HPA, "pathology.tsv.zip")

for url, path, desc in [
    ("https://v23.proteinatlas.org/download/normal_tissue.tsv.zip", normal_path, "normal tissue"),
    ("https://v23.proteinatlas.org/download/pathology.tsv.zip", pathology_path, "pathology"),
]:
    if not os.path.exists(path):
        print(f"  Downloading HPA {desc}...")
        try:
            resp = requests.get(url, timeout=120)
            resp.raise_for_status()
            with open(path, 'wb') as f:
                f.write(resp.content)
            print(f"  Saved: {os.path.basename(path)} ({len(resp.content) / 1e6:.1f} MB)")
        except Exception as e:
            print(f"  Failed to download {desc}: {e}")

# ══════════════════════════════════════════════════════════════════════════════
# 2. PARSE NORMAL TISSUE EXPRESSION
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("2. NORMAL TISSUE PROTEIN EXPRESSION")
print("=" * 70)

normal_results = []
staining_map = {"not detected": 0, "low": 1, "medium": 2, "high": 3}

try:
    if os.path.exists(normal_path):
        normal_df = pd.read_csv(normal_path, sep='\t', compression='zip')
    elif os.path.exists(hpa_path):
        # Parse from comprehensive file
        print("  Using comprehensive proteinatlas.tsv (extracting tissue columns)...")
        hpa_full = pd.read_csv(hpa_path, sep='\t', compression='zip')
        print(f"  Full HPA data: {len(hpa_full)} rows, {len(hpa_full.columns)} columns")
        # The comprehensive file has columns like "RNA tissue specificity", "Tissue expression cluster"
        # It also has IHC-related columns
        normal_df = hpa_full  # We'll filter by gene name below
    else:
        raise FileNotFoundError("No HPA data files found")
    print(f"  Normal tissue data: {len(normal_df)} rows")
    print(f"  Columns: {list(normal_df.columns)[:15]}...")

    # Check if this is the per-tissue file or the comprehensive file
    gene_col = "Gene name" if "Gene name" in normal_df.columns else (
        "Gene" if "Gene" in normal_df.columns else normal_df.columns[1])

    if "Tissue" in normal_df.columns and "Level" in normal_df.columns:
        # Standard normal_tissue.tsv format
        tissue_col = "Tissue"
        level_col = "Level"

        liver_normal = normal_df[
            (normal_df[gene_col].isin(sig_genes)) &
            (normal_df[tissue_col].str.lower().str.contains("liver", na=False))
        ]
        print(f"  Liver entries for signature genes: {len(liver_normal)}")

        for gene in sig_genes:
            gene_data = liver_normal[liver_normal[gene_col] == gene]
            if len(gene_data) > 0:
                levels = gene_data[level_col].str.lower().values
                scores = [staining_map.get(l.strip(), -1) for l in levels if isinstance(l, str)]
                max_score = max(scores) if scores else -1
                level_str = {0: "Not detected", 1: "Low", 2: "Medium", 3: "High"}.get(max_score, "N/A")
                cell_types = gene_data["Cell type"].unique() if "Cell type" in gene_data.columns else ["N/A"]
                normal_results.append({
                    "gene": gene, "tissue": "Liver (normal)",
                    "protein_level": level_str, "score": max_score,
                    "cell_types": "; ".join(str(c) for c in cell_types),
                    "n_entries": len(gene_data),
                })
                print(f"  {gene:<10} Normal liver: {level_str} (in {', '.join(str(c) for c in cell_types[:3])})")
            else:
                normal_results.append({"gene": gene, "tissue": "Liver (normal)",
                                        "protein_level": "No data", "score": -1,
                                        "cell_types": "N/A", "n_entries": 0})
                print(f"  {gene:<10} Normal liver: No data")
    else:
        # Comprehensive proteinatlas.tsv format
        print("  Using comprehensive format — extracting RNA tissue data")
        # Look for RNA expression columns related to liver
        liver_cols = [c for c in normal_df.columns if 'liver' in c.lower()]
        rna_liver_col = [c for c in liver_cols if 'rna' in c.lower() or 'nTPM' in c.lower()
                         or 'tissue' in c.lower()]
        print(f"  Liver-related columns: {liver_cols[:5]}")

        sig_data = normal_df[normal_df[gene_col].isin(sig_genes)]
        for gene in sig_genes:
            gene_data = sig_data[sig_data[gene_col] == gene]
            if len(gene_data) > 0:
                # Use RNA tissue expression if available
                if rna_liver_col:
                    val = gene_data[rna_liver_col[0]].values[0]
                    if pd.notna(val):
                        try:
                            val = float(val)
                            if val > 100: score = 3
                            elif val > 10: score = 2
                            elif val > 1: score = 1
                            else: score = 0
                        except:
                            score = -1
                    else:
                        score = -1
                else:
                    score = -1
                level_str = {0: "Not detected", 1: "Low", 2: "Medium", 3: "High"}.get(score, "N/A")
                normal_results.append({
                    "gene": gene, "tissue": "Liver (normal)",
                    "protein_level": level_str, "score": score,
                    "cell_types": "N/A (from RNA)", "n_entries": 1,
                })
                print(f"  {gene:<10} Normal liver: {level_str}")
            else:
                normal_results.append({"gene": gene, "tissue": "Liver (normal)",
                                        "protein_level": "No data", "score": -1,
                                        "cell_types": "N/A", "n_entries": 0})
                print(f"  {gene:<10} Normal liver: No data")

except Exception as e:
    print(f"  Error parsing normal tissue data: {e}")
    import traceback; traceback.print_exc()

# ══════════════════════════════════════════════════════════════════════════════
# 3. PARSE PATHOLOGY (CANCER) EXPRESSION
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("3. HCC PATHOLOGY PROTEIN EXPRESSION")
print("=" * 70)

pathology_results = []

try:
    if os.path.exists(pathology_path):
        path_df = pd.read_csv(pathology_path, sep='\t', compression='zip')
    elif os.path.exists(hpa_path):
        print("  Using comprehensive proteinatlas.tsv for pathology...")
        if 'hpa_full' not in dir():
            hpa_full = pd.read_csv(hpa_path, sep='\t', compression='zip')
        path_df = hpa_full
    else:
        raise FileNotFoundError("No pathology data")
    print(f"  Pathology data: {len(path_df)} rows")
    print(f"  Columns: {list(path_df.columns)[:15]}...")

    gene_col_p = "Gene name" if "Gene name" in path_df.columns else (
        "Gene" if "Gene" in path_df.columns else path_df.columns[1])

    if "Cancer" in path_df.columns:
        # Standard pathology.tsv format
        cancer_col = "Cancer"
        liver_cancer = path_df[
            (path_df[gene_col_p].isin(sig_genes)) &
            (path_df[cancer_col].str.lower().str.contains("liver", na=False))
        ]
        print(f"  Liver cancer entries for signature genes: {len(liver_cancer)}")

        for gene in sig_genes:
            gene_data = liver_cancer[liver_cancer[gene_col_p] == gene]
            if len(gene_data) > 0:
                high_col = [c for c in gene_data.columns if c.lower() == "high"]
                med_col = [c for c in gene_data.columns if c.lower() == "medium"]
                low_col = [c for c in gene_data.columns if c.lower() == "low"]
                nd_col = [c for c in gene_data.columns if c.lower() == "not detected"]

                n_high = gene_data[high_col[0]].sum() if high_col else 0
                n_med = gene_data[med_col[0]].sum() if med_col else 0
                n_low = gene_data[low_col[0]].sum() if low_col else 0
                n_nd = gene_data[nd_col[0]].sum() if nd_col else 0
                total = n_high + n_med + n_low + n_nd

                counts = {"High": n_high, "Medium": n_med, "Low": n_low, "Not detected": n_nd}
                dominant = max(counts, key=counts.get) if total > 0 else "N/A"

                pathology_results.append({
                    "gene": gene, "cancer": "Liver cancer (HCC)",
                    "n_high": int(n_high), "n_medium": int(n_med),
                    "n_low": int(n_low), "n_not_detected": int(n_nd),
                    "total_patients": int(total),
                    "dominant_staining": dominant,
                    "pct_positive": f"{(n_high + n_med + n_low) / total * 100:.0f}%" if total > 0 else "N/A",
                })
                coef = sig_coefs[gene]
                direction = "risk" if coef > 0 else "protective"
                print(f"  {gene:<10} HCC: H={int(n_high)} M={int(n_med)} L={int(n_low)} ND={int(n_nd)} "
                      f"({direction}, coef={coef:.3f})")
            else:
                pathology_results.append({
                    "gene": gene, "cancer": "Liver cancer (HCC)",
                    "n_high": 0, "n_medium": 0, "n_low": 0, "n_not_detected": 0,
                    "total_patients": 0, "dominant_staining": "No data",
                    "pct_positive": "N/A",
                })
                print(f"  {gene:<10} HCC: No data")
    else:
        # Comprehensive proteinatlas.tsv — extract pathology columns
        print("  Comprehensive format — looking for pathology columns...")
        pathology_cols = [c for c in path_df.columns if 'pathology' in c.lower()
                          or 'cancer' in c.lower()]
        liver_pathology_col = [c for c in path_df.columns
                                if 'liver' in c.lower() and ('cancer' in c.lower()
                                or 'pathology' in c.lower())]
        print(f"  Pathology columns: {pathology_cols[:8]}")
        print(f"  Liver pathology columns: {liver_pathology_col[:5]}")

        sig_data = path_df[path_df[gene_col_p].isin(sig_genes)]
        for gene in sig_genes:
            gene_data = sig_data[sig_data[gene_col_p] == gene]
            if len(gene_data) > 0 and liver_pathology_col:
                val = str(gene_data[liver_pathology_col[0]].values[0]).lower()
                # Map text to staining levels
                if 'high' in val: dominant = "High"
                elif 'medium' in val: dominant = "Medium"
                elif 'low' in val: dominant = "Low"
                elif 'not detected' in val: dominant = "Not detected"
                else: dominant = val.title() if val != 'nan' else "N/A"

                pathology_results.append({
                    "gene": gene, "cancer": "Liver cancer (HCC)",
                    "n_high": 0, "n_medium": 0, "n_low": 0, "n_not_detected": 0,
                    "total_patients": 0,
                    "dominant_staining": dominant,
                    "pct_positive": "N/A (from summary)",
                })
                coef = sig_coefs[gene]
                direction = "risk" if coef > 0 else "protective"
                print(f"  {gene:<10} HCC: {dominant} ({direction}, coef={coef:.3f})")
            else:
                pathology_results.append({
                    "gene": gene, "cancer": "Liver cancer (HCC)",
                    "n_high": 0, "n_medium": 0, "n_low": 0, "n_not_detected": 0,
                    "total_patients": 0, "dominant_staining": "No data",
                    "pct_positive": "N/A",
                })
                print(f"  {gene:<10} HCC: No data")

except Exception as e:
    print(f"  Error parsing pathology data: {e}")
    import traceback; traceback.print_exc()

# ══════════════════════════════════════════════════════════════════════════════
# 4. PROTEIN EXPRESSION VISUALIZATION
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("4. PROTEIN EXPRESSION FIGURES")
print("=" * 70)

if pathology_results:
    # Stacked bar plot: IHC staining distribution in HCC
    fig, ax = plt.subplots(figsize=(12, 6))
    genes_sorted = sorted(pathology_results, key=lambda x: sig_coefs.get(x["gene"], 0), reverse=True)
    gene_names = [r["gene"] for r in genes_sorted]
    n_h = [r["n_high"] for r in genes_sorted]
    n_m = [r["n_medium"] for r in genes_sorted]
    n_l = [r["n_low"] for r in genes_sorted]
    n_nd = [r["n_not_detected"] for r in genes_sorted]

    x = np.arange(len(gene_names))
    width = 0.6

    ax.bar(x, n_h, width, label='High', color='#d32f2f')
    ax.bar(x, n_m, width, bottom=n_h, label='Medium', color='#ff9800')
    ax.bar(x, n_l, width, bottom=[h + m for h, m in zip(n_h, n_m)], label='Low', color='#fdd835')
    ax.bar(x, n_nd, width, bottom=[h + m + l for h, m, l in zip(n_h, n_m, n_l)],
           label='Not detected', color='#e0e0e0')

    ax.set_xticks(x)
    ax.set_xticklabels(gene_names, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel("Number of HCC patients")
    ax.set_title("IHC Protein Staining in HCC — Signature Genes\n(Human Protein Atlas)",
                 fontsize=13, fontweight='bold')
    ax.legend(loc='upper right')

    # Add risk/protective annotation
    for i, gene in enumerate(gene_names):
        coef = sig_coefs[gene]
        color = 'red' if coef > 0 else 'blue'
        ax.annotate("+" if coef > 0 else "-", xy=(i, 0), xytext=(i, -2),
                     ha='center', fontsize=12, fontweight='bold', color=color)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "protein_ihc_staining.png"), dpi=200, bbox_inches='tight')
    print("Saved: protein_ihc_staining.png")

# Normal vs HCC comparison heatmap
if normal_results and pathology_results:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: Normal liver protein levels
    ax = axes[0]
    normal_scores = []
    for gene in sig_genes:
        nr = [r for r in normal_results if r["gene"] == gene]
        normal_scores.append(nr[0]["score"] if nr else -1)

    colors_normal = []
    for s in normal_scores:
        if s == 3: colors_normal.append('#d32f2f')
        elif s == 2: colors_normal.append('#ff9800')
        elif s == 1: colors_normal.append('#fdd835')
        elif s == 0: colors_normal.append('#e0e0e0')
        else: colors_normal.append('#bdbdbd')

    y = np.arange(len(sig_genes))
    ax.barh(y, [max(0, s) for s in normal_scores], color=colors_normal, edgecolor='black', linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(sig_genes, fontsize=10)
    ax.set_xlim(0, 3.5)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(["ND", "Low", "Med", "High"])
    ax.set_title("Normal Liver", fontsize=12, fontweight='bold')
    ax.set_xlabel("Protein Expression Level")

    # Right: HCC staining (% positive)
    ax = axes[1]
    hcc_pct = []
    for gene in sig_genes:
        pr = [r for r in pathology_results if r["gene"] == gene]
        if pr and pr[0]["total_patients"] > 0:
            total = pr[0]["total_patients"]
            pos = pr[0]["n_high"] + pr[0]["n_medium"] + pr[0]["n_low"]
            hcc_pct.append(pos / total * 100)
        else:
            hcc_pct.append(0)

    bar_colors = ['#d32f2f' if sig_coefs[g] > 0 else '#1565c0' for g in sig_genes]
    ax.barh(y, hcc_pct, color=bar_colors, edgecolor='black', linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(sig_genes, fontsize=10)
    ax.set_xlim(0, 105)
    ax.set_title("HCC Tumors", fontsize=12, fontweight='bold')
    ax.set_xlabel("% Patients with Detectable Protein")

    # Legend
    risk_patch = mpatches.Patch(color='#d32f2f', label='Risk gene (coef > 0)')
    prot_patch = mpatches.Patch(color='#1565c0', label='Protective gene (coef < 0)')
    ax.legend(handles=[risk_patch, prot_patch], loc='lower right', fontsize=9)

    plt.suptitle("Protein-Level Validation — Human Protein Atlas", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "protein_normal_vs_hcc.png"), dpi=200, bbox_inches='tight')
    print("Saved: protein_normal_vs_hcc.png")

# ══════════════════════════════════════════════════════════════════════════════
# 5. mRNA vs PROTEIN CONCORDANCE
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("5. mRNA-PROTEIN CONCORDANCE")
print("=" * 70)

# Load TCGA mRNA data
tcga_merged = pd.read_csv(os.path.join(DATA, "tcga", "tcga_ros_merged.csv"))

concordance_results = []
for gene in sig_genes:
    coef = sig_coefs[gene]
    direction_mrna = "upregulated" if coef > 0 else "downregulated"

    # mRNA in high vs low risk
    if "risk_score" in tcga_merged.columns and gene in tcga_merged.columns:
        med = tcga_merged["risk_score"].median()
        high_expr = tcga_merged.loc[tcga_merged["risk_score"] >= med, gene].mean()
        low_expr = tcga_merged.loc[tcga_merged["risk_score"] < med, gene].mean()
        fc_mrna = high_expr / (low_expr + 1e-10)
    else:
        fc_mrna = np.nan

    # Protein level in HCC
    pr = [r for r in pathology_results if r["gene"] == gene]
    if pr and pr[0]["total_patients"] > 0:
        prot_pct_pos = (pr[0]["n_high"] + pr[0]["n_medium"]) / pr[0]["total_patients"] * 100
        prot_dominant = pr[0]["dominant_staining"]
    else:
        prot_pct_pos = np.nan
        prot_dominant = "N/A"

    # Concordance: risk gene should have high protein in HCC
    if coef > 0:
        concordant = prot_pct_pos > 50 if not np.isnan(prot_pct_pos) else None
    else:
        concordant = prot_pct_pos < 50 if not np.isnan(prot_pct_pos) else None

    concordance_results.append({
        "gene": gene, "coef": coef,
        "mRNA_direction": direction_mrna,
        "mRNA_FC_high_vs_low": fc_mrna,
        "protein_pct_positive": prot_pct_pos,
        "protein_dominant": prot_dominant,
        "concordant": concordant,
    })
    status = "YES" if concordant else ("NO" if concordant is False else "N/A")
    print(f"  {gene:<10} mRNA: {direction_mrna:<14} FC={fc_mrna:.2f}  "
          f"Protein: {prot_dominant:<12} ({prot_pct_pos:.0f}% pos)  "
          f"Concordant: {status}")

conc_df = pd.DataFrame(concordance_results)
conc_df.to_csv(os.path.join(TABLES, "protein_concordance.csv"), index=False)
print("\nSaved: protein_concordance.csv")

n_concordant = sum(1 for r in concordance_results if r["concordant"] is True)
n_tested = sum(1 for r in concordance_results if r["concordant"] is not None)
print(f"\nConcordance: {n_concordant}/{n_tested} genes show mRNA-protein agreement")

# ══════════════════════════════════════════════════════════════════════════════
# 6. SAVE COMBINED TABLE
# ══════════════════════════════════════════════════════════════════════════════
all_protein = []
for gene in sig_genes:
    row = {"gene": gene, "coef": sig_coefs[gene]}
    nr = [r for r in normal_results if r["gene"] == gene]
    if nr:
        row["normal_liver_protein"] = nr[0]["protein_level"]
        row["normal_liver_cell_types"] = nr[0]["cell_types"]
    pr = [r for r in pathology_results if r["gene"] == gene]
    if pr:
        row["hcc_dominant_staining"] = pr[0]["dominant_staining"]
        row["hcc_pct_positive"] = pr[0]["pct_positive"]
        row["hcc_n_patients"] = pr[0]["total_patients"]
    cr = [r for r in concordance_results if r["gene"] == gene]
    if cr:
        row["mrna_protein_concordant"] = cr[0]["concordant"]
    all_protein.append(row)

pd.DataFrame(all_protein).to_csv(os.path.join(TABLES, "protein_validation_summary.csv"), index=False)
print("Saved: protein_validation_summary.csv")

print(f"\n✓ Protein-level validation complete.")
