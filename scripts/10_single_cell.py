"""
10_single_cell.py — Single-cell RNA-seq analysis of signature genes

1. Download HPA single-cell type data (pre-computed from HPA)
2. Show cell-type-specific expression of signature genes in liver
3. Create dotplot / heatmap of gene × cell-type expression
4. Determine if signature genes are tumor-intrinsic vs stromal
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
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
# 1. DOWNLOAD HPA SINGLE-CELL TYPE DATA
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("1. DOWNLOADING HPA SINGLE-CELL DATA")
print("=" * 70)

# HPA provides RNA single cell type data
sc_url = "https://v23.proteinatlas.org/download/rna_single_cell_type.tsv.zip"
sc_path = os.path.join(HPA, "rna_single_cell_type.tsv.zip")

if not os.path.exists(sc_path):
    print("  Downloading HPA single-cell type data...")
    try:
        resp = requests.get(sc_url, timeout=180)
        resp.raise_for_status()
        with open(sc_path, 'wb') as f:
            f.write(resp.content)
        print(f"  Saved ({len(resp.content) / 1e6:.1f} MB)")
    except Exception as e:
        print(f"  Failed: {e}")

# Also get RNA single cell type tissue data (tissue-specific)
sc_tissue_url = "https://v23.proteinatlas.org/download/rna_single_cell_type_tissue.tsv.zip"
sc_tissue_path = os.path.join(HPA, "rna_single_cell_type_tissue.tsv.zip")

if not os.path.exists(sc_tissue_path):
    print("  Downloading HPA tissue-specific single-cell data...")
    try:
        resp = requests.get(sc_tissue_url, timeout=180)
        resp.raise_for_status()
        with open(sc_tissue_path, 'wb') as f:
            f.write(resp.content)
        print(f"  Saved ({len(resp.content) / 1e6:.1f} MB)")
    except Exception as e:
        print(f"  Failed: {e}")
        # Fall back to general single cell data
        sc_tissue_path = None

# ══════════════════════════════════════════════════════════════════════════════
# 2. PARSE LIVER CELL-TYPE EXPRESSION
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("2. LIVER CELL-TYPE EXPRESSION OF SIGNATURE GENES")
print("=" * 70)

sc_results = []

# Try tissue-specific first (has liver cell types)
data_loaded = False
if sc_tissue_path and os.path.exists(sc_tissue_path):
    try:
        sc_df = pd.read_csv(sc_tissue_path, sep='\t', compression='zip')
        print(f"  Tissue-specific SC data: {len(sc_df)} rows")
        print(f"  Columns: {list(sc_df.columns)}")

        gene_col = "Gene name" if "Gene name" in sc_df.columns else sc_df.columns[1]
        tissue_col = [c for c in sc_df.columns if 'tissue' in c.lower()]
        cell_col = [c for c in sc_df.columns if 'cell' in c.lower() and 'type' in c.lower()]
        expr_cols = [c for c in sc_df.columns if c.lower() in ('ntpm', 'ptpm', 'nTPM')
                     or 'tpm' in c.lower() or 'read' in c.lower()]
        if not expr_cols:
            # Try numeric columns
            expr_cols = [c for c in sc_df.columns if sc_df[c].dtype in ['float64', 'int64']
                         and c not in [gene_col]]

        if tissue_col:
            tissue_col = tissue_col[0]
        else:
            tissue_col = None

        if cell_col:
            cell_col = cell_col[0]
        else:
            cell_col = None

        if expr_cols:
            expr_col = expr_cols[0]
        else:
            expr_col = None

        print(f"  Gene col: {gene_col}, Tissue col: {tissue_col}, Cell col: {cell_col}, Expr col: {expr_col}")

        if tissue_col and cell_col and expr_col:
            # Filter to liver tissue
            liver_sc = sc_df[
                (sc_df[gene_col].isin(sig_genes)) &
                (sc_df[tissue_col].str.lower().str.contains("liver", na=False))
            ]
            print(f"  Liver entries for signature genes: {len(liver_sc)}")

            if len(liver_sc) > 0:
                # Pivot: gene × cell_type
                sc_pivot = liver_sc.pivot_table(
                    index=gene_col, columns=cell_col, values=expr_col, aggfunc='mean'
                )
                # Filter to signature genes
                sc_pivot = sc_pivot.reindex(sig_genes).dropna(how='all')
                if len(sc_pivot) > 0:
                    data_loaded = True
                    sc_results = sc_pivot
                    print(f"  Cell types found: {list(sc_pivot.columns)}")
                    print(f"  Genes with data: {list(sc_pivot.index)}")
    except Exception as e:
        print(f"  Error parsing tissue-specific SC data: {e}")
        import traceback; traceback.print_exc()

# Fallback: use general single-cell type data
if not data_loaded and os.path.exists(sc_path):
    try:
        sc_df = pd.read_csv(sc_path, sep='\t', compression='zip')
        print(f"  General SC data: {len(sc_df)} rows")
        print(f"  Columns: {list(sc_df.columns)}")

        gene_col = "Gene name" if "Gene name" in sc_df.columns else sc_df.columns[1]
        cell_col = [c for c in sc_df.columns if 'cell' in c.lower() and 'type' in c.lower()]
        expr_cols = [c for c in sc_df.columns if 'tpm' in c.lower() or 'nTPM' in c.lower()]
        if not expr_cols:
            expr_cols = [c for c in sc_df.columns if sc_df[c].dtype in ['float64', 'int64']
                         and c not in [gene_col]]

        if cell_col:
            cell_col = cell_col[0]
        else:
            cell_col = sc_df.columns[2] if len(sc_df.columns) > 2 else None

        if expr_cols:
            expr_col = expr_cols[0]
        else:
            expr_col = sc_df.columns[-1]

        print(f"  Gene col: {gene_col}, Cell col: {cell_col}, Expr col: {expr_col}")

        if cell_col and expr_col:
            sig_sc = sc_df[sc_df[gene_col].isin(sig_genes)]
            print(f"  Entries for signature genes: {len(sig_sc)}")

            if len(sig_sc) > 0:
                sc_pivot = sig_sc.pivot_table(
                    index=gene_col, columns=cell_col, values=expr_col, aggfunc='mean'
                )
                sc_pivot = sc_pivot.reindex(sig_genes).dropna(how='all')

                # Keep most relevant cell types (liver-related + immune)
                liver_types = ['hepatocytes', 'cholangiocytes', 'kupffer cells',
                               'endothelial cells', 'stellate cells',
                               'hepatocyte', 'cholangiocyte', 'kupffer']
                immune_types = ['t-cells', 'b-cells', 'nk-cells', 'macrophages',
                                'monocytes', 'dendritic', 'neutrophils',
                                'T-cells', 'B-cells', 'NK-cells']
                relevant = []
                for ct in sc_pivot.columns:
                    ct_lower = ct.lower()
                    if any(lt in ct_lower for lt in liver_types + immune_types):
                        relevant.append(ct)

                if relevant:
                    sc_pivot = sc_pivot[relevant]
                elif len(sc_pivot.columns) > 20:
                    # Keep top-expressing cell types
                    top_cols = sc_pivot.mean().nlargest(15).index
                    sc_pivot = sc_pivot[top_cols]

                if len(sc_pivot) > 0:
                    data_loaded = True
                    sc_results = sc_pivot
                    print(f"  Cell types: {list(sc_pivot.columns)[:10]}...")
    except Exception as e:
        print(f"  Error parsing general SC data: {e}")
        import traceback; traceback.print_exc()

# ══════════════════════════════════════════════════════════════════════════════
# 3. VISUALIZATION
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("3. CELL-TYPE EXPRESSION VISUALIZATIONS")
print("=" * 70)

if isinstance(sc_results, pd.DataFrame) and len(sc_results) > 0:
    # Clean up cell type names
    sc_results.columns = [c.replace('_', ' ').title() if isinstance(c, str) else c
                          for c in sc_results.columns]

    # ── Dotplot (gene × cell type) ──
    fig, ax = plt.subplots(figsize=(max(10, len(sc_results.columns) * 0.8),
                                      max(5, len(sc_results) * 0.5 + 2)))

    # Normalize per gene (row-wise z-score)
    sc_norm = sc_results.copy()
    for gene in sc_norm.index:
        row = sc_norm.loc[gene]
        if row.std() > 0:
            sc_norm.loc[gene] = (row - row.mean()) / row.std()
        else:
            sc_norm.loc[gene] = 0

    # Create bubble plot
    cell_types = list(sc_results.columns)
    genes = list(sc_results.index)

    for i, gene in enumerate(genes):
        for j, ct in enumerate(cell_types):
            val = sc_results.loc[gene, ct] if pd.notna(sc_results.loc[gene, ct]) else 0
            z_val = sc_norm.loc[gene, ct] if pd.notna(sc_norm.loc[gene, ct]) else 0
            # Size proportional to expression
            size = max(10, min(300, val * 3))
            # Color by z-score
            color = plt.cm.RdBu_r((z_val + 3) / 6)  # Map z=-3..3 to 0..1
            ax.scatter(j, i, s=size, c=[color], edgecolors='black', linewidth=0.5, alpha=0.8)

    ax.set_xticks(range(len(cell_types)))
    ax.set_xticklabels(cell_types, rotation=60, ha='right', fontsize=9)
    ax.set_yticks(range(len(genes)))
    # Add risk/protective annotation to gene labels
    gene_labels = []
    for g in genes:
        coef = sig_coefs.get(g, 0)
        marker = " (+)" if coef > 0 else " (-)"
        gene_labels.append(g + marker)
    ax.set_yticklabels(gene_labels, fontsize=10)
    ax.set_title("Single-Cell Expression of Signature Genes\n(Human Protein Atlas)",
                 fontsize=13, fontweight='bold')

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=mcolors.Normalize(vmin=-3, vmax=3))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("Row Z-score", fontsize=10)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "single_cell_dotplot.png"), dpi=200, bbox_inches='tight')
    print("Saved: single_cell_dotplot.png")

    # ── Heatmap ──
    import seaborn as sns

    fig, ax = plt.subplots(figsize=(max(10, len(sc_results.columns) * 0.6),
                                      max(4, len(sc_results) * 0.45 + 2)))
    # Log-transform for better visualization
    sc_log = np.log2(sc_results + 1)

    sns.heatmap(sc_log, ax=ax, cmap="YlOrRd", xticklabels=True, yticklabels=gene_labels,
                cbar_kws={"label": "log2(nTPM + 1)", "shrink": 0.5},
                linewidths=0.5, linecolor='white')
    ax.set_title("Single-Cell Gene Expression — Liver Cell Types\n(Human Protein Atlas)",
                 fontsize=13, fontweight='bold')
    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Signature Gene")
    plt.xticks(rotation=60, ha='right', fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "single_cell_heatmap.png"), dpi=200, bbox_inches='tight')
    print("Saved: single_cell_heatmap.png")

    # ── Determine tumor-intrinsic vs stromal ──
    print("\n  Cell-type specificity analysis:")
    specificity_results = []
    hepatocyte_cols = [c for c in sc_results.columns if 'hepato' in c.lower()]
    immune_cols = [c for c in sc_results.columns if any(t in c.lower()
                   for t in ['t-cell', 'b-cell', 'nk', 'macro', 'mono', 'dendri', 'neutro',
                             'kupffer', 'immune'])]
    stromal_cols = [c for c in sc_results.columns if any(t in c.lower()
                    for t in ['stellat', 'fibro', 'endothe', 'smooth'])]

    for gene in genes:
        row = sc_results.loc[gene].dropna()
        if len(row) == 0:
            continue
        total = row.sum()
        if total == 0:
            continue

        hep_expr = row[hepatocyte_cols].mean() if hepatocyte_cols else 0
        imm_expr = row[immune_cols].mean() if immune_cols else 0
        str_expr = row[stromal_cols].mean() if stromal_cols else 0
        top_ct = row.idxmax()

        if hep_expr >= imm_expr and hep_expr >= str_expr:
            origin = "Tumor-intrinsic (hepatocyte)"
        elif imm_expr > hep_expr and imm_expr > str_expr:
            origin = "Immune-derived"
        elif str_expr > 0:
            origin = "Stromal"
        else:
            origin = "Ubiquitous"

        specificity_results.append({
            "gene": gene, "top_cell_type": top_ct,
            "hepatocyte_expr": hep_expr, "immune_expr": imm_expr,
            "stromal_expr": str_expr, "likely_origin": origin,
        })
        print(f"    {gene:<10} Top: {top_ct:<25} → {origin}")

    if specificity_results:
        spec_df = pd.DataFrame(specificity_results)
        spec_df.to_csv(os.path.join(TABLES, "single_cell_specificity.csv"), index=False)
        print("\n  Saved: single_cell_specificity.csv")

    # Save raw expression matrix
    sc_results.to_csv(os.path.join(TABLES, "single_cell_expression.csv"))
    print("  Saved: single_cell_expression.csv")

else:
    print("  No single-cell data available for visualization")
    print("  Creating mRNA-based cell-type inference instead...")

    # Fallback: Use TCGA mRNA data to correlate with known cell-type markers
    tcga = pd.read_csv(os.path.join(DATA, "tcga", "tcga_ros_merged.csv"))
    marker_genes = {
        "Hepatocyte": ["ALB", "HNF4A", "APOB", "CYP3A4"],
        "Cholangiocyte": ["KRT19", "KRT7", "EPCAM"],
        "Macrophage/Kupffer": ["CD68", "CD163", "CSF1R"],
        "T cell": ["CD3D", "CD3E", "CD8A", "CD4"],
        "B cell": ["CD19", "MS4A1", "CD79A"],
        "Endothelial": ["PECAM1", "CDH5", "VWF"],
        "Fibroblast": ["FAP", "ACTA2", "COL1A1"],
    }

    expr_full_path = os.path.join(DATA, "tcga", "tcga_lihc_expression_full.csv")
    if os.path.exists(expr_full_path):
        expr_full = pd.read_csv(expr_full_path, index_col=0)

        fig, ax = plt.subplots(figsize=(10, 8))
        corr_matrix = []
        ct_names = []

        for ct, markers in marker_genes.items():
            available_markers = [m for m in markers if m in expr_full.index]
            if available_markers:
                # Average marker expression per patient
                marker_expr = expr_full.loc[available_markers, tcga["patientId"]].mean(axis=0)
                ct_corrs = []
                for gene in sig_genes:
                    if gene in tcga.columns:
                        from scipy.stats import spearmanr
                        r, _ = spearmanr(tcga[gene].values, marker_expr.values)
                        ct_corrs.append(r)
                    else:
                        ct_corrs.append(0)
                corr_matrix.append(ct_corrs)
                ct_names.append(ct)

        if corr_matrix:
            import seaborn as sns
            corr_df = pd.DataFrame(corr_matrix, index=ct_names, columns=sig_genes)
            sns.heatmap(corr_df, ax=ax, cmap="RdBu_r", center=0,
                        annot=True, fmt=".2f", linewidths=0.5,
                        cbar_kws={"label": "Spearman r"})
            ax.set_title("Signature Gene Correlation with Cell-Type Markers\n(TCGA-LIHC mRNA)",
                         fontsize=13, fontweight='bold')
            plt.tight_layout()
            plt.savefig(os.path.join(FIGS, "single_cell_marker_corr.png"), dpi=200, bbox_inches='tight')
            print("  Saved: single_cell_marker_corr.png")

            corr_df.to_csv(os.path.join(TABLES, "cell_type_marker_correlations.csv"))
            print("  Saved: cell_type_marker_correlations.csv")

print(f"\n✓ Single-cell analysis complete.")
