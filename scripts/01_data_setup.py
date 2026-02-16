"""
01_data_setup.py — Build merged dataset with expanded ROS + ferroptosis gene set

Loads TCGA-LIHC expression + clinical data, filters to ~75 oxidative stress /
ferroptosis / iron metabolism genes, merges, and saves.
"""
import pandas as pd
import numpy as np
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data", "tcga")
os.makedirs(DATA, exist_ok=True)

# ── Expanded gene set (~75 genes) organized by pathway ──────────────────────
GENE_SETS = {
    "Nrf2-KEAP1 axis": ["NFE2L2", "KEAP1", "MAFG", "BACH1"],
    "Superoxide dismutases": ["SOD1", "SOD2", "SOD3"],
    "Catalase": ["CAT"],
    "Glutathione peroxidases": ["GPX1", "GPX2", "GPX3", "GPX4"],
    "Glutathione synthesis & recycling": ["GSR", "GCLC", "GCLM", "GSTP1"],
    "Thioredoxin system": ["TXNRD1", "TXNRD2", "TXN", "TXN2"],
    "Peroxiredoxins": ["PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5", "PRDX6"],
    "Nrf2 target genes": ["NQO1", "HMOX1", "HMOX2", "FTH1", "FTL", "SLC7A11", "SRXN1", "SQSTM1"],
    "Pentose phosphate (NADPH)": ["G6PD", "PGD", "ME1", "IDH1", "IDH2"],
    "Glutaredoxins": ["GLRX", "GLRX2", "MSRA"],
    "NADPH oxidases (ROS producers)": ["NOX1", "NOX4", "CYBB", "NCF1", "NCF2", "RAC1", "DUOX1", "DUOX2"],
    "Eicosanoid pathway": ["PTGS2", "ALOX5", "ALOX12", "ALOX15"],
    "Ferroptosis regulators": ["ACSL4", "LPCAT3", "TFRC", "SLC40A1", "STEAP3", "NCOA4", "HMGCR"],
    "Iron metabolism": ["HAMP", "TFR2", "SLC48A1", "IREB2", "ACO1"],
    "Heat shock proteins": ["HSPH1", "HSP90AA1", "HSPA1A", "HSPB1"],
    "Mitochondrial apoptosis/quality": ["BAX", "BAK1", "BCL2", "PINK1"],
}

all_genes = []
gene_to_pathway = {}
for pathway, genes in GENE_SETS.items():
    for g in genes:
        if g not in gene_to_pathway:
            gene_to_pathway[g] = pathway
            all_genes.append(g)

print(f"Total unique candidate genes: {len(all_genes)}")

# ── Load expression data ────────────────────────────────────────────────────
expr_path = os.path.join(DATA, "tcga_lihc_expression_full.csv")
print(f"Loading expression data from {expr_path} ...")
expr = pd.read_csv(expr_path, index_col=0)
print(f"  Expression matrix: {expr.shape[0]} genes × {expr.shape[1]} samples")

# ── Load clinical data ──────────────────────────────────────────────────────
clin_path = os.path.join(DATA, "tcga_lihc_clinical.csv")
clin = pd.read_csv(clin_path)
print(f"  Clinical data: {len(clin)} patients")

# ── Load stage data ─────────────────────────────────────────────────────────
stage_path = os.path.join(DATA, "tcga_lihc_stage.csv")
stage = pd.read_csv(stage_path)
# Fix column names if needed (stage CSV may have unnamed columns)
if "Unnamed: 0" in stage.columns and "0" in stage.columns:
    stage = stage.rename(columns={"Unnamed: 0": "patientId", "0": "tumor_stage"})
elif stage.columns[0] != "patientId":
    stage.columns = ["patientId", "tumor_stage"]
print(f"  Stage data: {len(stage)} patients")

# ── Filter to available genes ───────────────────────────────────────────────
available = [g for g in all_genes if g in expr.index]
missing = [g for g in all_genes if g not in expr.index]
print(f"\nGene availability: {len(available)}/{len(all_genes)}")
if missing:
    print(f"  Missing: {missing}")

# Report by pathway
print("\nAvailability by pathway:")
for pathway, genes in GENE_SETS.items():
    found = [g for g in genes if g in expr.index]
    print(f"  {pathway}: {len(found)}/{len(genes)}")
    if len(found) < len(genes):
        miss = [g for g in genes if g not in expr.index]
        print(f"    Missing: {miss}")

# ── Build merged dataset ────────────────────────────────────────────────────
# Extract expression for available genes, transpose to patient × gene
expr_sub = expr.loc[available].T
expr_sub.index.name = "patientId"
expr_sub = expr_sub.reset_index()

# Merge with clinical
merged = pd.merge(clin, expr_sub, on="patientId", how="inner")

# Merge stage data
if "tumor_stage" in stage.columns:
    stage_clean = stage[["patientId", "tumor_stage"]].dropna()
    if "tumor_stage" in merged.columns:
        # Replace empty stage column with actual data
        merged = merged.drop(columns=["tumor_stage"])
    merged = pd.merge(merged, stage_clean, on="patientId", how="left")

# Drop patients without survival data
n_before = len(merged)
merged = merged.dropna(subset=["OS_months", "OS_event"])
merged = merged[merged["OS_months"] > 0]
print(f"\nMerged dataset: {len(merged)} patients (dropped {n_before - len(merged)} without survival)")
print(f"  Events: {int(merged['OS_event'].sum())} deaths")
print(f"  Median follow-up: {merged['OS_months'].median():.1f} months")
print(f"  Gene columns: {len(available)}")

# Save
out_path = os.path.join(DATA, "tcga_ros_merged.csv")
merged.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")

# Save gene info
gene_info = pd.DataFrame([
    {"gene": g, "pathway": gene_to_pathway[g], "available": g in available}
    for g in all_genes
])
gene_info.to_csv(os.path.join(DATA, "gene_set_info.csv"), index=False)
print(f"Saved: gene_set_info.csv")
print("\nDone.")
