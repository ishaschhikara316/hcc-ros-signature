"""
05_external_validation.py — Validate on independent cohorts

GSE14520 (n=221, Chinese HBV-HCC, Affymetrix GPL3921) + multivariate Cox
GSE76427 (n=115, HBV-HCC, Illumina) + RFS analysis
ICGC LIRI-JP (n=232, Japanese HCC, RNA-seq via ICGC 25K) — fixed download
"""
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import requests
import gzip
import io
import json
import os
import xml.etree.ElementTree as ET
import warnings
warnings.filterwarnings('ignore')

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
GEO = os.path.join(DATA, "geo_cohorts")
MODEL = os.path.join(BASE, "results", "model")
TABLES = os.path.join(BASE, "results", "tables")
FIGS = os.path.join(BASE, "results", "figures")
os.makedirs(GEO, exist_ok=True)

# ── Load model ──────────────────────────────────────────────────────────────
with open(os.path.join(MODEL, "lasso_model.json")) as f:
    model = json.load(f)
selected_genes = model["genes"]
print(f"Signature genes ({len(selected_genes)}): {list(selected_genes.keys())}")

# ── Helper functions ────────────────────────────────────────────────────────

def download_file(url, local_path, desc=""):
    if os.path.exists(local_path):
        print(f"  Using cached: {os.path.basename(local_path)}")
        return True
    print(f"  Downloading {desc or os.path.basename(local_path)}...")
    try:
        resp = requests.get(url, timeout=180)
        resp.raise_for_status()
        with open(local_path, 'wb') as f:
            f.write(resp.content)
        return True
    except Exception as e:
        print(f"  FAILED: {e}")
        return False


def parse_geo_matrix(filepath):
    """Parse GEO series matrix → expression DataFrame + clinical DataFrame."""
    if filepath.endswith('.gz'):
        content = gzip.open(filepath, 'rt', errors='replace').read()
    else:
        content = open(filepath, 'r', errors='replace').read()
    lines = content.split('\n')

    sample_ids = []
    characteristics = {}
    platform = None
    data_lines = []
    in_data = False

    for line in lines:
        if line.startswith('!Sample_geo_accession'):
            sample_ids = [s.strip('"') for s in line.split('\t')[1:]]
        elif line.startswith('!Sample_platform_id'):
            parts = [s.strip('"') for s in line.split('\t')[1:]]
            if parts:
                platform = parts[0]
        elif line.startswith('!Sample_characteristics_ch1'):
            parts = [p.strip('"') for p in line.split('\t')[1:]]
            if parts and ':' in parts[0]:
                key = parts[0].split(':')[0].strip()
                vals = [p.split(':', 1)[-1].strip() if ':' in p else p for p in parts]
                orig = key
                idx = 1
                while key in characteristics:
                    key = f"{orig}_{idx}"
                    idx += 1
                characteristics[key] = vals
        elif line.startswith('!series_matrix_table_begin'):
            in_data = True
            continue
        elif line.startswith('!series_matrix_table_end'):
            in_data = False
        elif in_data and line.strip():
            data_lines.append(line)

    if not data_lines:
        return None, None, None, platform

    header = data_lines[0].split('\t')
    expr_samples = [s.strip('"') for s in header[1:]]
    rows = []
    row_ids = []
    for dl in data_lines[1:]:
        parts = dl.split('\t')
        if len(parts) <= 1:
            continue
        row_ids.append(parts[0].strip('"'))
        row_vals = []
        for x in parts[1:]:
            x = x.strip().strip('"')
            try:
                row_vals.append(float(x))
            except:
                row_vals.append(np.nan)
        rows.append(row_vals)

    expr = pd.DataFrame(rows, index=row_ids, columns=expr_samples)

    clin = pd.DataFrame(index=sample_ids)
    for key, vals in characteristics.items():
        if len(vals) == len(sample_ids):
            clin[key] = vals

    return expr, clin, sample_ids, platform


def get_platform_mapping(gpl):
    """Download GPL annotation, return probe→gene dict."""
    num = gpl.replace("GPL", "")
    prefix = "GPL" + num[:-3] + "nnn" if len(num) > 3 else "GPLnnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{prefix}/{gpl}/annot/{gpl}.annot.gz"
    cache_path = os.path.join(GEO, f"{gpl}.annot.gz")

    if not download_file(url, cache_path, f"{gpl} annotation"):
        return {}

    content = gzip.open(cache_path, 'rt', errors='replace').read()
    lines = content.split('\n')

    data_lines = []
    header = None
    for line in lines:
        if line.startswith('#') or line.startswith('^') or line.startswith('!'):
            continue
        if not line.strip():
            continue
        if header is None:
            header = line.split('\t')
            continue
        data_lines.append(line.split('\t'))

    if header is None or not data_lines:
        print(f"  Could not parse {gpl} annotation")
        return {}

    gene_col = None
    for i, col in enumerate(header):
        cl = col.lower().strip().strip('"')
        if cl in ('gene symbol', 'gene_symbol', 'symbol', 'genesymbol'):
            gene_col = i
            break
    if gene_col is None:
        for i, col in enumerate(header):
            cl = col.lower().strip().strip('"')
            if 'gene' in cl and 'symbol' in cl:
                gene_col = i
                break

    if gene_col is None:
        print(f"  No gene symbol column in: {[h.strip('\"') for h in header[:10]]}")
        return get_platform_mapping_soft(gpl)

    mapping = {}
    for row in data_lines:
        if len(row) > gene_col:
            probe = row[0].strip().strip('"')
            symbol = row[gene_col].strip().strip('"')
            symbol = symbol.split('///')[0].strip()
            if symbol and symbol.lower() not in ('', 'nan', '---', 'na', 'n/a'):
                mapping[probe] = symbol

    print(f"  {gpl} mapping: {len(mapping)} probes → genes")
    return mapping


def get_platform_mapping_soft(gpl):
    """Fallback: download SOFT format platform file."""
    num = gpl.replace("GPL", "")
    prefix = "GPL" + num[:-3] + "nnn" if len(num) > 3 else "GPLnnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{prefix}/{gpl}/soft/{gpl}_family.soft.gz"
    cache_path = os.path.join(GEO, f"{gpl}_family.soft.gz")

    if not download_file(url, cache_path, f"{gpl} SOFT"):
        return {}

    content = gzip.open(cache_path, 'rt', errors='replace').read()
    lines = content.split('\n')

    in_table = False
    header = None
    gene_col = None
    mapping = {}

    for line in lines:
        if line.startswith('!platform_table_begin'):
            in_table = True
            continue
        elif line.startswith('!platform_table_end'):
            break
        elif in_table:
            if header is None:
                header = line.split('\t')
                for i, col in enumerate(header):
                    cl = col.lower().strip()
                    if cl in ('gene symbol', 'gene_symbol', 'symbol', 'gene_assignment'):
                        gene_col = i
                        break
                if gene_col is None:
                    for i, col in enumerate(header):
                        if 'gene' in col.lower() and ('symbol' in col.lower() or 'name' in col.lower()):
                            gene_col = i
                            break
                continue
            if gene_col is not None:
                parts = line.split('\t')
                if len(parts) > gene_col:
                    probe = parts[0].strip()
                    symbol = parts[gene_col].strip()
                    if 'gene_assignment' in header[gene_col].lower() and '//' in symbol:
                        symbol = symbol.split('//')[1].strip()
                    symbol = symbol.split('///')[0].strip()
                    if symbol and symbol.lower() not in ('', 'nan', '---', 'na'):
                        mapping[probe] = symbol

    print(f"  {gpl} SOFT mapping: {len(mapping)} probes → genes")
    return mapping


def expr_to_genes(expr, mapping):
    """Convert probe-level to gene-level expression."""
    expr = expr.copy()
    expr["gene"] = expr.index.map(lambda x: mapping.get(str(x), np.nan))
    expr = expr[expr["gene"].notna() & (expr["gene"].astype(str) != "nan")]
    expr["m"] = expr.drop("gene", axis=1).mean(axis=1)
    expr = expr.sort_values("m", ascending=False).drop_duplicates(subset="gene", keep="first")
    expr = expr.set_index("gene").drop("m", axis=1)
    return expr


def compute_risk_score(expr, samples, coefs):
    """Compute risk scores on validation cohort."""
    found = [g for g in coefs if g in expr.index]
    missing = [g for g in coefs if g not in expr.index]
    if missing:
        print(f"  Found {len(found)}/{len(coefs)} genes, missing: {missing}")
    if len(found) < 3:
        return None, found, missing

    scores = np.zeros(len(samples))
    for gene in found:
        vals = expr.loc[gene, samples].values.astype(float)
        z = (vals - np.nanmean(vals)) / (np.nanstd(vals) + 1e-10)
        scores += coefs[gene] * z

    return scores, found, missing


def validate_cohort(name, os_months, os_event, risk_score, ax=None):
    """Standard validation: C-index, bootstrap CI, KM, Cox HR."""
    valid = pd.DataFrame({
        "os_months": os_months, "os_event": os_event, "risk_score": risk_score
    }).dropna()
    valid = valid[(valid["os_months"] > 0) & (valid["os_event"].isin([0, 1]))]

    if len(valid) < 20 or valid["os_event"].sum() < 5:
        print(f"  {name}: Too few patients/events ({len(valid)}, {int(valid['os_event'].sum())} events)")
        return None

    ci = concordance_index(valid["os_months"], -valid["risk_score"], valid["os_event"])
    boot = []
    for _ in range(1000):
        idx = np.random.choice(len(valid), len(valid), replace=True)
        bd = valid.iloc[idx]
        try:
            boot.append(concordance_index(bd["os_months"], -bd["risk_score"], bd["os_event"]))
        except:
            pass
    ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5]) if boot else (np.nan, np.nan)

    med = valid["risk_score"].median()
    high = valid[valid["risk_score"] >= med]
    low = valid[valid["risk_score"] < med]
    lr = logrank_test(high["os_months"], low["os_months"],
                      event_observed_A=high["os_event"], event_observed_B=low["os_event"])

    try:
        cph = CoxPHFitter()
        cph.fit(valid[["os_months", "os_event", "risk_score"]], duration_col="os_months", event_col="os_event")
        hr = np.exp(cph.params_["risk_score"])
        hr_ci = np.exp(cph.confidence_intervals_.values[0])
    except:
        hr, hr_ci = np.nan, [np.nan, np.nan]

    result = {
        "cohort": name, "n": len(valid), "events": int(valid["os_event"].sum()),
        "c_index": ci, "ci_lo": ci_lo, "ci_hi": ci_hi,
        "HR": hr, "HR_lower": hr_ci[0], "HR_upper": hr_ci[1],
        "logrank_p": lr.p_value,
    }
    print(f"  {name}: n={len(valid)}, events={int(valid['os_event'].sum())}, "
          f"C-index={ci:.3f} ({ci_lo:.3f}-{ci_hi:.3f}), HR={hr:.2f}, p={lr.p_value:.2e}")

    # Flag underpowered cohorts
    if valid["os_event"].sum() < 30:
        print(f"  ⚠ NOTE: Only {int(valid['os_event'].sum())} events — results may be underpowered")

    if ax is not None:
        kmf = KaplanMeierFitter()
        kmf.fit(high["os_months"], high["os_event"], label=f"High (n={len(high)})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color="red")
        kmf.fit(low["os_months"], low["os_event"], label=f"Low (n={len(low)})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color="blue")
        title = f"{name}\nC={ci:.3f}, HR={hr:.2f}, p={lr.p_value:.2e}"
        if valid["os_event"].sum() < 30:
            title += f"\n(underpowered: {int(valid['os_event'].sum())} events)"
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.set_xlabel("Months")
        ax.set_ylabel("Survival Probability")
        ax.legend(fontsize=8)
        ax.set_xlim(0, max(80, valid["os_months"].quantile(0.95)))

    return result


# ══════════════════════════════════════════════════════════════════════════════
# COHORT 1: GSE14520
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("COHORT 1: GSE14520 (Chinese HBV-HCC)")
print("=" * 70)

gse14520_results = None
gse14520_mv_results = []
# Keep references for signature comparison script
gse14520_expr_genes = None
gse14520_suppl = None
gse14520_matched = None
gse14520_surv_time_col = None
gse14520_surv_event_col = None

matrix_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/GSE14520-GPL3921_series_matrix.txt.gz"
matrix_path = os.path.join(GEO, "GSE14520-GPL3921_series_matrix.txt.gz")

if download_file(matrix_url, matrix_path, "GSE14520 matrix"):
    expr14, clin14, samples14, platform14 = parse_geo_matrix(matrix_path)
    if expr14 is not None:
        print(f"  Expression: {expr14.shape[0]} probes × {expr14.shape[1]} samples")
        if not platform14:
            platform14 = "GPL3921"
        print(f"  Platform: {platform14}")

        mapping14 = get_platform_mapping(platform14)
        if mapping14:
            expr14_genes = expr_to_genes(expr14, mapping14)
            gse14520_expr_genes = expr14_genes
            print(f"  Gene-level: {len(expr14_genes)} genes")

            suppl_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/suppl/GSE14520_Extra_Supplement.txt.gz"
            suppl_path = os.path.join(GEO, "GSE14520_Extra_Supplement.txt.gz")

            if download_file(suppl_url, suppl_path, "GSE14520 supplement"):
                suppl = pd.read_csv(suppl_path, sep='\t', compression='gzip')
                print(f"  Supplement: {len(suppl)} rows, columns: {list(suppl.columns)[:10]}")

                surv_time_col = None
                surv_event_col = None
                gsm_col = None
                for col in suppl.columns:
                    cl = col.lower()
                    if 'survival' in cl and 'month' in cl:
                        surv_time_col = col
                    elif 'survival' in cl and 'status' in cl:
                        surv_event_col = col
                    elif cl == 'affy_gsm':
                        gsm_col = col
                if gsm_col is None:
                    for col in suppl.columns:
                        if 'gsm' in col.lower():
                            gsm_col = col
                            break

                if surv_time_col and surv_event_col and gsm_col:
                    suppl = suppl.dropna(subset=[surv_time_col, surv_event_col, gsm_col])
                    suppl[gsm_col] = suppl[gsm_col].astype(str)
                    suppl = suppl.set_index(gsm_col)
                    gse14520_suppl = suppl
                    gse14520_surv_time_col = surv_time_col
                    gse14520_surv_event_col = surv_event_col

                    matched = [s for s in suppl.index if s in expr14_genes.columns]
                    gse14520_matched = matched
                    print(f"  Matched samples: {len(matched)}")

                    if len(matched) >= 20:
                        scores14, found14, miss14 = compute_risk_score(expr14_genes, matched, selected_genes)
                        if scores14 is not None:
                            gse14520_results = validate_cohort(
                                "GSE14520",
                                suppl.loc[matched, surv_time_col].values.astype(float),
                                suppl.loc[matched, surv_event_col].values.astype(float),
                                scores14
                            )

                            # ── GSE14520 Multivariate Cox ──────────────────────
                            print("\n  --- GSE14520 Multivariate Cox ---")
                            mv_df = pd.DataFrame({
                                "os_months": suppl.loc[matched, surv_time_col].values.astype(float),
                                "os_event": suppl.loc[matched, surv_event_col].values.astype(float),
                                "risk_score": scores14,
                            }, index=matched)

                            # Extract covariates from supplement
                            for col_name, mv_name in [("Age", "age"), ("Gender", "sex"),
                                                       ("TNM staging", "tnm_stage")]:
                                if col_name in suppl.columns:
                                    mv_df[mv_name] = suppl.loc[matched, col_name].values

                            # Encode sex
                            if "sex" in mv_df.columns:
                                mv_df["sex"] = (mv_df["sex"].astype(str).str.upper().str.strip() == "M").astype(float)
                                mv_df.loc[mv_df["sex"].isna(), "sex"] = np.nan

                            # Encode age
                            if "age" in mv_df.columns:
                                mv_df["age"] = pd.to_numeric(mv_df["age"], errors="coerce")

                            # Encode TNM stage
                            if "tnm_stage" in mv_df.columns:
                                def parse_tnm(s):
                                    s = str(s).strip().upper()
                                    for k, v in [("IV", 4), ("III", 3), ("II", 2), ("I", 1)]:
                                        if k in s:
                                            return v
                                    return np.nan
                                mv_df["tnm_stage"] = mv_df["tnm_stage"].apply(parse_tnm)

                            # Build multivariate model
                            mv_cols = ["os_months", "os_event", "risk_score"]
                            for cov in ["age", "sex", "tnm_stage"]:
                                if cov in mv_df.columns and mv_df[cov].notna().sum() > 30:
                                    mv_cols.append(cov)

                            mv_fit = mv_df[mv_cols].dropna()
                            mv_fit = mv_fit[mv_fit["os_months"] > 0]
                            print(f"  Multivariate model variables: {mv_cols[2:]}")
                            print(f"  Samples for multivariate: {len(mv_fit)}")

                            if len(mv_fit) >= 30 and len(mv_cols) > 3:
                                cph_mv = CoxPHFitter()
                                cph_mv.fit(mv_fit, duration_col="os_months", event_col="os_event")
                                print(f"\n  Multivariate Cox (GSE14520):")
                                print(cph_mv.summary[["coef", "exp(coef)", "exp(coef) lower 95%",
                                                       "exp(coef) upper 95%", "p"]].to_string())

                                for var in mv_cols[2:]:
                                    gse14520_mv_results.append({
                                        "Cohort": "GSE14520", "Analysis": "Multivariate",
                                        "Variable": var,
                                        "HR": cph_mv.summary.loc[var, "exp(coef)"],
                                        "CI_low": cph_mv.summary.loc[var, "exp(coef) lower 95%"],
                                        "CI_high": cph_mv.summary.loc[var, "exp(coef) upper 95%"],
                                        "p_value": cph_mv.summary.loc[var, "p"],
                                        "N": len(mv_fit),
                                    })
                            else:
                                print("  Insufficient data for multivariate model")
                else:
                    print(f"  Could not find survival columns in supplement: {list(suppl.columns)}")


# ══════════════════════════════════════════════════════════════════════════════
# COHORT 2: GSE76427 (OS + RFS)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("COHORT 2: GSE76427 (HBV-HCC, Illumina)")
print("=" * 70)

gse76427_results = None
gse76427_rfs_results = None

gse76_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76427/matrix/GSE76427_series_matrix.txt.gz"
gse76_path = os.path.join(GEO, "GSE76427_series_matrix.txt.gz")

if download_file(gse76_url, gse76_path, "GSE76427"):
    expr76, clin76, samples76, platform76 = parse_geo_matrix(gse76_path)
    if expr76 is not None:
        print(f"  Expression: {expr76.shape[0]} probes × {expr76.shape[1]} samples")
        print(f"  Platform: {platform76}")
        print(f"  Clinical columns: {list(clin76.columns)}")

        if platform76:
            mapping76 = get_platform_mapping(platform76)
        else:
            mapping76 = {}

        if mapping76:
            expr76_genes = expr_to_genes(expr76, mapping76)
        else:
            expr76_genes = expr76
        print(f"  Gene-level: {len(expr76_genes)} genes")

        # Filter tumor samples
        tumor_mask = pd.Series(True, index=clin76.index)
        for col in clin76.columns:
            if clin76[col].astype(str).str.contains('tumor', case=False, na=False).any():
                is_tumor = clin76[col].astype(str).str.contains('tumor', case=False, na=False)
                is_nontumor = clin76[col].astype(str).str.contains('non-tumor|adjacent|normal', case=False, na=False)
                mask = is_tumor & ~is_nontumor
                if mask.sum() > 10:
                    tumor_mask = mask
                    print(f"  Filtered to {mask.sum()} tumor samples")
                    break

        clin76_tumor = clin76[tumor_mask]

        # ── OS analysis ──
        tc = [c for c in clin76.columns if 'duryears_os' in c.lower()]
        ec = [c for c in clin76.columns if 'event_os' in c.lower()]

        if tc and ec:
            clin76_tumor["os_months"] = pd.to_numeric(clin76_tumor[tc[0]], errors='coerce') * 12
            clin76_tumor["os_event"] = pd.to_numeric(clin76_tumor[ec[0]], errors='coerce')

            valid76 = [s for s in clin76_tumor.index if s in expr76_genes.columns
                       and pd.notna(clin76_tumor.loc[s, "os_months"])
                       and clin76_tumor.loc[s, "os_months"] > 0]
            print(f"  Valid samples with OS: {len(valid76)}")

            scores76, found76, miss76 = compute_risk_score(expr76_genes, valid76, selected_genes)
            if scores76 is not None:
                gse76427_results = validate_cohort(
                    "GSE76427 (OS)",
                    clin76_tumor.loc[valid76, "os_months"].values.astype(float),
                    clin76_tumor.loc[valid76, "os_event"].values.astype(float),
                    scores76
                )
        else:
            print(f"  No OS columns found. Available: {list(clin76.columns)}")

        # ── RFS analysis ──
        rfs_tc = [c for c in clin76.columns if 'duryears_rfs' in c.lower()]
        rfs_ec = [c for c in clin76.columns if 'event_rfs' in c.lower()]

        if rfs_tc and rfs_ec:
            print("\n  --- GSE76427 Recurrence-Free Survival (RFS) ---")
            clin76_tumor["rfs_months"] = pd.to_numeric(clin76_tumor[rfs_tc[0]], errors='coerce') * 12
            clin76_tumor["rfs_event"] = pd.to_numeric(clin76_tumor[rfs_ec[0]], errors='coerce')

            valid76_rfs = [s for s in clin76_tumor.index if s in expr76_genes.columns
                           and pd.notna(clin76_tumor.loc[s, "rfs_months"])
                           and clin76_tumor.loc[s, "rfs_months"] > 0]
            print(f"  Valid samples with RFS: {len(valid76_rfs)}")

            if len(valid76_rfs) >= 20:
                scores76_rfs, _, _ = compute_risk_score(expr76_genes, valid76_rfs, selected_genes)
                if scores76_rfs is not None:
                    gse76427_rfs_results = validate_cohort(
                        "GSE76427 (RFS)",
                        clin76_tumor.loc[valid76_rfs, "rfs_months"].values.astype(float),
                        clin76_tumor.loc[valid76_rfs, "rfs_event"].values.astype(float),
                        scores76_rfs
                    )
        else:
            print(f"  No RFS columns found.")


# ══════════════════════════════════════════════════════════════════════════════
# COHORT 3: ICGC LIRI-JP (fixed download with fallback part files)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("COHORT 3: ICGC LIRI-JP (Japanese HCC, RNA-seq)")
print("=" * 70)

icgc_results = None
icgc_mv_results = []

ICGC_BASE = "https://object.genomeinformatics.org/icgc25k-open/release_28"
icgc_donor_cache = os.path.join(GEO, "liri_jp_donors.csv")
icgc_expr_cache = os.path.join(GEO, "liri_jp_expression_merged.csv")

if os.path.exists(icgc_donor_cache) and os.path.exists(icgc_expr_cache):
    print("  Using cached ICGC LIRI-JP data")
    icgc_donors = pd.read_csv(icgc_donor_cache, index_col=0)
    icgc_expr = pd.read_csv(icgc_expr_cache, index_col=0)
else:
    try:
        # Step 1: List all LIRI-JP donor IDs (proper XML namespace parsing)
        print("  Listing LIRI-JP donors from ICGC 25K...")
        list_url = (f"{ICGC_BASE.replace('/release_28', '')}?"
                    f"list-type=2&prefix=release_28/data/LIRI-JP/&delimiter=/&max-keys=500")
        resp = requests.get(list_url, timeout=60)
        if resp.status_code != 200:
            raise RuntimeError(f"Cannot list LIRI-JP donors: HTTP {resp.status_code}")

        root = ET.fromstring(resp.content)
        ns = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}
        prefixes = root.findall(".//s3:CommonPrefixes/s3:Prefix", ns)
        if not prefixes:
            prefixes = root.findall(".//CommonPrefixes/Prefix")
        donor_ids = []
        for p in prefixes:
            parts = p.text.strip("/").split("/")
            if len(parts) >= 3:
                donor_ids.append(parts[-1])
        print(f"  Found {len(donor_ids)} donors")
        if len(donor_ids) < 50:
            raise RuntimeError(f"Expected ~260 donors, found {len(donor_ids)}")

        # Step 2: Download headers
        print("  Downloading header files...")
        resp_hdr = requests.get(f"{ICGC_BASE}/headers/donor.tsv.gz", timeout=30)
        if resp_hdr.status_code != 200:
            raise RuntimeError(f"Cannot download donor header: HTTP {resp_hdr.status_code}")
        donor_header = gzip.decompress(resp_hdr.content).decode("utf-8").strip().split("\t")

        resp_hdr_exp = requests.get(f"{ICGC_BASE}/headers/exp_seq.tsv.gz", timeout=30)
        if resp_hdr_exp.status_code != 200:
            raise RuntimeError(f"Cannot download exp_seq header: HTTP {resp_hdr_exp.status_code}")
        exp_header = gzip.decompress(resp_hdr_exp.content).decode("utf-8").strip().split("\t")

        # Determine column indices
        gene_id_idx = exp_header.index("gene_id") if "gene_id" in exp_header else 7
        donor_id_idx = exp_header.index("icgc_donor_id") if "icgc_donor_id" in exp_header else 0
        value_idx = exp_header.index("normalized_read_count") if "normalized_read_count" in exp_header else 8
        print(f"  Donor header ({len(donor_header)}): {donor_header[:8]}")
        print(f"  Exp header ({len(exp_header)}): gene_id@{gene_id_idx}, donor@{donor_id_idx}, value@{value_idx}")

        # Step 3: Download donor clinical data with fallback
        print("  Downloading donor clinical data...")
        donor_rows = []
        failed_donors = 0
        sig_gene_names = set(selected_genes.keys())

        for i, did in enumerate(donor_ids):
            url = f"{ICGC_BASE}/data/LIRI-JP/{did}/donor/part-00002.gz"
            try:
                r = requests.get(url, timeout=15)
                if r.status_code == 200:
                    content = gzip.decompress(r.content).decode("utf-8").strip()
                    for line in content.split("\n"):
                        if line.strip():
                            donor_rows.append(line.split("\t"))
                elif r.status_code == 404:
                    # Fallback: try part-00000 through part-00004
                    for part_num in range(5):
                        url2 = f"{ICGC_BASE}/data/LIRI-JP/{did}/donor/part-{part_num:05d}.gz"
                        r2 = requests.get(url2, timeout=10)
                        if r2.status_code == 200:
                            content = gzip.decompress(r2.content).decode("utf-8").strip()
                            for line in content.split("\n"):
                                if line.strip():
                                    donor_rows.append(line.split("\t"))
                            break
                    else:
                        failed_donors += 1
                else:
                    failed_donors += 1
            except Exception:
                failed_donors += 1
            if (i + 1) % 50 == 0:
                print(f"    {i+1}/{len(donor_ids)} donors processed...")

        print(f"  Donor data: {len(donor_rows)} rows, {failed_donors} failed")

        # Step 4: Download expression data with fallback
        print(f"  Downloading expression data (filtering to {len(sig_gene_names)} signature genes)...")
        exp_rows = []
        failed_exp = 0
        donors_with_data = 0

        for i, did in enumerate(donor_ids):
            url = f"{ICGC_BASE}/data/LIRI-JP/{did}/exp_seq/part-00000.gz"
            try:
                r = requests.get(url, timeout=30)
                if r.status_code == 200:
                    content = gzip.decompress(r.content).decode("utf-8").strip()
                    found_any = False
                    for line in content.split("\n"):
                        if not line.strip():
                            continue
                        fields = line.split("\t")
                        if len(fields) > max(gene_id_idx, value_idx, donor_id_idx):
                            gene_name = fields[gene_id_idx]
                            if gene_name in sig_gene_names:
                                exp_rows.append(fields)
                                found_any = True
                    if found_any:
                        donors_with_data += 1
                elif r.status_code == 404:
                    for part_num in range(5):
                        url2 = f"{ICGC_BASE}/data/LIRI-JP/{did}/exp_seq/part-{part_num:05d}.gz"
                        r2 = requests.get(url2, timeout=20)
                        if r2.status_code == 200:
                            content = gzip.decompress(r2.content).decode("utf-8").strip()
                            found_any = False
                            for line in content.split("\n"):
                                if not line.strip():
                                    continue
                                fields = line.split("\t")
                                if len(fields) > max(gene_id_idx, value_idx, donor_id_idx):
                                    gene_name = fields[gene_id_idx]
                                    if gene_name in sig_gene_names:
                                        exp_rows.append(fields)
                                        found_any = True
                            if found_any:
                                donors_with_data += 1
                            break
                    else:
                        failed_exp += 1
                else:
                    failed_exp += 1
            except Exception:
                failed_exp += 1
            if (i + 1) % 50 == 0:
                print(f"    {i+1}/{len(donor_ids)} donors processed "
                      f"({len(exp_rows)} gene rows, {donors_with_data} donors with data)...")

        print(f"  Expression data: {len(exp_rows)} rows from {donors_with_data} donors, {failed_exp} failed")

        # Build donor DataFrame
        if len(donor_rows) > 0:
            icgc_donors = pd.DataFrame(donor_rows,
                                        columns=donor_header[:len(donor_rows[0])])
            # Deduplicate
            id_col = "icgc_donor_id" if "icgc_donor_id" in icgc_donors.columns else icgc_donors.columns[0]
            icgc_donors = icgc_donors.drop_duplicates(subset=[id_col], keep="first")
            icgc_donors.index = icgc_donors[id_col]
            print(f"  Donor DF: {icgc_donors.shape}")
            print(f"  Columns: {list(icgc_donors.columns)}")
        else:
            raise RuntimeError("No donor data downloaded")

        # Parse survival
        surv_time_col_icgc = None
        surv_event_col_icgc = None
        for col in icgc_donors.columns:
            if 'survival_time' in col.lower():
                surv_time_col_icgc = col
            if 'vital_status' in col.lower():
                surv_event_col_icgc = col

        if surv_time_col_icgc and surv_event_col_icgc:
            icgc_donors['os_time'] = pd.to_numeric(icgc_donors[surv_time_col_icgc], errors='coerce')
            # Convert days to months if median > 100
            if icgc_donors['os_time'].median() > 100:
                print(f"  Time in days (median={icgc_donors['os_time'].median():.0f}), converting to months")
                icgc_donors['os_months'] = icgc_donors['os_time'] / 30.44
            else:
                icgc_donors['os_months'] = icgc_donors['os_time']

            event_vals = icgc_donors[surv_event_col_icgc].astype(str).str.lower().str.strip()
            icgc_donors['os_event'] = event_vals.isin(['deceased', 'dead']).astype(int)
            print(f"  Vital status values: {icgc_donors[surv_event_col_icgc].unique()}")
        else:
            print(f"  WARNING: Missing survival columns. Available: {list(icgc_donors.columns)}")

        # Parse covariates for multivariate analysis
        # Age
        age_col = None
        for col in icgc_donors.columns:
            if 'age' in col.lower() and 'diagnosis' in col.lower():
                age_col = col
                break
        if age_col is None:
            for col in icgc_donors.columns:
                if 'age' in col.lower():
                    age_col = col
                    break
        if age_col:
            icgc_donors['age_years'] = pd.to_numeric(icgc_donors[age_col], errors='coerce')
            if icgc_donors['age_years'].median() > 200:
                icgc_donors['age_years'] = icgc_donors['age_years'] / 365.25
            print(f"  Age column: {age_col}, median={icgc_donors['age_years'].median():.1f} years")

        # Sex
        sex_col = None
        for col in icgc_donors.columns:
            if 'sex' in col.lower() or 'gender' in col.lower():
                sex_col = col
                break
        if sex_col:
            sex_vals = icgc_donors[sex_col].astype(str).str.lower().str.strip()
            icgc_donors['sex_binary'] = sex_vals.map({"male": 1, "female": 0, "m": 1, "f": 0}).astype(float)
            print(f"  Sex column: {sex_col}, values: {icgc_donors[sex_col].unique()}")

        # Stage
        stage_col = None
        for col in icgc_donors.columns:
            if 'stage' in col.lower() and 'tumour' in col.lower():
                stage_col = col
                break
        if stage_col is None:
            for col in icgc_donors.columns:
                if 'stage' in col.lower():
                    stage_col = col
                    break
        if stage_col:
            def parse_stage(s):
                s = str(s).strip().lower()
                for k, v in [("iv", 4), ("iii", 3), ("ii", 2), ("i", 1)]:
                    if k in s:
                        return v
                try:
                    val = int(float(s))
                    if 1 <= val <= 4:
                        return val
                except (ValueError, TypeError):
                    pass
                return np.nan
            icgc_donors['stage_numeric'] = icgc_donors[stage_col].apply(parse_stage)
            print(f"  Stage column: {stage_col}")
            print(f"  Stage parsed: {icgc_donors['stage_numeric'].value_counts().to_dict()}")

        # Build expression matrix
        if len(exp_rows) > 0:
            exp_df = pd.DataFrame(exp_rows, columns=exp_header[:len(exp_rows[0])])
            gene_col_name = "gene_id"
            value_col_name = "normalized_read_count"
            donor_id_col_name = "icgc_donor_id"

            found_direct = [g for g in selected_genes if g in set(exp_df[gene_col_name].unique())]
            print(f"  Signature genes in expression data: {len(found_direct)}/{len(selected_genes)} — {found_direct}")

            if len(found_direct) >= 3:
                exp_lasso = exp_df[exp_df[gene_col_name].isin(selected_genes.keys())].copy()
                exp_lasso[value_col_name] = pd.to_numeric(exp_lasso[value_col_name], errors='coerce')
                icgc_expr = exp_lasso.pivot_table(
                    index=gene_col_name, columns=donor_id_col_name,
                    values=value_col_name, aggfunc='mean'
                )
                print(f"  Expression matrix: {icgc_expr.shape}")
                # Log2 transform if raw counts
                if icgc_expr.median().median() > 50:
                    icgc_expr = np.log2(icgc_expr + 1)
                    print("  Applied log2 transformation")
            else:
                unique_genes = exp_df[gene_col_name].unique()
                print(f"  Gene IDs look like: {unique_genes[:5]}")
                print(f"  WARNING: Only {len(found_direct)} signature genes found — may use Ensembl IDs")
                icgc_expr = pd.DataFrame()
        else:
            icgc_expr = pd.DataFrame()

        # Cache
        if len(icgc_donors) > 0:
            icgc_donors.to_csv(icgc_donor_cache)
        if len(icgc_expr) > 0:
            icgc_expr.to_csv(icgc_expr_cache)

    except Exception as e:
        print(f"  Error downloading ICGC data: {e}")
        import traceback
        traceback.print_exc()
        icgc_donors = pd.DataFrame()
        icgc_expr = pd.DataFrame()

# Validate ICGC
if len(icgc_expr) > 0 and 'os_months' in icgc_donors.columns:
    common = [d for d in icgc_donors.index if d in icgc_expr.columns
              and pd.notna(icgc_donors.loc[d, 'os_months'])
              and icgc_donors.loc[d, 'os_months'] > 0]
    print(f"  Matched donors with survival: {len(common)}")

    if len(common) >= 20:
        scores_icgc, found_icgc, miss_icgc = compute_risk_score(icgc_expr, common, selected_genes)
        if scores_icgc is not None:
            icgc_results = validate_cohort(
                "ICGC LIRI-JP",
                icgc_donors.loc[common, "os_months"].values.astype(float),
                icgc_donors.loc[common, "os_event"].values.astype(float),
                scores_icgc
            )

            # ── ICGC Multivariate Cox ──
            print("\n  --- ICGC LIRI-JP Multivariate Cox ---")
            mv_icgc_df = pd.DataFrame({
                "os_months": icgc_donors.loc[common, "os_months"].values.astype(float),
                "os_event": icgc_donors.loc[common, "os_event"].values.astype(float),
                "risk_score": scores_icgc,
            }, index=common)

            # Add covariates
            mv_icgc_cols = ["os_months", "os_event", "risk_score"]
            for cov_col, cov_name in [("age_years", "age"), ("sex_binary", "sex"),
                                       ("stage_numeric", "stage")]:
                if cov_col in icgc_donors.columns:
                    mv_icgc_df[cov_name] = icgc_donors.loc[common, cov_col].values.astype(float)
                    if mv_icgc_df[cov_name].notna().sum() > 20:
                        mv_icgc_cols.append(cov_name)

            mv_icgc_fit = mv_icgc_df[mv_icgc_cols].dropna()
            mv_icgc_fit = mv_icgc_fit[mv_icgc_fit["os_months"] > 0]
            print(f"  Multivariate variables: {mv_icgc_cols[2:]}")
            print(f"  Samples for multivariate: {len(mv_icgc_fit)}")

            if len(mv_icgc_fit) >= 20 and len(mv_icgc_cols) > 3:
                cph_mv_icgc = CoxPHFitter()
                cph_mv_icgc.fit(mv_icgc_fit, duration_col="os_months", event_col="os_event")
                print(f"\n  Multivariate Cox (ICGC LIRI-JP):")
                print(cph_mv_icgc.summary[["coef", "exp(coef)", "exp(coef) lower 95%",
                                             "exp(coef) upper 95%", "p"]].to_string())

                for var in mv_icgc_cols[2:]:
                    icgc_mv_results.append({
                        "Cohort": "ICGC LIRI-JP", "Analysis": "Multivariate",
                        "Variable": var,
                        "HR": cph_mv_icgc.summary.loc[var, "exp(coef)"],
                        "CI_low": cph_mv_icgc.summary.loc[var, "exp(coef) lower 95%"],
                        "CI_high": cph_mv_icgc.summary.loc[var, "exp(coef) upper 95%"],
                        "p_value": cph_mv_icgc.summary.loc[var, "p"],
                        "N": len(mv_icgc_fit),
                    })
            else:
                print("  Insufficient data for multivariate model")


# ══════════════════════════════════════════════════════════════════════════════
# SAVE MULTIVARIATE RESULTS
# ══════════════════════════════════════════════════════════════════════════════
all_mv = gse14520_mv_results + icgc_mv_results
if all_mv:
    mv_out = pd.DataFrame(all_mv)
    mv_out.to_csv(os.path.join(TABLES, "validation_multivariate.csv"), index=False)
    print(f"\nSaved: validation_multivariate.csv ({len(mv_out)} rows)")
    print(mv_out.to_string(index=False))


# ══════════════════════════════════════════════════════════════════════════════
# COMBINED RESULTS
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("VALIDATION SUMMARY")
print("=" * 70)

all_results = []
if gse14520_results: all_results.append(gse14520_results)
if gse76427_results: all_results.append(gse76427_results)
if gse76427_rfs_results: all_results.append(gse76427_rfs_results)
if icgc_results: all_results.append(icgc_results)

if all_results:
    val_df = pd.DataFrame(all_results)
    val_df.to_csv(os.path.join(TABLES, "validation_results.csv"), index=False)
    print(val_df[["cohort", "n", "events", "c_index", "HR", "logrank_p"]].to_string(index=False))

    # ── KM plots for all cohorts ────────────────────────────────────────────
    n_c = len(all_results)
    fig, axes = plt.subplots(1, n_c, figsize=(6 * n_c, 5))
    if n_c == 1:
        axes = [axes]

    for i, res in enumerate(all_results):
        name = res["cohort"]
        if name == "GSE14520" and gse14520_results and gse14520_matched:
            scores_plot, _, _ = compute_risk_score(gse14520_expr_genes, gse14520_matched, selected_genes)
            if scores_plot is not None:
                validate_cohort(name,
                    gse14520_suppl.loc[gse14520_matched, gse14520_surv_time_col].values.astype(float),
                    gse14520_suppl.loc[gse14520_matched, gse14520_surv_event_col].values.astype(float),
                    scores_plot, ax=axes[i])
        elif "GSE76427" in name and "RFS" not in name:
            scores_plot76, _, _ = compute_risk_score(expr76_genes, valid76, selected_genes)
            if scores_plot76 is not None:
                validate_cohort(name,
                    clin76_tumor.loc[valid76, "os_months"].values.astype(float),
                    clin76_tumor.loc[valid76, "os_event"].values.astype(float),
                    scores_plot76, ax=axes[i])
        elif "GSE76427" in name and "RFS" in name:
            scores_rfs_plot, _, _ = compute_risk_score(expr76_genes, valid76_rfs, selected_genes)
            if scores_rfs_plot is not None:
                validate_cohort(name,
                    clin76_tumor.loc[valid76_rfs, "rfs_months"].values.astype(float),
                    clin76_tumor.loc[valid76_rfs, "rfs_event"].values.astype(float),
                    scores_rfs_plot, ax=axes[i])
        elif name == "ICGC LIRI-JP" and icgc_results:
            common_icgc = [d for d in icgc_donors.index if d in icgc_expr.columns
                           and pd.notna(icgc_donors.loc[d, 'os_months'])
                           and icgc_donors.loc[d, 'os_months'] > 0]
            scores_plot_icgc, _, _ = compute_risk_score(icgc_expr, common_icgc, selected_genes)
            if scores_plot_icgc is not None:
                validate_cohort(name,
                    icgc_donors.loc[common_icgc, "os_months"].values.astype(float),
                    icgc_donors.loc[common_icgc, "os_event"].values.astype(float),
                    scores_plot_icgc, ax=axes[i])

    plt.suptitle("External Validation — ROS/Ferroptosis Prognostic Signature", fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "validation_km.png"), dpi=200, bbox_inches='tight')
    print("\nSaved: validation_km.png")

    # ── Forest plot across cohorts ──────────────────────────────────────────
    train_perf = pd.read_csv(os.path.join(TABLES, "training_performance.csv"))
    all_for_forest = [{
        "cohort": "TCGA-LIHC (training)", "HR": train_perf["HR"].values[0],
        "HR_lower": train_perf["HR_ci_lo"].values[0], "HR_upper": train_perf["HR_ci_hi"].values[0],
        "c_index": train_perf["c_index"].values[0],
    }] + all_results

    fig, ax = plt.subplots(figsize=(8, max(3, len(all_for_forest) * 0.8 + 1)))
    for i, r in enumerate(all_for_forest):
        color = 'steelblue' if 'training' in r['cohort'] else 'darkorange'
        hr_lo = r.get('HR_lower', r.get('ci_lo', r['HR'] * 0.7))
        hr_hi = r.get('HR_upper', r.get('ci_hi', r['HR'] * 1.3))
        ax.plot([hr_lo, hr_hi], [i, i], color=color, linewidth=2.5)
        ax.scatter(r['HR'], i, color=color, s=80, zorder=5, edgecolors='black')
    ax.axvline(1.0, color='black', linestyle='--', alpha=0.5)
    ax.set_yticks(range(len(all_for_forest)))
    ax.set_yticklabels([f"{r['cohort']} (C={r['c_index']:.3f})" for r in all_for_forest])
    ax.set_xlabel("Hazard Ratio (95% CI)")
    ax.set_title("Cross-Cohort Validation — Forest Plot", fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "forest_plot_all.png"), dpi=200, bbox_inches='tight')
    print("Saved: forest_plot_all.png")

    # ── Fisher's combined p-value ───────────────────────────────────────────
    from scipy import stats as sp_stats
    pvals = [r["logrank_p"] for r in all_results if r["logrank_p"] > 0]
    if len(pvals) >= 2:
        chi2 = -2 * sum(np.log(p) for p in pvals)
        fisher_p = 1 - sp_stats.chi2.cdf(chi2, df=2 * len(pvals))
        print(f"\nFisher's combined p-value: {fisher_p:.2e} (from {len(pvals)} cohorts)")
else:
    print("No validation results obtained.")

# Save GSE14520 expression data for signature comparison
if gse14520_expr_genes is not None and gse14520_matched is not None:
    gse14520_cache_path = os.path.join(GEO, "gse14520_gene_expr_cache.csv")
    if not os.path.exists(gse14520_cache_path):
        gse14520_expr_genes.to_csv(gse14520_cache_path)
        print("Saved: gse14520_gene_expr_cache.csv (for signature comparison)")

print(f"\n✓ External validation complete.")
