# An 11-Gene ROS/Ferroptosis Prognostic Signature for Hepatocellular Carcinoma

Hepatocellular carcinoma (HCC) is the most common form of primary liver cancer and one of the leading causes of cancer-related deaths worldwide. Outcomes vary enormously between patients — some survive years after diagnosis, others deteriorate within months — and current clinical staging systems do a limited job of predicting who falls where. This project asks whether the molecular biology of oxidative stress and ferroptosis (a form of iron-dependent cell death driven by lipid peroxidation) can help us do better.

Reactive oxygen species (ROS) are a double-edged sword in cancer. At moderate levels, they promote tumor growth and survival signaling. At high levels, they trigger ferroptosis and kill cancer cells. Many HCC tumors rewire their antioxidant defenses — upregulating genes like thioredoxin reductase (TXNRD1), glutathione reductase (GSR), and the cystine transporter SLC7A11 — to keep ROS in a "Goldilocks zone" that supports proliferation without tipping into cell death. We hypothesized that the expression pattern of ROS/ferroptosis genes could serve as a molecular fingerprint that predicts patient survival.

## What We Did

### Data and Gene Selection

We started with 302 HCC patients from TCGA-LIHC (The Cancer Genome Atlas) and a curated set of ROS- and ferroptosis-related genes drawn from MSigDB, KEGG, and published literature. Using univariate Cox regression, we identified which of these genes had a statistically significant association with overall survival on their own.

### Building the Signature

We then used LASSO-penalized Cox regression — a method that simultaneously selects the most informative genes and estimates their prognostic weights — to distill the candidates down to an **11-gene signature**:

| Gene | Role | Direction |
|------|------|-----------|
| **TXNRD1** | Thioredoxin reductase; regenerates thioredoxin | Risk |
| **MAFG** | Transcription factor; partner of NRF2 in antioxidant response | Risk |
| **G6PD** | Glucose-6-phosphate dehydrogenase; NADPH production | Risk |
| **SQSTM1** | Autophagy receptor (p62); NRF2 pathway activator | Risk |
| **SLC7A11** | Cystine/glutamate antiporter; ferroptosis suppressor | Risk |
| **GSR** | Glutathione reductase; maintains reduced glutathione pool | Risk |
| **NCF2** | NADPH oxidase subunit; superoxide production | Risk |
| **HMOX1** | Heme oxygenase-1; heme catabolism, iron release | Risk |
| **GLRX2** | Glutaredoxin-2; mitochondrial redox homeostasis | Risk |
| **BACH1** | Transcriptional repressor of HMOX1 and antioxidant genes | Risk |
| **MSRA** | Methionine sulfoxide reductase; repairs oxidized proteins | Protective |

Ten of the eleven genes are upregulated in high-risk patients, reflecting a tumor that is actively ramping up antioxidant defenses. The one protective gene, MSRA, is downregulated in aggressive tumors — its loss removes a layer of oxidative damage repair.

Each patient receives a risk score: a weighted sum of these 11 gene expression values. Patients above the median score are classified as "high-risk."

### Validation

We validated the signature in three independent cohorts that were not used during model training:

| Cohort | N | Events | C-index | HR | Log-rank p |
|--------|---|--------|---------|-----|------------|
| **TCGA-LIHC** (training) | 302 | 129 | 0.700 | 3.40 | 3.9 x 10⁻⁶ |
| **GSE14520** (Chinese HBV cohort) | 221 | 85 | 0.596 | 1.79 | 0.004 |
| **ICGC LIRI-JP** (Japanese cohort) | 231 | 43 | 0.662 | 1.78 | 0.002 |
| **GSE76427** (European cohort, OS) | 115 | 23 | — | — | NS |
| **GSE76427** (RFS) | 108 | 48 | 0.571 | 1.46 | 0.125 |

**What these numbers mean:**
- **C-index** (concordance index): How well the model ranks patients — 0.5 is random guessing, 1.0 is perfect. Our values of 0.60–0.70 indicate moderate-to-good discrimination, which is typical for gene expression signatures in cancer.
- **HR** (hazard ratio): How much higher the death rate is in the high-risk group versus the low-risk group. An HR of 1.79 means high-risk patients face roughly 80% greater hazard.
- **Log-rank p**: Whether the survival difference between high- and low-risk groups is statistically significant. Below 0.05 is considered significant.
- GSE76427 OS was underpowered with only 23 death events, making it difficult to detect a real effect. When we looked at recurrence-free survival instead (48 events), the trend went in the expected direction.

### Multivariate Analysis

A prognostic signature is only useful if it tells you something beyond what you already know from standard clinical variables (age, sex, tumor stage). In multivariate Cox regression:

- **ICGC LIRI-JP**: Risk score remained independently significant (HR=1.79, p=0.004) after adjusting for age, sex, and stage.
- **GSE14520**: Risk score was borderline (HR=1.57, p=0.056), with TNM stage dominating (p=1.5 x 10⁻⁷) — not surprising given that staging captures a lot of prognostic information in this well-characterized cohort.

### Nomogram and Decision Curve Analysis

We combined the risk score with clinical variables (age, sex, tumor stage) into a points-based **nomogram** that predicts 1-, 3-, and 5-year survival probabilities for individual patients.

- The full nomogram achieved a C-index of **0.726** (95% CI: 0.673–0.774), significantly better than the risk score alone (C=0.703, p=0.019) and substantially better than stage alone (C=0.598).
- **Decision curve analysis** confirmed that the nomogram provides positive net benefit across a wide range of clinically relevant risk thresholds, meaning it would improve clinical decision-making compared to treating everyone or treating no one.

### Comparison With Existing Signatures

We compared our signature against three published prognostic signatures on the TCGA training set:

| Signature | Genes | C-index |
|-----------|-------|---------|
| Our ROS/Ferroptosis | 11 | 0.700 |
| Hong 8-gene OS | 8 | 0.702 |
| Buffa Hypoxia | 15 | 0.610 |
| MKI67 Proliferation | 4 | 0.637 |

Our signature performs comparably to the best existing signature (Hong et al.) while being rooted in a specific and therapeutically targetable biology — oxidative stress and ferroptosis — rather than generic prognostic genes.

## Biological Characterization

### Pathway Enrichment

The signature genes are strongly enriched in exactly the pathways we expected:

- **KEAP1-NRF2 pathway** (9/11 genes, Reactome, p = 1.2 x 10⁻¹⁶) — the master regulator of cellular antioxidant defense
- **Reactive Oxygen Species Pathway** (9 genes, MSigDB Hallmark, p = 1.3 x 10⁻¹⁹)
- **Ferroptosis** (3 genes, KEGG, p = 7.6 x 10⁻⁶)
- **Glutathione metabolism** (4 genes, KEGG, p = 2.2 x 10⁻⁷)

This confirms the signature is capturing genuine ROS/ferroptosis biology, not arbitrary prognostic noise.

### Immune Microenvironment

High-risk scores correlated with:
- **Decreased NK cell infiltration** (r = −0.17, p = 0.003) — suggesting immune evasion
- **Increased M1 macrophages** (r = 0.12, p = 0.044) — consistent with a pro-inflammatory, oxidatively stressed microenvironment

### Drug Sensitivity

Using cancer cell line pharmacogenomic data:
- **Erastin** (ferroptosis inducer via SLC7A11 inhibition): Strong positive correlation with risk score (r = 0.64, p < 10⁻³⁵), suggesting high-risk tumors may be particularly sensitive to ferroptosis-inducing therapy.
- **Doxorubicin**: Moderate positive correlation (r = 0.36), consistent with ROS-mediated mechanisms of action.
- **Sorafenib** (standard-of-care for advanced HCC): No significant correlation — the signature captures biology orthogonal to sorafenib targets.

### Single-Cell Resolution

Using Human Protein Atlas single-cell RNA-seq data from liver tissue, we determined the cellular origin of each signature gene:
- **5 genes** are primarily expressed in **hepatocytes** (TXNRD1, SQSTM1, GSR, MSRA, GLRX2) — tumor-intrinsic
- **3 genes** are primarily expressed in **immune cells** (G6PD in T cells, NCF2 and HMOX1 in Kupffer cells) — reflecting the tumor microenvironment
- **2 genes** are primarily expressed in **stromal cells** (MAFG, SLC7A11)

This mixed origin suggests the signature captures both tumor-autonomous redox rewiring and the oxidative stress response of the surrounding microenvironment.

### Protein-Level Validation

Using immunohistochemistry data from the Human Protein Atlas, 4 out of 10 evaluable genes showed concordance between mRNA upregulation and protein detection in HCC tissue (TXNRD1, MAFG, GSR, MSRA). The remaining discordances are expected — mRNA-protein correlation is notoriously imperfect, especially for enzymes with high catalytic turnover.

### Molecular Subtypes

Consensus clustering on the 11 signature genes identified two molecular subtypes:
- **Cluster 1** (n=41): High-risk, 70.7% mortality — an oxidatively stressed subtype
- **Cluster 2** (n=261): Low-risk, 38.3% mortality — a redox-balanced subtype

High-risk patients also had significantly higher tumor mutation burden (median 76 vs 66 mutations, p=0.012), consistent with the mutagenic effects of oxidative DNA damage.

## Repository Structure

```
hcc-ros-signature/
├── scripts/
│   ├── 01_data_setup.py            # Download TCGA-LIHC data, curate ROS/ferroptosis gene set
│   ├── 02_univariate_screening.py  # Univariate Cox regression for each gene
│   ├── 03_lasso_model.py           # LASSO-Cox model selection with cross-validation
│   ├── 04_training_evaluation.py   # Training set performance (KM, ROC, multivariate)
│   ├── 05_external_validation.py   # GSE14520, ICGC LIRI-JP, GSE76427 validation
│   ├── 06_signature_comparison.py  # Compare against Hong, Buffa, MKI67 signatures
│   ├── 07_biological_analyses.py   # Enrichment, immune, drug sensitivity
│   ├── 08_address_limitations.py   # Subgroup analyses, Schoenfeld residuals, demographics
│   ├── 09_protein_validation.py    # HPA immunohistochemistry protein-level evidence
│   ├── 10_single_cell.py           # Single-cell RNA-seq cell-type expression
│   ├── 11_nomogram_dca.py          # Nomogram, calibration, and decision curve analysis
│   └── 12_molecular_subtypes.py    # Consensus clustering and mutation landscape
├── data/                           # Downloaded and cached datasets (not tracked in git)
├── results/
│   ├── figures/                    # All publication-ready figures (39 PNGs)
│   └── tables/                     # All result tables (21 CSVs)
└── README.md
```

## How to Reproduce

```bash
# Install dependencies
pip install pandas numpy scipy lifelines scikit-learn matplotlib seaborn requests statsmodels

# Run the full pipeline (scripts are numbered in order)
for i in 01 02 03 04 05 06 07 08 09 10 11 12; do
    python3 -u scripts/${i}_*.py
done
```

Scripts download data automatically from TCGA (via GDC API), GEO, ICGC, Human Protein Atlas, and cBioPortal. An internet connection is required for the first run; subsequent runs use cached data.

## Key Takeaways

1. ROS and ferroptosis gene expression patterns carry real prognostic information in HCC, validated across ethnically and geographically diverse cohorts (American, Chinese, Japanese).
2. The signature is biologically coherent — it captures NRF2-mediated antioxidant rewiring, not statistical artifacts.
3. High-risk patients show features (NK cell depletion, high TMB, erastin sensitivity) that suggest potential therapeutic strategies: ferroptosis induction and immune checkpoint therapy.
4. The nomogram integrating molecular and clinical features outperforms either alone.

## Databases and Data Sources

This project relies on the following publicly available databases and resources:

| Database | Usage | Reference |
|----------|-------|-----------|
| [TCGA-LIHC](https://portal.gdc.cancer.gov/) | Training cohort (n=302): RNA-seq gene expression, clinical outcomes, mutation data | The Cancer Genome Atlas Research Network. *Cell*, 2017. DOI: [10.1016/j.cell.2017.05.046](https://doi.org/10.1016/j.cell.2017.05.046) |
| [GEO GSE14520](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520) | External validation (n=221): Chinese HBV-related HCC cohort, Affymetrix microarray | Roessler S, et al. *Cancer Res*, 2010. DOI: [10.1158/0008-5472.CAN-09-2453](https://doi.org/10.1158/0008-5472.CAN-09-2453) |
| [ICGC LIRI-JP](https://dcc.icgc.org/projects/LIRI-JP) | External validation (n=231): Japanese HCC cohort, RNA-seq | Fujimoto A, et al. *Nat Genet*, 2016. DOI: [10.1038/ng.3547](https://doi.org/10.1038/ng.3547) |
| [GEO GSE76427](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76427) | External validation (n=115): European HCC cohort, Illumina microarray | Grinchuk OV, et al. *Genome Med*, 2017. DOI: [10.1186/s13073-017-0438-y](https://doi.org/10.1186/s13073-017-0438-y) |
| [Human Protein Atlas](https://www.proteinatlas.org/) | Protein-level IHC validation, single-cell RNA-seq cell-type expression (v23) | Uhlen M, et al. *Science*, 2015. DOI: [10.1126/science.1260419](https://doi.org/10.1126/science.1260419) |
| [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) | Hallmark gene sets, GO terms for ROS/ferroptosis gene curation | Liberzon A, et al. *Cell Syst*, 2015. DOI: [10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004) |
| [KEGG](https://www.genome.jp/kegg/) | Pathway definitions (ferroptosis, glutathione metabolism) | Kanehisa M, Goto S. *Nucleic Acids Res*, 2000. DOI: [10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27) |
| [Reactome](https://reactome.org/) | Pathway enrichment analysis (KEAP1-NRF2 pathway) | Jassal B, et al. *Nucleic Acids Res*, 2020. DOI: [10.1093/nar/gkz1031](https://doi.org/10.1093/nar/gkz1031) |
| [Enrichr](https://maayanlab.cloud/Enrichr/) | Gene set enrichment analysis engine | Kuleshov MV, et al. *Nucleic Acids Res*, 2016. DOI: [10.1093/nar/gkw377](https://doi.org/10.1093/nar/gkw377) |
| [cBioPortal](https://www.cbioportal.org/) | Somatic mutation and copy number alteration data for TCGA-LIHC | Cerami E, et al. *Cancer Discov*, 2012. DOI: [10.1158/2159-8290.CD-12-0095](https://doi.org/10.1158/2159-8290.CD-12-0095) |
| [DepMap / GDSC2](https://depmap.org/portal/) | Cancer cell line drug sensitivity (IC50) data | Yang W, et al. *Nucleic Acids Res*, 2013. DOI: [10.1093/nar/gks1111](https://doi.org/10.1093/nar/gks1111) |

### Methods References

| Method | Library/Tool | Reference |
|--------|-------------|-----------|
| LASSO-Cox regression | [lifelines](https://lifelines.readthedocs.io/) | Davidson-Pilon C. *JOSS*, 2019. DOI: [10.21105/joss.01317](https://doi.org/10.21105/joss.01317) |
| Immune cell deconvolution | CIBERSORT (LM22 signatures) | Newman AM, et al. *Nat Methods*, 2015. DOI: [10.1038/nmeth.3337](https://doi.org/10.1038/nmeth.3337) |
| Decision curve analysis | Custom implementation | Vickers AJ, Elkin EB. *Med Decis Making*, 2006. DOI: [10.1177/0272989X06295361](https://doi.org/10.1177/0272989X06295361) |
| Consensus clustering | scikit-learn KMeans + AgglomerativeClustering | Monti S, et al. *Mach Learn*, 2003. DOI: [10.1023/A:1023949509487](https://doi.org/10.1023/A:1023949509487) |

### Competing Signatures Compared

| Signature | Reference |
|-----------|-----------|
| Hong 8-gene OS signature | Hong W, et al. *Ann Transl Med*, 2020. DOI: [10.21037/atm-20-3585](https://doi.org/10.21037/atm-20-3585) |
| Buffa hypoxia 15-gene signature | Buffa FM, et al. *Br J Cancer*, 2010. DOI: [10.1038/sj.bjc.6605612](https://doi.org/10.1038/sj.bjc.6605612) |
| MKI67 proliferation 4-gene signature | Whitfield ML, et al. *Mol Biol Cell*, 2002. DOI: [10.1091/mbc.01-12-0567](https://doi.org/10.1091/mbc.01-12-0567) |

## License

This project is for academic research purposes.
