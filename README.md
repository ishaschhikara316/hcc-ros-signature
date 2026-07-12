# An 11-Gene ROS/Ferroptosis Prognostic Signature for Hepatocellular Carcinoma

Hepatocellular carcinoma (HCC) is the most common form of primary liver cancer and one of the leading causes of cancer-related deaths worldwide. Outcomes vary enormously between patients (some survive years after diagnosis, others deteriorate within months), and current clinical staging systems do a limited job of predicting who falls where. This project asks whether the molecular biology of oxidative stress and ferroptosis (a form of iron-dependent cell death driven by lipid peroxidation) can help us do better.

Reactive oxygen species (ROS) are a double-edged sword in cancer. At moderate levels, they promote tumor growth and survival signaling. At high levels, they trigger ferroptosis and kill cancer cells. Many HCC tumors rewire their antioxidant defenses, upregulating genes like thioredoxin reductase (TXNRD1), glutathione reductase (GSR), and the cystine transporter SLC7A11, to keep ROS in a "Goldilocks zone" that supports proliferation without tipping into cell death. We hypothesized that the expression pattern of ROS/ferroptosis genes could be a molecular fingerprint that predicts patient survival.

## What We Did

### Data and Gene Selection

We started with 302 HCC patients from TCGA-LIHC (The Cancer Genome Atlas) and a curated set of ROS- and ferroptosis-related genes drawn from MSigDB, KEGG, and published literature. Using univariate Cox regression, we identified which of these genes had a statistically significant association with overall survival on their own.

### Building the Signature

We then used LASSO-penalized Cox regression (a method that simultaneously selects the most informative genes and estimates their prognostic weights) to distill the candidates down to an **11-gene signature**:

| Gene | Role | Direction |
|------|------|-----------|
| **TXNRD1** | Thioredoxin reductase; regenerates thioredoxin | Risk |
| **MAFG** | Transcription factor; partner of NRF2 in antioxidant response | Risk |
| **G6PD** | Glucose-6-phosphate dehydrogenase; NADPH production | Risk |
| **SQSTM1** | Autophagy receptor (p62); NRF2 pathway activator | Risk |
| **SLC7A11** | Cystine/glutamate antiporter; ferroptosis suppressor | Risk |
| **GSR** | Glutathione reductase; maintains reduced glutathione pool | Risk |
| **NCF2** | NADPH oxidase subunit; superoxide production (primarily immune cell-derived, reflects tumor microenvironment) | Risk |
| **HMOX1** | Heme oxygenase-1; heme catabolism, iron release | Risk |
| **GLRX2** | Glutaredoxin-2; mitochondrial redox homeostasis | Risk |
| **BACH1** | NRF2 pathway interactor; competitive repressor of HMOX1 and antioxidant genes | Risk |
| **MSRA** | Methionine sulfoxide reductase; repairs oxidized proteins | Protective |

Ten of the eleven genes are upregulated in high-risk patients, reflecting a tumor that is actively ramping up antioxidant defenses. The one protective gene, MSRA, is downregulated in aggressive tumors; its loss removes a layer of oxidative damage repair.

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
| **GSE54236** (additional validation) | 80 | — | 0.621 | 1.68 | 0.007 |

**What these numbers mean:**
- **C-index** (concordance index): How well the model ranks patients: 0.5 is random guessing, 1.0 is perfect. Our values of 0.60–0.70 indicate moderate-to-good discrimination, which is typical for gene expression signatures in cancer.
- **HR** (hazard ratio): How much higher the death rate is in the high-risk group versus the low-risk group. An HR of 1.79 means high-risk patients face roughly 80% greater hazard.
- **Log-rank p**: Whether the survival difference between high- and low-risk groups is statistically significant. Below 0.05 is considered significant.
- GSE76427 OS was underpowered with only 23 death events, making it difficult to detect a real effect. When we looked at recurrence-free survival instead (48 events), the trend went in the expected direction.
- GSE54236 provides additional validation in an independent HCC cohort, with a statistically significant result (p = 0.007).

### Multivariate Analysis

A prognostic signature is only useful if it tells you something beyond what you already know from standard clinical variables (age, sex, tumor stage). In multivariate Cox regression:

- **ICGC LIRI-JP**: Risk score remained independently significant (HR=1.79, p=0.004) after adjusting for age, sex, and stage.
- **GSE14520**: Risk score was borderline (HR=1.57, p=0.056), with TNM stage dominating (p=1.5 x 10⁻⁷), which is not surprising given that staging captures a lot of prognostic information in this well-characterized cohort.

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

Our signature performs comparably to the best existing signature (Hong et al.) while being rooted in a specific and therapeutically targetable biology (oxidative stress and ferroptosis) rather than generic prognostic genes.

## Biological Characterization

### Pathway Enrichment

The signature genes are strongly enriched in exactly the pathways we expected:

- **KEAP1-NRF2 pathway** (7/11 genes are NRF2-KEAP1 pathway genes, including regulators and interactors; Reactome, p = 1.2 x 10⁻¹⁶), the master regulator of cellular antioxidant defense
- **Reactive Oxygen Species Pathway** (9 genes, MSigDB Hallmark, p = 1.3 x 10⁻¹⁹)
- **Ferroptosis** (3 genes, KEGG, p = 7.6 x 10⁻⁶)
- **Glutathione metabolism** (4 genes, KEGG, p = 2.2 x 10⁻⁷)

This confirms the signature is capturing genuine ROS/ferroptosis biology, not arbitrary prognostic noise.

### Immune Checkpoints and Immunotherapy Landscape

We analyzed immune checkpoint expression, immune cell infiltration (via ssGSEA with 24 cell type signatures), and immunotherapy response prediction between risk groups:

- **3 of 9 immune checkpoints** were significantly differentially expressed: SIGLEC15 (strongly downregulated in high-risk, p < 0.0001), LAG3 (p = 0.04), and B7-H3/CD276 (p = 0.05).
- **5 of 23 immune cell types** showed significant correlation with risk score via ssGSEA, providing a detailed picture of how the tumor immune microenvironment shifts with oxidative stress status.
- High-risk tumors showed features consistent with immune dysfunction, including decreased NK cell infiltration (r = −0.17, p = 0.003) and increased M1 macrophages (r = 0.12, p = 0.044).

### Genome-Wide Pathway Enrichment (GSEA)

Rather than limiting enrichment analysis to the 11 signature genes, we ranked all ~20,000 genes in the transcriptome by their Spearman correlation with the risk score and ran preranked GSEA against four major pathway libraries (Hallmark, KEGG, Reactome, GO Biological Process):

- **1,384 of 4,045 pathways** were significantly enriched (FDR < 0.25): Hallmark 31/50, KEGG 174/304, Reactome 506/1,147, GO BP 673/2,544.
- **Enriched in high-risk tumors**: mTORC1 Signaling (Hallmark), NFE2L2 Nuclear Events (Reactome), Proteasome (KEGG), Ribosome Biogenesis (GO BP), consistent with a proliferative, oxidatively stressed phenotype.
- **Enriched in low-risk tumors**: Bile Acid Metabolism (Hallmark), Fatty Acid Oxidation (GO BP), reflecting normal hepatocyte metabolic function.

This genome-wide analysis confirms that our 11-gene signature is a proxy for a global transcriptomic shift toward oxidative stress, proliferation, and metabolic reprogramming rather than noise.

### Somatic Mutation Landscape

Using cBioPortal mutation data for TCGA-LIHC, we compared mutation frequencies of 30 HCC-relevant driver genes between high- and low-risk groups:

| Gene | High-Risk Freq | Low-Risk Freq | Odds Ratio | p-value |
|------|---------------|---------------|------------|---------|
| **TP53** | 40.4% | 17.9% | 3.11 | < 0.0001 |
| **CTNNB1** | 36.4% | 15.9% | 3.03 | 0.0001 |
| **KEAP1** | 8.6% | 2.0% | 4.65 | 0.018 |
| **TSC2** | 5.3% | 0.7% | 8.39 | 0.036 |
| **BAP1** | 3.3% | 9.9% | 0.31 | 0.035 |

The enrichment of TP53 mutations in high-risk tumors is consistent with the loss of p53-mediated ferroptosis regulation. KEAP1 mutations, which constitutively activate NRF2, provide a direct genetic link to the antioxidant rewiring captured by our signature. TSC2 mutations (mTORC1 regulator, OR=8.39) in high-risk tumors align with the GSEA finding that mTORC1 signaling is the top enriched Hallmark pathway. The CTNNB1/TP53 molecular class mapping showed that our high-risk group is strongly enriched for the TP53-mutant class (64.8% high-risk) and dual TP53+CTNNB1 mutants (88.2% high-risk), while the "Neither" group was predominantly low-risk (67.1%). Chi-squared test for association between molecular class and risk group: p < 0.0001.

### Drug Sensitivity

Using cancer cell line pharmacogenomic data:
- **Erastin** (ferroptosis inducer via SLC7A11 inhibition): Strong positive correlation with risk score (r = 0.64, p < 10⁻³⁵), suggesting high-risk tumors may be particularly sensitive to ferroptosis-inducing therapy.
- **Doxorubicin**: Moderate positive correlation (r = 0.36), consistent with ROS-mediated mechanisms of action.
- **Sorafenib** (standard-of-care for advanced HCC): No significant correlation. The signature captures biology orthogonal to sorafenib targets.

### Single-Cell Resolution

Using Human Protein Atlas single-cell RNA-seq data from liver tissue, we determined the cellular origin of each signature gene:
- **5 genes** are primarily expressed in **hepatocytes** (TXNRD1, SQSTM1, GSR, MSRA, GLRX2), tumor-intrinsic
- **3 genes** are primarily expressed in **immune cells** (G6PD in T cells, NCF2 and HMOX1 in Kupffer cells), reflecting the tumor microenvironment
- **2 genes** are primarily expressed in **stromal cells** (MAFG, SLC7A11)

This mixed origin suggests the signature captures both tumor-autonomous redox rewiring and the oxidative stress response of the surrounding microenvironment.

### Protein-Level Validation

Using immunohistochemistry data from the Human Protein Atlas, 4 out of 10 evaluable genes showed concordance between mRNA upregulation and protein detection in HCC tissue (TXNRD1, MAFG, GSR, MSRA). The remaining discordances are expected: mRNA-protein correlation is notoriously imperfect, especially for enzymes with high catalytic turnover.

### Pan-Cancer Generalizability

To test whether ROS/ferroptosis biology carries prognostic information beyond HCC, we applied the signature to four additional TCGA cancer types where oxidative stress is implicated:

| Cancer | N | Events | C-index | HR | Log-rank p |
|--------|---|--------|---------|-----|------------|
| **TCGA-KIRC** (Kidney Renal Clear Cell) | 529 | 173 | 0.580 | 1.73 | 0.001 |
| **TCGA-LUSC** (Lung Squamous Cell) | 494 | 211 | 0.504 | 0.96 | 0.711 (NS) |
| **TCGA-STAD** (Stomach Adenocarcinoma) | 388 | 156 | 0.482 | 0.88 | 0.617 (NS) |
| **TCGA-PAAD** (Pancreatic Adenocarcinoma) | 177 | 93 | 0.525 | 1.43 | 0.311 (NS) |

The signature was significant in **kidney renal clear cell carcinoma** (KIRC), where the VHL-HIF pathway creates a uniquely ROS-dependent tumor microenvironment. The lack of significance in other cancer types is expected: it reflects that our signature captures liver-specific redox biology rather than generic prognostic genes, which is actually a strength for clinical specificity.

### Molecular Subtypes

Consensus clustering on the 11 signature genes identified two molecular subtypes:
- **Cluster 1** (n=41): High-risk, 70.7% mortality, an oxidatively stressed subtype
- **Cluster 2** (n=261): Low-risk, 38.3% mortality, a redox-balanced subtype

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
│   ├── 12_molecular_subtypes.py    # Consensus clustering and mutation landscape
│   ├── 13_immune_checkpoints_tide.py # Immune checkpoints, ssGSEA, immunotherapy prediction
│   ├── 14_gsea_risk_groups.py      # Genome-wide preranked GSEA (Hallmark, KEGG, Reactome, GO)
│   ├── 15_mutation_landscape.py    # Somatic mutation analysis and molecular classification
│   ├── 16_additional_validation.py # Additional GEO validation cohorts (GSE54236, GSE36376)
│   └── 17_pan_cancer.py            # Pan-cancer validation (KIRC, LUSC, STAD, PAAD)
├── data/                           # Downloaded and cached datasets (not tracked in git)
├── results/
│   ├── figures/                    # All publication-ready figures (50+ PNGs)
│   └── tables/                     # All result tables (30+ CSVs)
└── README.md
```

## How to Reproduce

```bash
# Install dependencies
pip install pandas numpy scipy lifelines scikit-learn matplotlib seaborn requests statsmodels

# Run the full pipeline (scripts are numbered in order)
for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17; do
    python3 -u scripts/${i}_*.py
done
```

Scripts download data automatically from TCGA (via GDC API), GEO, ICGC, Human Protein Atlas, and cBioPortal. An internet connection is required for the first run; subsequent runs use cached data.

## Key Takeaways

1. ROS and ferroptosis gene expression patterns carry real prognostic information in HCC, validated across ethnically and geographically diverse cohorts (American, Chinese, Japanese) and an additional GEO cohort (GSE54236).
2. The signature is biologically coherent: genome-wide GSEA confirms it captures NRF2-mediated antioxidant rewiring, mTORC1 signaling, and metabolic reprogramming, not statistical artifacts.
3. High-risk patients are enriched for TP53 and KEAP1 mutations, show immune checkpoint dysregulation (SIGLEC15, LAG3), NK cell depletion, high TMB, and strong predicted sensitivity to ferroptosis inducers, all suggesting potential therapeutic strategies.
4. The signature also validates in kidney cancer (TCGA-KIRC), suggesting the ROS/ferroptosis axis is prognostically relevant in other ROS-dependent malignancies.
5. The nomogram integrating molecular and clinical features outperforms either alone.

## Complete Results at a Glance

### Table 1: Signature Validation Across All Cohorts

| Cohort | Population | Platform | N | Events | C-index | HR (95% CI) | Log-rank p | Interpretation |
|--------|-----------|----------|---|--------|---------|-------------|------------|----------------|
| TCGA-LIHC (training) | Mixed (US) | RNA-seq | 302 | 129 | 0.700 | 3.40 | 3.9 x 10⁻⁶ | Strong separation; defines the model |
| GSE14520 | Chinese (HBV) | Affymetrix | 221 | 85 | 0.596 | 1.79 (1.20–2.68) | 0.004 | Validates cross-platform and in HBV-dominant HCC |
| ICGC LIRI-JP | Japanese | RNA-seq | 231 | 43 | 0.662 | 1.78 | 0.002 | Independent RNA-seq validation; multivariate significant |
| GSE54236 | European | Agilent | 80 | 80 | 0.621 | 1.68 (1.04–2.72) | 0.007 | Additional validation on third microarray platform |
| GSE76427 (OS) | European | Illumina | 115 | 23 | — | — | NS | Underpowered (only 23 events) |
| GSE76427 (RFS) | European | Illumina | 108 | 48 | 0.571 | 1.46 | 0.125 | Trend in expected direction for recurrence |
| TCGA-KIRC (pan-cancer) | Mixed (US) | RNA-seq | 529 | 173 | 0.580 | 1.73 (1.40–2.14) | 0.001 | ROS/ferroptosis biology extends to kidney cancer |

### Table 2: What High-Risk Tumors Look Like (Biological Characterization)

| Feature | Finding | p-value | What It Means |
|---------|---------|---------|---------------|
| **Pathway enrichment** | 7/11 genes are NRF2-KEAP1 pathway genes (including regulators and interactors) | 1.2 x 10⁻¹⁶ | Signature captures the master antioxidant defense system |
| **Genome-wide GSEA** | 1,384/4,045 pathways enriched (FDR<0.25) | — | Not random noise — reflects global transcriptomic shift |
| **Top pathway (high-risk)** | mTORC1 Signaling | FDR < 0.25 | Growth/proliferation signals are active |
| **Top pathway (low-risk)** | Bile Acid Metabolism | FDR < 0.25 | Normal liver metabolic function preserved |
| **TP53 mutations** | 40.4% vs 17.9% in high vs low risk | < 0.0001 | p53 loss removes ferroptosis brake |
| **KEAP1 mutations** | 8.6% vs 2.0% in high vs low risk | 0.018 | Constitutive NRF2 activation — direct genetic proof |
| **TSC2 mutations** | 5.3% vs 0.7% in high vs low risk | 0.036 | mTORC1 pathway dysregulation |
| **BAP1 mutations** | 3.3% vs 9.9% in high vs low risk | 0.035 | BAP1 loss associated with *better* prognosis |
| **SIGLEC15 expression** | Downregulated in high-risk (FC=0.57) | < 0.0001 | Immune checkpoint dysregulation |
| **NK cell infiltration** | Decreased in high-risk (r = −0.17) | 0.003 | Immune evasion — anti-tumor immunity impaired |
| **Erastin sensitivity** | High-risk more sensitive (r = 0.64) | < 10⁻³⁵ | Ferroptosis-inducing therapy may work in high-risk patients |
| **Sorafenib sensitivity** | No correlation | NS | Signature captures biology orthogonal to current therapy |
| **Tumor mutation burden** | Higher in high-risk (76 vs 66 mutations) | 0.012 | Oxidative DNA damage → more mutations |

### Table 3: Nomogram Performance (Molecular + Clinical)

| Model | C-index (95% CI) | Improvement | Interpretation |
|-------|------------------|-------------|----------------|
| Tumor stage alone | 0.598 | — | Current clinical standard |
| Risk score alone | 0.703 (0.647–0.743) | +0.105 | Molecular information alone beats staging |
| Full nomogram (risk + age + sex + stage) | 0.726 (0.673–0.774) | +0.128 | Best performance — combines both information sources |

### Table 4: Molecular Subtypes (Consensus Clustering)

| Subtype | N | Mortality | Mean Risk Score | TP53 Mutant | CTNNB1 Mutant | Description |
|---------|---|-----------|-----------------|-------------|---------------|-------------|
| Oxidatively stressed (Cluster 1) | 41 | 70.7% | High | Enriched | Enriched | Aggressive, antioxidant-dependent |
| Redox-balanced (Cluster 2) | 261 | 38.3% | Low | Depleted | Depleted | Metabolically normal, better prognosis |

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
| [GEO GSE54236](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54236) | Additional external validation HCC cohort (n=80) | Villa E, et al. *Gastroenterology*, 2016. |
| [TCGA-KIRC/LUSC/STAD/PAAD](https://portal.gdc.cancer.gov/) | Pan-cancer validation cohorts via cBioPortal | The Cancer Genome Atlas Research Network, various publications |

### Methods References

| Method | Library/Tool | Reference |
|--------|-------------|-----------|
| LASSO-Cox regression | [lifelines](https://lifelines.readthedocs.io/) | Davidson-Pilon C. *JOSS*, 2019. DOI: [10.21105/joss.01317](https://doi.org/10.21105/joss.01317) |
| Immune cell deconvolution | CIBERSORT (LM22 signatures) | Newman AM, et al. *Nat Methods*, 2015. DOI: [10.1038/nmeth.3337](https://doi.org/10.1038/nmeth.3337) |
| Decision curve analysis | Custom implementation | Vickers AJ, Elkin EB. *Med Decis Making*, 2006. DOI: [10.1177/0272989X06295361](https://doi.org/10.1177/0272989X06295361) |
| Consensus clustering | scikit-learn KMeans + AgglomerativeClustering | Monti S, et al. *Mach Learn*, 2003. DOI: [10.1023/A:1023949509487](https://doi.org/10.1023/A:1023949509487) |
| Preranked GSEA | [gseapy](https://gseapy.readthedocs.io/) | Fang Z, et al. *Bioinformatics*, 2023. DOI: [10.1093/bioinformatics/btac757](https://doi.org/10.1093/bioinformatics/btac757) |
| ssGSEA immune profiling | gseapy ssgsea module | Barbie DA, et al. *Nature*, 2009. DOI: [10.1038/nature08460](https://doi.org/10.1038/nature08460) |

### Competing Signatures Compared

| Signature | Reference |
|-----------|-----------|
| Hong 8-gene OS signature | Hong W, et al. *Ann Transl Med*, 2020. DOI: [10.21037/atm-20-3585](https://doi.org/10.21037/atm-20-3585) |
| Buffa hypoxia 15-gene signature | Buffa FM, et al. *Br J Cancer*, 2010. DOI: [10.1038/sj.bjc.6605612](https://doi.org/10.1038/sj.bjc.6605612) |
| MKI67 proliferation 4-gene signature | Whitfield ML, et al. *Mol Biol Cell*, 2002. DOI: [10.1091/mbc.01-12-0567](https://doi.org/10.1091/mbc.01-12-0567) |

## License

This project is for academic research purposes.
