# Research Paper Proposal

## Redox Tightrope: An 11-Gene ROS/Ferroptosis Signature Predicts Survival and Ferroptosis Sensitivity in Hepatocellular Carcinoma

---

### Background

Hepatocellular carcinoma is the sixth most frequently diagnosed cancer in the world and the third leading cause of cancer-related death. What makes it particularly difficult to manage is how unpredictably it behaves. Two patients who look similar on paper — same tumor size, same stage, same liver function — can have very different outcomes. One may live five years or more; the other may succumb within months. The clinical staging systems we rely on today, such as BCLC and TNM, capture anatomy and liver function reasonably well, but they miss the molecular heterogeneity that actually drives these divergent trajectories.

Over the past decade, there has been growing interest in using gene expression patterns to fill this gap. The idea is straightforward: if a tumor's molecular state carries information about how aggressively it will behave, then measuring the right set of genes could help us identify which patients need more intensive treatment and which can be spared unnecessary toxicity. Several gene signatures for HCC have been proposed, but most are built from generic prognostic genes — genes that happen to correlate with survival without a clear connection to a specific biological process. This limits both their interpretability and their potential to guide therapy.

We take a different approach. Rather than casting a wide net across the entire transcriptome, we focus on a specific and well-characterized area of cancer biology: the interplay between reactive oxygen species (ROS) and ferroptosis.

### Why ROS and Ferroptosis?

The connection between oxidative stress and liver cancer is not incidental — it is central to how the disease develops and progresses. The liver is the body's primary metabolic organ. It processes toxins, drugs, and dietary compounds, all of which generate ROS as byproducts. Chronic liver diseases — hepatitis B and C, alcohol-related liver disease, metabolic dysfunction-associated steatohepatitis — produce sustained oxidative stress that damages DNA, promotes inflammation, and eventually drives malignant transformation. By the time HCC emerges, the tumor has been marinating in ROS for years.

But here is the paradox. While ROS contributed to the cancer's origin, the cancer now needs to keep ROS under control to survive. Too much oxidative stress triggers ferroptosis — a recently discovered form of regulated cell death driven by iron-dependent lipid peroxidation. Ferroptosis is distinct from apoptosis; it does not require caspases, and it is governed by its own set of molecular regulators, chief among them the cystine transporter SLC7A11 and the lipid repair enzyme GPX4. Tumors that cannot manage their ROS levels die by ferroptosis. Tumors that survive have learned to walk a tightrope: they tolerate enough ROS to fuel growth signaling while upregulating antioxidant defenses — the NRF2-KEAP1 pathway, the thioredoxin system, glutathione metabolism — to prevent ferroptotic collapse.

This redox balancing act is not a peripheral detail of tumor biology. It is a survival strategy, and it varies from patient to patient. We hypothesized that the pattern of ROS and ferroptosis gene expression in a given tumor reflects where it sits on this spectrum — and that this pattern, quantified as a gene signature, could predict clinical outcome.

There is also a therapeutic reason to care about this specific biology. Ferroptosis inducers like erastin and RSL3 are under active preclinical investigation. If we can identify which patients have tumors that are most dependent on their antioxidant defenses — the high-risk patients in our model — these may be exactly the patients who would benefit most from ferroptosis-targeted therapy. A signature rooted in generic survival genes cannot make this connection. A signature rooted in ROS/ferroptosis biology can.

### Aims

**Aim 1: Construct a parsimonious ROS/ferroptosis gene signature for HCC prognosis.**

We will curate a comprehensive set of ROS- and ferroptosis-related genes from established pathway databases (MSigDB Hallmark ROS pathway, KEGG ferroptosis and glutathione metabolism, Reactome KEAP1-NRF2 signaling) and published ferroptosis literature. Using TCGA-LIHC RNA-seq data from 302 patients as the training cohort, we will first screen genes by univariate Cox regression to retain those with a significant individual association with overall survival. We will then apply LASSO-penalized Cox regression with 10-fold cross-validation to select the smallest subset of genes that jointly predicts survival. The target is a model with 8 to 12 genes — large enough to capture the biology, small enough to be clinically translatable.

**Aim 2: Validate the signature across independent, geographically diverse cohorts.**

A model that only works on the data it was trained on is not useful. We will evaluate the signature in three external cohorts spanning different populations, etiologies, and profiling platforms:

- **GSE14520** (n=221): A Chinese cohort enriched for hepatitis B-related HCC, profiled on Affymetrix microarrays. This tests cross-platform generalizability and performance in an HBV-dominant population.
- **ICGC LIRI-JP** (n=231): A Japanese cohort profiled by RNA-seq through the International Cancer Genome Consortium. This tests generalizability in a second Asian population with mixed viral etiologies.
- **GSE76427** (n=115): A European cohort profiled on Illumina microarrays. This tests performance across a different ethnic background and viral landscape.

For each cohort, we will report Kaplan-Meier survival curves, log-rank tests, hazard ratios, concordance indices, and time-dependent ROC curves. We will also perform multivariate Cox regression to determine whether the risk score provides prognostic information independent of age, sex, and tumor stage.

**Aim 3: Build a clinically oriented nomogram integrating molecular and clinical features.**

Gene expression alone does not replace clinical judgement — it supplements it. We will combine the ROS/ferroptosis risk score with standard clinical variables (age, sex, tumor stage) into a points-based nomogram that predicts individualized 1-, 3-, and 5-year survival probabilities. We will evaluate the nomogram with calibration plots, formal C-index comparison against each component alone, and decision curve analysis to determine whether it provides genuine clinical net benefit across a range of decision thresholds.

**Aim 4: Characterize the biological and clinical landscape of the signature.**

To move beyond statistical association and toward mechanistic understanding, we will:

- Perform pathway enrichment analysis to confirm the signature captures ROS/ferroptosis biology rather than generic proliferation signals.
- Estimate immune cell composition using established deconvolution approaches to determine how the tumor immune microenvironment differs between risk groups.
- Correlate risk scores with drug sensitivity data from cancer cell line pharmacogenomic databases, with particular attention to ferroptosis inducers (erastin, RSL3), standard-of-care agents (sorafenib, lenvatinib), and cytotoxic chemotherapy (doxorubicin).
- Use single-cell RNA-seq reference data from the Human Protein Atlas to determine the cellular origin (hepatocyte, immune, stromal) of each signature gene.
- Validate protein-level expression using immunohistochemistry data to assess mRNA-protein concordance.
- Apply consensus clustering on the signature genes to identify molecular subtypes and examine their clinical and genomic characteristics, including tumor mutation burden.

### Preliminary Results

We have completed the full analysis pipeline and the results support the proposed study.

**The signature.** LASSO-Cox regression selected 11 genes: TXNRD1, MAFG, G6PD, SQSTM1, SLC7A11, GSR, NCF2, MSRA, HMOX1, GLRX2, and BACH1. Ten are risk genes (higher expression, worse outcome); one, MSRA, is protective. Nine of the eleven are direct transcriptional targets of the NRF2-KEAP1 pathway, the master regulator of cellular antioxidant defense. The model achieved a concordance index of 0.700 (95% CI: 0.647–0.743) in the training cohort, with a hazard ratio of 3.40 between high- and low-risk groups (p = 3.9 x 10⁻⁶).

**External validation.** The signature validated in GSE14520 (C-index = 0.596, HR = 1.79, p = 0.004) and ICGC LIRI-JP (C-index = 0.662, HR = 1.78, p = 0.002). In ICGC, the risk score remained independently significant in multivariate analysis after adjusting for age, sex, and stage (HR = 1.79, p = 0.004). GSE76427 overall survival was underpowered with only 23 events; recurrence-free survival (48 events) showed a trend in the expected direction (HR = 1.46).

**Nomogram.** The full nomogram (risk score + age + sex + stage) achieved a C-index of 0.726 (95% CI: 0.673–0.774), significantly outperforming the risk score alone (C = 0.703, p = 0.019) and tumor stage alone (C = 0.598). Decision curve analysis confirmed positive net benefit across clinically relevant thresholds at 1, 3, and 5 years.

**Biology.** Pathway enrichment confirmed strong enrichment for the KEAP1-NRF2 pathway (9/11 genes, p = 1.2 x 10⁻¹⁶), reactive oxygen species response (p = 1.3 x 10⁻¹⁹), and ferroptosis (p = 7.6 x 10⁻⁶). High-risk tumors showed decreased NK cell infiltration (r = −0.17, p = 0.003) and strong predicted sensitivity to the ferroptosis inducer erastin (r = 0.64, p < 10⁻³⁵). Single-cell analysis revealed that 5 signature genes are primarily expressed in hepatocytes (tumor-intrinsic), 3 in immune cells (microenvironment), and 2 in stromal cells, indicating the signature captures both tumor-autonomous and microenvironmental contributions to redox status. Consensus clustering identified two molecular subtypes — a high-risk oxidatively stressed subtype (n=41, 70.7% mortality) and a low-risk redox-balanced subtype (n=261, 38.3% mortality) — with the high-risk subtype carrying a significantly higher tumor mutation burden (p = 0.012).

**Subgroup robustness.** The signature was prognostic across early-stage (I–II, C = 0.688, p = 0.001) and late-stage (III–IV, C = 0.711, p < 0.001) disease, across younger (≤60, C = 0.739) and older (>60, C = 0.665) patients, and across both White/European (C = 0.622, p = 0.016) and Asian (C = 0.774, p < 0.001) subpopulations within TCGA.

### Why These Findings Matter

The core question behind this work is practical: can we do a better job of telling HCC patients what to expect, and can that knowledge change what we do for them?

**For prognosis**, the answer appears to be yes. Current staging systems classify patients into broad categories, but within any given stage, outcomes still vary widely. Our signature stratifies patients within the same stage — among Stage I–II patients, those with high risk scores had a hazard ratio of 3.37 compared to low-risk patients. This is the kind of information that could influence decisions about surveillance intensity, transplant eligibility, and adjuvant therapy. A patient classified as low-risk by both stage and molecular signature might be safely observed, while a patient who is early-stage but molecularly high-risk might benefit from closer follow-up or earlier intervention.

**For treatment selection**, the connection to ferroptosis biology opens a specific therapeutic angle. The finding that high-risk tumors are strongly predicted to be sensitive to erastin — a compound that kills cells by inhibiting the very transporter (SLC7A11) that is upregulated in our high-risk group — is not a coincidence. It follows directly from the biology: these tumors have become dependent on SLC7A11-mediated cystine import to feed their glutathione pool and suppress ferroptosis. Block that import, and the house of cards collapses. While ferroptosis inducers are not yet in clinical use for HCC, several are in preclinical and early clinical development. Our signature could serve as a companion biomarker — a way to select the patients most likely to respond.

**For understanding the disease**, the signature illuminates something fundamental about how HCC tumors survive. The fact that nine of eleven genes fall under the NRF2-KEAP1 pathway is striking. NRF2 activation is one of the most common molecular events in HCC, yet it is not routinely measured or incorporated into clinical decision-making. Our signature provides a transcriptomic readout of NRF2 pathway activity and shows that this readout predicts survival, immune evasion (via NK cell depletion), and mutational burden. The connection between high ROS defense, immune escape, and genomic instability paints a coherent picture: tumors that invest heavily in antioxidant defense create a microenvironment that is both mutagenic and immunosuppressive.

**For the field**, this work demonstrates that biology-driven gene selection — starting from a specific hypothesis about which pathways matter, rather than mining the transcriptome agnostically — produces signatures that are both prognostic and interpretable. The signature does not merely predict survival; it tells a biological story that connects to druggable targets and testable hypotheses. This is the difference between a statistical tool and a scientific one.

None of this means the signature is ready for the clinic tomorrow. It needs prospective validation, ideally in the context of a clinical trial. The protein-level concordance, while encouraging for some genes, was incomplete — a common limitation of mRNA-based biomarkers. And the European validation cohort was underpowered, leaving a gap in the evidence for Western patients. These are honest limitations, and we address them directly in the study.

But the foundation is solid: an 11-gene signature, validated across three continents, rooted in well-understood biology, connected to emerging therapies, and integrated into a nomogram that demonstrably outperforms stage alone. That is a useful starting point for both clinical translation and further biological investigation.

---

### Proposed Journal Targets

1. *Hepatology* — leading liver disease journal, strong fit for HCC prognostic biomarker studies
2. *Journal of Hepatology* — European counterpart, emphasis on translational work
3. *Cancer Research* — broader audience, appropriate for the biological characterization depth
4. *Briefings in Bioinformatics* — if positioned as a computational/bioinformatics contribution

### Proposed Timeline

| Phase | Duration |
|-------|----------|
| Manuscript writing and figure preparation | 4–6 weeks |
| Internal review and revision | 2 weeks |
| Submission | Week 8 |
