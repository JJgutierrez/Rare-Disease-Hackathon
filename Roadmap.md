```markdown
# Track 1 Engineering Roadmap: Causal Variant Prioritization (SDD Framework)
**Project:** SageBio "Rare Disease, Real Kid" MVA Hackathon 2026 (Track 1)[cite: 1]  
**Target:** Mosaic Variegated Aneuploidy (MVA) Proband Genomic Identification[cite: 1]  
**Methodology:** Spec-Driven Development (SDD) & Automated Verification Harness[cite: 1]  
**Date:** September 2026[cite: 1]  
**Status:** In Active Execution[cite: 1]  

---

## 1. System Specification & Goal Definition

### 1.1 Objective
Engineer a deterministic, reproducible, spec-driven computational pipeline to analyze raw Whole-Genome Sequencing (WGS) callsets (`WGS_EX2312012_HGWCNDSX7.vcf.gz`) alongside proband clinical presentation (`Challenge_Clinical_Phenotype_1.docx`) to isolate, rank, and phase causal pathogenic/compound-heterozygous variants responsible for Mosaic Variegated Aneuploidy (MVA)[cite: 1]. Output a schema-compliant, verified top-10 candidate CSV scored by the automated SageBio evaluation harness (`evaluation.py`)[cite: 1].

### 1.2 Fixed Constraints & Environmental Context
* **Platform / Compute:** macOS (Apple M-series unified memory architecture) running Antigravity IDE[cite: 1].
* **Storage Architecture & Single Source Policy:** All raw dataset files (~85 GB FASTQ callsets, VCFs, clinical docs), intermediate extractions, sliced VCFs, mini-BAM alignments, and work data MUST reside strictly on the external volume `/Volumes/WD_BLACK/HuggingFace Contest` as the single source of work data. Local repository paths under `data/` are symlinked to this volume and ignored via `.gitignore`.
* **Reference Genome:** GRCh38 / hg38 (`GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta`)[cite: 1]. Contig naming: **bare numbers** (`1`, `2`, ..., `15`, `X`, `Y`), strictly **no `chr` prefix**[cite: 1].
* **Coordinate Invariant:** Absolutely no hardcoded genomic boundaries in specs or source code. Coordinates are resolved dynamically at runtime via Ensembl REST API and cached locally to `specs/target_panel_regions.json`.
* **Submission Window:** August 24, 2026 – October 24, 2026[cite: 1].
* **Submission Quota:** Maximum 6 submissions total[cite: 1]. Every live submission must pass 100% local deterministic unit tests and schema verification before upload[cite: 1].
* **Post-Hackathon Compliance:** Gated dataset requires complete deletion of raw files post-competition and email confirmation to `MVAHackathon2026@synapse.org`[cite: 1].

---

## 2. SDD Verification Architecture & Success Metrics

### 2.1 Quantitative Acceptance Gates
Every pipeline iteration must satisfy the following invariant criteria before proceeding:

| Metric Gate | Threshold | Verification Method |
| :--- | :--- | :--- |
| **G1: Schema Compliance** | 100% adherence to 11-column specification[cite: 1] | `pytest tests/test_schema.py`[cite: 1] |
| **G2: Coordinate Sanity** | 0 coordinate syntax errors; bare contigs; POS within chromosome bounds[cite: 1] | `pytest tests/test_coordinates.py`[cite: 1] |
| **G3: Rank & EPCR Monotonicity** | $\le 10$ rows; $1.0 \ge \text{EPCR} \ge 0.0$; strictly monotonically decreasing[cite: 1] | `pytest tests/test_ranking.py`[cite: 1] |
| **G4: Allele Parity** | REF/ALT alleles match reference GRCh38 base at POS[cite: 1] | Biopython / local reference lookup script[cite: 1] |
| **G5: Phasing & CNV Integrity** | Disqualify *in-cis* variants; require dual hits or hemizygous CNV drop | Targeted alignment or phase-block inspection |

### 2.2 Submission Schema Invariant (`evaluation.py` Target)
```csv
proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type

```

* `proband_id`: Extracted exact proband identifier matching challenge ground-truth.


* `chrom_1`, `pos_1`, `ref_1`, `alt_1`: Primary allele coordinates and bases.


* `chrom_2`, `pos_2`, `ref_2`, `alt_2`: Paired secondary allele for compound heterozygotes (leave empty for single dominant/de novo proposals).


* `epcr`: Estimated Pathogenic Confidence Ratio / Probability $[0.0, 1.0]$.


* `finding_type`: Literal string `primary` or `secondary`.



### 2.3 Formal EPCR Scoring Model

Every candidate pair or single locus is scored via the deterministic formula:

$$\text{EPCR} = \min\left(1.0, \, \frac{S_{\text{allele}_1} + S_{\text{allele}_2}}{2} \times M_{\text{pheno}} \times M_{\text{phase}}\right)$$

* **Allele Deleteriousness Score ($S \in [0.0, 1.0]$):**
* Frameshift / Stop-gained / Canonical splice ($\pm 1, 2$ bp): **1.0**
* SpliceAI $\Delta \text{Score} \ge 0.50$: **0.90**
* AlphaMissense "likely pathogenic" or CADD $\ge 25$: **0.75**
* Missense (ambiguous/moderate) or SpliceAI $0.20 \le \Delta < 0.50$: **0.40**
* In-frame indel / synonymous: **0.10**
* Single allele dominant/de novo proposal: Set $S_{\text{allele}_2} = 1.0$ only if clear hemizygous/CNV deletion is demonstrated; otherwise, single variants default to $S_{\text{allele}_2} = 0.0$.


* **Phenotype Match Multiplier ($M_{\text{pheno}}$):**
* Core MVA loci (*BUB1B*, *TRIP13*, *CEP57*): **1.0**
* Extended SAC / Centrosomal loci (*MAD2L1*, *CDC20*, *BUB1*, *BUB3*, *TTK*, *CENPE*, *PLK1*): **0.70**


* **Phasing & State Multiplier ($M_{\text{phase}}$):**
* Confirmed *in-trans* (physical read-pair overlap, PID phase block, or confirmed heterozygous deletion): **1.0**
* Unphased (distant $>1\text{ kb}$, no read overlap, both variants rare/deleterious): **0.80**
* Confirmed *in-cis* (variants present on the same physical read or phase block): **0.0** (Disqualified)



---

## 3. Technical Roadmap & Phase Breakdown

```
[Phase 0: Spec & Test Harness]
       │
       ▼
[Phase 1: Clinical Parsing & Dynamic API Panel Coordinate Fetch]
       │
       ▼
[Phase 2: Streamed VCF Slicing & Quality Filtering]
       │
       ▼
[Phase 3: Targeted Mini-BAM Alignment, Phasing & CNV Exon Depth Audit]
       │
       ▼
[Phase 4: Multi-Model Pathogenicity & Splice Annotation]
       │
       ▼
[Phase 5: Biallelic Pairing, Local Validation & Submission Build]
       │
       ▼
[Phase 6: Leaderboard Submission Ledger & Controlled Iteration]

```

---

### Phase 0: SDD Test Harness Setup (Prerequisite)

**Goal:** Establish fail-fast test suites verifying submission artifacts before any model runs.

* **Spec File:** `specs/submission_schema.json`

* **Test Suites:**
* `tests/test_schema.py`: Verifies column count (11), names, casing, non-null fields, and row bounds ($1 \le N \le 10$).


* `tests/test_formatting.py`: Asserts no `chr` prefix, checks numeric POS, valid IUPAC nucleotides, and monotonic descending `epcr` values.


* `tests/test_mock_eval.py`: Simulates `evaluation.py` harness scoring logic against a mock truth set.




* **Exit Criteria:** `pytest tests/` runs and fails on invalid CSV fixtures, passes on compliant mock submission.



---

### Phase 1: Clinical Parsing & Dynamic Gene Panel Fetch

**Goal:** Parse proband phenotype and resolve exact GRCh38 gene coordinates via Ensembl REST API.

* **Tasks:**
1. Execute `scripts/01_parse_phenotype.py` on `Challenge_Clinical_Phenotype_1.docx`.


2. Extract: Proband ID string, age of onset, dysmorphisms, tumor history, and aneuploidy percentages.


3. Map clinical terms to Human Phenotype Ontology (HPO) codes (e.g., `HP:0000252` Microcephaly, `HP:0001511` IUGR, `HP:0002659` Wilms tumor susceptibility).


4. Run `scripts/utils/fetch_coords.py` against Ensembl REST API (`https://rest.ensembl.org/lookup/symbol/homo_sapiens/`) with uniform $\pm 50\text{ kb}$ padding for:
* Core MVA loci: `BUB1B`, `TRIP13`, `CEP57`

* Extended SAC/centrosomal loci: `MAD2L1`, `CDC20`, `BUB1`, `BUB3`, `TTK`, `CENPE`, `PLK1`



5. Cache verified coordinates to `specs/target_panel_regions.json`.


* **Deliverables:** `output/clinical_summary.json` and `specs/target_panel_regions.json`.


* **Exit Criteria:** `specs/target_panel_regions.json` exists, contains bare contigs (`15`, `5`, etc.), and validates against Ensembl GRCh38 release boundaries.

---

### Phase 2: Streamed VCF Ingestion & Quality Filtration

**Goal:** Stream the indexed WGS VCF and isolate candidate variants falling inside padded panel regions.

* **Tasks:**
1. Execute `scripts/02_slice_vcf.py` using `cyvcf2` streaming region queries.
2. Apply primary quality gates:
* `FILTER == "PASS"`

* Genotype Quality `GQ >= 20`

* Read Depth `DP >= 10`

* Exclude homozygous reference (`GT == 0/0`).




3. Check presence of Sentieon `PID` and `PGT` phase tags across panel variants.


* **Deliverables:** `output/panel_candidates_raw.parquet` and `output/panel_candidates_raw.csv`.
* **Exit Criteria:** Raw candidate set reduced from ~4.5M genomic variants to $< 500$ panel-specific variants.



---

### Phase 3: Targeted Mini-BAM Alignment, Phasing & CNV Exon Depth Audit

**Goal:** Resolve physical phase and detect heterozygous deletions without full-genome realignment.

* **Tasks:**
1. Build a mini-reference FASTA (`mini_panel_ref.fasta`) containing only the padded panel genomic regions (~2–3 MB total).
2. Stream raw FASTQs through `bwa-mem2 mem` / `minimap2` against `mini_panel_ref.fasta` to produce a sorted, indexed `mini_panel.bam` (~50–100 MB). Runtime: $<15$ minutes on Apple M-series.
3. **CNV / Exon Deletion Screen:** Run `samtools depth` across all known exons of *BUB1B*, *TRIP13*, and *CEP57*. Flag any continuous drop to ~50% baseline coverage indicating a heterozygous structural deletion.
4. **Physical Read-Pair Phasing:** For candidate variants separated by $< 600$ bp, inspect read pairs to verify whether mutations are *in-cis* or *in-trans*.


* **Deliverables:** `output/mini_panel.bam`, `output/exon_coverage_profiles.csv`, and `output/phasing_evidence.json`.
* **Exit Criteria:** Every candidate variant has documented physical phase context or is marked unphased ($M_{\text{phase}} = 0.80$).

---

### Phase 4: Multi-Model Pathogenicity & Splice Annotation

**Goal:** Annotate candidate variants with population frequencies and multi-algorithm deleteriousness scores.

* **Tasks:**
1. **Population Frequency Filter:** Query gnomAD v3/v4 via local tabix or Ensembl VEP API:
* Primary filter: Allele Frequency (AF) $< 0.01$ (1%).
* Secondary consideration: $0.005 \le \text{AF} < 0.01$.


2. **Coding Impact:** Query local precomputed AlphaMissense GRCh38 TSV (~90 MB) and CADD/REVEL databases.
3. **Splice Disruption:** Query the Broad Institute SpliceAI lookup service / REST API for all candidate variants. Record $\Delta \text{Scores}$ for acceptor gain/loss and donor gain/loss.


* **Deliverables:** `output/panel_candidates_annotated.parquet`.
* **Exit Criteria:** All retained variants have populated values for `gnomad_af`, `alphamissense_class`, `cadd_phred`, and `spliceai_max_ds`.

---

### Phase 5: Biallelic Pairing, Neutral Phasing EPCR & Submission Build

**Goal:** Synthesize same-gene candidate pairs, compute bounded additive EPCR scores with stable tie-breakers, exclude cis pairs, deploy advisory AI adjudication, and generate strict 11-column evaluator-compliant upload CSV.

* **Tasks:**
1. Group annotated variants strictly by `gene_symbol`.
2. Apply Phasing Rules:
   * `trans_supported`: $M_{\text{phase}} = 1.00$
   * `phase_ambiguous`: $M_{\text{phase}} = 0.50$ (`neutral_uncertainty`)
   * `cis_supported`: $M_{\text{phase}} = 0.00$ (**strictly excluded** from compound heterozygous shortlist)
3. Compute bounded additive EPCR using weights defined in `specs/pipeline_config.json`:
   $$\text{EPCR}_{\text{raw}} = \sum_i w_i M_i, \quad \text{EPCR}_{\text{final}} = \min(1.0, \max(0.0, \text{EPCR}_{\text{raw}}))$$
4. Perform deterministic tie-breaker sorting: `trans_supported` > core gene tier > consequence severity > lower AF > SpliceAI > genomic POS.
5. Generate deliverables:
   * `output/submissions/submission_track1_v1_enriched.csv` (internal audit file with raw scores and debug features).
   * `output/submissions/submission_track1_v1_final.csv` (strict 11-column evaluator upload CSV, top 10 rows).
   * `output/submissions/submission_track1_v1_audit.json` (strict step-by-step audit funnel chain with input/output SHA-256 hashes).
   * `output/submissions/submission_track1_v1_adjudication.json` & `validation.json` (advisory AI review cards).
6. Run local validation suite:
   ```bash
   .venv/bin/pytest -v
   ```

* **Deliverables:** `submission_track1_v1_final.csv`, `submission_track1_v1_enriched.csv`, `submission_track1_v1_audit.json`, `submission_track1_v1_adjudication.json`.
* **Exit Criteria:** `pytest tests/` exits 0 (13/13 passing); strict 11-column header match; 0 cis pairs in final top 10.

---

## Phase 6: Leaderboard Submission Ledger & Controlled Iteration

**Goal:** Optimize evaluation feedback across the 6 allowed submission attempts.

| Submission # | Strategy Focus | Rationale |
| --- | --- | --- |
| **Submission 1** | **Primary Best Candidate Pair (V1)** | Immediate, high-confidence submission (no throwaways). Benchmarks pipeline baseline. |
| **Submission 2** | **Alternative Splicing Hypothesis** | Prioritize deep intronic/cryptic splice variants flagged by SpliceAI $\ge 0.50$. |
| **Submission 3** | **Structural / Deletion Pairing** | Pair top SNV with confirmed exon coverage drop if CNV screen indicates hemizygosity. |
| **Submission 4** | **Secondary Gene Hypotheses** | Explore high-scoring candidates in *TRIP13*, *CEP57*, or *CDC20*. |
| **Submission 5** | **Expanded Panel Ensemble** | Integrate candidate alleles across secondary mitotic checkpoint regulators. |
| **Submission 6** | **Final Tuned Consensus** | Best composite configuration incorporating all feedback and leaderboard signals. |

* **Tracking Ledger:** Maintain `output/submissions/submission_log.md` recording submission timestamp, file SHA-256 hash, git commit, and returned evaluation score.

---

## 4. Antigravity Execution Commands

```bash
# 1. Activate venv
source .venv/bin/activate

# 2. Extract clinical data & proband ID
python scripts/01_parse_phenotype.py

# 3. Dynamically resolve gene panel coordinates from Ensembl REST
python scripts/utils/fetch_coords.py

# 4. Stream and slice VCF to candidate panel regions
python scripts/02_slice_vcf.py

# 5. Execute test harness
pytest tests/

```

```

```