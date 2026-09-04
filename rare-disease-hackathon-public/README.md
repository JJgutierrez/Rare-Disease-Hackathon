# Rare Disease, Real Kid: MVA Hackathon 2026 — Track 1 Model 2

**Team Name:** EventArgs  
**Submission Model:** Track 1 Model 2 — Genome-Wide Phenotype-First Candidate Recovery  
**License:** CC BY 4.0  

---

## Non-Clinical Research Disclaimer

This repository and pipeline contain research tools developed for the Sage Bionetworks "Rare Disease, Real Kid: MVA Hackathon 2026". **This work does not constitute medical care, diagnostic testing, or professional medical advice.** All candidate variant prioritizations are research hypotheses produced for competition benchmarking.

---

## Verbatim Required Organizer Acknowledgment

> "This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible."

---

## Data Privacy & Governance Statement

Controlled challenge data (VCF, BAM, FASTQ, clinical documents) and derived individual-level datasets (candidate parquets, candidate matrices, alignment logs) are **NOT** stored in this public repository. 

- This repository operates strictly on **synthetic mock data** located in `synthetic_examples/`.
- All controlled data and individual-level derivatives are processed in a local protected workspace (`rare-disease-hackathon-private/`) and will be deleted within **30 days post-hackathon**, followed by written email confirmation to `RarediseaserealkidMVAhackathon2026@synapse.org`.

---

## Overview & Methodology

Model 2 establishes a rule-aligned, genome-wide candidate recovery workflow following the zero-score baseline of Model 1 (V2).

### Architecture Highlights:
1. **Genome-Wide Unbiased Inventory**: Extracts all variant alleles directly from the raw callset without panel-restricted pre-filtering.
2. **Traceable Exclusion Audit**: Logs every variant filter decision into an auditable exclusion dataset (`model2_exclusions.parquet`).
3. **Phenotype-First Prioritization**: Scores candidates using HPO terms (`HP:0000252`, `HP:0001511`, `HP:0002671`, `HP:0001438`) extracted from the clinical case narrative.
4. **Implementation-Safe Evidence Formula & Phase Gate**:
   $$\text{RawEvidence} = w_p P + w_f F + w_i I + w_q Q + w_r R + w_g G + w_s S + w_h H$$
   - `confirmed_cis` pairs are excluded via a binary exclusion gate.
   - Unresolved phase ($H=0$) is treated as uncertainty, not evidence of a `trans` configuration.

---

## Reproducibility & Synthetic Mock Execution

To run the deterministic validation harness using synthetic mock data:

```bash
# 1. Activate Python virtual environment
source .venv/bin/activate

# 2. Run automated Pytest validation suite
pytest tests/ -v

# 3. Test local official evaluator harness
python src/evaluation.py synthetic_examples/mock_submission.csv
```
