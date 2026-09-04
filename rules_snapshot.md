# Rules Snapshot — Rare Disease MVA Hackathon 2026

**Date of Capture:** September 4, 2026  
**Source Platform:** Hugging Face Space `SageBio/rare-disease-real-kid-mva-hackathon-2026`  
**Purpose:** Authoritative rules ledger for Track 1 Model 2 engineering, validation, and submission.

---

## 1. Confirmed Contest Rules & Requirements

| Item | Confirmed Requirement | Implementation Rule |
| :--- | :--- | :--- |
| **Track 1 Objective** | Identify causal disease variant(s) in `PROBAND01` evaluated automatically against **NHS-validated causal variant(s)**. | Target true variant recovery via genome-wide candidate extraction and phenotype prioritization. |
| **Track 1 Metric** | Scored automatically using Rank Points (100 max) and F-max threshold metric. | Compute EPCR monotonically from raw evidence scores; prioritize candidate precision. |
| **`PROBAND01` Identifier** | Every submitted CSV row must specify `proband_id` as `PROBAND01`. | Enforce automated pre-flight assertion (`set(proband_id) == {"PROBAND01"}`). |
| **Contig Format** | Bare numeric/letter contig syntax accepted by portal (`4`, `15`, `X`, no `chr` prefix). | Strip `chr` prefixes during GRCh38 normalization. |
| **Singleton Convention** | Single-variant rows populate Locus 1 (`chrom_1`–`alt_1`) and leave Locus 2 (`chrom_2`–`alt_2`) **completely blank**. | Enforce 4 blank cells for Locus 2 on all singletons. |
| **Submission Deliverables** | 1. Track 1 Submission CSV (`submission_track1_model2_genomewide.csv`) <br>2. Track 1 Methods Form Workbook (`methods_description_form_eventargs_model2.xlsx`) <br>3. Written Submission Report <br>4. Public GitHub Repository (Code & synthetic test data only) <br>5. 3-Minute Video Pitch | Complete all 5 submission package components. |
| **LLM / GenAI Disclosure** | Mandatory disclosure of commercial Generative AI tool usage in the methods description form. | Disclose tool providers, models, tiers, and advisory-only usage accurately in the workbook. |
| **Data Privacy & Sharing** | **Strict Prohibition**: Participants may not release, grant access to, or reshare raw or derived individual-level data. | Store raw data & derived Parquet/CSV candidate matrices ONLY in local protected workspace `rare-disease-hackathon-private/`. |
| **30-Day Data Retention** | **Mandatory Deletion**: Delete all raw files AND derived individual-level datasets within 30 days post-hackathon. | Execute 30-day deletion checklist and email written confirmation to `RarediseaserealkidMVAhackathon2026@synapse.org`. |
| **Submission Licensing** | Submission deliverables released under CC BY 4.0 as required. Public repo will include a CC BY 4.0 `LICENSE`. | Include CC BY 4.0 `LICENSE` in clean public repo directory `rare-disease-hackathon-public/`. |
| **Publication Embargo** | No peer-reviewed manuscript submission using dataset prior to official organizer summary report/preprint. | Comply with embargo timeline for academic publications. |
| **Presentation Approval** | Conference abstracts or posters require prior written organizer approval. | Submit abstract/poster drafts for written organizer signoff. |
| **Required Attribution** | Mandatory verbatim acknowledgment text in README, report, pitch credits, and manuscripts. | Include exact verbatim organizer acknowledgment. |
| **Dataset Citation** | Cite Synapse dataset reference at time of publication. | Retrieve and use official Synapse citation reference. |

---

## 2. Verified Evaluator Contracts (Source: `scripts/utils/evaluation.py`)

| Parameter | Evaluator Contract | Verification Status |
| :--- | :--- | :--- |
| **Max Rows per Proband** | Up to **10 rows** per proband (`len(rows) <= 10`). | **VERIFIED IN EVALUATOR CODE** |
| **EPCR Range & Sorting** | `0 < epcr <= 1.0`. Evaluator sorts rows by `epcr` descending automatically. | **VERIFIED IN EVALUATOR CODE** |
| **`finding_type` Permitted Values**| Must be `"primary"` or `"secondary"` (case-insensitive, default `"primary"`). | **VERIFIED IN EVALUATOR CODE** |
| **Required Columns** | `proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr[,finding_type]` (`notes` is optional). | **VERIFIED IN EVALUATOR CODE** |
| **Singleton Locus 2 Format** | Blank Locus 2 (`chrom_2,pos_2,ref_2,alt_2` missing or empty strings) parses as 1-variant hypothesis (`frozenset([v1])`). | **VERIFIED IN EVALUATOR CODE** |
| **Paired Variant Format** | Populated Locus 2 parses as compound-heterozygous pair (`frozenset([v1, v2])`). | **VERIFIED IN EVALUATOR CODE** |
| **Partial Match Credit** | If compound het truth has 2 variants and full pair is not matched, a partial row matching 1 truth variant receives **50% of rank points** at its rank. | **VERIFIED IN EVALUATOR CODE** |
| **Rank Point Tiers** | Rank 1 = **100 pts**; Rank 2–3 = **50 pts**; Rank 4–5 = **25 pts**; Rank 6–10 = **10 pts**. | **VERIFIED IN EVALUATOR CODE** |

---

## 3. Pending Verification Parameters (Live Portal Operational Checks)

- [ ] **Max Track 1 Submissions Quota**: Total allowed submission attempts per team (default assumption: 6 total, 5 remaining).
- [ ] **Track 1 Judging Policy**: Confirm whether Track 1 uses the best, latest, or all submissions for final leaderboard ranking.
- [ ] **Eligible Variant Classes**: Confirm eligibility of noncoding, mosaic, CNV, and SV candidates.
- [ ] **Cross-Gene Pair Eligibility**: Confirm whether candidate pairs must belong to the same gene (default policy: restrict to same gene).

---

## 3. Required Verbatim Organizer Acknowledgment

> "This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible."
