# Track 1 Model 2 Roadmap — Rule-Aligned Genome-Wide Phenotype-First Candidate Recovery

## Purpose

Build a second, materially distinct Track 1 submission after Model 1 / V2 received a valid portal score of 0.0 with the feedback: **“No match — true variant(s) not found in submission.”**

This roadmap treats Model 2 as an evidence-backed **candidate-recovery** workflow. It does not assume that a plausible MVA gene or a locally simulated evaluator truth key is the hidden answer. The goal is to broaden coverage while preserving correct biology, traceability, and the challenge’s submission requirements.

## Non-negotiable contest rules

- Use `PROBAND01` in `proband_id` on **every** submission row. Do not use the internal VCF sample identifier `WGS_EX2312012`.
- Submit a Track 1 CSV in the template-supported structure:

  ```text
  proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type[,notes]
  ```

- Use locus 1 only and leave all locus-2 columns blank for a single-variant hypothesis.
- Use both locus-1 and locus-2 columns only for a genuine paired-variant hypothesis.
- Preserve exact GRCh38 chromosome, position, REF, and ALT values from the source call set after normalization checks.
- Keep chromosome formatting consistent with the accepted submission convention: bare numeric contigs, such as `4`, `5`, `6`, and `15`, not `chr4` or `chr15`.
- Submit the completed Track 1 methods workbook with the CSV, using the latest organizer-provided template.
- Complete the required commercial Generative-AI disclosure truthfully, listing only providers, models, plans/tiers, and uses that actually occurred.
- Do not publish or redistribute challenge-provided raw genomic data, clinical records, BAM/CRAM files, FASTQ files, or unredacted VCFs unless the challenge data-use terms explicitly permit it.
- Keep the project repository private if permitted during the hackathon, but make it public when the organizers’ final-evaluation requirements begin. Include the latest methods description in the linked repository if the organizer instructions require it.
- Treat every new model/approach as a separately documented submission. Do not overwrite or relabel the prior scored V2 baseline.

## Artifact and version policy

### Preserve Model 1

Archive the accepted first submission exactly as scored:

```text
releases/track1/model1_v2/
  submission_track1_v2_original_scored_zero.csv
  methods_description_form_eventargs_v2.xlsx
  manifest.json
  checksum.sha256
  portal_result.md
```

Record:

```text
Portal submission: 1
Proband: PROBAND01
Result: Rank points 0.0 / 100; F-max 0.000
Portal feedback: No match — true variant(s) not found in submission
Interpretation: The Model 1 submitted allele set did not overlap the evaluator's hidden true variant set.
```

### Create Model 2 separately

Use new names for all Model 2 work:

```text
runs/track1_model2_genomewide/
output/submissions/submission_track1_model2_genomewide.csv
output/submissions/submission_track1_model2_genomewide_audit.json
output/submissions/submission_track1_model2_genomewide_enriched.csv
methods_description_form_eventargs_model2.xlsx
```

Do not write a Model 2 output to `submission_track1_v2_final.csv`. Give Model 2 a fresh SHA-256, audit log, release manifest, evaluator output, and methods workbook.

## Model 2 hypothesis

**Model 2 / Genome-wide single-variant and compound-heterozygous recovery portfolio**

## Gate 0 — Authoritative Rule Verification

Do not generate a final Model 2 candidate portfolio or spend another submission attempt until the following items are confirmed from the current organizer app, current FAQ/rules, current Track 1 template, and official evaluator implementation:

- Maximum Track 1 submissions allowed per team.
- Maximum rows allowed in one Track 1 CSV.
- Exact required/optional CSV columns, including whether `notes` is optional.
- Allowed `finding_type` values and whether multiple `primary` rows are permitted.
- Exact scoring treatment for singleton rows and paired-variant rows.
- Whether partial pair matches receive credit.
- Whether duplicate variants or pair-plus-singleton duplication is allowed or penalized.
- Whether the evaluator uses the best, latest, or all Track 1 submissions.
- Eligibility of SNVs, indels, noncoding variants, mosaic variants, CNVs, SVs, and mitochondrial variants.
- Whether a pair must be same-gene.
- Exact data-use, data-retention, deletion, and public-repository requirements.
- Required repository license, if any.

Save the evidence in `rules_snapshot.md` with date, source URL/path, exact quoted wording, and any screenshots.

## Gate 0A — Data Governance, Privacy, and Release Compliance

Controlled challenge data and all individual-level derived data may be processed only in the participant’s protected working environment for hackathon purposes. They must not be committed to any repository, uploaded to third-party public services, included in a presentation or report, shared with unauthorized people, or retained after the allowed 30-day post-hackathon deletion window.

Before any public repository, report, video, presentation, or social post:

- [ ] Confirm no VCF, BAM, CRAM, FASTQ, phenotype document, participant screenshot, raw alignment image, or individual-level intermediate/derived dataset is committed or uploaded.
- [ ] Confirm no raw data, derivative parquet/CSV, Jupyter notebook output, cache, or test fixture contains participant-level genomic information.
- [ ] Review Git history, Git LFS objects, GitHub Actions artifacts, release attachments, issue attachments, and any connected cloud-sync folders—not only the current working tree—for controlled data or derived individual-level outputs.
- [ ] Run a repository secret/data scan and manually inspect Git history before public release.
- [ ] Use synthetic examples or schema-only mock data for reproducibility demonstrations.
- [ ] Keep source data and derived datasets only in approved local/private environments during the hackathon.
- [ ] Do not attempt to contact the data subject, family members, or MVA Society contacts.
- [ ] Report actual or suspected unauthorized disclosure through Sage Bionetworks' Privacy and Compliance process via Sage Help Center.

Within 30 days after hackathon close:

- [ ] Delete raw VCF, BAM/CRAM, FASTQ, phenotype files, and access downloads.
- [ ] Delete all derived individual-level data: Parquet, CSV, annotated outputs, evidence matrices, alignment artifacts, logs containing variant-level data, caches, notebooks, and temporary directories.
- [ ] Delete copies from local disks, external drives, cloud instances, notebook environments, private repos, backups, recycle bins/trash, and synchronized folders.
- [ ] Preserve only code, non-sensitive documentation, synthetic examples, and organizer-permitted anonymized outputs.
- [ ] Email `RarediseaserealkidMVAhackathon2026@synapse.org` confirming deletion is complete.

## Required organizer acknowledgment

For all public communications, reports, READMEs, pitch credits, and publications arising from this work, use this exact verbatim text:

> "This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible."

## Model 2 candidate classes

Model 1 concentrated on a narrow MVA/spindle-checkpoint panel and paired recessive hypotheses. Model 2 will evaluate multiple plausible inheritance and variant-class paths in parallel:

1. High-impact loss-of-function candidates, including stop-gained, frameshift, essential-splice, and start-lost variants.
2. Moderate-impact missense variants with high phenotype/functional support.
3. Phenotype-supported gene candidates: variants in genes prioritized by the phenotype-aware ranker. Noncoding or intronic variants may advance only when they also meet the predeclared splice/functional-evidence policy.
4. Recessive homozygous candidates.
5. Same-gene compound-heterozygous candidate pairs (confirmed cis excluded via binary exclusion gate).
6. Canonical splice variants as a distinct category; near-splice, intronic, and noncoding candidates under a predeclared versioned policy. Missing splice scores are recorded as unknown and do not contribute positive splice evidence. Missingness is not interpreted as evidence of either benignity or pathogenicity.
7. Low-allele-fraction/mosaic candidates, if present in the source call set.
8. Structural/CNV candidates only if supported by an eligible source call set or a properly documented independent analysis.

The Model 2 objective is broader **truth recovery**, not filling every possible submission row with speculative candidates.

## Phase 2 — Rebuild the candidate universe

Do not rely only on `panel_candidates_annotated.parquet` if that file was created after restricting to MVA-panel retrieval regions. Build a full genome-wide candidate inventory directly from the original source VCF/call set.

Create these distinct outputs:

```text
runs/track1_model2_genomewide/
  model2_all_variants.parquet
  model2_exclusions.parquet
  model2_phenotype_input.json
  model2_evidence_matrix.csv
```

For every retained or excluded candidate, preserve:

```text
source_variant_id, chrom, pos, ref, alt, gene, transcript, consequence, impact,
filter, gt, dp, gq, allele_balance, variant_allele_fraction, population_frequency,
inheritance_bucket, phase_status, included, exclusion_reason
```

**Gate:** Confirm that the candidate universe includes genes and variants outside the original MVA panel. A panel-only rerun is not a materially independent Model 2 approach.

## Phase 3 — Candidate evidence gates

Every final-row candidate must pass all applicable gates.

| Gate | Requirement |
|---|---|
| Source presence | Candidate exactly matches a record in the original source call set. |
| Reference correctness | GRCh38 chromosome, position, REF, ALT, and normalization are verified. |
| QC | Genotype, depth, GQ, allele balance, mapping/read context, and filter status are recorded and acceptable under a declared policy. |
| Mechanism | The functional consequence is compatible with the proposed model. |
| Gene/disease relevance | Phenotype-aware or curated evidence supports inclusion; plausible gene membership alone is insufficient. |
| Population evidence | Frequency is evaluated with build- and allele-matched data where available; absence from a database is not pathogenicity evidence. |
| Pair validity | Both alleles meet eligibility; same-gene/inheritance rationale is documented; confirmed cis pairs are excluded. |
| Phase | Phase is recorded as confirmed trans, confirmed cis, or unresolved. Unresolved phase must not be called trans. |
| Intronic/noncoding evidence | Include only with documented splice prediction, functional evidence, or a clearly declared exploratory policy. |
| Distinctness | The row adds coverage of a meaningfully different hypothesis and is not a redundant permutation. |

## Phase 4 — Phenotype-first independent ranking

Run a genome-wide phenotype-aware prioritizer such as Exomiser or LIRICAL using only HPO terms supported by the challenge clinical material.

### Required records

- Input HPO list/Phenopacket.
- HPO ontology version.
- Tool name, version, container/image or environment lock file.
- Disease and phenotype database release/version.
- Genome build: GRCh38.
- Command/configuration file.
- Inheritance-model settings.
- Full tool output and output hash.

### Required comparison

Create a `model2_evidence_matrix.csv` for the top candidates, with at least:

```text
candidate_id, gene, chrom, pos, ref, alt, candidate_type,
VEP consequence, genotype, DP, GQ, allele balance, filter status,
phenotype rank, phenotype score, disease match, population frequency,
ClinVar/ClinGen status where used, phase status, cis exclusion,
evidence supporting, evidence contradicting, missing evidence,
final inclusion decision, final EPCR
```

Use phenotype ranking as an independent signal. Do not tune HPO inputs or score weights after seeing candidate results solely to promote a desired gene.

## Phase 5 — Pair and singleton policy

### Single-variant rows

Use a single-variant row only when the candidate supports a credible single-allele model, such as an appropriate heterozygous, dominant, de novo, mosaic, or other documented mechanism.

For a valid singleton row:

```csv
PROBAND01,chrom,pos,ref,alt,,,,,epcr,finding_type
```

All four locus-2 cells must be blank.

### Paired rows

Use a pair only when both alleles plausibly form one causal model. Required documentation:

- Both variants are present in the source data.
- Both have compatible genotype and functional evidence.
- Gene relationship and inheritance rationale are explicit.
- Confirmed cis pairs are excluded from a recessive compound-heterozygous model.
- Unresolved phase remains unresolved.

### Redundancy control

Do not submit a pair plus its component singleton variants merely to hedge scoring unless the official rules explicitly permit and make that strategy meaningful. If allowed, each row must still be biologically defensible and must have a documented reason for inclusion.

Prefer a smaller, high-evidence, diverse portfolio over a row-maximized list of speculative variants.

## Phase 6 — Ranking and EPCR policy

Do not assign EPCR values merely as descending labels such as 0.95, 0.94, 0.93.

Implement and freeze a documented score policy before selecting final candidates. Example components:

```text
phenotype support
+ functional consequence
+ genotype/inheritance compatibility
+ source-call quality
+ population rarity evidence
+ gene-disease validity
+ phase support or phase uncertainty penalty
+ splice/SV evidence where applicable
```

Map the resulting normalized evidence score to a 0–1 EPCR value consistently. Document:

- Component definitions.
- Weights or decision rules.
- Missing-data treatment.
- How ties are resolved.
- Why every final candidate’s EPCR is greater than or equal to the next row.

**Required labeling policy:** Use one `primary` row for the best-supported principal hypothesis. Mark alternative hypotheses as `secondary` unless the official rules explicitly define a different policy.

## Phase 7 — Build the Model 2 submission

Create a fresh CSV only after all gates pass:

```text
output/submissions/submission_track1_model2_genomewide.csv
```

Minimum validation checks:

```text
[ ] Header matches the current official template/evaluator requirement.
[ ] Every row has proband_id = PROBAND01.
[ ] All chromosomes use accepted formatting.
[ ] All positions are positive integers.
[ ] REF and ALT are valid and exactly source-matched.
[ ] Singleton rows have all locus-2 fields blank.
[ ] Pair rows have all locus-2 fields completed.
[ ] No malformed partial second locus exists.
[ ] EPCR values are numeric, in the allowed range, and follow the frozen ranking policy.
[ ] finding_type values follow official permitted values.
[ ] Row count is within the official maximum.
[ ] No duplicate or redundant rows without explicit documented justification.
[ ] All final variants are present in the original call set.
[ ] Phase and inheritance claims remain scientifically accurate.
```

Run the official/local evaluator exactly as supplied by the challenge. Then create:

```text
submission_track1_model2_genomewide.csv.sha256
submission_track1_model2_genomewide_evaluator.log
submission_track1_model2_genomewide_audit.json
submission_track1_model2_genomewide_manifest.json
```

## Phase 8 — Update the Track 1 methods form

Use the current `methods_description_form.xlsx` and save a separate Model 2 copy:

```text
methods_description_form_eventargs_model2.xlsx
```

Required content:

- Team name: `EventArgs`.
- Model number: `Model 2 / Genome-wide single-variant and compound-heterozygous recovery portfolio`.
- Clear statement that Model 2 broadens the candidate universe beyond Model 1’s MVA-focused pair model after a no-match benchmark result.
- Accurate model, public-data, compute/runtime, manual-review, and generative-AI disclosures.
- Explicit statement that AI was advisory only and did not autonomously generate/alter final variant calls or unsupported biological claims.
- Accurate single-variant and paired-candidate handling.
- Clear secondary-finding policy.
- Limitations: hidden truth unknown, phase constraints, variant-class sensitivity, no clinical diagnostic claim, and any unavailability of parental/long-read validation.

Do not copy Model 1’s exact methods description if Model 2 has a different search universe, ranking method, candidate policy, or evidence sources.

## Phase 9 — Submission and audit

Upload only these Model 2 files:

```text
submission_track1_model2_genomewide.csv
methods_description_form_eventargs_model2.xlsx
```

Before final portal submission:

```text
[ ] Confirm the portal recognizes PROBAND01.
[ ] Confirm the displayed row count is expected.
[ ] Confirm no schema or parsing error appears.
[ ] Confirm the selected file is Model 2, not the Model 1 archive.
[ ] Screenshot validation and submission confirmation.
```

After submission, log:

```text
Portal submission number
Submission timestamp in local time and UTC
Portal submission ID
CSV SHA-256
Methods workbook filename/version
Portal rank points
F-max
F-max threshold
Portal feedback verbatim
```

Never overwrite Model 1 artifacts after Model 2 is submitted.

## Implementation Sequence

1. Complete `rules_snapshot.md` using live portal, current evaluator, template, and official rules text.
2. Archive accepted zero-score Model 1/V2 submission outside active output directory.
3. Create protected/private Model 2 workspace and verify exclusion from Git, cloud sync, LFS, and CI.
4. Generate complete genome-wide normalized variant inventory.
5. Build full exclusion audit before applying aggressive filters.
6. Create HPO input from clinical material and run phenotype-first prioritization.
7. Apply evidence formula and inheritance/phase gates.
8. Select small, nonredundant candidate portfolio within confirmed row limit.
9. Generate Model 2 CSV, methods workbook, audit manifest, evaluator log, and SHA-256.
10. Run local tests and current organizer evaluator.
11. Validate portal-recognized `PROBAND01` and row count before submission.
12. Prepare report, clean repository, and 3-minute pitch without exposing controlled data.

## Final Go/No-Go Checklist

Do **NOT** upload Model 2 if any of the following remain unresolved:
- [ ] Maximum rows/attempts are still unverified.
- [ ] Candidate portfolio was generated from speculative hidden-answer coverage.
- [ ] A final variant is absent from the source call set.
- [ ] A singleton lacks a coherent single-allele mechanism.
- [ ] A recessive pair is confirmed cis or lacks a same-gene/inheritance rationale.
- [ ] An intronic/noncoding candidate lacks positive functional/splice support.
- [ ] EPCR values are arbitrary.
- [ ] Final methods workbook still describes Model 1.
- [ ] Controlled data or derived individual-level outputs are present in the public repo.

## Official-Evaluator Test Gate

- [ ] Run the current evaluator implementation supplied by the organizers against the final CSV.
- [ ] Save the exact command, evaluator source version/commit if available, stdout/stderr, and exit code in the Model 2 manifest.
- [ ] Confirm the live portal sees the intended number of rows and `PROBAND01` before final submission.

## Manual Final Review

- [ ] Confirm each submitted coordinate and REF/ALT allele against the original GRCh38 source VCF.
- [ ] Confirm every singleton has an evidence-supported single-allele mechanism.
- [ ] Confirm every pair has a documented same-gene/inheritance rationale and that no confirmed cis pair is included.
- [ ] Review the top candidate evidence matrix for contradictory population, phenotype, QC, or inheritance evidence.
- [ ] Confirm the Model 2 CSV uses a new filename (`submission_track1_model2_genomewide.csv`) and has not overwritten any Model 1/V2 artifact.
- [ ] Confirm that the final methods workbook describes Model 2—not Model 1—and accurately discloses all commercial GenAI use.
- [ ] Confirm that the portal-selected file is the Model 2 artifact and that the portal recognizes `PROBAND01` before submission.

## Deliverables

A Model 2 submission is ready only when these are complete:

```text
[ ] rules_snapshot.md
[ ] model2 candidate-universe outputs
[ ] model2 exclusion audit
[ ] phenotype-aware ranking configuration and results
[ ] evidence matrix for finalists
[ ] frozen ranking/EPCR policy
[ ] submission_track1_model2_genomewide.csv
[ ] evaluator output showing success
[ ] SHA-256 checksum
[ ] manifest and audit JSON
[ ] methods_description_form_eventargs_model2.xlsx
[ ] portal submission log and screenshots
```

## Success criteria

Model 2 is successful as an engineering and competition submission only if it is:

- Valid under the current portal rules.
- Materially distinct from Model 1.
- Broader than the original MVA-panel pair-only approach.
- Evidence-backed rather than slot-maximized speculation.
- Reproducible from versioned inputs, configuration, code, and output hashes.
- Scientifically cautious about phase, causality, and clinical interpretation.
