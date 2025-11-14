# vcf-triage-cli

[![CI](https://github.com/munaberhe/vcf-triage-cli/actions/workflows/ci.yml/badge.svg)](https://github.com/munaberhe/vcf-triage-cli/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](#license)

A minimal, test-backed **CLI for germline VCF triage**. It filters single-sample VCFs by PASS/QUAL/DP/AF and heterozygous allele balance, extracts key VEP annotations (gene, consequence, HGVS), supports a **gene allowlist** (`--genes`), and can write a short **summary report** (`--summary`). Built to demonstrate practical GLH/NHS-style workflows.

---

## Features

- **Filters**: `FILTER=PASS` (optional), `QUAL ≥ min-qual`, `DP ≥ min-dp`, `AF < max-af` (if present), **AB 0.3–0.7** for hets (if `AD` present).
- **Annotation extract**: from VEP `ANN` → gene, consequence, HGVS c./p.
- **Gene allowlist**: keep only genes in a provided file (`--genes`).
- **Summary**: quick counts + consequence breakdown (`--summary`).
- **Whitespace-robust**: accepts tab- or space-delimited toy VCFs.
- **No heavy deps**: pure Python; optional `conda` env for reproducibility.
- **Tests + CI**: `pytest` with GitHub Actions workflow.

---

## Why this matters (clinical rationale)

- **Depth (DP)** ensures enough reads to trust a call.
- **QUAL** proxies variant caller confidence.
- **Population frequency (AF)** keeps rare variants for rare disease contexts.
- **Allele balance (AB)** for heterozygotes guards against technical bias.
- Clear **CSV output** helps move from raw VCF to clinician-ready notes.

---

## Installation

**Option A: System Python **
```bash
git clone https://github.com/munaberhe/vcf-triage-cli.git
cd vcf-triage-cli
python -m pip install -U pip
python -m pip install pytest
```

**Option B: Conda environment 
```bash
git clone https://github.com/munaberhe/vcf-triage-cli.git
cd vcf-triage-cli
conda create -n vcf-triage python=3.11 -y
conda activate vcf-triage
python -m pip install pytest
```

> Optional speed-up: you can install `cyvcf2` later for bigger files. This minimal CLI doesn’t require it.

---

## Quick start

Run the tool on the included sample:

```bash
python -m src.vcf_triage --vcf data/sample.vcf --out out.csv
```

Open `out.csv` to see the kept variants.

---

## Examples

**Basic:**
```bash
python -m src.vcf_triage --vcf data/sample.vcf --out out.csv
```

**Keep only a gene subset & write a summary:**
```bash
printf "GENE1\n" > data/keep_genes.txt
python -m src.vcf_triage --vcf data/sample.vcf --out out.csv \
  --genes data/keep_genes.txt --summary summary.txt
cat summary.txt
```

**Include non-PASS sites (exploratory):**
```bash
python -m src.vcf_triage --vcf data/sample.vcf --out out_all.csv --include-nonpass
```

**Tighter thresholds:**
```bash
python -m src.vcf_triage --vcf data/sample.vcf --out out_20x.csv --min-dp 20 --min-qual 35
```

---

## Command-line reference

```
usage: vcf_triage.py --vcf <path> --out <path> [options]

required:
  --vcf PATH            Input .vcf or .vcf.gz (single-sample)
  --out PATH            Output CSV path

filters (defaults shown):
  --min-dp INT          Minimum depth to keep (default: 10)
  --min-qual FLOAT      Minimum QUAL to keep (default: 30.0)
  --max-af FLOAT        Maximum AF if present (default: 0.01)
  --include-nonpass     Include sites where FILTER != PASS

extras:
  --genes PATH          File of genes to keep (one per line)
  --summary PATH        Write a tiny text summary to this path
```

**Output columns:**
```
chrom,pos,ref,alt,gene,consequence,hgvs_c,hgvs_p,af,gt,dp,ab,filters
```

---

## Data & formats

- **VCF**: expects standard fields; works with `INFO/AF`, `FORMAT/DP`, `FORMAT/AD`, and VEP `INFO/ANN` if present.
- **Allele Balance (AB)**: computed as `alt / (ref + alt)` from `AD` when available.
- **HGVS**: pulled from `ANN` fields when present (VEP layout).

---

## Project structure

```
vcf-triage-cli/
  data/
    sample.vcf
  src/
    vcf_triage.py
  tests/
    conftest.py
    test_basic.py
    test_genes.py
    test_summary.py
  .github/
    workflows/
      ci.yml
```

---

## Testing

Run all tests:
```bash
pytest -q
```

If `pytest` can’t import `src/`, `tests/conftest.py` adds the repo root to `PYTHONPATH`.

---

## CI

A minimal GitHub Actions workflow runs on every push/PR:
```
.github/workflows/ci.yml
```
The badge at the top of this README reflects the current status.

---

## Roadmap 

- [ ] Add `--vcf-index` support and `.vcf.gz` indexing tips.
- [ ] Canonical transcript selection for `HGVS` (VEP/Ensembl API notes).
- [ ] Optional `--vep-cache` / `--dbnsfp` annotations (docs only).
- [ ] Streamlit mini-dashboard for CSV review.
- [ ] IGV batch generator from CSV (links to BAM/CRAM).

---

## Information Governance (NHS context)

- No PHI in paths, filenames, commit messages, or screenshots.
- Use pseudonymised IDs; restricted access; approved secure transfer only.
- Document retention & deletion policies; keep an audit trail.

---

## Screenshots

_Add a small screenshot of `out.csv` here for recruiters._
```markdown
![triage-output](docs/out_csv_screenshot.png)
```

---
## IGV batch generation

```bash
# single BAM for all rows
python -m src.igv_batch --csv out.csv \
  --bam "/absolute/path/to/sample.bam" \
  --out igv_batch.txt --genome GRCh38 --flank 100 \
  --snapshot-dir igv_snapshots --snapshot-prefix triage

## Notes
- VCF format: HTS/SAM spec  
- VEP `ANN` convention  
- Inspired by routine triage steps in germline diagnostics (GLH/NHS)

