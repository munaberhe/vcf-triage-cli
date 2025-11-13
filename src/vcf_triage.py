"""
VCF Triage CLI (whitespace-robust) + --genes allowlist + --summary report

Filters a single-sample VCF by:
  - PASS (unless --include-nonpass)
  - QUAL (>= --min-qual)
  - Depth DP (>= --min-dp)
  - Population AF (< --max-af) when AF exists in INFO
  - Allele balance for hets (0.3–0.7) when AD exists
  - Optional gene allowlist via --genes (keep only listed genes)

Also extracts basic VEP ANN fields (gene, consequence, HGVS c./p.).

Usage:
  python -m src.vcf_triage --vcf data/sample.vcf --out out.csv
  python -m src.vcf_triage --vcf data/sample.vcf --out out.csv --genes data/keep_genes.txt --summary summary.txt
"""

import argparse
import csv
import gzip
import io
import re
from typing import Dict, Iterator, Optional, Set, Tuple


def open_text_auto(path: str):
    """Open plain text or gzipped text as a text file-handle."""
    return io.TextIOWrapper(gzip.open(path, "rb")) if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def parse_info(s: str) -> Dict[str, str]:
    """Parse the INFO field (key=value;key2;...)."""
    out: Dict[str, str] = {}
    for field in s.split(";"):
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
            out[k] = v
        else:
            out[field] = "True"
    return out


def parse_ann(info: Dict[str, str]) -> Tuple[str, str, str, str]:
    """
    Parse the first VEP ANN record.
    Typical order: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|...|HGVSc|HGVSp
    """
    ann = info.get("ANN", "")
    if not ann:
        return "", "", "", ""
    first = ann.split(",")[0]
    parts = first.split("|")
    gene = parts[3] if len(parts) > 3 else ""
    consequence = parts[1] if len(parts) > 1 else ""
    hgvsc = parts[9] if len(parts) > 9 else ""
    hgvsp = parts[10] if len(parts) > 10 else ""
    return gene, consequence, hgvsc, hgvsp


def allele_balance(fmt_map: Dict[str, str]) -> Optional[float]:
    """Compute AB = alt/(ref+alt) from AD if present."""
    ad = fmt_map.get("AD")
    if not ad:
        return None
    try:
        ref, alt = [int(x) for x in ad.split(",")[:2]]
        tot = ref + alt
        return round(alt / tot, 3) if tot else None
    except Exception:
        return None


def simple_vcf_iter(path: str) -> Iterator[Dict[str, object]]:
    """Tiny single-sample VCF iterator. Accepts tab OR space-delimited toy files."""
    with open_text_auto(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                # header split on any whitespace to be robust to toy files
                _header = re.split(r"\s+", line.strip())
                continue
            if not line.strip():
                continue
            # data line
            fields = re.split(r"\s+", line.strip())
            if len(fields) < 8:
                # skip malformed lines
                continue
            chrom, pos, _id, ref, alt, qual, flt, info = fields[:8]
            fmt = fields[8] if len(fields) > 8 else ""
            sample = fields[9] if len(fields) > 9 else ""
            fmt_map: Dict[str, str] = {}
            if fmt and sample:
                keys = fmt.split(":")
                vals = sample.split(":")
                fmt_map = {k: (vals[i] if i < len(vals) else "") for i, k in enumerate(keys)}
            yield {
                "CHROM": chrom,
                "POS": int(pos),
                "REF": ref,
                "ALT": alt.split(","),  # list
                "QUAL": None if qual in (".", "") else float(qual),
                "FILTER": "PASS" if flt == "." else flt,
                "INFO": parse_info(info),
                "FORMAT": fmt_map,
            }


def triage(
    vcf_path: str,
    out_csv: str,
    min_dp: int = 10,
    min_qual: float = 30.0,
    max_af: float = 0.01,
    include_nonpass: bool = False,
    genes_set: Optional[Set[str]] = None,
    summary_path: Optional[str] = None,
) -> bool:
    """Core filtering + CSV writer (with optional gene allowlist and summary)."""
    cols = ["chrom", "pos", "ref", "alt", "gene", "consequence", "hgvs_c", "hgvs_p", "af", "gt", "dp", "ab", "filters"]
    kept = 0
    by_consequence: Dict[str, int] = {}

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)

        for rec in simple_vcf_iter(vcf_path):
            flt = str(rec["FILTER"])
            if not include_nonpass and flt != "PASS":
                continue

            qual = float(rec["QUAL"] or 0.0)
            if qual < min_qual:
                continue

            info = rec["INFO"]  # type: ignore
            af_raw = info.get("AF")
            af = float(af_raw) if af_raw not in (None, "") else None
            if af is not None and af >= max_af:
                continue

            fmt = rec["FORMAT"]  # type: ignore
            dp = int(fmt.get("DP", 0) or 0)
            if dp < min_dp:
                continue

            gt = fmt.get("GT", "")
            ab = allele_balance(fmt)
            if gt in ("0/1", "1/0") and ab is not None and not (0.3 <= ab <= 0.7):
                continue

            gene, consequence, hgvsc, hgvsp = parse_ann(info)  # type: ignore

            # Gene allowlist—keep only if gene in provided set
            if genes_set is not None:
                if not gene or gene not in genes_set:
                    continue

            for alt in rec["ALT"]:  # type: ignore
                w.writerow(
                    [
                        rec["CHROM"],
                        rec["POS"],
                        rec["REF"],
                        alt,
                        gene,
                        consequence,
                        hgvsc,
                        hgvsp,
                        af if af is not None else "",
                        gt,
                        dp,
                        ab if ab is not None else "",
                        flt,
                    ]
                )
                kept += 1
                if consequence:
                    by_consequence[consequence] = by_consequence.get(consequence, 0) + 1

    if summary_path:
        with open(summary_path, "w") as s:
            s.write(f"Kept variants: {kept}\n")
            if by_consequence:
                s.write("By consequence:\n")
                for k, v in sorted(by_consequence.items(), key=lambda x: (-x[1], x[0])):
                    s.write(f"  {k}: {v}\n")

    return True


def main():
    ap = argparse.ArgumentParser(description="Filter and summarise a VCF into CSV")
    ap.add_argument("--vcf", required=True, help="Path to input .vcf or .vcf.gz")
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument("--min-dp", type=int, default=10, help="Minimum depth (DP) to keep a site")
    ap.add_argument("--min-qual", type=float, default=30.0, help="Minimum QUAL to keep a site")
    ap.add_argument("--max-af", type=float, default=0.01, help="Maximum allele frequency (AF) if present in INFO")
    ap.add_argument("--include-nonpass", action="store_true", help="Include sites where FILTER != PASS")
    ap.add_argument("--genes", help="File of genes to keep (one per line)")
    ap.add_argument("--summary", help="Write a small text summary to this path")
    args = ap.parse_args()

    genes_set: Optional[Set[str]] = None
    if args.genes:
        with open(args.genes) as gf:
            genes_set = {line.strip() for line in gf if line.strip()}

    triage(
        args.vcf,
        args.out,
        args.min_dp,
        args.min_qual,
        args.max_af,
        args.include_nonpass,
        genes_set,
        args.summary,
    )


if __name__ == "__main__":
    main()

