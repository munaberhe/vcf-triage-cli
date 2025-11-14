"""
IGV batch generator from triage CSV.

Usage examples:
  # Single BAM for all rows
  python -m src.igv_batch --csv out.csv --bam /path/to/sample.bam --out igv_batch.txt

  # Per-row BAM/CRAM path stored in a CSV column (e.g., "bam_path")
  python -m src.igv_batch --csv out.csv --bam-col bam_path --out igv_batch.txt

Then in IGV:
  1) Open IGV (desktop)
  2) Tools -> Run Batch Script -> select igv_batch.txt
IGV will load the file(s), navigate to each locus, and (if snapshot commands exist) save images.

Notes:
- This script just writes an IGV batch file; it does not run IGV itself.
- Your triage CSV must include at least: chrom,pos. Optional: gene, consequence, etc.
"""

from __future__ import annotations
import argparse
import csv
import os
from pathlib import Path

def sanitize_filename(s: str) -> str:
    # Keep it conservative for filesystem safety
    bad = '<>:"/\\|?* \t'
    out = "".join(ch if ch not in bad else "_" for ch in s)
    # limit length a bit
    return out[:80] if len(out) > 80 else out

def make_igv_batch(
    csv_path: str,
    out_path: str = "igv_batch.txt",
    genome: str = "GRCh38",
    bam: str | None = None,
    bam_col: str | None = None,
    flank: int = 100,
    snapshot_dir: str | None = None,
    snapshot_prefix: str | None = None,
) -> int:
    """
    Returns number of loci written.
    """
    out_lines: list[str] = []
    out_lines.append(f"new")
    out_lines.append(f"genome {genome}")

    rows = list(csv.DictReader(open(csv_path, newline="")))
    if not rows:
        Path(out_path).write_text("\n".join(out_lines) + "\n")
        return 0

    # If a single BAM is provided, load once
    if bam:
        out_lines.append(f"load {bam}")

    # If snapshots requested, set directory
    if snapshot_dir:
        out_lines.append(f"snapshotDirectory {snapshot_dir}")

    count = 0
    for r in rows:
        chrom = (r.get("chrom") or r.get("CHROM") or "").strip()
        pos_str = (r.get("pos") or r.get("POS") or "").strip()
        if not chrom or not pos_str.isdigit():
            # skip malformed line
            continue
        pos = int(pos_str)
        start = max(1, pos - flank)
        end = pos + flank

        # If BAM per row: load it before goto (IGV tolerates multiple loads)
        row_bam_path = None
        if bam_col:
            row_bam_path = (r.get(bam_col) or "").strip()
            if row_bam_path:
                out_lines.append(f"load {row_bam_path}")

        out_lines.append(f"goto {chrom}:{start}-{end}")

        # Create a snapshot line if requested
        if snapshot_dir is not None:
            gene = (r.get("gene") or r.get("GENE") or "").strip()
            consequence = (r.get("consequence") or r.get("Consequences") or r.get("CONSEQUENCE") or "").strip()
            label_parts = [chrom, str(pos)]
            if gene:
                label_parts.append(gene)
            if consequence:
                # keep short
                label_parts.append(consequence.split("&")[0][:24])
            base_name = "_".join(label_parts)
            if snapshot_prefix:
                base_name = f"{snapshot_prefix}_{base_name}"
            fname = sanitize_filename(base_name) + ".png"
            out_lines.append(f"snapshot {fname}")

        count += 1

    # Optionally sort and collapse tracks to a tidy view
    # (You can adjust or remove these lines to taste)
    # out_lines.insert(2, "sort base")
    # out_lines.insert(3, "collapse")

    Path(out_path).write_text("\n".join(out_lines) + "\n")
    return count

def main():
    ap = argparse.ArgumentParser(description="Generate an IGV batch script from a triage CSV.")
    ap.add_argument("--csv", required=True, help="Path to triage CSV (must have chrom,pos)")
    ap.add_argument("--out", default="igv_batch.txt", help="Output IGV batch file")
    ap.add_argument("--genome", default="GRCh38", help="Genome ID to use in IGV (e.g., GRCh37, GRCh38)")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--bam", help="Single BAM/CRAM to load for all loci")
    group.add_argument("--bam-col", help="CSV column name that contains BAM/CRAM path per row")
    ap.add_argument("--flank", type=int, default=100, help="bp padding around POS")
    ap.add_argument("--snapshot-dir", help="Directory for IGV snapshots (IGV will create and write)")
    ap.add_argument("--snapshot-prefix", help="Prefix for snapshot filenames")

    args = ap.parse_args()

    # Create snapshot dir if specified (IGV will also handle this typically, but safe to create)
    if args.snapshot_dir and not os.path.exists(args.snapshot_dir):
        os.makedirs(args.snapshot_dir, exist_ok=True)

    n = make_igv_batch(
        csv_path=args.csv,
        out_path=args.out,
        genome=args.genome,
        bam=args.bam,
        bam_col=args.bam_col,
        flank=args.flank,
        snapshot_dir=args.snapshot_dir,
        snapshot_prefix=args.snapshot_prefix,
    )
    print(f"Wrote {n} loci to {args.out}")

if __name__ == "__main__":
    main()

