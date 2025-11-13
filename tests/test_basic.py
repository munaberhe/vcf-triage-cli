import csv, pathlib
from src.vcf_triage import triage

def test_filters(tmp_path: pathlib.Path):
    out = tmp_path/"out.csv"
    triage("data/sample.vcf", str(out), min_dp=10, min_qual=30)
    rows = list(csv.DictReader(open(out)))
    poses = {r["pos"] for r in rows}
    assert "1000" in poses      # kept
    assert "2000" not in poses  # low DP/QUAL
    assert "3000" not in poses  # not PASS
    assert "4000" not in poses  # AF too high

