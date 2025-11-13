import pathlib
from src.vcf_triage import triage

def test_summary_written(tmp_path: pathlib.Path):
    out = tmp_path/"out.csv"
    summary = tmp_path/"summary.txt"
    triage("data/sample.vcf", str(out), min_dp=10, min_qual=30, summary_path=str(summary))
    txt = summary.read_text()
    assert "Kept variants: 1" in txt
    assert "missense_variant" in txt  # from sample data

