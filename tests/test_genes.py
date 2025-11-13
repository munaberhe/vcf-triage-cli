import csv, pathlib
from src.vcf_triage import triage

def test_gene_filter(tmp_path: pathlib.Path):
    out = tmp_path/"out.csv"
    genes = {"GENE1"}
    triage("data/sample.vcf", str(out), min_dp=10, min_qual=30, genes_set=genes)
    rows = list(csv.DictReader(open(out)))
    poses = {r["pos"] for r in rows}
    assert poses == {"1000"}  # only the GENE1 missense survives

