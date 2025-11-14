import csv
from pathlib import Path
from src.igv_batch import make_igv_batch

def test_make_batch_single_bam(tmp_path: Path):
    # create tiny CSV
    csv_path = tmp_path / "triage.csv"
    csv_path.write_text("chrom,pos,gene,consequence\n1,1000,GENE1,missense_variant\n")
    out = tmp_path / "igv_batch.txt"
    n = make_igv_batch(str(csv_path), str(out), genome="GRCh38", bam="/x/sample.bam", flank=50, snapshot_dir=str(tmp_path/"shots"), snapshot_prefix="demo")
    assert n == 1
    txt = out.read_text()
    assert "genome GRCh38" in txt
    assert "load /x/sample.bam" in txt
    assert "goto 1:950-1050" in txt
    assert "snapshotDirectory" in txt
    assert "snapshot demo_1_1000_GENE1" in txt

def test_make_batch_bam_col(tmp_path: Path):
    # CSV with per-row bam path
    csv_path = tmp_path / "triage.csv"
    csv_path.write_text("chrom,pos,gene,consequence,bam_path\n1,2000,GENE2,synonymous_variant,/x/sample2.cram\n")
    out = tmp_path / "igv_batch.txt"
    n = make_igv_batch(str(csv_path), str(out), genome="GRCh38", bam=None, bam_col="bam_path", flank=25, snapshot_dir=None)
    assert n == 1
    txt = out.read_text()
    assert "load /x/sample2.cram" in txt
    assert "goto 1:1975-2025" in txt
    # No snapshot lines because snapshot_dir not set
    assert "snapshot " not in txt

