import pandas as pd
import pytest

import sequana.repeats.gquad as gquad_module
from sequana.repeats.gquad import GQuadruplex

from . import test_dir

MEASLES = f"{test_dir}/data/fasta/measles.fa"


def _write_fasta(tmpdir, seq, name="planted.fa"):
    path = tmpdir.join(name)
    path.write(f">s\n{seq}\n")
    return str(path)


def test_gq_columns_and_invariants():
    g = GQuadruplex(MEASLES)
    g.run(progress=False)
    assert list(g.df.columns) == ["seqid", "start", "end", "strand", "islands", "runs", "max", "sequence"]
    for row in g.df.itertuples(index=False):
        assert row.end - row.start == len(row.sequence)
        assert row.runs >= 4
        assert row.max >= g.min_rep
        assert row.islands >= 1
        assert row.strand in ("+", "-")
        # reported sequence is always G-rich (rev-comp taken for minus strand)
        assert row.sequence.upper().count("G") >= row.sequence.upper().count("C")


def test_gq_planted_plus(tmpdir):
    fasta = _write_fasta(tmpdir, "N" * 4 + "GGGTTGGGTTGGGTTGGG" + "N" * 4)
    g = GQuadruplex(fasta)
    g.run(progress=False)
    assert len(g.df) == 1
    r = g.df.iloc[0]
    assert (r["start"], r["end"], r["strand"]) == (4, 22, "+")
    assert (r["islands"], r["runs"], r["max"]) == (4, 4, 3)
    assert r["sequence"] == "GGGTTGGGTTGGGTTGGG"


def test_gq_planted_minus(tmpdir):
    # C-rich on forward strand -> minus-strand GQ, sequence is the rev-comp (G-rich)
    fasta = _write_fasta(tmpdir, "N" * 4 + "CCCAACCCAACCCAACCC" + "N" * 4)
    g = GQuadruplex(fasta)
    g.run(progress=False)
    assert len(g.df) == 1
    r = g.df.iloc[0]
    assert r["strand"] == "-"
    assert r["sequence"] == "GGGTTGGGTTGGGTTGGG"
    assert (r["islands"], r["runs"], r["max"]) == (4, 4, 3)


def test_gq_needs_four_runs(tmpdir):
    # only 3 islands -> npos = 3 < 4 -> nothing
    fasta = _write_fasta(tmpdir, "N" * 4 + "GGGTTGGGTTGGG" + "N" * 4)
    g = GQuadruplex(fasta)
    g.run(progress=False)
    assert g.df.empty


def test_gq_max_above_min(tmpdir):
    # two len-9 islands: 4 runs of size 4 fit across them -> max = 4
    fasta = _write_fasta(tmpdir, "N" * 4 + "G" * 9 + "TT" + "G" * 9 + "N" * 4)
    g = GQuadruplex(fasta)
    g.run(progress=False)
    assert len(g.df) == 1
    r = g.df.iloc[0]
    assert (r["islands"], r["runs"], r["max"]) == (2, 4, 4)


def test_gq_spacer_too_large(tmpdir):
    # islands separated by 8 (> max_spacer 7) are not merged -> no single 4-island motif
    fasta = _write_fasta(tmpdir, "N" * 4 + "GGG" + "T" * 8 + "GGG" + "T" * 8 + "GGG" + "T" * 8 + "GGG" + "N" * 4)
    g = GQuadruplex(fasta)
    g.run(progress=False)
    assert g.df.empty


def test_gq_numba_fallback_parity():
    g1 = GQuadruplex(MEASLES)
    g1.run(progress=False)
    orig = gquad_module._HAS_NUMBA
    try:
        gquad_module._HAS_NUMBA = False
        g2 = GQuadruplex(MEASLES)
        g2.run(progress=False)
    finally:
        gquad_module._HAS_NUMBA = orig
    pd.testing.assert_frame_equal(g1.df, g2.df)


def test_gq_to_gff_and_bed(tmpdir):
    fasta = _write_fasta(tmpdir, "N" * 4 + "GGGTTGGGTTGGGTTGGG" + "N" * 4)
    g = GQuadruplex(fasta)
    g.run(progress=False)
    gff = tmpdir.join("gq.gff")
    g.to_gff(gff)
    lines = [l for l in open(gff) if not l.startswith("#")]
    assert lines
    fields = lines[0].split("\t")
    assert fields[2] == "G_Quadruplex_Motif"
    assert int(fields[3]) == int(g.df.iloc[0].start) + 1
    assert ";islands=" in fields[8] and ";runs=" in fields[8] and ";max=" in fields[8]
    bed = tmpdir.join("gq.bed")
    g.to_bed(bed)
    assert bed.check()


def test_gq_to_bed_empty_raises(tmpdir):
    g = GQuadruplex(_write_fasta(tmpdir, "N" * 50))
    g.run(progress=False)
    with pytest.raises(ValueError):
        g.to_bed(str(tmpdir.join("empty.bed")))
