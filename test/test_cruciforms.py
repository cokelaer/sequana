import pandas as pd

import sequana.repeats.cruciforms as cruciforms_module
from sequana.repeats.cruciforms import Cruciforms
from sequana.tools import reverse_complement

from . import test_dir

MEASLES = f"{test_dir}/data/fasta/measles.fa"


def _write_fasta(tmpdir, seq, name="planted.fa"):
    path = tmpdir.join(name)
    path.write(f">s\n{seq}\n")
    return str(path)


def test_crux_bed(tmpdir):
    # smoke test: run + BED export on a real genome
    c = Cruciforms(MEASLES)
    c.run(progress=False)
    out = tmpdir.join("test.bed")
    c.to_bed(out)
    assert out.check()


def test_cruciforms_columns_and_invariants():
    c = Cruciforms(MEASLES)
    c.run(progress=False)
    assert list(c.df.columns) == ["seqid", "start", "end", "length", "stem", "spacer", "subset", "sequence"]
    seq = open(MEASLES).read().split("\n", 1)[1].replace("\n", "").upper()
    for row in c.df.itertuples(index=False):
        stem, spacer = row.stem, row.spacer
        # geometry is consistent
        assert row.end - row.start == 2 * stem + spacer
        assert stem >= c.min_stem_len
        assert c.min_spacer <= spacer <= c.max_spacer
        # the two arms are reverse complements of each other (inverted repeat)
        left = seq[row.start : row.start + stem]
        right = seq[row.end - stem : row.end]
        assert left == reverse_complement(right)
        # short-loop / long-loop coupling
        assert stem >= (c.min_stem_len if spacer <= c.short_spacer_max else c.min_stem_len_long)


def test_cruciforms_planted(tmpdir):
    # an inverted repeat isolated by N flanks (N never pairs) => detected
    arm = "ACGGTTAC"
    motif = arm + "AT" + reverse_complement(arm)
    fasta = _write_fasta(tmpdir, "N" * 6 + motif + "N" * 6)
    c = Cruciforms(fasta)
    c.run(progress=False)
    assert not c.df.empty
    # a hit spanning exactly the planted motif exists
    hit = c.df[(c.df.start == 6) & (c.df.end == 6 + len(motif))]
    assert not hit.empty
    assert hit["stem"].max() >= len(arm)


def test_cruciforms_ignores_pure_mirror(tmpdir):
    # a homopurine mirror repeat (arm + reverse(arm)) is NOT an inverted repeat
    arm = "AGGAAGGAGA"
    motif = arm + "AA" + arm[::-1]
    fasta = _write_fasta(tmpdir, "N" * 6 + motif + "N" * 6)
    c = Cruciforms(fasta)
    c.run(progress=False)
    assert c.df.empty


def test_cruciforms_numba_fallback_parity(tmpdir):
    fasta = _write_fasta(tmpdir, "N" * 6 + "ACGGTTAC" + "AT" + reverse_complement("ACGGTTAC") + "N" * 6)
    c1 = Cruciforms(fasta)
    c1.run(progress=False)
    orig = cruciforms_module._HAS_NUMBA
    try:
        cruciforms_module._HAS_NUMBA = False
        c2 = Cruciforms(fasta)
        c2.run(progress=False)
    finally:
        cruciforms_module._HAS_NUMBA = orig
    pd.testing.assert_frame_equal(c1.df, c2.df)


def test_cruciforms_to_gff(tmpdir):
    fasta = _write_fasta(tmpdir, "N" * 6 + "ACGGTTAC" + "AT" + reverse_complement("ACGGTTAC") + "N" * 6)
    c = Cruciforms(fasta)
    c.run(progress=False)
    out = tmpdir.join("ir.gff")
    c.to_gff(out)
    lines = [l for l in open(out) if not l.startswith("#")]
    assert lines
    fields = lines[0].split("\t")
    assert fields[2] == "Inverted_Repeat"
    # gff is 1-based: start column == df start + 1
    assert int(fields[3]) == int(c.df.iloc[0].start) + 1
