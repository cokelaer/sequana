import pandas as pd

import sequana.repeats.mirror as mirror_module
from sequana.repeats.mirror import MirrorRepeats
from sequana.tools import reverse_complement

from . import test_dir

MEASLES = f"{test_dir}/data/fasta/measles.fa"


def _write_fasta(tmpdir, seq, name="planted.fa"):
    path = tmpdir.join(name)
    path.write(f">s\n{seq}\n")
    return str(path)


def test_mirror_columns_and_invariants():
    m = MirrorRepeats(MEASLES)
    m.run(progress=False)
    assert list(m.df.columns) == ["seqid", "start", "end", "length", "repeat", "spacer", "subset", "sequence"]
    seq = open(MEASLES).read().split("\n", 1)[1].replace("\n", "").upper()
    for row in m.df.itertuples(index=False):
        rep, spacer = row.repeat, row.spacer
        assert row.end - row.start == 2 * rep + spacer
        assert rep >= m.min_repeat
        assert m.min_spacer <= spacer <= m.max_spacer
        # right arm is the reverse (NOT complement) of the left arm
        left = seq[row.start : row.start + rep]
        right = seq[row.end - rep : row.end]
        assert left == right[::-1]


def test_mirror_planted(tmpdir):
    arm = "ACGTTGCAAC"  # mixed bases -> not a triplex, still a mirror repeat
    motif = arm + "GG" + arm[::-1]
    fasta = _write_fasta(tmpdir, "N" * 6 + motif + "N" * 6)
    m = MirrorRepeats(fasta)
    m.run(progress=False)
    assert not m.df.empty
    hit = m.df[(m.df.start == 6) & (m.df.end == 6 + len(motif))]
    assert not hit.empty
    assert hit["repeat"].max() >= len(arm)


def test_mirror_ignores_inverted_repeat(tmpdir):
    # a (non-palindromic) inverted repeat is reverse-complement, not a mirror
    arm = "ACGTTGCAAC"
    motif = arm + "GG" + reverse_complement(arm)
    fasta = _write_fasta(tmpdir, "N" * 6 + motif + "N" * 6)
    m = MirrorRepeats(fasta)
    m.run(progress=False)
    assert m.df.empty


def test_mirror_includes_homopurine(tmpdir):
    # H-DNA triplex motifs are a subset of mirror repeats
    arm = "AGGAAGGAGA"
    motif = arm + "AA" + arm[::-1]
    fasta = _write_fasta(tmpdir, "N" * 6 + motif + "N" * 6)
    m = MirrorRepeats(fasta)
    m.run(progress=False)
    assert not m.df.empty


def test_mirror_n_does_not_extend(tmpdir):
    fasta = _write_fasta(tmpdir, "N" * 60)
    m = MirrorRepeats(fasta)
    m.run(progress=False)
    assert m.df.empty


def test_mirror_numba_fallback_parity(tmpdir):
    arm = "ACGTTGCAAC"
    fasta = _write_fasta(tmpdir, "N" * 6 + arm + "GG" + arm[::-1] + "N" * 6)
    m1 = MirrorRepeats(fasta)
    m1.run(progress=False)
    orig = mirror_module._HAS_NUMBA
    try:
        mirror_module._HAS_NUMBA = False
        m2 = MirrorRepeats(fasta)
        m2.run(progress=False)
    finally:
        mirror_module._HAS_NUMBA = orig
    pd.testing.assert_frame_equal(m1.df, m2.df)


def test_mirror_to_gff(tmpdir):
    arm = "ACGTTGCAAC"
    fasta = _write_fasta(tmpdir, "N" * 6 + arm + "GG" + arm[::-1] + "N" * 6)
    m = MirrorRepeats(fasta)
    m.run(progress=False)
    out = tmpdir.join("mr.gff")
    m.to_gff(out)
    lines = [l for l in open(out) if not l.startswith("#")]
    assert lines
    fields = lines[0].split("\t")
    assert fields[2] == "Mirror_Repeat"
    assert int(fields[3]) == int(m.df.iloc[0].start) + 1
